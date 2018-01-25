"""
HCV Drug Resistance Rule Parser definition
"""

from functools import reduce, total_ordering
from pyparsing import (Literal, nums, Word, Forward, Optional, Regex,
                       infixNotation, delimitedList, opAssoc, ParseException)

from pyvdrm.drm import MissingPositionError
from pyvdrm.drm import AsiExpr, AsiBinaryExpr, DRMParser
from pyvdrm.vcf import MutationSet


def update_flags(fst, snd):
    for k in snd:
        if k in fst:
            fst[k].append(snd[k])
        else:
            fst[k] = snd[k]  # this could be achieved with a defaultdict
    return fst


@total_ordering
class Score(object):
    """Encapsulate a score and the residues that support it"""

    def __init__(self, score, residues, flags=None):
        """ Initialize.

        :param bool|float score: value of the score
        :param residues: sequence of Mutations
        :param flags: dictionary of user defined strings and supporting Mutations
        """
        self.score = score
        self.residues = set(residues)

        if flags is None:
            self.flags = {}
        else:
            self.flags = flags

    def __add__(self, other):
        flags = update_flags(self.flags, other.flags)
        return Score(self.score + other.score,
                     self.residues | other.residues,
                     flags)

    def __sub__(self, other):
        flags = update_flags(self.flags, other.flags)
        return Score(self.score - other.score,
                     self.residues | other.residues,
                     flags)

    def __repr__(self):
        return "Score({!r}, {!r})".format(self.score, self.residues)

    def __eq__(self, other):
        return self.score == other.score

    def __lt__(self, other):
        # the total_ordering decorator populates the other 5 comparison
        # operations. Implement them explicitly if this causes performance
        # issues
        return self.score < other.score

    def __bool__(self):
        return self.score


class BoolTrue(AsiExpr):
    """Boolean True constant"""
    def __call__(self, *args):
        return Score(True, [])


class BoolFalse(AsiExpr):
    """Boolean False constant"""
    def __call__(self, *args):
        return Score(False, [])


class AndExpr(AsiExpr):
    """Fold boolean AND on children"""

    def __call__(self, mutations):
        scores = map(lambda f: f(mutations), self.children[0])
        scores = [Score(False, []) if s is None else s for s in scores]
        if not scores:
            raise ValueError

        residues = set()

        for s in scores:
            residues |= s.residues
            if not s.score:
                return Score(False, [])

        return Score(True, residues)


class OrExpr(AsiBinaryExpr):
    """Boolean OR on children (binary only)"""

    def __call__(self, mutations):
        arg1, arg2 = self.children

        score1 = arg1(mutations)
        score2 = arg2(mutations)

        if score1 is None:
            score1 = Score(False, [])
        if score2 is None:
            score2 = Score(False, [])

        return Score(score1.score or score2.score,
                     score1.residues | score2.residues)


class EqualityExpr(AsiExpr):
    """ASI2 style inequality expressions"""

    def __init__(self, label, pos, children):
        super().__init__(label, pos, children)
        self.operation, limit = children
        self.limit = int(limit)

    def __call__(self, x):
        if self.operation == 'ATLEAST':
            return x >= self.limit
        elif self.operation == 'EXACTLY':
            return x == self.limit
        elif self.operation == 'NOMORETHAN':
            return x <= self.limit

        raise NotImplementedError


class ScoreExpr(AsiExpr):
    """Score expressions propagate DRM scores"""

    def __call__(self, mutations):

        flags = {}
        if len(self.children) == 4:
            operation, _, flag, _ = self.children
            flags[flag] = []
            score = 0  # should be None

        elif len(self.children) == 3:
            operation, minus, score = self.children
            if minus != '-':  # this is parsing the expression twice, refactor
                raise ValueError
            score = -1 * int(score)

        elif len(self.children) == 2:
            operation, score = self.children
            score = int(score)

        else:
            raise ValueError

        # evaluate operation and return score
        result = operation(mutations)
        if result is None:
            return None

        if result.score is False:
            return Score(0, [])
        return Score(score, result.residues, flags=flags)


class ScoreList(AsiExpr):
    """Lists of scores are either summed or maxed"""

    def __call__(self, mutations):
        operation, *rest = self.children
        if operation == 'MAX':
            terms = rest
            func = max
        else:
            # the default operation is sum
            terms = self.children
            func = sum
        scores = [f(mutations) for f in terms]
        matched_scores = [score.score for score in scores if score.score]
        residues = reduce(lambda x, y: x | y,
                          (score.residues for score in scores))
        flags = {}
        for score in scores:
            flags.update(score.flags)
        return Score(bool(matched_scores) and func(matched_scores),
                     residues,
                     flags)


class SelectFrom(AsiExpr):
    """Return True if some number of mutations match"""

    def typecheck(self, tokens):
        # if type(tokens[0]) != EqualityExpr:
        #     raise TypeError()
        pass

    def __call__(self, mutations):
        operation, *rest = self.children
        # the head of the arg list must be an equality expression
       
        scored = [f(mutations) for f in rest]
        passing = sum(bool(score.score) for score in scored)

        return Score(operation(passing),
                     reduce(lambda x, y: x | y,
                            (item.residues for item in scored)))


class AsiScoreCond(AsiExpr):
    """Score condition"""

    label = "ScoreCond"

    def __call__(self, args):
        """Score conditions evaluate a list of expressions and sum scores"""
        return sum((f(args) for f in self.children), Score(False, set()))


class AsiMutations(object):
    """List of mutations given an ambiguous pattern"""

    def __init__(self, _label=None, _pos=None, args=None):
        """Initialize set of mutations from a potentially ambiguous residue
        """
        self.mutations = MutationSet(''.join(args))

    def __repr__(self):
        if self.mutations is None:
            return "AsiMutations()"
        return "AsiMutations(args={!r})".format(str(self.mutations))

    def __call__(self, env):
        is_found = False
        for mutation_set in env:
            is_found |= mutation_set.pos == self.mutations.pos
            intersection = self.mutations.mutations & mutation_set.mutations
            if len(intersection) > 0:
                return Score(True, intersection)

        if not is_found:
            # Some required positions were not found in the environment.
            raise MissingPositionError('Missing position {}.'.format(
                self.mutations.pos))
        return Score(False, set())


class HCVR(DRMParser):
    """HCV Resistance Syntax definition"""

    def parser(self, rule):

        select = Literal('SELECT').suppress()
        except_ = Literal('EXCEPT')
        exactly = Literal('EXACTLY')
        atleast = Literal('ATLEAST')

        from_ = Literal('FROM').suppress()

        max_ = Literal('MAX')

        and_ = Literal('AND').suppress()
        or_ = Literal('OR').suppress()
        # min_ = Literal('MIN')

        notmorethan = Literal('NOTMORETHAN')
        l_par = Literal('(').suppress()
        r_par = Literal(')').suppress()

        quote = Literal('"')

        mapper = Literal('=>').suppress()
        integer = Word(nums)

        residue = Optional(Regex(r'[A-Z]')) + integer + Regex(r'\!?[diA-Z]+')
        residue.setParseAction(AsiMutations)

        # Syntax of expressions
        excludestatement = except_ + residue

        quantifier = exactly | atleast | notmorethan
        inequality = quantifier + integer
        inequality.setParseAction(EqualityExpr)

        select_quantifier = infixNotation(inequality,
                                          [(and_, 2, opAssoc.LEFT, AndExpr),
                                           (or_, 2, opAssoc.LEFT, OrExpr)])

        residue_list = l_par + delimitedList(residue) + r_par

        # so selectstatement.eval :: [Mutation] -> Maybe Bool
        selectstatement = select + select_quantifier + from_ + residue_list
        selectstatement.setParseAction(SelectFrom)

        bool_ = (Literal('TRUE').suppress().setParseAction(BoolTrue) |
                 Literal('FALSE').suppress().setParseAction(BoolFalse))

        booleancondition = Forward()
        condition = residue | excludestatement | selectstatement | bool_

        booleancondition << infixNotation(condition,
                                          [(and_, 2, opAssoc.LEFT, AndExpr),
                                           (or_, 2, opAssoc.LEFT, OrExpr)]) | condition

        score = Optional(Literal('-')) + integer | quote + Regex(r'[a-zA-Z0-9 _]+') + quote
        scoreitem = booleancondition + mapper + score
        scoreitem.setParseAction(ScoreExpr)
        scorelist = max_ + l_par + delimitedList(scoreitem) + r_par |\
            delimitedList(scoreitem)
        scorelist.setParseAction(ScoreList)

        scorecondition = Literal('SCORE FROM').suppress() +\
            l_par + delimitedList(scorelist) + r_par

        scorecondition.setParseAction(AsiScoreCond)

        statement = booleancondition | scorecondition

        try:
            return statement.parseString(rule)
        except ParseException as ex:
            ex.msg = 'Error in HCVR: ' + ex.markInputline()
            raise
