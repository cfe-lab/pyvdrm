"""
HCV Drug Resistance Rule Parser definition
"""

from functools import reduce, total_ordering
from pyparsing import (Literal, nums, Word, Forward, Optional, Regex,
                       infixNotation, delimitedList, opAssoc, alphas)
from pyvdrm.drm import AsiExpr, AsiBinaryExpr, AsiUnaryExpr, DRMParser
from pyvdrm.vcf import MutationSet

def update_flags(fst, snd):
    for k in snd:
        if k in fst:
            fst[k].append(snd[k])
        else:
            fst[k] = snd[k] # this chould be achieved with a defaultdict
    return fst


def maybe_foldl(func, noneable):
    """Safely fold a function over a potentially empty list of
    potentially null values"""
    if noneable is None:
        return None
    clean = [x for x in noneable if x is not None]
    if not clean:
        return None
    return reduce(func, clean)


def maybe_map(func, noneable):
    if noneable is None:
        return None
    r_list = []
    for x in noneable:
        if x is None:
            continue
        result = func(x)
        if result is None:
            continue
        r_list.append(result)
    if not r_list:
        return None
    return r_list


@total_ordering
class Score(object):
    """Encapsulate a score and the residues that support it"""

    residues = set([])
    score = None
    flags = {} # allow a score expression to raise a user defined string

    def __init__(self, score, residues, flags={}):
        """ Initialize.

        :param bool|float score: value of the score
        :param residues: sequence of Mutations
        :param flags: dictionary of user defined strings and supporting Mutations
        """
        self.score = score
        self.residues = set(residues)
        self.flags = flags

    def __add__(self, other):
        flags = update_flags(self.flags, other.flags)
        return Score(self.score + other.score, self.residues | other.residues,
                flags)

    def __sub__(self, other):
        flags = update_flags(self.flags, other.flags)
        return Score(self.score - other.score, self.residues | other.residues,
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


class Negate(AsiExpr):
    """Unary negation of boolean child"""
    def __call__(self, mutations):
        child_score = self.children[0](mutations)
        if child_score is None:
            return Score(True, []) # TODO: propagate negative residues
        return Score(not child_score.score, child_score.residues)


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

        residues = set([])
        for s in scores:
            if not s.score:
                return Score(False, [])
            residues = residues | s.residues

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
            score = 0 # should be None

        elif len(self.children) == 3:
            operation, minus, score = self.children
            if minus != '-': # this is parsing the expression twice, refactor
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
            return maybe_foldl(max, [f(mutations) for f in rest])

        # the default operation is sum
        return maybe_foldl(lambda x, y: x+y, [f(mutations) for f in self.children])


class SelectFrom(AsiExpr):
    """Return True if some number of mutations match"""

    def typecheck(self, tokens):
        # if type(tokens[0]) != EqualityExpr:
        #     raise TypeError()
        pass

    def __call__(self, mutations):
        operation, *rest = self.children
        # the head of the arg list must be an equality expression
       
        scored = list(maybe_map(lambda f: f(mutations), rest))
        passing = len(scored) 

        if operation(passing):
            return Score(True, maybe_foldl(
                lambda x, y: x.residues.union(y.residues), scored))
        else:
            return None


class AsiScoreCond(AsiExpr):
    """Score condition"""

    label = "ScoreCond"

    def __call__(self, args):
        """Score conditions evaluate a list of expressions and sum scores"""
        return maybe_foldl(lambda x, y: x+y, map(lambda x: x(args), self.children))


class AsiMutations(object):
    """List of mutations given an ambiguous pattern"""

    def __init__(self, _label=None, _pos=None, args=None):
        """Initialize set of mutations from a potentially ambiguous residue
        """
        self.mutations = args and MutationSet(''.join(args))

    def __repr__(self):
        if self.mutations is None:
            return "AsiMutations()"
        return "AsiMutations(args={!r})".format(str(self.mutations))

    def __call__(self, env):
        for mutation_set in env:
            intersection = self.mutations.mutations & mutation_set.mutations
            if len(intersection) > 0:
                return Score(True, intersection)
        return None


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

        mutation = Optional(Regex(r'[A-Z]')) + integer + Regex(r'[diA-Z]+')
        mutation.setParseAction(AsiMutations)

        not_ = Literal('NOT').suppress() + mutation
        not_.setParseAction(Negate)

        residue = mutation | not_
        # integer + l_par + not_ + Regex(r'[A-Z]+') + r_par
        # roll this next rule into the mutation object

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

        bool_ = Literal('TRUE').suppress().setParseAction(BoolTrue) |\
                Literal('FALSE').suppress().setParseAction(BoolFalse)

        booleancondition = Forward()
        condition = residue | excludestatement | selectstatement | bool_

        booleancondition << infixNotation(condition,
                                          [(and_, 2, opAssoc.LEFT, AndExpr),
                                           (or_, 2, opAssoc.LEFT, OrExpr)]) | condition

        score = Optional(Literal('-')) + integer | quote + Word(alphas) + quote
        scoreitem = booleancondition + mapper + score
        scoreitem.setParseAction(ScoreExpr)
        scorelist = max_ + l_par + delimitedList(scoreitem) + r_par |\
            delimitedList(scoreitem)
        scorelist.setParseAction(ScoreList)

        scorecondition = Literal('SCORE FROM').suppress() +\
            l_par + delimitedList(scorelist) + r_par

        scorecondition.setParseAction(AsiScoreCond)

        statement = booleancondition | scorecondition

        return statement.parseString(rule)
