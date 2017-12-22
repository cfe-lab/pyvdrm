"""
ASI2 Parser definition
"""

from functools import reduce, total_ordering
from pyparsing import (Literal, nums, Word, Forward, Optional, Regex,
                       infixNotation, delimitedList, opAssoc)
from pyvdrm.drm import AsiExpr, AsiBinaryExpr, AsiUnaryExpr, DRMParser
from pyvdrm.vcf import MutationSet


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

    def __init__(self, score, residues):
        """ Initialize.

        :param bool|float score: value of the score
        :param residues: sequence of Mutations
        """
        self.score = score
        self.residues = set(residues)

    def __add__(self, other):
        return Score(self.score + other.score, self.residues | other.residues)

    def __sub__(self, other):
        return Score(self.score - other.score, self.residues | other.residues)

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


class Negate(AsiUnaryExpr):
    """Unary negation of boolean child"""
    def __call__(self, mutations):
        arg = self.children(mutations)
        return Score(not arg.score, arg.residues)


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
    """ASI2 inequality expressions"""

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
        if len(self.children) == 3:
            operation, minus, score = self.children
            if minus != '-':
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
        return Score(score, result.residues)


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
        intersection = set(env) & self.mutations.mutations
        if len(intersection) > 0:
            return Score(True, intersection)
        return None


class ASI2(DRMParser):
    """ASI2 Syntax definition"""

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
        mapper = Literal('=>').suppress()
        integer = Word(nums)

        mutation = Optional(Regex(r'[A-Z]')) + integer + Regex(r'[diA-Z]+')
        mutation.setParseAction(AsiMutations)

        not_ = Literal('NOT').suppress() + mutation
        not_.setParseAction(Negate)

        residue = mutation | not_
        # integer + l_par + not_ + Regex(r'[A-Z]+') + r_par
        # roll this next rule into the mutation object

        # Syntax of ASI expressions
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

        booleancondition = Forward()
        condition = residue | excludestatement | selectstatement

        booleancondition << infixNotation(condition,
                                          [(and_, 2, opAssoc.LEFT, AndExpr),
                                           (or_, 2, opAssoc.LEFT, OrExpr)]) | condition

        scoreitem = booleancondition + mapper + Optional(Literal('-')) + integer
        scoreitem.setParseAction(ScoreExpr)
        scorelist = max_ + l_par + delimitedList(scoreitem) + r_par |\
            delimitedList(scoreitem)
        scorelist.setParseAction(ScoreList)

        scorecondition = Literal('SCORE FROM').suppress() +\
            l_par + delimitedList(scorelist) + r_par

        scorecondition.setParseAction(AsiScoreCond)

        statement = booleancondition | scorecondition

        return statement.parseString(rule)
