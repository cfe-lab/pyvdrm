"""
HCV Drug Resistance Rule Parser definition, version 2
"""

from functools import reduce, total_ordering
from pyparsing import (Literal, nums, Word, Forward, Optional, Regex, Combine,
                       infixNotation, delimitedList, opAssoc, ParseException)

from pyvdrm.drm import MissingPositionError
from pyvdrm.drm import AsiExpr, AsiBinaryExpr, DRMParser
from pyvdrm.vcf import MutationSet


class BoolTrue(AsiExpr):
    """Boolean True constant"""
    def __call__(self, *args):
        return True


class BoolFalse(AsiExpr):
    """Boolean False constant"""
    def __call__(self, *args):
        return False


class AndExpr(AsiBinaryExpr):
    """Boolean AND on children"""

    def __call__(self, mutations):
        larg, rarg = self.children
        return larg(mutations) and rarg(mutations)


class OrExpr(AsiBinaryExpr):
    """Boolean OR on children (binary only)"""

    def __call__(self, mutations):
        larg, rarg = self.children
        return larg or rarg


class ScoreExpr(AsiExpr):
    """Score expressions propagate DRM scores"""

    def __call__(self, mutations):
        child = self.children[0]
        if type(child) == ExprList:
            return child(mutations)
        else:
            return float(child)


class ExprList(AsiExpr):

    def __call__(self, mutations):
        operation, *rest = self.children
        if operation == 'MAX':
            terms = rest
            func = max
        elif operation == 'MIN':
            terms = rest
            func = min
        elif operation == 'MEAN':
            terms = rest
            func = lambda x: sum(x) / len(x)
        else:
            # the default operation is sum
            terms = self.children
            func = sum

        scores = [f(mutations) for f in terms]

        return func(scores)


class Condition(AsiExpr):
    def __call__(self, mutations):
        if len(self.children) == 1:
            child = self.children[0]
            return child(mutations)


class Expr(AsiExpr):
    def __call__(self, mutations):
        if len(self.children) == 1:
            return float(self.children[0])

        condition, scoreExpr = self.children
        if condition(mutations):
            return scoreExpr(mutations)
        else:
            return 0.0

class AsiMutations(object):
    """List of mutations given an ambiguous pattern"""

    def __init__(self, _label=None, _pos=None, args=None):
        """Initialize set of mutations from a potentially ambiguous residue
        """
        self.mutations = MutationSet(''.join(args))

    def __repr__(self):
        return "AsiMutations(args={!r})".format(str(self.mutations))

    def __call__(self, env):
        is_found = False
        for mutation_set in env:
            is_found |= mutation_set.pos == self.mutations.pos
            intersection = self.mutations.mutations & mutation_set.mutations
            if len(intersection) > 0:
                return True

        if not is_found:
            # Some required positions were not found in the environment.
            raise MissingPositionError('Missing position {}.'.format(
                self.mutations.pos))
        return False


class HCVR2(DRMParser):
    """HCV Resistance Syntax definition, version 2"""

    def parser(self, rule):
        # literals
        max_ = Literal('MAX')
        min_ = Literal('MIN')
        mean = Literal('MEAN')

        and_ = Literal('AND').suppress()
        or_ = Literal('OR').suppress()

        l_par = Literal('(').suppress()
        r_par = Literal(')').suppress()

        mapper = Literal('=>').suppress()

        # complex atoms
        integer = Word(nums)
        residue = Optional(Regex(r'[A-Z]')) + integer + Regex(r'\!?[diA-Z]+')
        residue.setParseAction(AsiMutations)
        
        float_ = Combine(Optional(Literal('-')) + integer +\
                         Optional(Literal('.') + Optional(integer)))

        accumulator = max_ | min_ | mean

        bool_ = (Literal('TRUE').suppress() |
                 Literal('FALSE').suppress())

        # compound expressions
        booleancondition = Forward()

        condition = residue | bool_ | l_par + booleancondition + r_par
        condition.setParseAction(Condition)

        booleancondition << infixNotation(condition,
                                          [(and_, 2, opAssoc.LEFT, AndExpr),
                                           (or_, 2, opAssoc.LEFT, OrExpr)])

        expr_list = Forward()
        score = float_ | expr_list
        score.setParseAction(ScoreExpr)

        expr = condition + mapper + score | float_
        expr.setParseAction(Expr)
        expr_list << Optional(accumulator) + l_par + delimitedList(expr) + r_par
        expr_list.setParseAction(ExprList)

        try:
            return score.parseString(rule)
        except ParseException as ex:
            ex.msg = 'Error in HCVR2: ' + ex.markInputline()
            raise

    def __call__(self, mutations):
        return self.dtree(mutations)
