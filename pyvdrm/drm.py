"""Asi2 specification
"""

from functools import total_ordering
from abc import ABCMeta, abstractmethod

class Score(object):
    """Encapsulate a score and the residues that support it"""

    def __init__(self, score, residues, flags=None):
        """ Initialize.

        :param score: value of the score
        :param residues: sequence of Mutations
        :param flags: dictionary of user defined strings and supporting Mutations
        """
        self.score = score
        self.residues = set(residues)

        if flags is None:
            self.flags = {}
        else:
            self.flags = flags
    
    @staticmethod
    def update_flags(fst, snd):
        for k in snd:
            if k in fst:
                fst[k].append(snd[k])
            else:
                fst[k] = snd[k]  # this could be achieved with a defaultdict
        return fst

    def __repr__(self):
        return "Score({!r}, {!r})".format(#self.__class__.__name__,
                                         self.score,
                                         self.residues)


@total_ordering
class IntScore(Score):
    def __add__(self, other):
        flags = self.update_flags(self.flags, other.flags)

        # Here is where we define the arithmetic behaviour of flag terms
        result = self.score 
        if not other.score is None:
            result = self.score + other.score

        return IntScore(result,
                        self.residues | other.residues,
                        flags)

    def __sub__(self, other):
        flags = self.update_flags(self.flags, other.flags)

        result = self.score
        if not other.score is None:
            result = self.score - other.score

        return IntScore(result,
                        self.residues | other.residues,
                        flags)

    def __eq__(self, other):
        return self.score == other.score

    def __lt__(self, other):
        # the total_ordering decorator populates the other 5 comparison
        # operations. Implement them explicitly if this causes performance
        # issues

        if other.score is None:
            raise TypeError

        return self.score < other.score

    def __bool__(self):
        raise TypeError

    def __int__(self):
        return self.score


class BoolScore(Score):
    def __eq__(self, other):
        # is this used?
        return self.score == other.score

    def __and__(self, other):
        flags = self.update_flags(self.flags, other.flags)
        if self.score and other.score:
            return BoolScore(True, self.residues | other.residues, flags)
        return BoolScore(False, set())

    def __or__(self, other):
        flags = self.update_flags(self.flags, other.flags)
        if self.score or other.score:
            return BoolScore(True, self.residues | other.residues, flags)
        return BoolScore(False, set())

    def __int__(self):
        raise TypeError


class AsiParseError(Exception):
    pass


class MissingPositionError(Exception):
    pass


class DRMParser(metaclass=ABCMeta):
    """abstract class for DRM rule parsers/evaluators"""

    def __init__(self, rule):
        """drug resistance mutation callers are initialized with rule strings,
            the initialized parser has a callable decision tree object
        """
        self.rule = rule
        self.dtree, *rest = self.parser(rule)

    @abstractmethod
    def parser(self, rule_string):
        """The parser returns a decision tree based on the rule string"""
        pass

    def __call__(self, mutations):
        score = self.dtree(mutations)
        if score is None:
            return False
        return score.score

    def __repr__(self):
        return self.rule


class DrmExpr(object):
    """A callable ASI2 expression"""

    children = []
    label = None

    def __init__(self, _label, _pos, terms):
        """By default we assume the head of the arg list is the operation"""

        self.typecheck(terms.asList())
        self.children = terms

        if not self.label:
            self.label = str(type(self))

    def typecheck(self, terms):
        """Override typecheck method to define runtime errors"""
        pass

    def __call__(self, args):
        """Evaluate child tokens with args"""
        return self.children(args)


class DrmBinaryExpr(DrmExpr):
    """Subclass with syntactic sugar for boolean ops"""

    def __init__(self, label, pos, tokens):
        super().__init__(label, pos, tokens)
        self.children = tokens[0]

    def typecheck(self, tokens):
        if len(tokens[0]) != 2:
            raise AsiParseError

    def __repr__(self):
        arg1, arg2 = self.children
        return "{} {} {}".format(arg1, type(self), arg2)


class DrmUnaryExpr(DrmExpr):
    """Subclass for atoms and unary ops"""

    def typecheck(self, tokens):
        if isinstance(tokens[0], list):
            raise AsiParseError

    def __repr__(self):
        return str(self.children)
