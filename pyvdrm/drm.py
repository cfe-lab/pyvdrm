"""Asi2 specification
"""

from abc import ABCMeta, abstractmethod


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


class AsiExpr(object):
    """A callable ASI2 expression"""

    children = []
    label = None

    def __init__(self, _label, _pos, tokens):
        """By default we assume the head of the arg list is the operation"""

        self.typecheck(tokens.asList())
        self.children = tokens

        if not self.label:
            self.label = str(type(self))

    def typecheck(self, tokens):
        """Override typecheck method to define runtime errors"""
        pass

    def __call__(self, args):
        """Evaluate child tokens with args"""
        return self.children(args)


class AsiBinaryExpr(AsiExpr):
    """Subclass with syntactic sugar for boolean ops"""

    def __init__(self, label, pos, tokens):
        super().__init__(label, pos, tokens)
        self.children = tokens[0]

    def typecheck(self, tokens):
        if len(tokens[0]) != 2:
            pass
            #raise AsiParseError

    def __repr__(self):
        arg1, arg2 = self.children
        return "{} {} {}".format(arg1, type(self), arg2)


class AsiUnaryExpr(AsiExpr):
    """Subclass for atoms and unary ops"""

    def typecheck(self, tokens):
        if isinstance(tokens[0], list):
            raise AsiParseError

    def __repr__(self):
        return str(self.children)
