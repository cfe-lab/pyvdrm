import os
import unittest

from pyparsing import ParseException

from pyvdrm.drm import MissingPositionError
from pyvdrm.hcvr2 import HCVR2, AsiMutations, Score
from pyvdrm.vcf import Mutation, MutationSet, VariantCalls

from pyvdrm.tests.test_vcf import add_mutations


# noinspection SqlNoDataSourceInspection,SqlDialectInspection
class TestRuleParser(unittest.TestCase):

    def test_anything(self):
        HCVR2("( 10N => -1.0 )")

    def test_expr(self):
        HCVR2("( 151M => -1.0, 130M => 2.0 )")

    def test_boolean_or(self):
        HCVR2("( ( 151M AND 69i ) => 1.0)")

    def test_nested_condition(self):
        rule = HCVR2("(L41L => 3, 67N => (2N => 2)))")
