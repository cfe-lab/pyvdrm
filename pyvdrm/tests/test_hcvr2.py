import os
import unittest

from pyparsing import ParseException

from pyvdrm.drm import MissingPositionError
from pyvdrm.hcvr2 import HCVR2, AsiMutations
from pyvdrm.vcf import Mutation, MutationSet, VariantCalls

from pyvdrm.tests.test_vcf import add_mutations

# noinspection SqlNoDataSourceInspection,SqlDialectInspection
class TestRuleParser(unittest.TestCase):

    def test_anything(self):
        r = HCVR2("( 10N => -1.0 )")
        self.assertEqual(r(VariantCalls("10N")), -1.0)
        self.assertEqual(r(VariantCalls("10X")), 0.0)

    def test_expr(self):
        r = HCVR2("( 151M => -1.0, 130M => 2.0 )")
        self.assertEqual(r(VariantCalls("151M 130M")), 1.0)
        self.assertEqual(r(VariantCalls("151M 130T")), -1.0)

    def test_boolean_or(self):
        r = HCVR2("( ( 151M AND 69i ) => 1.0)")
        self.assertEqual(r(VariantCalls("151M 69i")), 1.0)
        self.assertEqual(r(VariantCalls("151S 69i")), 0.0)

    def test_nested_condition(self):
        r = HCVR2("(L41L => 3, 67N => (2N => 2)))")
        self.assertEqual(r(VariantCalls("67N 2N L41S")), 2.0)

    def test_mean_expr(self):
        r = HCVR2("MEAN (L41L => 3, 67N => (2N => 2)))")
        self.assertEqual(r(VariantCalls("67N 2N L41R")), 1.0)
        self.assertEqual(r(VariantCalls("67N 2N L41L")), 2.5)

    def test_sum_expr(self):
        r = HCVR2("(0.5, L41L => 3, 67N => (2N => 2))")
        self.assertEqual(r(VariantCalls("67N 2N L41R")), 2.5)
        self.assertEqual(r(VariantCalls("67N 2N L41L")), 5.5)

    def test_bools_expr(self):
        r = HCVR2("(TRUE => (0.5, 0.25), FALSE => 5, L41L => 3, 67N => (2N => 2))")
        self.assertEqual(r(VariantCalls("67N 2N L41R")), 2.75)
        self.assertEqual(r(VariantCalls("67N 2N L41L")), 5.75)
