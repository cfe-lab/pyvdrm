import unittest
from pyvdrm.hcvr import HCVR, AsiMutations, Score
from pyvdrm.vcf import Mutation, MutationSet, VariantCalls


# noinspection SqlNoDataSourceInspection,SqlDialectInspection
class TestRuleParser(unittest.TestCase):

    def test_stanford_ex1(self):
        HCVR("151M OR 69i")

    def test_stanford_ex2(self):
        rule = HCVR("SELECT ATLEAST 2 FROM (41L, 67N, 70R, 210W, 215F, 219Q)")
        m1 = MutationSet('41L')
        m2 = MutationSet('67N')
        m3 = MutationSet('70N')
        self.assertTrue(rule([m1, m2]))
        self.assertFalse(rule([m1, m3]))

    def test_stanford_ex3(self):
        HCVR("SELECT ATLEAST 2 AND NOTMORETHAN 2 FROM (41L, 67N, 70R, 210W, 215FY, 219QE)")

    def test_stanford_ex4(self):
        HCVR("215FY AND NOT 184VI")

    def test_stanford_rest(self):
        examples = ["SCORE FROM (65R => 20, 74V => 20, 184VI => 20)",
                    "151M AND EXCLUDE 69i",
                    # "69(NOT TDN)",
                    "215F OR 215Y",
                    "SCORE FROM (101P => 40, 101E => 30, 101HN => 15, 101Q => 5 )",
                    "SCORE FROM ( MAX  (101P => 40, 101E => 30, 101HN => 15, 101Q => 5 ))",
                    "(184V AND 115F) => 20"
                    "3N AND 9N",
                    "2N OR 9N AND 2N",
                    "3N AND (2N AND (4N OR 2N))"]

        for ex in examples:
            x = HCVR(ex)
            self.assertEqual(ex, x.rule)

    def test_asi2_compat(self):
        q = "SCORE FROM ( 98G => 10, 100I => 40,\
                          MAX (101P => 40, 101E => 30, 101HN => 15, 101Q => 5) )"
        HCVR(q)


# noinspection SqlNoDataSourceInspection,SqlDialectInspection
class TestRuleSemantics(unittest.TestCase):
    def test_score_from(self):
        rule = HCVR("SCORE FROM ( 100G => 10, 101D => 20 )")
        self.assertEqual(rule(VariantCalls("100G 102G")), 10)

    def test_score_negate(self):
        rule = HCVR("SCORE FROM ( NOT 100G => 10, NOT 101SD => 20 )")
        self.assertEqual(rule(VariantCalls("100G 102G")), 20)
        self.assertEqual(rule(VariantCalls("100S 101S")), 10)

    def test_score_residues(self):
        rule = HCVR("SCORE FROM ( 100G => 10, 101D => 20 )")
        expected_residue = repr({Mutation('S100G')})

        result = rule.dtree(VariantCalls("S100G R102G"))

        self.assertEqual(expected_residue, repr(result.residues))

    def test_score_from_max(self):
        rule = HCVR("SCORE FROM (MAX (100G => 10, 101D => 20, 102D => 30))")
        self.assertEqual(rule(VariantCalls("100G 101D")), 20)
        self.assertEqual(rule(VariantCalls("10G 11D")), False)

    def test_score_from_max_neg(self):
        rule = HCVR("SCORE FROM (MAX (100G => -10, 101D => -20, 102D => 30))")
        self.assertEqual(rule(VariantCalls("100G 101D")), -10)
        self.assertEqual(rule(VariantCalls("10G 11D")), False)

    def test_bool_and(self):
        rule = HCVR("1G AND (2T AND 7Y)")
        self.assertEqual(rule(VariantCalls("2T 7Y 1G")), True)
        self.assertEqual(rule(VariantCalls("2T 3Y 1G")), False)
        self.assertEqual(rule(VariantCalls("7Y 1G 2T")), True)
        self.assertEqual(rule([]), False)

    def test_bool_constants(self):
        rule = HCVR("TRUE OR 1G")
        self.assertEqual(rule(VariantCalls("2G")), True)
        rule = HCVR("FALSE AND 1G")
        self.assertEqual(rule(VariantCalls("1G")), False)
        rule = HCVR("TRUE OR (FALSE AND TRUE)")
        self.assertEqual(rule(VariantCalls("1G")), True)

    def test_bool_or(self):
        rule = HCVR("1G OR (2T OR 7Y)")
        self.assertTrue(rule(VariantCalls("2T")))
        self.assertFalse(rule(VariantCalls("3T")))
        self.assertTrue(rule(VariantCalls("1G")))
        self.assertFalse(rule([]))

    def test_select_from_atleast(self):
        rule = HCVR("SELECT ATLEAST 2 FROM (2T, 7Y, 3G)")
        self.assertTrue(rule(VariantCalls("2T 7Y 1G")))
        self.assertFalse(rule(VariantCalls("2T 4Y 5G")))
        self.assertTrue(rule(VariantCalls("3G 9Y 2T")))

    def test_score_from_exactly(self):
        rule = HCVR("SELECT EXACTLY 1 FROM (2T, 7Y)")
        score = rule(VariantCalls("2T 7Y 1G"))
        self.assertEqual(0, score)

    def test_score_comment(self):
        rule = HCVR("SCORE FROM (100G => 10, 200T => 3, 100S => \"comment\")")
        self.assertEqual(rule(VariantCalls("100G")), 10)
        result = rule.dtree(VariantCalls("100S 200T"))
        self.assertEqual(result.score, 3)
        self.assertTrue("comment" in result.flags)

class TestActualRules(unittest.TestCase):
    def test_hivdb_rules_parse(self):
        for line in open("pyvdrm/tests/HIVDB.rules"):
            r = HCVR(line)
            self.assertEqual(line, r.rule)

    def test_chained_and(self):
        rule = HCVR("""
        SCORE FROM(41L => 5, 62V => 5, MAX ( 65E => 10, 65N =>
        30, 65R => 45 ), MAX ( 67E => 5, 67G => 5, 67H => 5, 67N => 5, 67S =>
        5, 67T => 5, 67d => 30 ), 68d => 15, MAX ( 69G => 10, 69i => 60, 69d =>
        15 ), MAX ( 70E => 15, 70G => 15, 70N => 15, 70Q => 15, 70R => 5, 70S
        => 15, 70T => 15, 70d => 15 ), MAX ( 74I => 30, 74V => 30 ), 75I => 5,
        77L => 5, 115F => 60, 116Y => 10, MAX ( 151L => 30, 151M => 60 ), MAX(
        184I => 15, 184V => 15 ), 210W => 5, MAX ( 215A => 5, 215C => 5, 215D
        => 5, 215E => 5, 215F => 10, 215I => 5, 215L => 5, 215N => 5, 215S =>
        5, 215V => 5, 215Y => 10 ), MAX ( 219E => 5, 219N => 5, 219Q => 5, 219R
        => 5 ), (40F AND 41L AND 210W AND 215FY) => 5, (41L AND 210W) => 10,
        (41L AND 210W AND 215FY) => 5, (41L AND 44AD AND 210W AND 215FY) => 5,
        (41L AND 67EGN AND 215FY) => 5, (67EGN AND 215FY AND 219ENQR) => 5,
        (67EGN AND 70R AND 184IV AND 219ENQR) => 20, (67EGN AND 70R AND
        219ENQR) => 10, (70R AND 215FY) => 5, (74IV AND 184IV) => 15, (77L AND
        116Y AND 151M) => 10, MAX ((210W AND 215ACDEILNSV) => 5, (210W AND
        215FY) => 10), MAX ((41L AND 215ACDEILNSV) => 5, (41L AND 215FY) =>
        15))
        """)
        self.assertEqual(rule(VariantCalls("40F 41L 210W 215Y")), 65)
        self.assertEqual(rule(VariantCalls("41L 210W 215F")), 60)
        self.assertEqual(rule(VariantCalls("40F 210W 215Y")), 25)
        self.assertEqual(rule(VariantCalls("40F 67G 215Y")), 15)


class TestAsiMutations(unittest.TestCase):
    def test_init_args(self):
        expected_mutation_set = MutationSet('Q80KR')
        m = AsiMutations(args='Q80KR')

        self.assertEqual(expected_mutation_set, m.mutations)
        self.assertEqual(expected_mutation_set.wildtype, m.mutations.wildtype)

    def test_init_none(self):
        m = AsiMutations()

        self.assertIsNone(m.mutations)

    def test_repr(self):
        expected_repr = "AsiMutations(args='Q80KR')"
        m = AsiMutations(args='Q80KR')

        r = repr(m)

        self.assertEqual(expected_repr, r)

    def test_repr_none(self):
        expected_repr = "AsiMutations()"
        m = AsiMutations()

        r = repr(m)

        self.assertEqual(expected_repr, r)


class TestScore(unittest.TestCase):
    def test_init(self):
        expected_value = 10
        expected_mutations = {Mutation('A23R')}

        score = Score(expected_value, expected_mutations)

        self.assertEqual(expected_value, score.score)
        self.assertEqual(expected_mutations, score.residues)

    def test_repr(self):
        expected_repr = "Score(10, {Mutation('A23R')})"
        score = Score(10, {Mutation('A23R')})

        r = repr(score)

        self.assertEqual(expected_repr, r)


if __name__ == '__main__':
    unittest.main()
