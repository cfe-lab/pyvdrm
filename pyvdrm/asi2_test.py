import unittest
from functools import reduce
from pyvdrm.hivdb import ASI2
from pyvdrm.vcf import Mutation, MutationSet

def mus(mutations):
    if mutations == '':
        return set([])
    return reduce(lambda x, y: x.union(y),
            map(lambda x: set(MutationSet.from_string(x)), mutations.split()))

class TestRuleParser(unittest.TestCase):

    def test_stanford_ex1(self):
        ASI2("151M OR 69i")
    def test_stanford_ex2(self):
        rule = ASI2("SELECT ATLEAST 2 FROM (41L, 67N, 70R, 210W, 215F, 219Q)")
        m1 = Mutation(41, 'L')
        m2 = Mutation(67, 'N')
        m3 = Mutation(70 ,'N')
        self.assertTrue(rule([m1, m2]))
        self.assertFalse(rule([m1, m3]))

    def test_stanford_ex3(self):
        ASI2("SELECT ATLEAST 2 AND NOTMORETHAN 2 FROM (41L, 67N, 70R, 210W, 215FY, 219QE)")
    def test_stanford_ex4(self):
        ASI2("215FY AND NOT 184VI")


    def test_stanford_rest(self):
        examples = [
        "SCORE FROM (65R => 20, 74V => 20, 184VI => 20)",
        "151M AND EXCLUDE 69i",
#        "69(NOT TDN)",
        "215F OR 215Y",
        "SCORE FROM (101P => 40, 101E => 30, 101HN => 15, 101Q => 5 )",
        "SCORE FROM ( MAX  (101P => 40, 101E => 30, 101HN => 15, 101Q => 5 ))",
        "(184V AND 115F) => 20"
        "3N AND 9N",
        "2N OR 9N AND 2N",
        "3N AND (2N AND (4N OR 2N))",
        ]

        for ex in examples:
            print(ex)
            x = ASI2(ex)

    def test_asi2_compat(self):
        q = "SCORE FROM ( 98G => 10, 100I => 40,\
                          MAX (101P => 40, 101E => 30, 101HN => 15, 101Q => 5) )"
        ASI2(q)

class TestRuleSemantics(unittest.TestCase):
    
    def test_score_from(self):
        rule = ASI2("SCORE FROM ( 100G => 10, 101D => 20 )")
        self.assertEqual(rule(mus("100G 102G")), 10)

    def test_score_from_max(self):
        rule = ASI2("SCORE FROM (MAX (100G => 10, 101D => 20, 102D => 30))")
        self.assertEqual(rule(mus("100G 101D")), 20)
        self.assertEqual(rule(mus("10G 11D")), False)

    def test_score_from_max_neg(self):
        rule = ASI2("SCORE FROM (MAX (100G => -10, 101D => -20, 102D => 30))")
        self.assertEqual(rule(mus("100G 101D")), -10)
        self.assertEqual(rule(mus("10G 11D")), False) 


    def test_bool_and(self):
        rule = ASI2("1G AND (2T AND 7Y)")
        self.assertEqual(rule(mus("2T 7Y 1G")), True)
        self.assertEqual(rule(mus("2T 1Y 1G")), False)
        self.assertEqual(rule(mus("7Y 1G 2T")), True)
        self.assertEqual(rule([]), False)

    def test_bool_or(self):
        rule = ASI2("1G OR (2T OR 7Y)")
        self.assertTrue(rule(mus("2T")))
        self.assertFalse(rule(mus("3T")))
        self.assertTrue(rule(mus("1G")))
        self.assertFalse(rule([]))

    def test_select_from_atleast(self):
        rule = ASI2("SELECT ATLEAST 2 FROM (2T, 7Y, 3G)")
        self.assertTrue(rule(mus("2T 7Y 1G")))
        self.assertFalse(rule(mus("2T 0Y 0G")))
        self.assertTrue(rule(mus("3G 9Y 2T")))

    def test_score_from_exactly(self):
        rule = ASI2("SELECT EXACTLY 1 FROM (2T, 7Y)")
        score = rule(mus("2T 7Y 1G"))

class TestActualRules(unittest.TestCase):
    def test_hivdb_rules_parse(self):
        for line in open("HIVDB.rules"):
            print(line)
            r = ASI2(line.strip())

    def test_chained_and(self):
        rule = ASI2("""
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
        self.assertEqual(rule(mus("40F 41L 210W 215Y")), 65)
        self.assertEqual(rule(mus("41L 210W 215F")), 60)
        self.assertEqual(rule(mus("40F 210W 215Y")), 25)
        self.assertEqual(rule(mus("40F 67G 215Y")), 15) 
if __name__ == '__main__':
    unittest.main()