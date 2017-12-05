import unittest
from pyvdrm.vcf import Mutation

class TestMutationParser(unittest.TestCase):

    def test_noWt(self):
        map(Mutation.from_string, "20N 10NNN 1ASDF 0X".split())

    def test_eq(self):
        r1 = Mutation.from_string("20A")
        r2 = Mutation.from_string("20A")
        self.assertEqual(r2, r1)

    def test_ineq_pos(self):
        r1 = Mutation.from_string("20A")
        r2 = Mutation.from_string("21A")
        self.assertNotEqual(r2, r1)

    def test_ineq_aa(self):
        r1 = Mutation.from_string("20A")
        r2 = Mutation.from_string("20C")
        self.assertNotEqual(r2, r1)


if __name__ == '__main__':
    unittest.main()
