import unittest
from pyvdrm.vcf import Mutation, MutationSet, call_mutations


class TestMutationParser(unittest.TestCase):

    def test_noWt(self):
        for residue in "20N 10NNN 1ASDF 0X".split():
            mutations = MutationSet.from_string(residue)
            self.assertIsNone(mutations.wildtype)

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


class TestCallMutations(unittest.TestCase):
    def test(self):
        reference = 'ACHE'
        sample = 'ICRE'
        expected_mutations = [Mutation.from_string(residue)
                              for residue in ('A1I', 'H3R')]

        mutations = call_mutations(reference, sample)

        self.assertEqual(expected_mutations, mutations)


if __name__ == '__main__':
    unittest.main()
