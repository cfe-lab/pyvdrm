import unittest
from pyvdrm.vcf import Mutation, MutationSet, VariantCalls


class TestMutation(unittest.TestCase):
    def test_init_text(self):
        expected_wildtype = 'Q'
        expected_position = 20
        expected_variant = 'A'
        m = Mutation('Q20A')

        self.assertEqual(expected_wildtype, m.wildtype)
        self.assertEqual(expected_position, m.pos)
        self.assertEqual(expected_variant, m.variant)

    def test_init_details(self):
        expected_wildtype = 'Q'
        expected_position = 20
        expected_variant = 'A'
        m = Mutation(wildtype=expected_wildtype,
                     pos=expected_position,
                     variant=expected_variant)

        self.assertEqual(expected_wildtype, m.wildtype)
        self.assertEqual(expected_position, m.pos)
        self.assertEqual(expected_variant, m.variant)

    def test_init_bad_text(self):
        expected_message = \
            r'Mutation text expects wild type \(optional\), position, and one variant\.'
        with self.assertRaisesRegex(ValueError, expected_message):
            Mutation('!20A')

    def test_init_multiple_variants(self):
        expected_message = \
            r'Mutation text only allows one variant\.'
        with self.assertRaisesRegex(ValueError, expected_message):
            Mutation('Q20AE')

    def test_repr(self):
        m = Mutation('Q20A')
        r = repr(m)
        self.assertEqual("Mutation('Q20A')", r)

    def test_str(self):
        m = Mutation('Q20A')
        s = str(m)
        self.assertEqual("Q20A", s)

    def test_no_wildtype(self):
        m = Mutation('20A')
        r = repr(m)
        self.assertEqual("Mutation('20A')", r)
        self.assertIsNone(m.wildtype)

    def test_eq(self):
        r1 = Mutation("20A")
        r2 = Mutation("20A")
        self.assertEqual(r2, r1)

    def test_ineq_pos(self):
        r1 = Mutation("20A")
        r2 = Mutation("21A")
        self.assertNotEqual(r2, r1)

    def test_ineq_aa(self):
        r1 = Mutation("20A")
        r2 = Mutation("20C")
        self.assertNotEqual(r2, r1)

    def test_equal_with_wildtype_mismatch(self):
        set1 = Mutation('Q1A')
        set2 = Mutation('R1A')

        expected_message = 'Wild type mismatch between Q1A and R1A'
        with self.assertRaisesRegex(ValueError, expected_message):
            if set1 == set2:
                pass

    def test_immutable(self):
        m = Mutation('Q1A')

        with self.assertRaises(AttributeError):
            m.pos = 2


class TestMutationSet(unittest.TestCase):
    def test_init_text(self):
        expected_wildtype = 'Q'
        expected_position = 1
        expected_mutations = {Mutation('Q1A'), Mutation('Q1C')}

        ms = MutationSet('Q1AC')

        self.assertEqual(expected_wildtype, ms.wildtype)
        self.assertEqual(expected_position, ms.pos)
        self.assertEqual(expected_mutations, ms.mutations)

    def test_init_variants(self):
        expected_wildtype = 'Q'
        expected_position = 1
        expected_mutations = {Mutation('Q1A'), Mutation('Q1C')}

        ms = MutationSet(wildtype=expected_wildtype,
                         pos=expected_position,
                         variants='AC')

        self.assertEqual(expected_wildtype, ms.wildtype)
        self.assertEqual(expected_position, ms.pos)
        self.assertEqual(expected_mutations, ms.mutations)

    def test_init_mutations(self):
        expected_wildtype = 'Q'
        expected_position = 1
        expected_mutations = {Mutation('Q1A'), Mutation('Q1C')}

        ms = MutationSet(mutations=expected_mutations)

        self.assertEqual(expected_wildtype, ms.wildtype)
        self.assertEqual(expected_position, ms.pos)
        self.assertEqual(expected_mutations, ms.mutations)

    def test_init_mutations_some_without_wildtype(self):
        expected_wildtype = 'Q'
        expected_position = 1
        expected_mutations = {Mutation('Q1A'), Mutation('1C')}

        ms = MutationSet(mutations=expected_mutations)

        self.assertEqual(expected_wildtype, ms.wildtype)
        self.assertEqual(expected_position, ms.pos)
        self.assertEqual(expected_mutations, ms.mutations)

    def test_init_no_mutations(self):
        expected_wildtype = 'Q'
        expected_position = 1
        expected_mutations = set()

        ms = MutationSet(wildtype=expected_wildtype,
                         pos=expected_position,
                         mutations=expected_mutations)

        self.assertEqual(expected_wildtype, ms.wildtype)
        self.assertEqual(expected_position, ms.pos)
        self.assertEqual(expected_mutations, ms.mutations)

    def test_no_wildtypes(self):
        expected_position = 10
        expected_mutations = {Mutation('10C')}

        ms = MutationSet(pos=expected_position,
                         mutations=expected_mutations)

        self.assertIsNone(ms.wildtype)
        self.assertEqual(expected_position, ms.pos)
        self.assertEqual(expected_mutations, ms.mutations)

    def test_init_no_mutations_without_wildtype(self):
        position = 1
        mutations = set()

        with self.assertRaisesRegex(ValueError,
                                    r'No wildtype and no variants\.'):
            MutationSet(pos=position, mutations=mutations)

    def test_init_no_mutations_without_position(self):
        wildtype = 'Q'
        mutations = set()

        with self.assertRaisesRegex(ValueError,
                                    r'No position and no variants\.'):
            MutationSet(wildtype=wildtype, mutations=mutations)

    def test_init_position_mismatch(self):
        wildtype = 'Q'
        position = 2
        mutations = {Mutation('Q1A'), Mutation('Q1C')}

        with self.assertRaisesRegex(ValueError,
                                    r'Multiple positions found: 1, 2\.'):
            MutationSet(wildtype=wildtype, pos=position, mutations=mutations)

    def test_init_multiple_positions(self):
        mutations = {Mutation('Q3A'), Mutation('Q2C')}

        with self.assertRaisesRegex(ValueError,
                                    r'Multiple positions found: 2, 3\.'):
            MutationSet(mutations=mutations)

    def test_init_wildtype_mismatch(self):
        wildtype = 'R'
        position = 1
        mutations = {Mutation('Q1A'), Mutation('Q1C')}

        with self.assertRaisesRegex(ValueError,
                                    r'Multiple wildtypes found: Q, R\.'):
            MutationSet(wildtype=wildtype, pos=position, mutations=mutations)

    def test_init_multiple_wildtypes(self):
        mutations = {Mutation('R1A'), Mutation('Q1C')}

        with self.assertRaisesRegex(ValueError,
                                    r'Multiple wildtypes found: Q, R\.'):
            MutationSet(mutations=mutations)

    def test_repr(self):
        expected_repr = "MutationSet('Q1AC')"
        ms = MutationSet('Q1AC')

        r = repr(ms)

        self.assertEqual(expected_repr, r)

    def test_str(self):
        expected_str = 'Q1AC'
        ms = MutationSet(expected_str)

        s = str(ms)

        self.assertEqual(expected_str, s)

    def test_equal(self):
        set1 = MutationSet('Q1AC')
        set2 = MutationSet('Q1CA')
        set3 = MutationSet('Q2AC')
        set4 = MutationSet('R2AC')
        set5 = MutationSet('Q1A')

        self.assertEqual(set1, set2)
        self.assertNotEqual(set1, set3)
        self.assertNotEqual(set1, set4)
        self.assertNotEqual(set1, set5)

    def test_equal_with_wildtype_mismatch(self):
        set1 = MutationSet('Q1AC')
        set2 = MutationSet('R1AC')

        expected_message = 'Wild type mismatch between Q1AC and R1AC'
        with self.assertRaisesRegex(ValueError, expected_message):
            if set1 == set2:
                pass

    def test_hash(self):
        hash1 = hash(MutationSet('Q1AC'))
        hash2 = hash(MutationSet('Q1CA'))
        hash3 = hash(MutationSet('Q2AC'))
        hash4 = hash(MutationSet('Q1A'))
        hash5 = hash(MutationSet('R1AC'))

        self.assertEqual(hash1, hash2)
        self.assertNotEqual(hash1, hash3)
        self.assertNotEqual(hash1, hash4)
        self.assertEqual(hash1, hash5)

    def test_no_wildtype(self):
        for text in "20N 10NNN 1ASDF 0X".split():
            mutations = MutationSet(text)
            self.assertIsNone(mutations.wildtype)

    def test_immutable(self):
        ms = MutationSet('Q80KR')

        with self.assertRaises(AttributeError):
            ms.wildtype = 'D'


class TestVariantCalls(unittest.TestCase):
    def test_init_text(self):
        expected_mutation_sets = {MutationSet('A1IL'), MutationSet('H3R')}

        calls = VariantCalls('A1IL H3R')

        self.assertIsNone(calls.reference)
        self.assertEqual(expected_mutation_sets, calls.mutation_sets)

    def test_init_single_sequence(self):
        reference = 'ACHE'
        sample = 'ICRE'
        expected_calls = VariantCalls('A1I C2C H3R E4E')

        calls = VariantCalls(reference=reference, sample=sample)

        self.assertEqual(reference, calls.reference)
        self.assertEqual(expected_calls, calls)

    def test_init_multiple_sequences(self):
        reference = 'ACHE'
        sample = ['IN', 'C', 'HR', 'E']
        expected_calls = VariantCalls('A1IN C2C H3HR E4E')

        calls = VariantCalls(reference=reference, sample=sample)

        self.assertEqual(expected_calls, calls)

    def test_init_bad_length(self):
        reference = 'ACHE'
        sample = 'ICREL'

        with self.assertRaisesRegex(
                ValueError, r'Reference length was 4 and sample length was 5\.'):
            VariantCalls(reference=reference, sample=sample)

    def test_init_with_reference(self):
        expected_reference = 'ASH'
        expected_repr = "VariantCalls('A1IL H3R')"

        calls = VariantCalls('1IL 3R', reference=expected_reference)
        r = repr(calls)

        self.assertEqual(expected_reference, calls.reference)
        self.assertEqual(expected_repr, r)

    def test_repr(self):
        expected_repr = "VariantCalls('A1IL H3R')"
        calls = VariantCalls('A1IL H3R')

        r = repr(calls)

        self.assertEqual(expected_repr, r)

    def test_str(self):
        expected_str = 'A1IL H3R'
        calls = VariantCalls(expected_str)

        s = str(calls)

        self.assertEqual(expected_str, s)

    def test_eq(self):
        calls1 = VariantCalls('A1IL H3R')
        calls2 = VariantCalls('H3R A1IL')
        calls3 = VariantCalls('A1IL H3Q')
        calls4 = VariantCalls('A1IL')

        self.assertEqual(calls1, calls2)
        self.assertNotEqual(calls1, calls3)
        self.assertNotEqual(calls1, calls4)

    def test_hash(self):
        hash1 = hash(VariantCalls('A1IL H3R'))
        hash2 = hash(VariantCalls('H3R A1IL'))
        hash3 = hash(VariantCalls('A1IL H3Q'))
        hash4 = hash(VariantCalls('A1IL'))

        self.assertEqual(hash1, hash2)
        self.assertNotEqual(hash1, hash3)
        self.assertNotEqual(hash1, hash4)

    def test_iter(self):
        calls = VariantCalls('A1IL H3R')
        expected_mutation_sets = {MutationSet('A1IL'), MutationSet('H3R')}

        mutation_sets = set(calls)

        self.assertEqual(expected_mutation_sets, mutation_sets)

    def test_in(self):
        calls = VariantCalls('A1IL H3R')
        mutation_set1 = MutationSet('H3R')
        mutation_set2 = MutationSet('H4R')

        self.assertIn(mutation_set1, calls)
        self.assertNotIn(mutation_set2, calls)

    def test_immutable(self):
        calls = VariantCalls('A1IL H3R')

        with self.assertRaises(AttributeError):
            calls.reference = 'ASH'


if __name__ == '__main__':
    unittest.main()
