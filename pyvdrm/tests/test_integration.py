"""
Integration tests for pyvdrm - testing end-to-end usage examples
"""
import unittest

from pyvdrm.asi2 import ASI2
from pyvdrm.hcvr import HCVR
from pyvdrm.vcf import Mutation, VariantCalls


class TestASI2Integration(unittest.TestCase):
    """Test ASI2 algorithm with real-world examples"""

    def test_readme_example(self):
        """Test the example similar to README"""
        # Define a rule
        rule = ASI2("SCORE FROM (MAX (100T => 20, 282N => 15))")

        # Evaluate against mutations
        score = rule(VariantCalls("100d 282N"))

        self.assertEqual(15, score)

    def test_basic_scoring_rule(self):
        """Test a basic scoring rule"""
        rule = ASI2("SCORE FROM (65R => 20, 74V => 20, 184VI => 20)")

        # Test with one mutation
        self.assertEqual(20, rule(VariantCalls("65R 74d 184d")))

        # Test with multiple mutations
        self.assertEqual(40, rule(VariantCalls("65R 74V 184d")))

        # Test with all mutations
        self.assertEqual(60, rule(VariantCalls("65R 74V 184V")))

    def test_boolean_and_rule(self):
        """Test boolean AND logic"""
        rule = ASI2("100G AND 200T")

        self.assertTrue(rule(VariantCalls("100G 200T")))
        self.assertFalse(rule(VariantCalls("100G 200d")))
        self.assertFalse(rule(VariantCalls("100d 200T")))

    def test_boolean_or_rule(self):
        """Test boolean OR logic"""
        rule = ASI2("100G OR 200T")

        self.assertTrue(rule(VariantCalls("100G 200d")))
        self.assertTrue(rule(VariantCalls("100d 200T")))
        self.assertTrue(rule(VariantCalls("100G 200T")))
        self.assertFalse(rule(VariantCalls("100d 200d")))

    def test_select_atleast(self):
        """Test SELECT ATLEAST operator"""
        rule = ASI2("SELECT ATLEAST 2 FROM (41L, 67N, 70R)")

        self.assertTrue(rule(VariantCalls("41L 67N 70d")))
        self.assertTrue(rule(VariantCalls("41L 67N 70R")))
        self.assertFalse(rule(VariantCalls("41L 67d 70d")))

    def test_max_operator(self):
        """Test MAX operator in scoring"""
        rule = ASI2("SCORE FROM (MAX (100P => 40, 100E => 30, 100H => 15))")

        # Should take the maximum score
        self.assertEqual(30, rule(VariantCalls("100E")))
        self.assertEqual(40, rule(VariantCalls("100P")))
        self.assertEqual(15, rule(VariantCalls("100H")))

    def test_mutation_from_sequence(self):
        """Test creating mutations from aligned sequences"""
        reference = "ACHE"
        sample = "ICRE"

        calls = VariantCalls(reference=reference, sample=sample)

        # Should have mutations at positions 1 and 3
        self.assertEqual(4, len(calls))  # All positions present

        # Create a rule to check specific mutations
        rule = ASI2("1I AND 3R")
        self.assertTrue(rule(calls))


class TestHCVRIntegration(unittest.TestCase):
    """Test HCVR algorithm with real-world examples"""

    def test_basic_scoring(self):
        """Test basic HCVR scoring"""
        rule = HCVR("SCORE FROM (100G => 10, 200T => 20)")

        self.assertEqual(10, rule(VariantCalls("100G 200d")))
        self.assertEqual(20, rule(VariantCalls("100d 200T")))
        self.assertEqual(30, rule(VariantCalls("100G 200T")))

    def test_negative_mutations(self):
        """Test NOT operator (negative mutations)"""
        rule = HCVR("SCORE FROM (100!G => 10)")

        # Score when NOT G
        self.assertEqual(10, rule(VariantCalls("100T 200d")))

        # No score when it IS G
        self.assertEqual(0, rule(VariantCalls("100G 200d")))

    def test_min_operator(self):
        """Test MIN operator"""
        rule = HCVR("SCORE FROM (MIN (100G => 40, 100E => 30, 100H => 15))")

        # Should take the minimum score
        self.assertEqual(15, rule(VariantCalls("100H")))

    def test_multiple_variants_at_position(self):
        """Test mutations with multiple variants at a position"""
        rule = HCVR("SCORE FROM (100GE => 20)")

        # Should match if either G or E
        self.assertEqual(20, rule(VariantCalls("100G")))
        self.assertEqual(20, rule(VariantCalls("100E")))
        self.assertEqual(0, rule(VariantCalls("100T")))


class TestMutationAPI(unittest.TestCase):
    """Test the Mutation and MutationSet API"""

    def test_mutation_creation(self):
        """Test creating mutations in different ways"""
        # From string
        m1 = Mutation("Q80K")
        self.assertEqual("Q", m1.wildtype)
        self.assertEqual(80, m1.pos)
        self.assertEqual("K", m1.variant)

        # From parameters
        m2 = Mutation(wildtype="Q", pos=80, variant="K")
        self.assertEqual(m1, m2)

        # Without wildtype
        m3 = Mutation("80K")
        self.assertIsNone(m3.wildtype)
        self.assertEqual(80, m3.pos)

    def test_variant_calls_from_text(self):
        """Test creating VariantCalls from text"""
        calls = VariantCalls("A1I H3R E4D")

        self.assertEqual(3, len(calls))
        self.assertIn("A1I", str(calls))
        self.assertIn("H3R", str(calls))
        self.assertIn("E4D", str(calls))

    def test_variant_calls_from_sequences(self):
        """Test creating VariantCalls from sequences"""
        reference = "AHEC"
        sample = "IRDC"

        calls = VariantCalls(reference=reference, sample=sample)

        # All positions should be present
        self.assertEqual(4, len(calls))

        # Convert to string and check
        calls_str = str(calls)
        self.assertIn("A1I", calls_str)
        self.assertIn("H2R", calls_str)
        self.assertIn("E3D", calls_str)
        self.assertIn("C4C", calls_str)


class TestErrorHandling(unittest.TestCase):
    """Test error handling and validation"""

    def test_invalid_mutation_format(self):
        """Test that invalid mutation formats raise errors"""
        with self.assertRaises(ValueError):
            Mutation("!20A")

        with self.assertRaises(ValueError):
            Mutation("Q20")  # Missing variant

    def test_invalid_rule_syntax(self):
        """Test that invalid rule syntax raises ParseException"""
        from pyparsing import ParseException

        with self.assertRaises(ParseException):
            ASI2("SCORE FROM ( 10R => 2;0 )")

    def test_sequence_length_mismatch(self):
        """Test that mismatched sequence lengths raise error"""
        with self.assertRaises(ValueError):
            VariantCalls(reference="ACE", sample="ACHED")


if __name__ == '__main__':
    unittest.main()
