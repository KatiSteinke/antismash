# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from antismash.common import pattern_search


class TestParsePattern(unittest.TestCase):
    # def test_pattern_validity(self):
    #     with self.assertRaisesRegex(ValueError, 'Invalid format'):
    #         pattern_search.check_pattern('Z-L')
    #     with self.assertRaisesRegex(ValueError, 'Invalid format: incorrect protein start in pattern'):
    #         pattern_search.check_pattern('A-L-<')
    #     with self.assertRaisesRegex(ValueError, 'Invalid format: incorrect protein end in pattern'):
    #         pattern_search.check_pattern('A->-L')
    #     with self.assertRaisesRegex(ValueError, 'Invalid format: multiple sets of repeats'):
    #         pattern_search.check_pattern('A(2,3)(3,4)-L')
    #     with self.assertRaisesRegex(ValueError, 'Invalid format: multiple sets of brackets in one element'):
    #         pattern_search.check_pattern('[AL][YG]-L')
    #


    def test_simple_amino(self):
        simple_result = pattern_search.SimpleAmino('A')
        assert simple_result.element == 'A'
        with self.assertRaisesRegex(ValueError, 'Invalid amino acid'):
            pattern_search.SimpleAmino('B')
        assert simple_result.match('A')
        assert not simple_result.match('C')

    def test_any_amino(self):
        any_result = pattern_search.AnyAmino('x')
        with self.assertRaisesRegex(ValueError, 'Attempting to use defined amino acid as AnyAmino'):
            pattern_search.AnyAmino('A')

    def test_multiple_amino(self):
        multiple_result = pattern_search.MultipleAmino('[AC]')
        assert multiple_result.element == {'A', 'C'}
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.MultipleAmino('[AC')
        with self.assertRaisesRegex(ValueError, 'Invalid amino acid'):
            pattern_search.MultipleAmino('[AB]')
        with self.assertRaisesRegex(ValueError, 'Invalid string for MultipleAmino'):
            pattern_search.MultipleAmino('[AC]A')

    def test_negated_amino(self):
        result = pattern_search.NegatedAmino('{AC}')
        assert result.element == {'A', 'C'}
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.NegatedAmino('{AC')
        with self.assertRaisesRegex(ValueError, 'Invalid amino acid'):
            pattern_search.NegatedAmino('{AB}')
        with self.assertRaisesRegex(ValueError, 'Invalid string for NegatedAmino'):
            pattern_search.NegatedAmino('{AC}A')

    def test_pattern(self):
        result = pattern_search.Pattern('A')
        self.assertIsInstance(result.elements[0], pattern_search.SimpleAmino)

        result = pattern_search.Pattern('A-C')
        assert len(result.elements) == 2

        result = pattern_search.Pattern('x')
        self.assertIsInstance(result.elements[0], pattern_search.AnyAmino)

        result = pattern_search.Pattern('[AC]')
        self.assertIsInstance(result.elements[0], pattern_search.MultipleAmino)

        result = pattern_search.Pattern('{AC}')
        self.assertIsInstance(result.elements[0], pattern_search.NegatedAmino)

    def test_invalid_pattern(self):
        with self.assertRaises(ValueError):
            pattern_search.Pattern('AA')
        with self.assertRaises(ValueError):
            pattern_search.Pattern('*')

    def test_find(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('Y')
        assert pattern.find(sequence) == -1

        pattern = pattern_search.Pattern('A')
        assert pattern.find(sequence) == 1

        pattern = pattern_search.Pattern('A-G')
        assert pattern.find(sequence) == 1

        pattern = pattern_search.Pattern('A-Y')
        assert pattern.find(sequence) == -1

        pattern = pattern_search.Pattern('A-T')
        assert pattern.find(sequence) == 6

    def test_find_anchor(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('Y')
        assert pattern.find_anchor(sequence) == -1

        pattern = pattern_search.Pattern('A')
        assert pattern.find_anchor(sequence) == 1
        assert pattern.find_anchor(sequence, 2) == 6
        assert pattern.find_anchor(sequence, 7) == -1
        # test out of bounds value
        assert pattern.find_anchor(sequence, 12) == -1
        # test searching from the end
        assert pattern.find_anchor(sequence, -3) == 6
        assert pattern.find_anchor(sequence, -1) == -1