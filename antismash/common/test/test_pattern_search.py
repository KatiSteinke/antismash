# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import pattern_search


class TestParsePattern(unittest.TestCase):

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
        assert multiple_result.options == {'A', 'C'}
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.MultipleAmino('[AC')
        with self.assertRaisesRegex(ValueError, 'Invalid amino acid'):
            pattern_search.MultipleAmino('[AB]')
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):  # TODO: better error message...
            pattern_search.MultipleAmino('[AC]A')
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.MultipleAmino('{AC]')
        with self.assertRaisesRegex(ValueError, 'No valid options provided'):
            pattern_search.MultipleAmino('[]')

    def test_negated_amino(self):
        result = pattern_search.NegatedAmino('{AC}')
        assert result.options == {'A', 'C'}
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.NegatedAmino('{AC')
        with self.assertRaisesRegex(ValueError, 'Invalid amino acid'):
            pattern_search.NegatedAmino('{AB}')
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.NegatedAmino('{AC}A')
        with self.assertRaisesRegex(ValueError, 'Brackets do not match'):
            pattern_search.NegatedAmino('[AC}')
        with self.assertRaisesRegex(ValueError, 'No valid options provided'):
            pattern_search.NegatedAmino('{}')

    def test_multiple_repeats(self):
        result = pattern_search.MultipleAmino('[AT](2)')
        assert result.repeats == 2
        result = pattern_search.MultipleAmino('[AT](24)')
        assert result.repeats == 24


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

        assert pattern.find('') == -1

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

    def test_find_any(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x')
        assert pattern.find(sequence) == 0

    def test_find_any_combos(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x-A')
        assert pattern.find(sequence) == 0
        assert pattern.find('A') == -1
        assert pattern.find('') == -1
        assert pattern.find('AC') == -1

        # three-letter combinations:
        pattern = pattern_search.Pattern('x-A-x')
        assert pattern.find(sequence) == 0
        assert pattern.find('A') == -1

    def test_find_multiple(self):
        sequence = 'MAGICHAT'
        # simplest but inefficient: one acceptable amino acid
        pattern = pattern_search.Pattern('[A]')
        assert pattern.find(sequence) == 1
        assert pattern.find('M') == -1
        assert pattern.find('') == -1

        # actual use case: multiple acceptable amino acids
        pattern = pattern_search.Pattern('[AC]')
        # This finds an A...
        assert pattern.find(sequence) == 1
        # ...or a C
        assert pattern.find('CAT') == 0
        assert pattern.find('KITTY') == -1

    def test_find_negated(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('{A}')
        assert pattern.find(sequence) == 0
        assert pattern.find('A') == -1
        assert pattern.find('') == -1

        pattern = pattern_search.Pattern('{MAGIC}')
        assert pattern.find(sequence) == 5
        assert pattern.find('MAGIC') == -1
        assert pattern.find('GAMMA') == -1

    def test_find_any_repeats(self):
        sequence = 'MAGICHAT'

        pattern = pattern_search.Pattern('M-x-x-x-C')
        assert pattern.find(sequence) == 0

        pattern = pattern_search.Pattern('M-x(3)-C')
        assert pattern.find(sequence) == 0

    def test_find_multiple_repeats(self):
        sequence = 'MAGICHAT'

        pattern = pattern_search.Pattern('[AT](2)')
        assert pattern.find(sequence) == 6
        assert pattern.find('A') == -1

    def test_find_negated_repeats(self):
        pattern = pattern_search.Pattern('{A}(3)')
        assert pattern.find('BB') == -1

