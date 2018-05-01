# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import pattern_search


class TestParsePattern(unittest.TestCase):

    def test_simple_amino(self):
        simple_result = pattern_search.SimpleAmino('A')
        assert simple_result.amino == 'A'
        with self.assertRaisesRegex(ValueError, 'Invalid amino acid'):
            pattern_search.SimpleAmino('B')
        assert simple_result.match('A')
        assert not simple_result.match('C')

    def test_any_amino(self):
        any_result = pattern_search.AnyAmino('x')
        assert any_result.max_repeats == 1
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
        assert result.min_repeats == 2
        result = pattern_search.MultipleAmino('[AT](24)')
        assert result.min_repeats == 24

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

    def test_find_any(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x')
        assert pattern.find(sequence) == 0

    def test_find_single_then_any(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('M-x')
        assert pattern.find('M') == -1
        assert pattern.find('') == -1
        assert pattern.find('AM') == -1
        assert pattern.find(sequence) == 0

    def test_find_any_combos(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x-A')
        assert pattern.find('A') == -1
        assert pattern.find('') == -1
        assert pattern.find('AC') == -1
        assert pattern.find(sequence) == 0


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

    def test_find_range(self):
        sequence = 'MYKITTYYYY'

        pattern = pattern_search.Pattern('Y(1,2)')
        assert pattern.find(sequence) == 1

        pattern = pattern_search.Pattern('Y(2,4)')
        assert pattern.find(sequence) == 6

        pattern = pattern_search.Pattern('T(2,3)-Y')
        assert pattern.find(sequence) == 4

        pattern = pattern_search.Pattern('[TY](2,3)')
        assert pattern.find(sequence) == 4

        pattern = pattern_search.Pattern('{MY}(2,3)')
        assert pattern.find(sequence) == 2

        pattern = pattern_search.Pattern('x(2,3)')
        assert pattern.find(sequence) == 0

        pattern = pattern_search.Pattern('T-Y(0,1)')
        assert pattern.find(sequence) == 4

        pattern = pattern_search.Pattern('K-Y(0,1)-I')
        assert pattern.find(sequence) == 2

        pattern = pattern_search.Pattern('K-x(0,1)-I')
        assert pattern.find(sequence) == 2

    def test_parse_repeats(self):
        assert pattern_search.parse_repeats('', 0) == (1, 1)
        assert pattern_search.parse_repeats('(5)', 0) == (5, 5)
        assert pattern_search.parse_repeats('(5,6)', 0) == (5, 6)

        assert pattern_search.parse_repeats('A(5,6)', -5) == (5, 6)

        with self.assertRaises(ValueError):
            pattern_search.parse_repeats('(5', 0)
        with self.assertRaises(ValueError):
            pattern_search.parse_repeats('(5-C(2)', 0)
        with self.assertRaisesRegex(ValueError, 'Invalid repeat'):
            pattern_search.parse_repeats('(5))', 0)

    def test_match_with_next(self):
        simple_amino = pattern_search.SimpleAmino('A')
        any_amino = pattern_search.AnyAmino('x', simple_amino)
        assert any_amino.match_including_following('MAGICHAT')
        assert any_amino.match_including_following('MAGICHAT').distance == 2

        assert any_amino.match_including_following('MAGICHAT', 5)

    def test_range_with_optional(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('H-A-T-S(0,1)')
        assert pattern.find(sequence) == 5

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

    def test_find_all(self):
        sequence = 'HEYKITTYKITTY'
        pattern = pattern_search.Pattern('Y')
        assert pattern.find_all(sequence) == [2, 7, 12]

        pattern = pattern_search.Pattern('K-I-T(2)-Y')
        assert pattern.find_all(sequence) == [3, 8]

        sequence = 'CATINTHEHAT'
        pattern = pattern_search.Pattern('x-A-T')
        assert pattern.find_all(sequence) == [0, 8]

        pattern = pattern_search.Pattern('[CH]-A-T')
        assert pattern.find_all(sequence) == [0, 8]

        pattern = pattern_search.Pattern('{M}-A-T')
        assert pattern.find_all(sequence) == [0, 8]
        # TODO: figure out sequence for checking overlapping matches

        sequence = 'HEEEEY'
        pattern = pattern_search.Pattern('E(1,4)-Y')
        assert pattern.find_all(sequence) == [1, 2, 3, 4]
