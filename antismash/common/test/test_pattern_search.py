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
            pattern_search.Pattern('AA.')
        with self.assertRaises(ValueError):
            pattern_search.Pattern('*.')

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
        with self.assertRaises(ValueError):
            pattern_search.parse_repeats('(5,6,8)', 0)


class TestFindSimple(unittest.TestCase):
    def test_find(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('Y.')
        assert pattern.find(sequence) == -1
        pattern = pattern_search.Pattern('A.')
        assert pattern.find(sequence) == 1

        pattern = pattern_search.Pattern('A-G.')
        assert pattern.find(sequence) == 1

        pattern = pattern_search.Pattern('A-Y.')
        assert pattern.find(sequence) == -1

        pattern = pattern_search.Pattern('A-T.')
        assert pattern.find(sequence) == 6

        assert pattern.find('') == -1

    def test_find_range_simple(self):
        sequence = 'MYKITTYYYY'

        pattern = pattern_search.Pattern('Y(1,2).')
        assert pattern.find(sequence) == 1

        pattern = pattern_search.Pattern('Y(2,4).')
        assert pattern.find(sequence) == 6

        pattern = pattern_search.Pattern('T(2,3)-Y.')
        assert pattern.find(sequence) == 4

    def test_find_optional_simple(self):
        sequence = 'MYKITTYYYY'
        pattern = pattern_search.Pattern('T-Y(0,1).')
        assert pattern.find(sequence) == 4

        pattern = pattern_search.Pattern('K-Y(0,1)-I.')
        assert pattern.find(sequence) == 2

    def test_optional_at_end(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('H-A-T-S(0,1).')
        assert pattern.find(sequence) == 5


class TestFindAny(unittest.TestCase):
    def test_find_any(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x.')
        assert pattern.find(sequence) == 0

    def test_find_single_then_any(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('M-x.')
        assert pattern.find('M') == -1
        assert pattern.find('') == -1
        assert pattern.find('AM') == -1
        assert pattern.find(sequence) == 0

    def test_find_any_combos(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x-A.')
        assert pattern.find('A') == -1
        assert pattern.find('') == -1
        assert pattern.find('AC') == -1
        assert pattern.find(sequence) == 0

        # three-letter combinations:
        pattern = pattern_search.Pattern('x-A-x.')
        assert pattern.find(sequence) == 0
        assert pattern.find('A') == -1

    def test_find_any_repeats(self):
        sequence = 'MAGICHAT'

        pattern = pattern_search.Pattern('M-x-x-x-C.')
        assert pattern.find(sequence) == 0

        pattern = pattern_search.Pattern('M-x(3)-C.')
        assert pattern.find(sequence) == 0

    def test_find_range_any(self):
        sequence = 'MYKITTYYYY'

        pattern = pattern_search.Pattern('x(2,3).')
        assert pattern.find(sequence) == 0

        pattern = pattern_search.Pattern('K-x(0,1)-I.')
        assert pattern.find(sequence) == 2


class TestFindMultiple(unittest.TestCase):
    def test_find_multiple(self):
        sequence = 'MAGICHAT'
        # simplest but inefficient: one acceptable amino acid
        pattern = pattern_search.Pattern('[A].')
        assert pattern.find(sequence) == 1
        assert pattern.find('M') == -1
        assert pattern.find('') == -1

        # actual use case: multiple acceptable amino acids
        pattern = pattern_search.Pattern('[AC].')
        # This finds an A...
        assert pattern.find(sequence) == 1
        # ...or a C
        assert pattern.find('CAT') == 0
        assert pattern.find('KITTY') == -1

    def test_find_multiple_repeats(self):
        sequence = 'MAGICHAT'

        pattern = pattern_search.Pattern('[AT](2).')
        assert pattern.find(sequence) == 6
        assert pattern.find('A') == -1

    def test_find_range_multiple(self):
        sequence = 'MYKITTYYYY'

        pattern = pattern_search.Pattern('[TY](2,3).')
        assert pattern.find(sequence) == 4


class TestFindNegated(unittest.TestCase):
    def test_find_negated(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('{A}.')
        assert pattern.find(sequence) == 0
        assert pattern.find('A') == -1
        assert pattern.find('') == -1

        pattern = pattern_search.Pattern('{MAGIC}.')
        assert pattern.find(sequence) == 5
        assert pattern.find('MAGIC') == -1
        assert pattern.find('GAMMA') == -1

    def test_find_negated_repeats(self):
        pattern = pattern_search.Pattern('{A}(3).')
        assert pattern.find('BB') == -1

    def test_find_range_negated(self):
        sequence = 'MYKITTYYYY'
        pattern = pattern_search.Pattern('{MY}(2,3).')
        assert pattern.find(sequence) == 2


class TestTermini(unittest.TestCase):
    def test_termini(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('<M.')
        assert pattern.find(sequence) == 0
        assert pattern.find('AMAGICHAT') == -1

        pattern = pattern_search.Pattern('T>.')
        assert pattern.find(sequence) == 7
        assert pattern.find('MAGICHATS') == -1

        pattern = pattern_search.Pattern('<M-x(6)-T>.')
        assert pattern.find(sequence) == 0
        assert pattern.find('MAGICHATS') == -1
        assert pattern.find('AMAGICHAT') == -1

        sequence = 'MMMCHEESE'
        pattern = pattern_search.Pattern('<M(3).')
        assert pattern.find(sequence) == 0

        pattern = pattern_search.Pattern('<M(3,4).')
        assert pattern.find(sequence) == 0

    def test_nterm_bad_offset(self):
        assert not pattern_search.NTerminalAmino('<M').match('MAGICHAT', -1)

    def test_cterm_bad_offset(self):
        assert not pattern_search.CTerminalAmino('T>').match('MAGICHAT', -1)

    def test_invalid_ctermini(self):
        with self.assertRaisesRegex(ValueError, 'Invalid pattern'):
            pattern_search.Pattern('A>-T.')

        with self.assertRaisesRegex(ValueError, 'Invalid pattern'):
            pattern_search.Pattern('A>-T>.')

    def test_invalid_ntermini(self):
        with self.assertRaisesRegex(ValueError, 'Invalid pattern'):
            pattern_search.Pattern('<M-<A.')

        with self.assertRaisesRegex(ValueError, 'Invalid pattern'):
            pattern_search.Pattern('M-<A.')

    def test_amino_or_cterminus_hit(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('T-[T>].')
        assert pattern.find(sequence) == 7

    def test_amino_or_cterminus_miss(self):
        pattern = pattern_search.Pattern('A-[M>].')
        assert pattern.find('MAGICHAT') == -1
        assert not pattern.find_all('MAGICHAT')

    def test_amino_or_cterminus_hit_multiple(self):
        pattern = pattern_search.Pattern('A-[GT>].')
        assert pattern.find('MAGICHAT') == 1
        assert len(pattern.find_all('MAGICHAT')) == 2

    def test_optional_cterminus_repeat(self):
        pattern = pattern_search.Pattern('A-[T>](1,2).')
        assert pattern.find('MAGICHAT') == 6
        hits = pattern.find_all('MAGICHAT')
        assert len(hits) == 1
        assert hits[0].start == 6 and hits[0].end == 8
        assert "MAGICHAT"[hits[0].start:hits[0].end] == "AT"

        pattern = pattern_search.Pattern('[T>](1,2).')
        assert pattern.find('MAGICHAT') == 7
        assert len(pattern.find_all('MAGICHAT')) == 1

    def test_invalid_optional_or_cterm(self):
        with self.assertRaises(ValueError):
            pattern_search.Pattern('A>-[T>].')

        with self.assertRaises(ValueError):
            pattern_search.Pattern('A>-C-[T>].')


class TestMatchWithFollowing(unittest.TestCase):
    def test_match_with_next(self):
        simple_amino = pattern_search.SimpleAmino('A')
        any_amino = pattern_search.AnyAmino('x', simple_amino)
        assert any_amino.match_including_following('MAGICHAT')

        assert any_amino.match_including_following('MAGICHAT').distance == 2

        assert any_amino.match_including_following('MAGICHAT', 5)


class TestMatchAll(unittest.TestCase):
    def test_match_all_lengths_12(self):
        sequence = 'MAAGGIC'
        multiple_as = pattern_search.SimpleAmino('A(1,2)', pattern_search.SimpleAmino('G(1,2)'))
        simple_m = pattern_search.SimpleAmino('M', multiple_as)
        assert len(simple_m.match_all_possible(sequence)) == 2

    def test_match_all_lengths_04(self):
        sequence = 'MAAAAGIC'
        multiple_as = pattern_search.SimpleAmino('A(0,4)')
        simple_m = pattern_search.SimpleAmino('M', multiple_as)
        assert len(simple_m.match_all_possible(sequence)) == 5

    def test_match_all_lengths_03(self):
        sequence = 'AAATAAAA'
        multi_a_1 = pattern_search.SimpleAmino('A(0,3)')
        simple_t = pattern_search.SimpleAmino('T', multi_a_1)
        multi_a_2 = pattern_search.SimpleAmino('A(0,3)', simple_t)
        assert len(multi_a_2.match_all_possible(sequence)) == 4

    def test_leading_repeat(self):
        sequence = 'MMAGIC'
        result = pattern_search.Pattern('M(1,2)-A.').head.match_all_possible(sequence)
        assert len(result) == 1

    def test_leading_multiple_repeats(self):
        sequence = 'MMAGIC'
        result = pattern_search.Pattern('M(1,2)-A(0,1).').head.match_all_possible(sequence)
        assert len(result) == 3


class TestFindAll(unittest.TestCase):
    def test_find_all_leading_repeat(self):
        sequence = 'MMAGIC'

        result = pattern_search.Pattern('M(1,2)-A.').find_all(sequence)
        assert len(result) == 2

    def test_find_all_leading_multiple_repeats(self):
        sequence = 'MMAGIC'

        result = pattern_search.Pattern('M(1,2)-A(0,1).').find_all(sequence)
        assert len(result) == 5

    def test_find_all_multiple_repeats(self):
        sequence = 'AAATAAAA'
        pattern = pattern_search.Pattern('A(0,3)-T-A(0,3).')
        assert len(pattern.find_all(sequence)) == 16

    def test_find_all_start(self):
        sequence = 'TAAATAAAATA'
        pattern = pattern_search.Pattern('<T-A(0,3).')
        assert len(pattern.find_all(sequence)) == 4

    def test_find_all_end(self):
        pattern = pattern_search.Pattern('A(0,3)-T-A>.')
        assert len(pattern.find_all('TAAATAAAATA')) == 4

    def test_find_all_basic_types(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('x(3).')
        assert len(pattern.find_all(sequence)) == 6

        pattern = pattern_search.Pattern('[MH]-[AGI](0,3).')
        assert len(pattern.find_all(sequence)) == 6

        pattern = pattern_search.Pattern('{AG}(1,3).')
        assert len(pattern.find_all(sequence)) == 8

    def test_distance(self):
        pattern = pattern_search.Pattern('M-A-G-I-C.')
        result = pattern.find_all('MAGICHAT')
        assert result[0].end == 5
        assert result[0].length == 5

        result = pattern.find_all('AMAGICHAT')
        assert result[0].end == 6
        assert result[0].length == 5
