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
        with self.assertRaisesRegex(ValueError, 'Invalid pattern'):
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
        with self.assertRaisesRegex(ValueError, 'Invalid pattern'):
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

    def test_parse_termini(self):
        assert pattern_search.parse_terminus('A') == (False, 'A', False)
        assert pattern_search.parse_terminus('<A') == (True, 'A', False)
        assert pattern_search.parse_terminus('A>') == (False, 'A', True)
        assert pattern_search.parse_terminus('<A>') == (True, 'A', True)


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

        pattern = pattern_search.Pattern('Y(2,3)>.')
        assert pattern.find('KITTYYY') == 4

    def test_any_end(self):
        pattern = pattern_search.Pattern('A-x(1,2)>.')
        assert pattern.find('MAGICHAT') == 6

    def test_any_start(self):
        pattern = pattern_search.Pattern('<M(1,2).')
        assert pattern.find('MAGICHAT') == 0

    def multiple_start(self):
        pattern = pattern_search.Pattern('<[MA].')
        assert pattern.find('MAGICHAT') == 0

    def test_start_and_end(self):
        pattern = pattern_search.Pattern('<A>.')
        assert pattern.head.nterm
        assert pattern.head.cterm
        assert pattern.find('A') == 0
        assert pattern.find('AHAT') == -1

    def test_invalid_amino_start_end(self):
        with self.assertRaises(ValueError):
            pattern = pattern_search.Pattern('<AA>.')
            print(pattern.elements)

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

    def test_invalid_cterm_multiple(self):
        with self.assertRaisesRegex(ValueError, 'MultipleAmino can only have one of optional and fixed C terminus'):
            pattern_search.Pattern('[A>]>.')

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


class TestOffsetCheck(unittest.TestCase):
    def test_offset_ok_element(self):
        sequence = 'MAGICHAT'

        simple_amino = pattern_search.SimpleAmino('A')
        assert simple_amino.offset_ok(sequence, 0)
        assert not simple_amino.offset_ok(sequence, 8)

        simple_with_repeats = pattern_search.SimpleAmino('A(2)')
        assert simple_with_repeats.offset_ok(sequence, 0)
        assert not simple_with_repeats.offset_ok(sequence, 8)

    def test_offset_ok_element_nterm(self):
        sequence = 'MAGICHAT'

        simple_nterm = pattern_search.SimpleAmino('M', nterm=True)
        assert simple_nterm.offset_ok(sequence, 0)
        assert not simple_nterm.offset_ok(sequence, 1)

        nterm_with_repeats = pattern_search.SimpleAmino('M(2)', nterm=True)
        assert nterm_with_repeats.offset_ok(sequence, 0)
        assert nterm_with_repeats.offset_ok(sequence, 1)
        assert not nterm_with_repeats.offset_ok(sequence, 2)

        nterm_range = pattern_search.SimpleAmino('M(1,4)', nterm=True)
        assert nterm_range.offset_ok(sequence, 0)
        assert nterm_range.offset_ok(sequence, 3)
        assert not nterm_range.offset_ok(sequence, 4)

    def test_offset_ok_element_cterm(self):
        sequence = 'MAGICHAT'

        simple_cterm = pattern_search.SimpleAmino('M', cterm=True)
        assert simple_cterm.offset_ok(sequence, 7)
        assert not simple_cterm.offset_ok(sequence, 6)

        cterm_with_repeats = pattern_search.SimpleAmino('M(2)', cterm=True)
        assert cterm_with_repeats.offset_ok(sequence, 7)
        assert cterm_with_repeats.offset_ok(sequence, 6)
        assert not cterm_with_repeats.offset_ok(sequence, 2)

        cterm_range = pattern_search.SimpleAmino('M(1,4)', cterm=True)
        assert cterm_range.offset_ok(sequence, 7)
        assert cterm_range.offset_ok(sequence, 4)
        assert not cterm_range.offset_ok(sequence, 0)

    def test_offset_ok_multiple(self):
        sequence = 'MAGICHAT'

        multiple_plain = pattern_search.MultipleAmino('[MAT]')
        assert multiple_plain.offset_ok(sequence, 0)
        assert not multiple_plain.offset_ok(sequence, 8)

        multiple_with_repeat = pattern_search.MultipleAmino('[MAT](2)')
        assert multiple_with_repeat.offset_ok(sequence, 0)
        assert not multiple_with_repeat.offset_ok(sequence, 8)

    def test_offset_ok_multiple_nterm(self):
        sequence = 'MAGICHAT'

        multiple_nterm = pattern_search.MultipleAmino('[MAT]', nterm=True)
        assert multiple_nterm.offset_ok(sequence, 0)
        assert not multiple_nterm.offset_ok(sequence, 1)

        multiple_nterm_with_repeat = pattern_search.MultipleAmino('[MAT](2)', nterm=True)
        assert multiple_nterm_with_repeat.offset_ok(sequence, 0)
        assert multiple_nterm_with_repeat.offset_ok(sequence, 1)
        assert not multiple_nterm_with_repeat.offset_ok(sequence, 2)

        multiple_nterm_range = pattern_search.MultipleAmino('[MAT](1,4)', nterm=True)
        assert multiple_nterm_range.offset_ok(sequence, 0)
        assert multiple_nterm_range.offset_ok(sequence, 3)
        assert not multiple_nterm_range.offset_ok(sequence, 4)

    def test_offset_ok_multiple_cterm(self):
        sequence = 'MAGICHAT'

        multiple_cterm = pattern_search.MultipleAmino('[MAT]', cterm=True)
        assert multiple_cterm.offset_ok(sequence, 7)
        assert not multiple_cterm.offset_ok(sequence, 6)

        multiple_cterm_with_repeat = pattern_search.MultipleAmino('[MAT](2)', cterm=True)
        assert multiple_cterm_with_repeat.offset_ok(sequence, 7)
        assert multiple_cterm_with_repeat.offset_ok(sequence, 6)
        assert not multiple_cterm_with_repeat.offset_ok(sequence, 2)

        multiple_cterm_range = pattern_search.MultipleAmino('[MAT](1,4)', cterm=True)
        assert multiple_cterm_range.offset_ok(sequence, 7)
        assert multiple_cterm_range.offset_ok(sequence, 4)
        assert not multiple_cterm_range.offset_ok(sequence, 0)

    def test_offset_ok_optional_cterm(self):
        sequence = 'MAGICHAT'

        optional_cterm = pattern_search.MultipleAmino('[T>]')
        assert optional_cterm.offset_ok(sequence, 7)
        assert optional_cterm.offset_ok(sequence, 8)


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
        pattern = pattern_search.Pattern('<T.')
        assert len(pattern.find_all(sequence)) == 1

        pattern = pattern_search.Pattern('<T-A(0,3).')
        assert len(pattern.find_all(sequence)) == 4

    def test_find_all_end(self):
        sequence = 'TAAATAAAATA'
        pattern = pattern_search.Pattern('A>.')
        assert len(pattern.find_all(sequence)) == 1

        pattern = pattern_search.Pattern('A(0,3)-T-A>.')
        assert len(pattern.find_all('TAAATAAAATA')) == 4

    def test_find_all_range_end(self):
        sequence = 'KITTYYYY'
        pattern = pattern_search.Pattern('Y(2,5)>.')
        assert len(pattern.find_all(sequence)) == 3

    def test_find_all_range_start(self):
        pattern = pattern_search.Pattern('<M(1,3).')
        assert len(pattern.find_all('MMAGIC')) == 2

    def test_multiple_end(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('[AT]>.')
        assert len(pattern.find_all(sequence)) == 1

    def test_multiple_start(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('<[MA].')
        assert len(pattern.find_all(sequence)) == 1

    def test_find_all_range_multi_start(self):
        assert len(pattern_search.Pattern('<[MA](1,3).').find_all('MAGIC')) == 2

    def test_find_all_range_multi_end(self):
        assert len(pattern_search.Pattern('[AT](1,3)>.').find_all('MAGICHAT')) == 2

    def test_negated_end(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('{D}>.')
        assert len(pattern.find_all(sequence)) == 1

    def test_negated_start(self):
        sequence = 'MAGICHAT'
        pattern = pattern_search.Pattern('<{N}.')
        assert len(pattern.find_all(sequence)) == 1

    def test_find_all_range_not_start(self):
        assert len(pattern_search.Pattern('<{G}(1,3).').find_all('MAGIC')) == 2

    def test_find_all_range_not_end(self):
        assert len(pattern_search.Pattern('{H}(1,3)>.').find_all('MAGICHAT')) == 2

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
        assert len(result[0]) == 5

        result = pattern.find_all('AMAGICHAT')
        assert result[0].end == 6
        assert len(result[0]) == 5
