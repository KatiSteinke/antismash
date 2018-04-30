# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""A collection of functions for parsing protein motif patterns of a
    defined format and searching for them in an amino acid sequence.
"""

from typing import List, Set, Tuple  # pylint:disable=unused-import

AMINOS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
          'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}


class Match:
    def __init__(self, hit: bool, distance: int = -1) -> None:
        self.hit = hit
        if hit:
            assert distance > -1
        self._distance = distance

    @property
    def distance(self) -> int:
        if not self.hit:
            raise ValueError('Cannot access distance without a match')
        return self._distance

    def __bool__(self) -> bool:
        return self.hit


class Element:
    def __init__(self, min_repeats: int, max_repeats: int, next_element: "Element" = None) -> None:
        self.min_repeats = min_repeats
        self.max_repeats = max_repeats
        self.next_element = next_element

    def match(self, sequence: str) -> Match:
        raise NotImplementedError('Missing match implementation')

    def match_including_following(self, sequence: str) -> Match:
        for repeat in range(self.max_repeats, self.min_repeats - 1, -1):
            if repeat:
                if len(sequence) < repeat:
                    continue
                matches = True
                for char in sequence[:repeat]:
                    match = self.match(char)
                    matches = matches and match
                if not matches:
                    continue
            if self.next_element:
                nextmatch = self.next_element.match_including_following(sequence[repeat:])
            else:
                return Match(True, repeat)
            if nextmatch:
                return Match(True, repeat + nextmatch.distance)
        return Match(False)


class SimpleAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        if element[0] not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.amino = element[0]
        min_repeats, max_repeats = parse_repeats(element[1:])
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str) -> Match:
        if not sequence:
            return Match(False)
        return Match(sequence[0] == self.amino, 1)

    def __str__(self) -> str:
        return "SimpleAmino(%s)(%d,%d)" % (self.amino, self.min_repeats, self.max_repeats)

class AnyAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        if element[0] != 'x':
            raise ValueError('Attempting to use defined amino acid as AnyAmino')
        min_repeats, max_repeats = parse_repeats(element[1:])
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str) -> Match:
        # TODO: should invalid amino be an error?
        return Match(len(sequence) > 0, 1)

    def __str__(self) -> str:
        return "AnyAmino(%d,%d)" % (self.min_repeats, self.max_repeats)


class MultipleAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        self.options, min_repeats, max_repeats = parse_options(element, "[", "]")
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str) -> Match:
        if not sequence:
            return Match(False)
        return Match(sequence[0] in self.options, 1)

    def __str__(self) -> str:
        return "MultipleAmino(%s)(%d,%d)" % (self.options, self.min_repeats, self.max_repeats)


class NegatedAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        self.options, min_repeats, max_repeats = parse_options(element, "{", "}")
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str) -> Match:
        if not sequence:
            return Match(False)
        return Match(sequence[0] not in self.options, 1)

    def __str__(self) -> str:
        return "MultipleAmino(%s)(%d,%d)" % (self.options, self.min_repeats, self.max_repeats)


class Pattern:
    def __init__(self, pattern_string: str) -> None:
        self.elements = [None]  # type: List[Element]
        for part in pattern_string.split('-')[::-1]:
            if part[0] in AMINOS:
                self.elements.append(SimpleAmino(part, self.elements[-1]))
            elif part[0] == 'x':
                self.elements.append(AnyAmino(part, self.elements[-1]))
            elif part[0] == '[':
                self.elements.append(MultipleAmino(part, self.elements[-1]))
            elif part[0] == '{':
                self.elements.append(NegatedAmino(part, self.elements[-1]))
            else:
                raise ValueError('Invalid pattern: %s' % pattern_string)
        self.head = self.elements[-1]
        assert self.head is not None

    def find(self, sequence: str, idx: int = 0) -> int:
        anchor_idx = self.find_anchor(sequence, idx)
        if anchor_idx < 0:
            return -1
        while anchor_idx >= 0:
            matches = self.head.match_including_following(sequence[anchor_idx:])
            if matches:
                return anchor_idx
            anchor_idx = self.find_anchor(sequence, anchor_idx + 1)
        return -1

    def find_anchor(self, sequence: str, idx: int = 0) -> int:
        length = len(sequence)
        if idx < 0:
            idx = length + idx
        while idx < length:
            if self.head.match(sequence[idx:]):
                return idx
            idx += 1
        return -1

    def find_all(self, sequence: str) -> List[int]:
        results = []
        match = self.find(sequence)
        while match >= 0:
            results.append(match)
            match = self.find(sequence, match + 1)
        return results


def parse_repeats(sequence: str) -> Tuple[int, int]:
    if not sequence:
        return 1, 1
    if not sequence.startswith("(") or not sequence.endswith(")"):
        raise ValueError("Brackets do not match")
    try:
        repeats = [int(number) for number in (sequence[1:-1]).split(",")]
    except ValueError:
        raise ValueError("Invalid repeat: %s" % sequence)
    if not 1 <= len(repeats) <= 2:
        raise ValueError("Invalid repeat")
    if len(repeats) == 1:
        return repeats[0], repeats[0]
    return repeats[0], repeats[1]


def parse_options(sequence: str, start: str, end: str) -> Tuple[Set[str], int, int]:
    if sequence[0] != start:
        raise ValueError("Brackets do not match")
    options = set()  # type: Set[str]
    idx = 1
    while idx < len(sequence) and sequence[idx] != end:
        if sequence[idx] not in AMINOS:
            raise ValueError("Invalid amino acid")
        options.add(sequence[idx])
        idx += 1
    if idx == len(sequence):
        raise ValueError("Brackets do not match")
    min_repeats, max_repeats = parse_repeats(sequence[idx+1:])
    if not options:
        raise ValueError("No valid options provided")
    return options, min_repeats, max_repeats
