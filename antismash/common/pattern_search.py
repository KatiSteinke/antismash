# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""A collection of functions for parsing protein motif patterns of a defined format and searching for them in an amino
    acid sequence.
"""

from typing import List, Set

AMINOS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}


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
    def match(self, sequence: str) -> Match:
        raise NotImplementedError('Missing match implementation')


class SimpleAmino(Element):
    def __init__(self, element: str) -> None:
        if element not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.element = element
        self.repeats = 1

    def match(self, sequence: str) -> Match:
        if not len(sequence) >= self.repeats:
            return Match(False)
        return Match(sequence[0] == self.element, self.repeats)


class AnyAmino(Element):
    def __init__(self, element: str) -> None:
        if element[0] != 'x':
            raise ValueError('Attempting to use defined amino acid as AnyAmino')
        self.repeats = parse_repeats(element[1:])

    def match(self, sequence: str) -> Match:
        # TODO: should invalid amino be an error?
        return Match(len(sequence) >= self.repeats, self.repeats)


class MultipleAmino(Element):
    def __init__(self, element: str) -> None:
        self.options = set()  # type: Set[str]
        length = len(element)

        if length < 3:
            raise ValueError('Invalid string for MultipleAmino')

        if element[0] != "[":
            raise ValueError('Invalid string for MultipleAmino')

        idx = 1
        while idx < length:
            char = element[idx]
            if char in AMINOS:
                self.options.add(char)
            elif char == "]":
                break
            else:
                raise ValueError('Invalid amino acid')
            idx += 1

        if idx == length:
            raise ValueError('Brackets do not match')
        self.repeats = parse_repeats(element[idx+1:])

    def match(self, sequence: str) -> Match:
        if not len(sequence) >= self.repeats:
            return Match(False)
        return Match(set(sequence[:self.repeats]).issubset(self.options), self.repeats)


class NegatedAmino(Element):
    def __init__(self, element: str) -> None:
        self.options = set()  # type: Set[str]
        length = len(element)

        if length < 3:
            raise ValueError('Invalid string for NegatedAmino')

        if element[0] != "{":
            raise ValueError('Invalid string for NegatedAmino')

        idx = 1
        while idx < length:
            char = element[idx]
            if char in AMINOS:
                self.options.add(char)
            elif char == "}":
                break
            else:
                raise ValueError('Invalid amino acid')
            idx += 1

        if idx == length:
            raise ValueError('Brackets do not match')
        self.repeats = parse_repeats(element[idx+1:])

    def match(self, sequence: str) -> Match:
        if not len(sequence) >= self.repeats:
            return Match(False)
        return Match(set(sequence[:self.repeats]).isdisjoint(self.options), self.repeats)


class Pattern:
    def __init__(self, pattern_string: str) -> None:
        self.elements = []  # type: List[Element]
        for part in pattern_string.split('-'):
            if part[0] in AMINOS:
                self.elements.append(SimpleAmino(part))
            elif part[0] == 'x':
                self.elements.append(AnyAmino(part))
            elif part[0] == '[':
                self.elements.append(MultipleAmino(part))
            elif part[0] == '{':
                self.elements.append(NegatedAmino(part))
            else:
                raise ValueError('Invalid pattern: %s' % pattern_string)

    def find(self, sequence: str) -> int:
        anchor_idx = self.find_anchor(sequence)
        if anchor_idx < 0:
            return -1
        while anchor_idx >= 0:
            matches = True
            idx = anchor_idx
            for element in self.elements:
                match = element.match(sequence[idx:])
                matches = matches and match
                if not match:
                    break
                idx += match.distance
            if matches:
                return anchor_idx
            anchor_idx = self.find_anchor(sequence, anchor_idx+1)
        return -1

    def find_anchor(self, sequence: str, idx: int = 0) -> int:
        anchor_element = self.elements[0]
        length = len(sequence)
        if idx < 0:
            idx = length + idx
        while idx < length:
            if anchor_element.match(sequence[idx:]):
                return idx
            idx += 1
        return -1


def parse_repeats(sequence: str) -> int:
    if not sequence:
        return 1
    if not sequence.startswith("(") or not sequence.endswith(")"):
        raise ValueError("Brackets do not match")
    return int(sequence[1:-1])