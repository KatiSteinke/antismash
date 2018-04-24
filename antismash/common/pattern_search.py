# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""A collection of functions for parsing protein motif patterns of a defined format and searching for them.
"""

from typing import List

# result classes: base class? Check validity of input in init! Then only check type of element later

AMINOS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}


class Element:
    pass


class SimpleAmino(Element):
    def __init__(self, element: str) -> None:
        if element not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.element = element

    def match(self, sequence: str) -> bool:
        return sequence[0] == self.element


class AnyAmino(Element):
    def __init__(self, element: str) -> None:
        if element != 'x':
            raise ValueError('Attempting to use defined amino acid as AnyAmino')


class MultipleAmino(Element):
    def __init__(self, element: str) -> None:
        self.element = set()
        length = len(element)

        if length < 3:
            raise ValueError('Invalid string for MultipleAmino')

        if element[0] != "[":
            raise ValueError('Invalid string for MultipleAmino')

        idx = 1
        while idx < length:
            char = element[idx]
            if char in AMINOS:
                self.element.add(char)
            elif char == "]":
                break
            else:
                raise ValueError('Invalid amino acid')
            idx += 1

        if idx == length:
            raise ValueError('Brackets do not match')
        if idx < length - 1:
            raise ValueError('Invalid string for MultipleAmino')


class NegatedAmino(Element):
    def __init__(self, element: str) -> None:
        self.element = set()
        length = len(element)

        if length < 3:
            raise ValueError('Invalid string for NegatedAmino')

        if element[0] != "{":
            raise ValueError('Invalid string for NegatedAmino')

        idx = 1
        while idx < length:
            char = element[idx]
            if char in AMINOS:
                self.element.add(char)
            elif char == "}":
                break
            else:
                raise ValueError('Invalid amino acid')
            idx += 1

        if idx == length:
            raise ValueError('Brackets do not match')
        if idx < length - 1:
            raise ValueError('Invalid string for NegatedAmino')


class Pattern:
    def __init__(self, pattern_string: str) -> None:
        self.elements = []
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
                idx += 1
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
