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

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if not self.hit:
            return "Match(False)"
        return "Match(%s, %s)" % (self.hit, self.distance)


class MatchLocation:
    def __init__(self, start: int, end: int) -> None:
        self.start = start
        self.end = end

    @property
    def length(self) -> int:
        return self.end - self.start

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "Match from {} to {} (length: {})".format(self.start, self.end, self.length)


class Element:
    def __init__(self, min_repeats: int, max_repeats: int, next_element: "Element" = None) -> None:
        self.min_repeats = min_repeats
        self.max_repeats = max_repeats
        self.next_element = next_element

    def match(self, sequence: str, offset: int = 0) -> Match:
        raise NotImplementedError('Missing match implementation')

    def match_including_following(self, sequence: str, offset: int = 0) -> Match:
        repeats = list(range(self.max_repeats, self.min_repeats - 1, -1))
        for repeat in repeats:  # start with highest number repeats: greedy
            if repeat:
                if (len(sequence) - offset) < repeat - 1:  # if repeats are longer than remaining sequence...
                    continue  # ...try the next lowest number of repeats
                matches = True
                for idx in range(offset, offset + repeat):
                    match = self.match(sequence, idx)
                    matches = matches and match
                if not matches:  # if no match was found for the current number of repeats...
                    continue  # ...try the next lowest number of repeats
            if self.next_element:
                nextmatch = self.next_element.match_including_following(sequence, offset + repeat)
            else:
                return Match(True, repeat)
            if nextmatch:  # since Match instances are evaluated as bool here
                return Match(True, repeat + nextmatch.distance)
        return Match(False)

    def match_all_possible(self, sequence: str, offset: int = 0) -> List[MatchLocation]:
        results = []
        repeats = list(range(self.max_repeats, self.min_repeats - 1, -1))
        for repeat in repeats:  # start with highest number repeats: greedy
            if repeat:
                if (len(sequence) - offset) < repeat - 1:  # if repeats are longer than remaining sequence...
                    continue  # ...try the next lowest number of repeats
                matches = True
                for idx in range(offset, offset + repeat):
                    match = self.match(sequence, idx)
                    matches = matches and match
                    if not matches:
                        break
                if not matches:
                    continue
            if self.next_element:
                nextmatch = self.next_element.match_all_possible(sequence, offset + repeat)
                if nextmatch:
                    newmatch = [MatchLocation(offset, location.end) for location in nextmatch]
                    results.extend(newmatch)
            else:
                if offset + repeat <= len(sequence):
                    results.append(MatchLocation(offset, offset + repeat))
                # otherwise covered since it's end terminal that would have the
                # same coordinates
        for match in results:
            assert match.end <= len(sequence), "overlong slice: %s" % match
        return results


class SimpleAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        if element[0] not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.amino = element[0]
        min_repeats, max_repeats = parse_repeats(element, 1)
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str, offset: int = 0) -> Match:
        if not 0 <= offset < len(sequence):
            return Match(False)
        return Match(sequence[offset] == self.amino, 1)

    def __str__(self) -> str:
        return "SimpleAmino(%s)(%d,%d)" % (self.amino, self.min_repeats, self.max_repeats)


class NTerminalAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        if element[1] not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.amino = element[1]
        min_repeats, max_repeats = parse_repeats(element, 2)
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str, offset: int = 0) -> Match:
        if not 0 <= offset < len(sequence):
            return Match(False)
        return Match(sequence.startswith(self.amino) and sequence[offset] == self.amino, 1)

    def __str__(self) -> str:
        return "NTerminalAmino(%s)(%d,%d)" % (self.amino, self.min_repeats, self.max_repeats)


class CTerminalAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        if element[0] not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.amino = element[0]
        min_repeats, max_repeats = parse_repeats(element.strip('>'), 1)
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str, offset: int = 0) -> Match:
        if not 0 <= offset < len(sequence):
            return Match(False)
        return Match(sequence.endswith(self.amino) and offset == (len(sequence) - 1), 1)

    def __str__(self) -> str:
        return "CTerminalAmino(%s)(%d,%d)" % (self.amino, self.min_repeats, self.max_repeats)


class AnyAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        if element[0] != 'x':
            raise ValueError('Attempting to use defined amino acid as AnyAmino')
        min_repeats, max_repeats = parse_repeats(element, 1)
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str, offset: int = 0) -> Match:
        # TODO: should invalid amino be an error?
        return Match(len(sequence) > offset, 1)

    def __str__(self) -> str:
        return "AnyAmino(%d,%d)" % (self.min_repeats, self.max_repeats)


class MultipleAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        self.options, min_repeats, max_repeats, self.optional_c_terminus = parse_options(element, "[", "]")
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str, offset: int = 0) -> Match:
        if not 0 <= offset < len(sequence):
            return Match(self.optional_c_terminus and 0 <= offset <= len(sequence), 0)
        return Match(sequence[offset] in self.options, 1)

    def __str__(self) -> str:
        return "MultipleAmino(%s%s)(%d,%d)" % (self.options, ">" if self.optional_c_terminus else "",
                                               self.min_repeats, self.max_repeats)


class NegatedAmino(Element):
    def __init__(self, element: str, next_element: Element = None) -> None:
        self.options, min_repeats, max_repeats, self.optional_c_terminus = parse_options(element, "{", "}")
        super().__init__(min_repeats, max_repeats, next_element)

    def match(self, sequence: str, offset: int = 0) -> Match:
        if not 0 <= offset < len(sequence):
            return Match(False)
        return Match(sequence[offset] not in self.options, 1)

    def __str__(self) -> str:
        return "MultipleAmino(%s)(%d,%d)" % (self.options, self.min_repeats, self.max_repeats)


class Pattern:
    def __init__(self, pattern_string: str) -> None:
        assert pattern_string.endswith('.')
        self.elements = [None]  # type: List[Element]
        self.n_terminus = None  # type: NTerminalAmino
        self.c_terminus = None  # type: CTerminalAmino
        for part in pattern_string.strip('.').split('-')[::-1]:
            if part[-1] == '>':
                if self.elements[-1] is None and not self.c_terminus:  # only if it's last element
                    self.c_terminus = CTerminalAmino(part, self.elements[-1])
                    self.elements.append(self.c_terminus)
                else:
                    raise ValueError('Invalid pattern: %s' % pattern_string)
            elif part[0] in AMINOS:
                self.elements.append(SimpleAmino(part, self.elements[-1]))
            elif part[0] == 'x':
                self.elements.append(AnyAmino(part, self.elements[-1]))
            elif part[0] == '[':
                self.elements.append(MultipleAmino(part, self.elements[-1]))
            elif part[0] == '{':
                self.elements.append(NegatedAmino(part, self.elements[-1]))
            elif part[0] == '<' and pattern_string.startswith('<') and not self.n_terminus:
                self.n_terminus = NTerminalAmino(part, self.elements[-1])
                self.elements.append(self.n_terminus)
            else:
                raise ValueError('Invalid pattern: %s' % pattern_string)
        self.head = self.elements[-1]
        assert self.head is not None

    def find(self, sequence: str) -> int:
        if self.n_terminus and not sequence.startswith(self.n_terminus.amino):
            return -1
        if self.c_terminus and not sequence.endswith(self.c_terminus.amino):
            return -1
        anchor_idx = 0
        while anchor_idx < len(sequence):
            matches = self.head.match_including_following(sequence, anchor_idx)
            if matches:
                return anchor_idx
            anchor_idx += 1
            if self.n_terminus:
                break
        return -1

    def find_all(self, sequence: str) -> List[MatchLocation]:
        results = []
        index = 0
        while index < len(sequence):
            results.extend(self.head.match_all_possible(sequence, index))
            index += 1
            if self.n_terminus:
                break
        return results


def parse_repeats(sequence: str, offset: int) -> Tuple[int, int]:
    if (len(sequence) - 1) < offset:
        return 1, 1
    if sequence[offset] != "(" or not sequence.endswith(")"):
        raise ValueError("Brackets do not match")
    try:
        repeats = [int(number) for number in (sequence[offset + 1:-1]).split(",")]
    except ValueError:
        raise ValueError("Invalid repeat: %s" % sequence)
    if not 1 <= len(repeats) <= 2:
        raise ValueError("Invalid repeat")
    if len(repeats) == 1:
        return repeats[0], repeats[0]
    return repeats[0], repeats[1]


def parse_options(sequence: str, start: str, end: str) -> Tuple[Set[str], int, int, bool]:  # TODO: make class!
    if sequence[0] != start:
        raise ValueError("Brackets do not match")
    options = set()  # type: Set[str]
    optional_c_terminus = False
    idx = 1
    while idx < len(sequence) and sequence[idx] != end:
        if sequence[idx] not in AMINOS:
            if sequence[idx] == ">" and not optional_c_terminus:
                optional_c_terminus = True
                idx += 1
                continue
            else:
                raise ValueError("Invalid amino acid")
        options.add(sequence[idx])
        idx += 1
    if idx == len(sequence):
        raise ValueError("Brackets do not match")
    min_repeats, max_repeats = parse_repeats(sequence, idx+1)
    if not options:
        raise ValueError("No valid options provided")
    return options, min_repeats, max_repeats, optional_c_terminus
