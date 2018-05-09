# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""A collection of functions for parsing protein motif patterns of a
    defined format and searching for them in an amino acid sequence.
"""

from typing import List, Set, Tuple  # pylint:disable=unused-import

AMINOS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
          'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}


class Match:
    """Holds information on whether or not an Element was found in a sequence
        and how long the match was.
        Evaluates to True if there was a hit, False otherwise.
        If there is a hit, a distance must always be provided.
    """
    def __init__(self, hit: bool, distance: int = -1) -> None:
        self.hit = hit
        if hit:
            assert distance > -1
        self._distance = distance

    @property
    def distance(self) -> int:
        """Get the distance of the Match."""
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
    """Holds information about the location (start and end) of an Element in a
        protein sequence.
    """
    def __init__(self, start: int, end: int) -> None:
        self.start = start
        self.end = end

    def __len__(self) -> int:
        return self.end - self.start

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "Match from {} to {} (length: {})".format(self.start, self.end, len(self))


class Element:
    """A base class for all elements that can form a Pattern."""
    def __init__(self, min_repeats: int, max_repeats: int, next_element: "Element" = None,
                 nterm: bool = False, cterm: bool = False) -> None:
        self.min_repeats = min_repeats
        self.max_repeats = max_repeats
        self.next_element = next_element
        if next_element and next_element.nterm:
            raise ValueError('Invalid pattern')
        if cterm and next_element:
            raise ValueError('Invalid pattern')
        self.nterm = nterm
        self.cterm = cterm

    def match(self, sequence: str, offset: int = 0) -> Match:
        raise NotImplementedError('Missing match implementation')

    def offset_ok(self, sequence: str, offset: int) -> bool:
        if not 0 <= offset < len(sequence):
            return False
        if self.nterm and offset >= self.max_repeats:
            return False
        if self.cterm and offset < len(sequence) - self.max_repeats:
            return False
        return True

    def match_including_following(self, sequence: str, offset: int = 0) -> Match:
        """Match the current Element and all subsequent Elements
            in a Pattern in the input sequence at the given offset.

            Arguments:
                  sequence: the sequence to be searched
                  offset: the position in the sequence at which to search

            Returns:
                  A Match: False if the Element couldn't be found,
                  True and the amount of repeats found if this is the last Element in the Pattern,
                  True and the distance of this and all following Elements otherwise.
        """
        for repeat in range(self.max_repeats, self.min_repeats - 1, -1):  # start with highest number repeats: greedy
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
        """Find the lengths of all possible versions of the Element and all subsquent
            Elements in the Pattern that match in the sequence at the given offset.

            Arguments:
                sequence: the sequence to be searched
                offset: the position in the sequence at which to search

            Returns:
                A list of all MatchLocations for the different versions of the Pattern that
                match the sequence. The start is always the offset, the end varies.
        """
        results = []
        for repeat in range(self.max_repeats, self.min_repeats - 1, -1):  # start with highest number repeats: greedy
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


class SimpleAmino(Element):  # TODO: logic for C terminus/N terminus parsing here! Can't have N as a next or a C with a next != None
    """An Element representing only one specific acceptable amino acid,
        as well as the times it must be repeated, and the next Element in
        the Pattern.
    """
    def __init__(self, element: str, next_element: Element = None, nterm: bool = False, cterm: bool = False) -> None:
        if element[0] not in AMINOS:
            raise ValueError('Invalid amino acid')
        self.amino = element[0]
        min_repeats, max_repeats = parse_repeats(element, 1)
        super().__init__(min_repeats, max_repeats, next_element, nterm, cterm)

    def match(self, sequence: str, offset: int = 0) -> Match:
        """Attempt to match the specified amino acid in an amino acid
            sequence at the given offset.

            Arguments:
                sequence: the amino acid sequence to match
                offset: the position at which to match

            Returns:
                A Match: True if the specified amino acid is at the position in the sequence,
                False if it isn't or the sequence is shorter than the specified offset.
        """
        return Match(self.offset_ok(sequence, offset) and sequence[offset] == self.amino, 1)

    def __str__(self) -> str:
        return "SimpleAmino(%s)(%d,%d)" % (self.amino, self.min_repeats, self.max_repeats)



class AnyAmino(Element):
    """An Element representing that any amino acid is acceptable at this position,
        as well as the number of characters that can be any amino acid, and the
        next Element in the Pattern.
    """
    def __init__(self, element: str, next_element: Element = None, nterm: bool = False, cterm: bool = False) -> None:
        if element[0] != 'x':
            raise ValueError('Attempting to use defined amino acid as AnyAmino')
        min_repeats, max_repeats = parse_repeats(element, 1)
        super().__init__(min_repeats, max_repeats, next_element, nterm, cterm)

    def match(self, sequence: str, offset: int = 0) -> Match:
        """Verify that any amino acid exists in the specified position in an
            amino acid sequence.

            Arguments:
                sequence: the amino acid sequence to match
                offset: the position at which to match

            Returns:
                  A Match - True if the offset is still within the sequence,
                  False otherwise.
        """
        # TODO: should invalid amino be an error?
        return Match(self.offset_ok(sequence, offset), 1)

    def __str__(self) -> str:
        return "AnyAmino(%d,%d)" % (self.min_repeats, self.max_repeats)


class MultipleAmino(Element):
    """An Element representing a set of amino acids that are acceptable at this
        position, as well as the number of characters that must satisfy this
        condition, and the next Element in the Pattern.
    """
    def __init__(self, element: str, next_element: Element = None, nterm: bool = False, cterm: bool = False) -> None:
        parsed_options = parse_options(element, "[", "]")
        self.options = parsed_options.options
        min_repeats = parsed_options.min_repeats
        max_repeats = parsed_options.max_repeats
        self.optional_c_terminus = parsed_options.optional_cterm
        super().__init__(min_repeats, max_repeats, next_element, nterm, cterm)

    def match(self, sequence: str, offset: int = 0) -> Match:
        """Attempt to match the acceptable amino acids in an amino acid
            sequence at the given offset.

            Arguments:
                sequence: the amino acid sequence to match
                offset: the position at which to match

            Returns:
                A Match: True if the amino acid at the given position is an acceptable amino acid,
                False if it isn't or the sequence is shorter than the specified offset.
        """
        if not 0 <= offset < len(sequence):
            return Match(self.optional_c_terminus and 0 <= offset <= len(sequence), 0)
        return Match(self.offset_ok(sequence, offset) and sequence[offset] in self.options, 1)

    def __str__(self) -> str:
        return "MultipleAmino(%s%s)(%d,%d)" % (self.options, ">" if self.optional_c_terminus else "",
                                               self.min_repeats, self.max_repeats)


class NegatedAmino(Element):
    """An Element representing a set of amino acids that are not acceptable at this
        position, as well as the number of characters that must satisfy this
        condition, and the next Element in the Pattern.
    """
    def __init__(self, element: str, next_element: Element = None, nterm: bool = False, cterm: bool = False) -> None:
        parsed_options = parse_options(element, "{", "}")
        self.options = parsed_options.options
        min_repeats = parsed_options.min_repeats
        max_repeats = parsed_options.max_repeats
        super().__init__(min_repeats, max_repeats, next_element, nterm, cterm)

    def match(self, sequence: str, offset: int = 0) -> Match:
        """Attempt to match the acceptable amino acids in an amino acid
            sequence at the given offset.

            Arguments:
                sequence: the amino acid sequence to match
                offset: the position at which to match

            Returns:
                A Match: True if the amino acid at the given position is not an unacceptable amino acid,
                False if it isn't or the sequence is shorter than the specified offset.
        """
        return Match(self.offset_ok(sequence, offset) and sequence[offset] not in self.options, 1)

    def __str__(self) -> str:
        return "MultipleAmino(%s)(%d,%d)" % (self.options, self.min_repeats, self.max_repeats)


class Pattern:
    """A pattern representing a conserved amino acid motif in a protein.
        The input string must end in a period and can only contain the following
        elements, separated by dashes:
        - simple IUPAC one-letter codes for amino acids when only this amino acid is
          acceptable in the position
        - lowercase x to represent any amino acid being acceptable in this position
        - multiple one-letter codes encased in square brackets to represent these
          amino acids all being acceptable in the position (e.g. [MA] matches M and A)
        - one or more one-letter codes encased in curly brackets to represent these
          amino acids not being acceptable in the position (e.g. {M} matches any but M)
        - < followed by an one-letter code in the beginning of the pattern to represent
          it must occur in the N-terminus of the sequence
        - a one-letter code followed by > at the end of the pattern to represent it must
          occur in the C-terminus of the sequence
        - alternatively, > inside square brackets to represent the end of the protein or
          any of the other amino acids specified inside the brackets
        Repeats may be specified for all of these by giving the number of repeats in
        parentheses, either as a single number or a range, e.g. K-I-T(2)-Y(1,3)>.
    """
    def __init__(self, pattern_string: str) -> None:
        assert pattern_string.endswith('.')
        self.elements = [None]  # type: List[Element]
        for part in pattern_string.strip('.').split('-')[::-1]:
            # preprocess part, remove <>
            nterm, part, cterm = parse_terminus(part)
            if part[0] in AMINOS:
                self.elements.append(SimpleAmino(part, self.elements[-1], nterm, cterm))
            elif part[0] == 'x':
                self.elements.append(AnyAmino(part, self.elements[-1], nterm, cterm))
            elif part[0] == '[':
                self.elements.append(MultipleAmino(part, self.elements[-1], nterm, cterm))
            elif part[0] == '{':
                self.elements.append(NegatedAmino(part, self.elements[-1], nterm, cterm))
            else:
                raise ValueError('Invalid pattern: %s' % pattern_string)
        self.head = self.elements[-1]
        assert self.head is not None

    def find(self, sequence: str) -> int:
        """Find the first occurrence of the Pattern in an amino acid sequence.

            Arguments:
                sequence: the amino acid sequence to search

            Returns:
                The position of the first occurrence of the Pattern,
                or -1 if there is no match.
        """
        anchor_idx = 0
        while anchor_idx < len(sequence):
            matches = self.head.match_including_following(sequence, anchor_idx)
            print(anchor_idx)
            if matches:
                print('Match found at', anchor_idx)
                return anchor_idx
            anchor_idx += 1
        return -1

    def find_all(self, sequence: str) -> List[MatchLocation]:
        """Find all occurrences of the Pattern in an amino acid sequence.

            Arguments:
                sequence: the amino acid sequence to search

            Returns:
                A list of MatchLocation representations of all locations of
                occurrences of the Pattern.
        """
        results = []
        index = 0
        while index < len(sequence):
            results.extend(self.head.match_all_possible(sequence, index))
            index += 1
        return results


class AmbiguityOptions:
    """A class holding results of parsing an element representing multiple acceptable
        or unacceptable amino acids.
        """
    def __init__(self, options: Set[str], min_repeats: int, max_repeats: int, optional_cterm: bool) -> None:
        assert options
        self.options = options
        self.min_repeats = min_repeats
        self.max_repeats = max_repeats
        self.optional_cterm = optional_cterm


def parse_repeats(sequence: str, offset: int) -> Tuple[int, int]:
    """Parse the number of repeats for an Element from a string.

        Arguments:
            sequence: the string representation of the Element
            offset: the position at which the information about the repeats starts.

        Returns:
            The minimum and maximum amount of repeats as a tuple;
            minimum and maximum are identical if there is only one value.
            If there are no repeats, returns 1, 1.
    """
    print(sequence, offset)
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


def parse_options(sequence: str, start: str, end: str) -> AmbiguityOptions:
    """Parse the amino acid options and repeats for a MultipleAmino or NegatedAmino
        Element from a string.

        Arguments:
            sequence: the string representation of the Element
            start: the opening bracket ([ or {)
            end: the closing bracket (] or })
            offset: the position at which the information about the options starts

        Returns:
            The set of amino acids specified by the string, as well as the number of repeats
            and whether or not the end of the protein can occur instead, summarized in an
            AmbiguityOptions instance.

    """
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
    return AmbiguityOptions(options, min_repeats, max_repeats, optional_c_terminus)


def parse_terminus(element: str) -> Tuple[bool, str, bool]:
    cterm = False
    nterm = False
    if element.startswith('<'):
        nterm = True
        element = element[1:]
    if element.endswith('>'):
        cterm = True
        element = element[:-1]
    return nterm, element, cterm
