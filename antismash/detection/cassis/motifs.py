# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Motif related functions and classes for CASSIS """

from collections import defaultdict
import csv
import logging
import os
from typing import Any, Dict, List, Optional, Set
from xml.etree import cElementTree as ElementTree

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from antismash.detection import cassis

from .pairings import PROMOTER_RANGE, Pairing
from .promoters import Promoter


class Motif(Pairing):
    """ An individual motif, tracks a location as a plus/minus pair, an evalue (score)
        and any hits.
    """
    def __init__(self, plus: int, minus: int, score: Optional[float] = None,
                 hits: Optional[Dict[str, int]] = None) -> None:
        super().__init__(plus, minus)
        self._score = None  # type: Optional[float]
        if score:
            self.score = score
        self.seqs = []  # type: List[str]
        self.hits = defaultdict(lambda: 0)  # type: Dict[str, int]
        if hits:
            for key, val in hits.items():
                self.hits[str(key)] = int(val)

    @property
    def score(self) -> Optional[float]:
        """ The score/evalue of a motif """
        return self._score

    @score.setter
    def score(self, score: float) -> None:
        self._score = float(score)

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, Motif)
                and self.plus == other.plus
                and self.minus == other.minus
                and self._score == other.score
                and self.seqs == other.seqs
                and self.hits == other.hits)

    def __repr__(self) -> str:
        return "Motif(%s, score=%s, hits=%s)" % (self.pairing_string, self.score, self.hits)


def generate_motifs(meme_dir: str, anchor_promoter: int, promoters: List[Promoter]) -> List[Motif]:
    """Prepare sets of promoter sequences and motif subdirectories"""
    promoter_sets = []

    if not os.path.exists(meme_dir):
        os.makedirs(meme_dir)

    # prepare sets of promoter sequences (MEME input)
    indices = set()  # type: Set[Tuple[int, int]] # to monitor unique start_index/end_index
    for pm in PROMOTER_RANGE:
        start_index = anchor_promoter - pm.minus
        end_index = anchor_promoter + pm.plus

        if start_index < 0:  # anchor promoter near beginning of record --> truncate
            if cassis.VERBOSE_DEBUG:
                logging.debug("Promoter set " + pm.pairing_string + " exceeds upstream record border")
            start_index = 0

        if end_index > len(promoters) - 1:  # anchor promoter near end of record --> truncate
            if cassis.VERBOSE_DEBUG:
                logging.debug("Promoter set " + pm.pairing_string + " exceeds downstream record border")
            end_index = len(promoters) - 1

        # discard promoter sets, which reappear due to truncation
        if (start_index, end_index) not in indices:
            indices.add((start_index, end_index))

            # check (again, compare init of PROMOTER_RANGE) if the promoter set has at least 4 promoters
            if end_index - start_index + 1 >= 4:
                promoter_sets.append(Motif(pm.plus, pm.minus))

                pm_dir = os.path.join(meme_dir, pm.pairing_string)
                if not os.path.exists(pm_dir):
                    os.makedirs(pm_dir)

                # write promoter sequences to fasta file, in respective "plus-minus" subdir
                with open(os.path.join(pm_dir, "promoters.fasta"), "w") as pm_handle:
                    for i in range(start_index, end_index + 1):
                        seq = SeqRecord(promoters[i].seq,
                                        id=promoters[i].get_id(),
                                        description="length={}bp".format(len(promoters[i].seq)))
                        if i == anchor_promoter:  # mark anchor gene
                            seq.id += "__ANCHOR"  # must be part of id, otherwise MEME woun't recognize it
                        SeqIO.write(seq, pm_handle, "fasta")
        else:
            if cassis.VERBOSE_DEBUG:
                logging.debug("Duplicate promoter set " + pm.pairing_string)

    return promoter_sets


def filter_meme_results(meme_dir: str, promoter_sets: List[Motif], anchor: str) -> List[Motif]:
    """Analyse and filter MEME results"""
    for motif in promoter_sets:
        xml_file = os.path.join(meme_dir, motif.pairing_string, "meme.xml")
        root = ElementTree.parse(xml_file).getroot()
        reason = root.find("model/reason_for_stopping").text
        anchor_seq_id = ""

        # no motif found for given e-value cutoff :-(
        if "Stopped because motif E-value > " in reason:
            if cassis.VERBOSE_DEBUG:
                logging.debug("MEME: motif %s; e-value exceeds cutoff", motif.pairing_string)

        # motif(s) found :-)
        elif "Stopped because requested number of motifs (1) found" in reason:
            # find anchor genes' sequence_id
            training_set = root.findall("training_set/sequence")  # all promoter sequences passed to MEME
            for element in training_set:
                if "__ANCHOR" in element.attrib["name"]:
                    anchor_seq_id = element.attrib["id"]  # e.g. id=sequence_1

            # only accept motifs which occur in the anchor genes promoter
            # sequences which contributed to the motif
            contributing_sites = root.findall("motifs/motif/contributing_sites/contributing_site")
            if anchor_seq_id in map(lambda site: site.attrib["sequence_id"], contributing_sites):
                # save motif score
                motif.score = float(root.find("motifs/motif").attrib["e_value"])  # one motif, didn't ask MEME for more

                # save sequence sites which represent the motif
                motif.seqs = ["".join(map(lambda letter: letter.attrib["letter_id"],
                                          site.findall("site/letter_ref")))
                              for site in contributing_sites]
                # write sites to fasta file
                with open(os.path.join(meme_dir, str(motif),
                                       "binding_sites.fasta"), "w") as handle:
                    handle.write(">{}__{}\n".format(anchor, str(motif)))
                    handle.write("\n".join(motif.seqs))
                if cassis.VERBOSE_DEBUG:
                    logging.debug("MEME: motif %s; e-value = %s", motif, motif.score)
            else:
                if cassis.VERBOSE_DEBUG:
                    logging.debug("MEME: motif %s; does not occur in anchor gene promoter", motif)

        # unexpected reason, don't know why MEME stopped :-$
        else:
            logging.error("MEME stopped unexpectedly (reason: " + reason + ")")

    return [motif for motif in promoter_sets if motif.score is not None]


def filter_fimo_results(motifs: List[Motif], fimo_dir: str, promoters: List[Promoter],
                        anchor_promoter: int) -> List[Motif]:
    """Analyse and filter FIMO results"""

    filtered = []

    for motif in motifs:
        assert isinstance(motif, Motif), type(motif)
        with open(os.path.join(fimo_dir, str(motif), "fimo.txt"), "r") as handle:
            table_reader = csv.reader(handle, delimiter="\t")
            for row in table_reader:
                # skip comment lines
                if row[0].startswith("#"):
                    continue
                seq_id = row[1]
                motif.hits[seq_id] += 1

        # write binding sites per promoter to file
        with open(os.path.join(fimo_dir, str(motif), "bs_per_promoter.csv"), "w") as handle:
            table_writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            table_writer.writerow(["#", "promoter", "binding sites"])  # table head
            for i, promoter in enumerate(promoters):
                table_writer.writerow([i+1, promoter.get_id(), motif.hits.get(promoter.get_id(), 0)])

        percentage = len(motif.hits) / len(promoters) * 100
        if percentage == 0.0:
            # too low
            if cassis.VERBOSE_DEBUG:
                logging.debug("FIMO: motif %s; occurs in %d promoters (no hits)",
                              motif, len(motif.hits))
            continue
        elif percentage > cassis.MAX_PERCENTAGE:
            # too high
            if cassis.VERBOSE_DEBUG:
                logging.debug("FIMO: %s; occurs in %d promoters; %.2f%% of all promoters (too many)",
                              motif, len(motif.hits), percentage)
            continue
        elif promoters[anchor_promoter].get_id() not in motif.hits:  # not in achor promoter
            # no site in anchor promoter
            if cassis.VERBOSE_DEBUG:
                logging.debug("FIMO: motif %s; no hits in the promoter of the anchor gene", motif)
            continue

        # everything ok
        if cassis.VERBOSE_DEBUG:
            logging.debug("FIMO: motif %s; occurs in %d promoters; %.2f%% of all promoters",
                          motif, len(motif.hits), percentage)
        filtered.append(motif)

    return filtered
