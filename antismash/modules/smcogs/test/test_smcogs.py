# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.feature import GeneFunction, CDSFeature, FeatureLocation
from antismash.common.test import helpers
from antismash.common.hmmscan_refinement import HMMResult
from antismash.modules.smcogs import classify, SMCOGResults


class TestSMCOGLoad(unittest.TestCase):
    def test_load(self):
        # this mostly just tests that the cog annotation file isn't corrupted
        annotations = classify.load_cog_annotations()
        assert len(annotations) == 301
        for key, function in annotations.items():
            assert isinstance(function, GeneFunction), "cog annotation %s has bad type" % key


class TestAddingToRecord(unittest.TestCase):
    def test_classification_with_colon(self):
        # since SMCOG id and description are stored in a string separated by :,
        # ensure that descriptions containing : are properly handled
        cds = helpers.DummyCDS(locus_tag="test")
        record = helpers.DummyRecord(features=[cds], seq="A"*100)
        record.add_cluster(helpers.DummyCluster(0, 100))
        results = SMCOGResults(record.id)
        results.best_hits[cds.get_name()] = HMMResult("SMCOG1212:sodium:dicarboxylate_symporter",
                                                      0, 100, 2.3e-126, 416)
        results.add_to_record(record)
        gene_functions = cds.gene_functions.get_by_tool("smcogs")
        assert len(gene_functions) == 1
        assert str(gene_functions[0]).startswith("transport (smcogs) SMCOG1212:sodium:dicarboxylate_symporter"
                                                 " (Score: 416; E-value: 2.3e-126)")
