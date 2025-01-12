"""Tests for the covplot module."""

import os
import tempfile
import unittest

import pandas as pd

from virmet.covplot import best_species


class TestCovplot(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.orgs_file = os.path.join(self.tmpdir, "orgs_list.tsv")
        df = pd.DataFrame(
            {
                "organism": [
                    "organism A",
                    "organism B",
                    "species A",
                    "species B",
                ],
                "reads": [40, 30, 20, 10],
            }
        )
        df.to_csv(self.orgs_file, sep="\t", header=True, index=False)

    def test_best_species(self):
        bsa = best_species(self.orgs_file, "orga")
        self.assertEqual(bsa, "organism A")

        bsa = best_species(self.orgs_file, "organism B")
        self.assertEqual(bsa, "organism B")

        bsa = best_species(self.orgs_file, "species ")
        self.assertEqual(bsa, "species A")
