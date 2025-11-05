"""Tests for the common module."""

import os
import tempfile
import unittest

from datetime import date, timedelta

from virmet.common import (
    viral_query,
    ftp_down,
    run_child,
)


def parse_file_line(line):
    ftl = line.strip().split(":")[1][1:]
    if ftl.startswith("gzip compressed") and "extra field" in ftl:
        return "bgzipped"
    if ftl.startswith("gzip compressed") and "extra field" not in ftl:
        return "gzipped"
    if ftl.startswith("ASCII"):
        return "ascii"


class TestFTPDown(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        # small file, 335 KB
        self.remote_1 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.2wayconspseudos.gtf.gz"
        # again small file
        self.remote_2 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/_README.TXT"

    def test_nodecompress(self):
        out_file = os.path.join(self.tmpdir, "gtf.txt.gz")
        ftp_down(self.remote_1, out_file)
        ftl = run_child(f"file {out_file}")
        os.remove(out_file)
        ft = parse_file_line(ftl)
        self.assertEqual(ft, "gzipped")

    def test_decompress(self):
        out_file = os.path.join(self.tmpdir, "gtf.txt")
        ftp_down(self.remote_1, out_file)
        ftl = run_child(f"file {out_file}")
        os.remove(out_file)
        ft = parse_file_line(ftl)
        self.assertEqual(ft, "ascii")

    def test_append(self):
        out_file = os.path.join(self.tmpdir, "README.TXT")
        open(out_file, "w").close()
        ftp_down(self.remote_2, out_file)
        with open(out_file) as f:
            n_lines_1 = sum(1 for _ in f)
        ftp_down(self.remote_2, out_file)
        with open(out_file) as f:
            n_lines_2 = sum(1 for _ in f)
        os.remove(out_file)
        self.assertEqual(n_lines_2, 2 * n_lines_1)


class TestMisc(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.fasta = open(os.path.join(self.tmpdir, "tmp.fasta"), "w")
        self.fasta.write(">gi|1234|xyz\nAGCTAGC\n>gi|ABCD\nATCG\n")
        self.fasta.close()
        self.remote_2 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/_README.TXT"

    def test_viral_query(self):
        all_accs = viral_query(
            DB_DIR_UPDATE=self.tmpdir,
            viral_db="n",
            update_min_date=(date.today() - timedelta(days=10)).strftime(
                "%Y/%m/%d"
            ),
        )

        vir_acc = len(all_accs)
        self.assertGreater(vir_acc, 2)
        accs_again = viral_query(
            DB_DIR_UPDATE=self.tmpdir,
            viral_db="n",
            update_min_date=(date.today() - timedelta(days=10)).strftime(
                "%Y/%m/%d"
            ),
        )
        vir_acc_again = len(accs_again)
        self.assertEqual(vir_acc, vir_acc_again)


if __name__ == "__main__":
    print("xyz")
    unittest.main()
