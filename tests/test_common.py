"""Tests for the common module."""

import os
import tempfile
import unittest

from virmet.common import (
    bact_fung_query,
    ftp_down,
    get_gids,
    multiple_download,
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
        # big file, 39 MB
        # self.remote_1 = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz'
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

    def test_gids(self):
        ids = get_gids(self.fasta.name)
        self.assertTrue("1234" in ids)
        self.assertTrue("ABCD" in ids)

    def test_bact_query(self):
        all_urls = bact_fung_query(query_type="bacteria", download=True, target_folder=self.tmpdir)
        bac_lines = len(all_urls)
        self.assertGreater(bac_lines, 100)
        self.assertTrue(os.path.exists(os.path.join(self.tmpdir,"bacteria_refseq_info.tsv")))
        os.rename(os.path.join(self.tmpdir, "bacteria_refseq_info.tsv"), os.path.join(self.tmpdir, "xyz.tsv"))
        urls_again = bact_fung_query(
            query_type="bacteria", download=False, info_file="xyz.tsv", target_folder=self.tmpdir
        )
        bac_lines_again = len(urls_again)
        self.assertEqual(bac_lines, bac_lines_again)
        os.remove(os.path.join(self.tmpdir, "xyz.tsv"))

    def test_fung_query(self):
        all_urls = bact_fung_query(query_type="fungi", download=True, target_folder=self.tmpdir)
        fung_lines = len(all_urls)
        self.assertGreater(fung_lines, 10)
        self.assertTrue(os.path.exists(os.path.join(self.tmpdir, "fungi_refseq_info.tsv")))
        os.rename(os.path.join(self.tmpdir,"fungi_refseq_info.tsv"), os.path.join(self.tmpdir, "xyz.tsv"))
        urls_again = bact_fung_query(
            query_type="fungi", download=False, info_file="xyz.tsv", target_folder=self.tmpdir
        )
        fung_lines_again = len(urls_again)
        self.assertEqual(fung_lines, fung_lines_again)
        os.remove(os.path.join(self.tmpdir, "xyz.tsv"))

    def test_multi_download(self):
        tmpf = os.path.join(self.tmpdir, "tmp_multi_down.txt")
        open(tmpf, "w").close()
        # download same file twice
        dl_pair = tmpf, [self.remote_2, self.remote_2]
        multiple_download(dl_pair)
        self.assertTrue(os.path.exists(tmpf))
        with open(tmpf) as f:
            lines = list(f)
        os.remove(tmpf)
        m = len(lines)
        # even number of lines
        self.assertEqual(m % 2, 0)
        halfm = int(m / 2)
        self.assertEqual(lines[0], lines[halfm])
        self.assertEqual(lines[1], lines[1 + halfm])


if __name__ == "__main__":
    print("xyz")
    unittest.main()
