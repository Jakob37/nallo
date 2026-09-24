import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "generate_gens_metadata.py"


class GensMetadataTest(unittest.TestCase):
    def run_metadata(self, with_upd=False):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        (root / "ref.fai").write_text("chr1\t100\t0\t0\t0\nchr2\t100\t0\t0\t0\n")
        (root / "sample.roh").write_text(
            "# ROH output\n"
            "RG\tchild\tchr1\t1\t20\t20\t5\t90\n"
            "RG\tchild\tchr1\t15\t30\t16\t5\t90\n"
            "RG\tchild\tchr2\t1\t10\t10\t5\t84\n"
            "RG\tother\tchr2\t1\t20\t20\t5\t90\n"
        )
        (root / "ratios.tsv").write_text(
            "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
            "chr1\t1\t10\t0\nchr1\t11\t20\t1\nchrX\t1\t10\t1\n"
        )
        args = [sys.executable, str(SCRIPT), "--sample", "child", "--sex", "M",
                "--fai", "ref.fai", "--roh", "sample.roh", "--ratios", "ratios.tsv"]
        if with_upd:
            (root / "regions.bed").write_text("chr1\t40\t50\tORIGIN=PATERNAL;TYPE=ISODISOMY\n")
            (root / "sites.bed").write_text(
                "chr1\t40\t41\tUPD_PATERNAL_ORIGIN\n"
                "chr1\t42\t43\tUPD_MATERNAL_ORIGIN\n"
                "chr1\t44\t45\tANTI_UPD\n"
                "chr1\t46\t47\tUNINFORMATIVE\n"
            )
            args += ["--upd-regions", "regions.bed", "--upd-sites", "sites.bed"]
        subprocess.run(args, cwd=root, check=True, capture_output=True, text=True)
        return root

    def test_coverage_roh_and_missing_upd(self):
        root = self.run_metadata()
        self.assertIn("%Autosomal LOH\t15.0", (root / "child.meta.tsv").read_text())
        chrom = (root / "child.chrom_meta.tsv").read_text()
        self.assertIn("1\tEstimated copy number\t2.83", chrom)
        self.assertIn("X\tEstimated copy number\t2.0", chrom)
        self.assertNotIn("2\tEstimated copy number", chrom)
        self.assertNotIn("Total SNPs", chrom)
        self.assertFalse((root / "child.gens_track.upd.bed").exists())
        self.assertIn("chr1\t0\t20\tLOH", (root / "child.gens_track.roh.bed").read_text())

    def test_trio_upd_follows_olwgs_non_informative_definition(self):
        root = self.run_metadata(with_upd=True)
        chrom = (root / "child.chrom_meta.tsv").read_text()
        self.assertIn("1\tTotal SNPs\t4", chrom)
        self.assertIn("1\tNon-informative\t3", chrom)
        self.assertIn("1\tNon-informative (%)\t75.0", chrom)
        self.assertIn("1\tMismatch mother\t1", chrom)
        self.assertIn("1\tMismatch father\t1", chrom)
        self.assertIn("1\tAnti-UPD\t1", chrom)
        self.assertIn("2\tTotal SNPs\t0", chrom)
        self.assertIn("chr1\t40\t50\tUniparental segment (PATERNAL)",
                      (root / "child.gens_track.upd.bed").read_text())

    def test_empty_upd_analysis_is_distinct_from_missing_upd(self):
        root = self.run_metadata(with_upd=True)
        (root / "sites.bed").write_text("")
        (root / "regions.bed").write_text("")
        subprocess.run(
            [sys.executable, str(SCRIPT), "--sample", "child", "--sex", "F",
             "--fai", "ref.fai", "--roh", "sample.roh", "--ratios", "ratios.tsv",
             "--upd-regions", "regions.bed", "--upd-sites", "sites.bed"],
            cwd=root, check=True, capture_output=True, text=True,
        )
        chrom = (root / "child.chrom_meta.tsv").read_text()
        self.assertIn("1\tTotal SNPs\t0", chrom)
        self.assertIn("X\tEstimated copy number\t4.0", chrom)
        self.assertTrue((root / "child.gens_track.upd.bed").exists())
        self.assertEqual((root / "child.gens_track.upd.bed").read_text(), "")

    def test_missing_autosomal_reference_length_fails(self):
        root = self.run_metadata()
        (root / "ref.fai").write_text("chr2\t100\t0\t0\t0\n")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--sample", "child", "--sex", "M",
             "--fai", "ref.fai", "--roh", "sample.roh", "--ratios", "ratios.tsv"],
            cwd=root, capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing reference length for ROH chromosome 1", result.stderr)


if __name__ == "__main__":
    unittest.main()
