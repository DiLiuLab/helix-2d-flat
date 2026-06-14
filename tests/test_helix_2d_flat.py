import os
import tempfile
import unittest

import helix_2D_flatV3 as helix


class SequenceTests(unittest.TestCase):
    def test_cleans_fasta_and_common_separators(self):
        text = ">example\nacg-t 12\n; ignored comment\n"
        self.assertEqual(helix.clean_sequence_text(text), "ACGT")

    def test_rejects_unsupported_letters(self):
        with self.assertRaisesRegex(ValueError, "Unsupported sequence"):
            helix.clean_sequence_text("ACGTZ")

    def test_detects_dna_and_rna(self):
        self.assertEqual(helix.detect_nucleic_acid("ACGT", "auto"), "dna")
        self.assertEqual(helix.detect_nucleic_acid("ACGU", "auto"), "rna")
        with self.assertRaisesRegex(ValueError, "both T and U"):
            helix.detect_nucleic_acid("AUT", "auto")

    def test_complements_iupac_bases(self):
        self.assertEqual(
            "".join(helix.complement_base(base, "dna") for base in "ACGTRYN"),
            "TGCAYRN",
        )
        self.assertEqual(helix.complement_base("A", "rna"), "U")


class DrawingTests(unittest.TestCase):
    def test_rows_are_antiparallel(self):
        options = helix.DrawingOptions(left_direction="up-to-bottom")
        rows = helix.make_rows("ACG", options)
        self.assertEqual([row.left_number for row in rows], [1, 2, 3])
        self.assertEqual([row.right_number for row in rows], [3, 2, 1])
        self.assertEqual([row.right_base for row in rows], ["T", "G", "C"])

    def test_svg_contains_expected_structure(self):
        options = helix.DrawingOptions(
            show_base_pair_lines=True,
            transparent=True,
            terminal_style="phosphate",
        )
        svg = helix.create_helix_svg("GC", options)
        self.assertIn("<svg ", svg)
        self.assertIn('<g id="bonds">', svg)
        self.assertIn('stroke-dasharray="5 5"', svg)
        self.assertNotIn('width="100%" height="100%"', svg)
        self.assertEqual(svg.count("<polygon "), 4)

    def test_write_svg_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = os.path.join(temp_dir, "nested", "duplex.svg")
            helix.write_svg_file("AT", output, helix.DrawingOptions())
            self.assertTrue(os.path.isfile(output))
            with open(output, encoding="utf-8") as handle:
                self.assertTrue(handle.read().startswith('<?xml version="1.0"'))

    def test_output_path_derivation(self):
        self.assertEqual(
            helix.derive_output_path(None, os.path.join("data", "seq.fasta")),
            os.path.join("data", "seq_helix2d.svg"),
        )
        self.assertEqual(helix.derive_output_path("figure", None), "figure.svg")


if __name__ == "__main__":
    unittest.main()
