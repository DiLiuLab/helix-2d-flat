import os
import unittest
import xml.etree.ElementTree as ET

import helix_2D_flatV3 as helix


SVG_NS = "http://www.w3.org/2000/svg"
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURE_PATH = os.path.join(PROJECT_ROOT, "assets", "parameter_definitions.svg")


def tag(name):
    return "{{{}}}{}".format(SVG_NS, name)


def normalized_element(element):
    return (
        element.tag,
        tuple(sorted(element.attrib.items())),
        (element.text or "").strip(),
        tuple(normalized_element(child) for child in element),
    )


class ParameterFigureTests(unittest.TestCase):
    def setUp(self):
        self.figure = ET.parse(FIGURE_PATH).getroot()

    def test_molecule_matches_script_generated_svg(self):
        options = helix.DrawingOptions(
            show_numbers=False,
            terminal_labels=False,
            show_base_pair_lines=False,
            transparent=True,
        )
        generated = ET.fromstring(helix.create_helix_svg("AG", options))
        figure_molecule = self.figure.find(
            "{}[@id='script-generated-duplex']".format(tag("g"))
        )
        self.assertIsNotNone(figure_molecule)

        for group_id in ("bonds", "nucleotides", "terminal-labels"):
            expected = generated.find("{}[@id='{}']".format(tag("g"), group_id))
            actual = figure_molecule.find(
                "{}[@id='{}']".format(tag("g"), group_id)
            )
            self.assertEqual(
                normalized_element(actual),
                normalized_element(expected),
                group_id,
            )

    def test_only_geometry_parameters_are_labeled(self):
        expected = {
            "base_bond_length",
            "base_height",
            "pair_gap",
            "phosphate_kink_length",
            "phosphate_radius",
            "purine_width",
            "pyrimidine_width",
            "row_spacing",
            "sugar_radius",
        }
        labels = {
            (element.text or "").strip()
            for element in self.figure.findall(".//{}".format(tag("text")))
            if element.attrib.get("class") == "measure-label"
        }
        self.assertEqual(labels, expected)

    def test_vertical_dimensions_have_vertical_labels(self):
        labels = {
            (element.text or "").strip(): element
            for element in self.figure.findall(".//{}".format(tag("text")))
        }
        for label in ("base_height", "row_spacing"):
            self.assertIn("rotate(-90", labels[label].attrib.get("transform", ""))


if __name__ == "__main__":
    unittest.main()
