#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import convert_to_division_model_averaging as converter


SAMPLE_XML = """<?xml version="1.0"?>
<beast>
    <cenancestorTreeLikelihood id="treeLikelihood" divisionModel="fission" heightRules="true">
        <patterns idref="alignment"/>
        <treeModel idref="treeModel"/>
        <siteModel idref="siteModel"/>
        <strictClockCenancestorBranchRates idref="branchRates"/>
    </cenancestorTreeLikelihood>
    <operators id="operators">
        <scaleOperator weight="0.1"><parameter idref="clock.rate"/></scaleOperator>
    </operators>
    <mcmc id="mcmc" operatorAnalysis="results/run.ops">
        <posterior id="posterior">
            <prior id="prior">
                <uniformPrior lower="0" upper="1"><parameter idref="clock.rate"/></uniformPrior>
            </prior>
            <likelihood id="likelihood">
                <cenancestorTreeLikelihood idref="treeLikelihood"/>
            </likelihood>
        </posterior>
        <operators idref="operators"/>
        <log id="screenLog"><posterior idref="posterior"/></log>
        <log id="filelog" fileName="results/run.multi.log"><posterior idref="posterior"/></log>
        <logTree id="treeFileLog" fileName="result"><treeModel idref="treeModel"/></logTree>
    </mcmc>
    <!-- Power-posterior sampling for MLE estimation -->
    <marginalLikelihoodEstimator id="MLE">
        <samplers><mcmc idref="mcmc"/></samplers>
        <log id="MLELog" fileName="results/run.MLE.log"><posterior idref="posterior"/></log>
    </marginalLikelihoodEstimator>
    <!-- "MLE" estimation using the harmonic mean estimator -->
    <harmonicMeanAnalysis fileName="results/run.multi.log">
        <likelihoodColumn name="likelihood"/>
    </harmonicMeanAnalysis>
</beast>
"""


class ConversionTests(unittest.TestCase):
    def test_default_conversion(self):
        converted = converter.convert_xml(SAMPLE_XML)
        ET.fromstring(converted)

        self.assertIn("<modelAveragingCenancestorTreeLikelihood", converted)
        self.assertNotIn('divisionModel="fission"', converted)
        self.assertIn('id="divisionModel" value="0" lower="0" upper="3"', converted)
        self.assertIn('uniformIntegerOperator lower="0" upper="3" weight="0.5"', converted)
        self.assertIn('divisionModelPrior weights="1.0 1.0 1.0 1.0"', converted)
        self.assertEqual(converted.count('<parameter idref="divisionModel"/>'), 4)
        self.assertNotIn("marginalLikelihoodEstimator", converted)
        self.assertNotIn("Power-posterior sampling for MLE estimation", converted)
        self.assertIn("harmonicMeanAnalysis", converted)
        self.assertIn('operatorAnalysis="results/run.divisionModelAveraged.ops"', converted)
        self.assertEqual(converted.count('fileName="results/run.multi.divisionModelAveraged.log"'), 2)
        self.assertIn('fileName="result.divisionModelAveraged"', converted)

    def test_custom_values_and_named_initial_models(self):
        for name, expected in converter.MODEL_VALUES.items():
            with self.subTest(name=name):
                converted = converter.convert_xml(
                    SAMPLE_XML,
                    operator_weight=2.25,
                    prior_weights=(0.0, 1.0, 2.0, 3.0),
                    initial_model=converter.parse_initial_model(name),
                    preserve_mle=True,
                )
                self.assertIn(f'id="divisionModel" value="{expected}"', converted)
                self.assertIn('weight="2.25"', converted)
                self.assertIn('weights="0.0 1.0 2.0 3.0"', converted)
                self.assertIn("marginalLikelihoodEstimator", converted)
                self.assertIn('fileName="results/run.MLE.divisionModelAveraged.log"', converted)
                self.assertIn("harmonicMeanAnalysis", converted)

    def test_output_target_renaming(self):
        cases = {
            "run.log": "run.divisionModelAveraged.log",
            "run.MLE.log": "run.MLE.divisionModelAveraged.log",
            "/tmp/a/run.trees": "/tmp/a/run.divisionModelAveraged.trees",
            "result": "result.divisionModelAveraged",
            ".hidden": ".hidden.divisionModelAveraged",
            "run.divisionModelAveraged.log": "run.divisionModelAveraged.log",
        }
        for source, expected in cases.items():
            with self.subTest(source=source):
                self.assertEqual(converter.rename_output_target(source), expected)

    def test_rejects_already_converted_or_ambiguous_xml(self):
        converted = converter.convert_xml(SAMPLE_XML)
        with self.assertRaisesRegex(converter.ConversionError, "already uses"):
            converter.convert_xml(converted)
        with self.assertRaisesRegex(converter.ConversionError, "found 0"):
            converter.convert_xml("<beast/>")

    def test_argument_validation(self):
        for value in ("identity", "budding", "fission", "split-fission", "0", "3"):
            converter.parse_initial_model(value)
        for value in ("-1", "4", "unknown"):
            with self.subTest(value=value), self.assertRaises(Exception):
                converter.parse_initial_model(value)
        for value in ("0", "-1", "nan", "inf"):
            with self.subTest(value=value), self.assertRaises(Exception):
                converter.parse_operator_weight(value)
        for value in ("1 2 3", "0 0 0 0", "1 -1 1 1", "1 1 inf 1"):
            with self.subTest(value=value), self.assertRaises(Exception):
                converter.parse_prior_weights(value)


class CommandLineTests(unittest.TestCase):
    def run_converter(self, extra_arguments):
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            input_path = directory_path / "input.xml"
            output_path = directory_path / "output.xml"
            input_path.write_text(SAMPLE_XML, encoding="utf-8")
            command = [
                sys.executable,
                str(Path(converter.__file__)),
                "-i",
                str(input_path),
                "-o",
                str(output_path),
                *extra_arguments,
            ]
            result = subprocess.run(command, text=True, capture_output=True, check=False)
            output = output_path.read_text(encoding="utf-8") if output_path.exists() else ""
            return result, output

    def test_short_arguments(self):
        result, output = self.run_converter(
            ["-w", "1.25", "-p", "4 3 2 1", "-m", "split-fission", "--MLE"]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('weight="1.25"', output)
        self.assertIn('weights="4.0 3.0 2.0 1.0"', output)
        self.assertIn('value="3"', output)
        self.assertIn("marginalLikelihoodEstimator", output)

    def test_long_argument_aliases(self):
        args = converter.build_argument_parser().parse_args(
            [
                "--input",
                "input.xml",
                "--output",
                "output.xml",
                "--operator-weight",
                "1.5",
                "--prior-weights",
                "1 2 3 4",
                "--initial-model",
                "budding",
            ]
        )
        self.assertEqual(args.input, Path("input.xml"))
        self.assertEqual(args.output, Path("output.xml"))
        self.assertEqual(args.operator_weight, 1.5)
        self.assertEqual(args.prior_weights, (1.0, 2.0, 3.0, 4.0))
        self.assertEqual(args.initial_model, 1)

    def test_long_arguments_and_mle_forms(self):
        forms = [
            ([], False),
            (["--noMLE"], False),
            (["--MLE", "0"], False),
            (["--MLE", "1"], True),
            (["--MLE"], True),
        ]
        for arguments, expected_mle in forms:
            with self.subTest(arguments=arguments):
                result, output = self.run_converter(arguments)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual("marginalLikelihoodEstimator" in output, expected_mle)
                self.assertIn("harmonicMeanAnalysis", output)

    def test_invalid_mle_value_and_mutually_exclusive_flags(self):
        for arguments in (["--MLE", "2"], ["--MLE", "--noMLE"]):
            with self.subTest(arguments=arguments):
                result, _ = self.run_converter(arguments)
                self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
