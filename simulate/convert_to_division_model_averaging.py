#!/usr/bin/env python3
"""Convert a fixed PHYFUM division-model XML to model averaging."""

from __future__ import annotations

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Match, Sequence


MODEL_VALUES = {
    "identity": 0,
    "budding": 1,
    "fission": 2,
    "split-fission": 3,
}
OUTPUT_MARKER = "divisionModelAveraged"


class ConversionError(ValueError):
    """Raised when an input XML cannot be converted unambiguously."""


def parse_initial_model(value: str) -> int:
    normalized = value.strip().lower().replace("_", "-")
    if normalized in MODEL_VALUES:
        return MODEL_VALUES[normalized]
    try:
        numeric_value = int(normalized)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "initial model must be identity, budding, fission, split-fission, or 0..3"
        ) from error
    if numeric_value not in range(4):
        raise argparse.ArgumentTypeError("initial model must be between 0 and 3")
    return numeric_value


def parse_operator_weight(value: str) -> float:
    try:
        weight = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("operator weight must be a number") from error
    if not math.isfinite(weight) or weight <= 0:
        raise argparse.ArgumentTypeError("operator weight must be finite and positive")
    return weight


def parse_prior_weights(value: str) -> tuple[float, float, float, float]:
    fields = value.replace(",", " ").split()
    if len(fields) != 4:
        raise argparse.ArgumentTypeError("prior weights must contain exactly four values")
    try:
        weights = tuple(float(field) for field in fields)
    except ValueError as error:
        raise argparse.ArgumentTypeError("prior weights must be numbers") from error
    if any(not math.isfinite(weight) or weight < 0 for weight in weights):
        raise argparse.ArgumentTypeError("prior weights must be finite and nonnegative")
    if not any(weight > 0 for weight in weights):
        raise argparse.ArgumentTypeError("at least one prior weight must be positive")
    return weights  # type: ignore[return-value]


def format_number(value: float) -> str:
    return str(value)


def child_indent(body: str, fallback: str) -> str:
    match = re.search(r"(?m)^([ \t]*)<", body)
    return match.group(1) if match else fallback + "\t"


def indentation_unit(parent_indent: str, nested_indent: str) -> str:
    if nested_indent.startswith(parent_indent) and len(nested_indent) > len(parent_indent):
        return nested_indent[len(parent_indent) :]
    return "\t"


def insert_before_closing(body: str, addition: str) -> str:
    trailing_whitespace = re.search(r"\s*$", body)
    if trailing_whitespace is None:
        return body + addition
    position = trailing_whitespace.start()
    prefix = body[:position]
    separator = "" if prefix.endswith(("\n", "\r")) else "\n"
    trailing = body[position:]
    if addition.endswith(("\n", "\r")):
        trailing = re.sub(r"^\r?\n", "", trailing, count=1)
    return prefix + separator + addition + trailing


def add_division_model_to_likelihood(xml_text: str, initial_model: int) -> str:
    pattern = re.compile(
        r"(?P<line_indent>^[ \t]*)<cenancestorTreeLikelihood\b(?P<attrs>[^>]*)>"
        r"(?P<body>.*?)"
        r"(?P=line_indent)</cenancestorTreeLikelihood>",
        re.MULTILINE | re.DOTALL,
    )
    matches = list(pattern.finditer(xml_text))
    if len(matches) != 1:
        raise ConversionError(
            f"expected exactly one fixed cenancestorTreeLikelihood definition; found {len(matches)}"
        )

    match = matches[0]
    attrs = re.sub(
        r"\s+divisionModel\s*=\s*(['\"]).*?\1",
        "",
        match.group("attrs"),
        flags=re.DOTALL,
    )
    body = match.group("body")
    line_indent = match.group("line_indent")
    indent = child_indent(body, line_indent)
    inner_indent = indent + indentation_unit(line_indent, indent)
    division_xml = (
        f"{indent}<divisionModel>\n"
        f'{inner_indent}<parameter id="divisionModel" value="{initial_model}" lower="0" upper="3"/>\n'
        f"{indent}</divisionModel>\n"
    )

    branch_rate_line = re.search(
        r"(?m)^[ \t]*<[^>]*CenancestorBranchRates\b[^>]*/>[ \t]*\n?", body
    )
    if branch_rate_line:
        body = body[: branch_rate_line.start()] + division_xml + body[branch_rate_line.start() :]
    else:
        body = insert_before_closing(body, division_xml)

    replacement = (
        f'{line_indent}<modelAveragingCenancestorTreeLikelihood{attrs}>'
        f"{body}{line_indent}</modelAveragingCenancestorTreeLikelihood>"
    )
    converted = xml_text[: match.start()] + replacement + xml_text[match.end() :]
    return re.sub(
        r"<cenancestorTreeLikelihood\b(?P<attrs>[^>]*)/>",
        r"<modelAveragingCenancestorTreeLikelihood\g<attrs>/>",
        converted,
    )


def add_operator(xml_text: str, operator_weight: float) -> str:
    pattern = re.compile(
        r"(?P<line_indent>^[ \t]*)<operators\b(?P<attrs>[^>]*)>"
        r"(?P<body>.*?)"
        r"(?P=line_indent)</operators>",
        re.MULTILINE | re.DOTALL,
    )
    matches = list(pattern.finditer(xml_text))
    if len(matches) != 1:
        raise ConversionError(f"expected exactly one operators definition; found {len(matches)}")
    match = matches[0]
    indent = child_indent(match.group("body"), match.group("line_indent"))
    inner_indent = indent + indentation_unit(match.group("line_indent"), indent)
    addition = (
        f"{indent}<uniformIntegerOperator lower=\"0\" upper=\"3\" "
        f"weight=\"{format_number(operator_weight)}\">\n"
        f'{inner_indent}<parameter idref="divisionModel"/>\n'
        f"{indent}</uniformIntegerOperator>\n"
    )
    body = insert_before_closing(match.group("body"), addition)
    replacement = (
        f'{match.group("line_indent")}<operators{match.group("attrs")}>'
        f'{body}{match.group("line_indent")}</operators>'
    )
    return xml_text[: match.start()] + replacement + xml_text[match.end() :]


def add_prior(xml_text: str, prior_weights: Sequence[float]) -> str:
    pattern = re.compile(
        r"(?P<line_indent>^[ \t]*)<prior\b(?P<attrs>[^>]*)>"
        r"(?P<body>.*?)"
        r"(?P=line_indent)</prior>",
        re.MULTILINE | re.DOTALL,
    )
    matches = [match for match in pattern.finditer(xml_text) if re.search(r'\bid\s*=\s*(["\'])prior\1', match.group("attrs"))]
    if len(matches) != 1:
        raise ConversionError(f"expected exactly one prior definition with id=\"prior\"; found {len(matches)}")
    match = matches[0]
    indent = child_indent(match.group("body"), match.group("line_indent"))
    inner_indent = indent + indentation_unit(match.group("line_indent"), indent)
    weights = " ".join(format_number(weight) for weight in prior_weights)
    addition = (
        f'{indent}<divisionModelPrior weights="{weights}">\n'
        f'{inner_indent}<parameter idref="divisionModel"/>\n'
        f"{indent}</divisionModelPrior>\n"
    )
    body = insert_before_closing(match.group("body"), addition)
    replacement = (
        f'{match.group("line_indent")}<prior{match.group("attrs")}>'
        f'{body}{match.group("line_indent")}</prior>'
    )
    return xml_text[: match.start()] + replacement + xml_text[match.end() :]


def add_division_model_logs(xml_text: str) -> str:
    pattern = re.compile(
        r"(?P<line_indent>^[ \t]*)<log\b(?P<attrs>[^>]*)>"
        r"(?P<body>.*?)"
        r"</log>",
        re.MULTILINE | re.DOTALL,
    )

    def add_parameter(match: Match[str]) -> str:
        attrs = match.group("attrs")
        id_match = re.search(r'\bid\s*=\s*(["\'])(?P<id>[^"\']+)\1', attrs)
        if not id_match or id_match.group("id").lower() not in {"screenlog", "filelog"}:
            return match.group(0)
        body = match.group("body")
        if re.search(r'<parameter\b[^>]*\bidref\s*=\s*(["\'])divisionModel\1', body):
            return match.group(0)
        indent = child_indent(body, match.group("line_indent"))
        body = insert_before_closing(body, f'{indent}<parameter idref="divisionModel"/>\n')
        return f'{match.group("line_indent")}<log{attrs}>{body}</log>'

    return pattern.sub(add_parameter, xml_text)


def remove_mle_block(xml_text: str) -> str:
    block_pattern = re.compile(
        r"(?P<indent>^[ \t]*)<marginalLikelihoodEstimator\b[^>]*>.*?"
        r"(?P=indent)</marginalLikelihoodEstimator>[ \t]*(?:\r?\n)?",
        re.MULTILINE | re.DOTALL,
    )
    matches = list(block_pattern.finditer(xml_text))
    if len(matches) > 1:
        raise ConversionError(
            f"expected at most one marginalLikelihoodEstimator block; found {len(matches)}"
        )
    if not matches:
        return xml_text

    match = matches[0]
    start = match.start()
    preceding = xml_text[:start]
    comment_match = re.search(
        r"(?P<space>^[ \t]*\r?\n)?^[ \t]*<!--\s*"
        r"Power-posterior sampling for MLE estimation\s*-->[ \t]*\r?\n$",
        preceding,
        re.MULTILINE | re.IGNORECASE,
    )
    if comment_match:
        start = comment_match.start()
    return xml_text[:start] + xml_text[match.end() :]


def rename_output_target(target: str) -> str:
    slash_position = max(target.rfind("/"), target.rfind("\\"))
    directory = target[: slash_position + 1]
    filename = target[slash_position + 1 :]
    if filename == OUTPUT_MARKER or filename.endswith(f".{OUTPUT_MARKER}"):
        return target
    marker_with_dots = f".{OUTPUT_MARKER}."
    if marker_with_dots in filename:
        return target
    dot_position = filename.rfind(".")
    if dot_position > 0 and dot_position < len(filename) - 1:
        filename = f"{filename[:dot_position]}.{OUTPUT_MARKER}{filename[dot_position:]}"
    else:
        filename = f"{filename}.{OUTPUT_MARKER}"
    return directory + filename


def rename_output_attributes(xml_text: str) -> str:
    pattern = re.compile(
        r"(?P<name>\b(?:fileName|operatorAnalysis))(?P<spacing>\s*=\s*)"
        r"(?P<quote>[\"'])(?P<value>.*?)(?P=quote)"
    )

    def rename(match: Match[str]) -> str:
        return (
            f'{match.group("name")}{match.group("spacing")}{match.group("quote")}'
            f'{rename_output_target(match.group("value"))}{match.group("quote")}'
        )

    return pattern.sub(rename, xml_text)


def validate_xml(xml_text: str, description: str) -> None:
    try:
        ET.fromstring(xml_text)
    except ET.ParseError as error:
        raise ConversionError(f"{description} is not well-formed XML: {error}") from error


def convert_xml(
    xml_text: str,
    operator_weight: float = 0.5,
    prior_weights: Sequence[float] = (1.0, 1.0, 1.0, 1.0),
    initial_model: int = 0,
    preserve_mle: bool = False,
) -> str:
    validate_xml(xml_text, "input")
    if "<modelAveragingCenancestorTreeLikelihood" in xml_text:
        raise ConversionError("input already uses model-averaged division likelihood")

    converted = xml_text
    if not preserve_mle:
        converted = remove_mle_block(converted)
    converted = add_division_model_to_likelihood(converted, initial_model)
    converted = add_operator(converted, operator_weight)
    converted = add_prior(converted, prior_weights)
    converted = add_division_model_logs(converted)
    converted = rename_output_attributes(converted)
    validate_xml(converted, "converted output")
    return converted


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, type=Path, help="fixed-model XML")
    parser.add_argument("-o", "--output", required=True, type=Path, help="converted XML")
    parser.add_argument(
        "-w", "--operator-weight", type=parse_operator_weight, default=0.5
    )
    parser.add_argument(
        "-p",
        "--prior-weights",
        type=parse_prior_weights,
        default=(1.0, 1.0, 1.0, 1.0),
        metavar='"W0 W1 W2 W3"',
    )
    parser.add_argument("-m", "--initial-model", type=parse_initial_model, default=0)
    mle_group = parser.add_mutually_exclusive_group()
    mle_group.add_argument(
        "--MLE",
        nargs="?",
        const=1,
        default=0,
        type=int,
        choices=(0, 1),
        metavar="{0,1}",
        help="preserve the power-posterior MLE block (bare --MLE means 1)",
    )
    mle_group.add_argument(
        "--noMLE",
        action="store_true",
        help="remove the power-posterior MLE block (default)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="allow input replacement or output overwrite"
    )
    return parser


def run(arguments: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(arguments)
    input_path = args.input.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    if not input_path.is_file():
        parser.error(f"input is not a readable file: {input_path}")
    if input_path == output_path and not args.overwrite:
        parser.error("input and output are the same; use --overwrite to replace the input")
    if output_path.exists() and not args.overwrite:
        parser.error(f"output already exists: {output_path}; use --overwrite to replace it")
    if not output_path.parent.is_dir():
        parser.error(f"output directory does not exist: {output_path.parent}")

    try:
        source = input_path.read_text(encoding="utf-8")
        converted = convert_xml(
            source,
            operator_weight=args.operator_weight,
            prior_weights=args.prior_weights,
            initial_model=args.initial_model,
            preserve_mle=bool(args.MLE) and not args.noMLE,
        )
        output_path.write_text(converted, encoding="utf-8")
    except (OSError, UnicodeError, ConversionError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(run())
