import pytest
from pathlib import Path


@pytest.fixture()
def example_config():
    return Path("tests/data/config.yaml")


@pytest.fixture()
def defaults():
    return Path("src/phyfum/workflow/config.schema.yml")
