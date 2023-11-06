import yaml
from os.path import abspath


def readDefaults(yaml_path: str) -> dict:
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    default_values = {key: value["default"] for key, value in config.get("properties", {}).items() if "default" in value}
    return default_values


def makeConfig(config: dict, defaults: dict) -> dict:
    """
    Remove Nones from config and add missing variables from the defaults
    """
    config = {key: value for key, value in config.items() if value is not None and value != "None"}

    # Merge 'config' with 'defaults' to add missing variables
    config = {**defaults, **config}
    if not config.get("mle_iterations"):
        config["mle_iterations"] = int(config.get("iterations") / config.get("mle_steps"))
    if not config.get("mle_sampling"):
        config["mle_sampling"] = int(config.get("iterations") / 1000)
    if config.get("blacklist"):
        config["blacklist"] = abspath(config.get("blacklist"))
    return config
