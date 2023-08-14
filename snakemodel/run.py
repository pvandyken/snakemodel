#!/usr/bin/env python3
import os
from snakebids.app import SnakeBidsApp
from snakebids.exceptions import ConfigError


def print_design(app: SnakeBidsApp):
    if app.config["print_design"]:
        app.config["targets_by_analysis_level"][app.config["analysis_level"]] = [
            "test_patsy_dm"
        ]

def pick_component(app: SnakeBidsApp):
    selected = app.config["component"]
    if selected not in app.config["pybids_inputs"]:
        choices = "\n    ".join(app.config["pybids_inputs"])
        raise ConfigError(f"--component must be set to one of\n    {choices}")

    param_map = app.config["pybids_inputs"]["param_map"]
    app.config["pybids_inputs"] = {
        selected: app.config["pybids_inputs"][selected]
    }
    # if "filter_param_map" in app.args.args_dict:
        # app.config["pybids_inputs"]["param_map"] = param_map

def main():
    app = SnakeBidsApp(
        os.path.abspath(os.path.dirname(__file__)),
        plugins=[
            pick_component,
            print_design,
        ],
    )
    app.run_snakemake()


if __name__ == "__main__":
    main()