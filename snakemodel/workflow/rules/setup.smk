import os
import re
import tempfile
import functools as ft
import shlex
import itertools as it

from bids.layout import parse_file_entities
import more_itertools as itx
from snakebids import bids, generate_inputs, filter_list
from snakebids.exceptions import ConfigError
import pandas as pd
import numpy as np

from pathlib import Path
from snakeboost import Tar, Pyscript, XvfbRun, PipEnv, Boost, Datalad, Env
import snakeboost.bash as sh
from templateflow import api as tflow
from lib.participants import filter_participants

from patsy import dmatrix
import sklearn.impute as impute
from sklearn import impute, preprocessing as preproc, pipeline, compose
# def get_labels(label):
#     cli = config.get(label, None)
#     vals = cli if isinstance(cli, int) else cli.split(",") if cli is not None else None
#     return vals

# def get_participants(participant_file):
#     subs = pd.read_csv(participant_file, sep='\t')
#     return subs['participant_id'].map(lambda s: s[4:])

# participant_label = get_labels("participant_label")
# exclude_participant_label = (
#     get_labels("exclude_participant_label") if participant_label is None else None
# )
# if participant_label is None and exclude_participant_label is None:
#     try:
#         participant_label = list(get_participants(
#             Path(config['bids_dir'], 'derivatives', 'snakedwi-0.1.0', 'participants.tsv')
#         ))
#     except FileNotFoundError:
#         pass



tmpdir = eval(workflow.default_resources._args.get("tmpdir"), {"system_tmpdir": tempfile.gettempdir()})

if workflow.run_local:
    workflow.shadow_prefix = os.environ.get("SLURM_TMPDIR")

###
# Input Globals
###
inputs = generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs_run'],
    derivatives=config.get("derivatives", False),
    participant_label=config.get('participant_label'),
    exclude_participant_label=config.get('exclude_participant_label'),
    pybidsdb_dir=config.get("pybidsdb_dir"),
)

outdir = config.get("outdir")
settings = dict(setting.split('=') for setting in config.get("model_settings", []) or [])
if config["model"] == "tbss":
    settings["parameter"] = settings.get("parameter", "FA")
fa = config["component"]
if "param_map" in inputs:
    comp = "param_map"
else:
    comp = fa



###
# Output Globals
###

work = Path(tmpdir) / "snakemodel"
source = Path(config['output_dir']) / "sourcedata"
output = Path(config['output_dir'])
qc = Path(output)/"qc"

# Unique ID for easy naming in temporary files

def _uid(comp = None, entities = None):
    parts = comp.wildcards.values() if comp is not None else []
    wrapped_entities = (f"{{{entity}}}" for entity in entities or [])
    return ".".join(it.chain(parts, wrapped_entities))

def tempout(rulename: str, comp = None, extension = "", *entities):
    uid = _uid(comp, entities)
    return temp(work.joinpath(uid, rulename).with_suffix(extension))

def log(rulename: str, comp = None, *entities):
    uid = _uid(comp, entities)
    return Path("code", "logs", uid, rulename).with_suffix(".log")

def benchmark(rulename: str, comp = None, *entities):
    uid = _uid(comp, entities)
    return Path("code", "benchmarks", uid,  rulename).with_suffix(".tsv")

def resource(path):
    return os.path.join(workflow.basedir, "..", "resources", path)


# ### Model variables
# class Contrast:
#     def __init__(self, spec: str):
#         lt = "<" in spec
#         gt = ">" in spec
#         if lt and gt:
#             raise Exception(f"Only one of '<' and '>' may be in a spec (got: {spec})")
#         if lt:
#             self.lower = 
        

# if not config['unpaired_ttest']:
#     raise Exception("--unpaired-ttest must be supplied a value")
# group_col = config['unpaired_ttest']
# groups = list({

# })

###
# Utility functions
###
boost = Boost(work, logger)
pyscript = Pyscript(workflow.basedir)


###
# Pipenvs
###

design_matrix_env = PipEnv(
    packages = [
        "more-itertools",
        "numpy",
        "pandas",
        "scikit-learn",
        "scipy",
        "/scratch/knavynde/snakeboost",
        "h5py",
    ],
    flags = config.get("pip-flags", ""),
    root = work,
)