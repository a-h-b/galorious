import os
import shutil
import gzip
import yaml
from re import match
from copy import deepcopy
import subprocess

from snakemake.utils import validate

def dict_merge(a, b):
    """
    Deep merge 2 dicts together
    """
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.items():
        if k in result and isinstance(result[k], dict):
            result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result

# default executable for snakemake
shell.executable("bash")

# custom configuration file
CUSTOM_CONFIG_PATH = os.environ.get("CONFIGFILE","../nothing/here")

# get parameters from the command line
OUTPUTDIR = os.path.expandvars(config['outputdir'])

EMAIL = config['email']
if EMAIL != "":
    if not re.fullmatch(r"^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$", EMAIL):
        EMAIL = ""
        print("Your email address is not valid, you will not receive notifications.")

validate(config, schema="../../schemas/schema.yaml")

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))

STEPS = sorted(os.environ.get("STEPS", config['steps']).split())

FILTER = os.environ.get("FILTER", config['illumina_filter']).split()

# temporary directory will be stored inside the OUTPUTDIR directory
# unless an absolute path is set
TMPDIR = os.path.expandvars(config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.abspath(os.path.join(OUTPUTDIR, TMPDIR))
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

