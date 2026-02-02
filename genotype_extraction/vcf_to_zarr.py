#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from tqdm import tqdm
import time
import os.path
import pandas as pd
import numpy as np
import scipy
import h5py
import allel
from allel import chunked
import pandas as pd
import zarr
import itertools
import argparse
parser = argparse.ArgumentParser(description="Example script")
parser.add_argument("--vcf", required=True, help="Input file vcf ")
parser.add_argument("--out_zarr", help="zarr output file ")

args = parser.parse_args()

print('scikit-allel', allel.__version__)

vcf = args.vcf
out_zarr = args.out_zarr
allel.vcf_to_zarr(vcf, out_zarr)
