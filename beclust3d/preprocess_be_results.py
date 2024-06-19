"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.2

"""

import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from DSSPparser import parser
from DSSPparser import parseDSSP
import time
from pathlib import Path
import numpy as np
import statistics as st
import re
import json
import os.path
from os import path
import subprocess
import shutil
from biopandas.pdb import PandasPdb
import csv
import math;
import seaborn as sns



