import multiprocessing as mp
from Bio import Entrez
import argparse as ap
import sys
import time
from pathlib import Path
import pickle
import pmidprocessor as pp


pmid = sys.argv[-1]
pp.write_pickle(pmid)


