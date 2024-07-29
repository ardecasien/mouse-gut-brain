# make pathway file

nano mimenet_pathways.sh

#!/bin/bash

module load cuDNN/8.2.1/CUDA-11.3 python/3.9

python MiMeNet_train.py -micro pathway-data-matched-to-expression.csv -metab exp-data-matched-to-pathways.csv -micro_norm None -metab_norm None -net_params network_parameters.txt -num_background 10 -num_run_cv 10 -output mimeout-exp-pathways-matched 

echo done

# submit job

sbatch --partition=gpu --gres=gpu:k80:1,lscratch:10 --mem=50g -c14 --time=12:00:00 mimenet_pathways.sh

# dependencies 

import pickle
import biom
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import scipy.cluster.hierarchy as shc
# from skbio.stats.composition import clr 
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import KFold
from src.models.MiMeNet import MiMeNet, tune_MiMeNet 
from src.models.MLPNN import MLPNN 
from scipy.cluster.hierarchy import cut_tree
from scipy.stats import mannwhitneyu
import tensorflow as tf



