###########
## sync files
###########

rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/MiMeNet_train.py decasienar@helix.nih.gov:/data/DNU/alex/microbiome
rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/fwmousedata/exp-data-matched-to-pathways.csv decasienar@helix.nih.gov:/data/DNU/alex/microbiome
rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/fwmousedata/pathway-data-matched-to-expression.csv decasienar@helix.nih.gov:/data/DNU/alex/microbiome
rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/fwmousedata/exp-data-matched-to-families.csv decasienar@helix.nih.gov:/data/DNU/alex/microbiome
rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/fwmousedata/family-data-matched-to-expression.csv decasienar@helix.nih.gov:/data/DNU/alex/microbiome
rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/network_parameters.txt decasienar@helix.nih.gov:/data/DNU/alex/microbiome

rsync -avzP  /Users/decasienar/Desktop/Papers/Microbiome/MiMeNet-master/ decasienar@helix.nih.gov:~/.local/lib/python3.9/site-packages

sinteractive --gres=gpu:p100:1 

cd /data/DNU/alex/microbiome

# pip install scikit-bio
# mkdir ~/.local/lib/python3.9
# mkdir ~/.local/lib/python3.9/site-packages

# make pathway file

nano mimenet_pathways.sh

#!/bin/bash

module load  cuDNN/8.2.1/CUDA-11.3 python/3.9

python MiMeNet_train.py -micro pathway-data-matched-to-expression.csv -metab exp-data-matched-to-pathways.csv -micro_norm None -metab_norm None -net_params network_parameters.txt -num_background 10 -num_run_cv 10 -output mimeout-exp-pathways-matched 

echo done

# make family file

nano mimenet_families.sh

#!/bin/bash

module load  cuDNN/8.2.1/CUDA-11.3 python/3.9

python MiMeNet_train.py -micro family-data-matched-to-expression.csv -metab exp-data-matched-to-families.csv -micro_norm None -metab_norm None -net_params network_parameters.txt -num_background 10 -num_run_cv 10 -output mimeout-exp-families-matched 

echo done

# submit jobs

sbatch --partition=gpu --gres=gpu:k80:1,lscratch:10 --mem=50g -c14 --time=12:00:00 mimenet_pathways.sh
sbatch --partition=gpu --gres=gpu:k80:1,lscratch:10 --mem=100g -c5 --time=12:00:00 mimenet_families.sh

# export

rsync -avzP decasienar@helix.nih.gov:"/data/DNU/alex/microbiome/results/mimeout-exp-pathways-matched/Images/*.png" /Users/decasienar/Desktop/Papers/Microbiome/results-exp-pathways-matched
rsync -avzP decasienar@helix.nih.gov:"/data/DNU/alex/microbiome/results/mimeout-exp-pathways-matched/*.csv" /Users/decasienar/Desktop/Papers/Microbiome/results-exp-pathways-matched
rsync -avzP decasienar@helix.nih.gov:"/data/DNU/alex/microbiome/results/mimeout-exp-pathways-matched/CV/interaction_score_matrix.csv" /Users/decasienar/Desktop/Papers/Microbiome/results-exp-pathways-matched

rsync -avzP decasienar@helix.nih.gov:"/data/DNU/alex/microbiome/results/mimeout-exp-families-matched/Images/*.png" /Users/decasienar/Desktop/Papers/Microbiome/results-exp-families-matched
rsync -avzP decasienar@helix.nih.gov:"/data/DNU/alex/microbiome/results/mimeout-exp-families-matched/*.csv" /Users/decasienar/Desktop/Papers/Microbiome/results-exp-families-matched
rsync -avzP decasienar@helix.nih.gov:"/data/DNU/alex/microbiome/results/mimeout-exp-families-matched/CV/interaction_score_matrix.csv" /Users/decasienar/Desktop/Papers/Microbiome/results-exp-families-matched

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



