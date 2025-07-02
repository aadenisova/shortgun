%load_ext autoreload
import os
#import numpy as np 
import glob
from skbio.stats.composition import clr
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from utils import project_pca, calculate_Q_metrics, classify_enterotypes

from warnings import simplefilter
simplefilter('ignore')


data_orig = {}
y_orig = {}

blacklist = []

for tax in ['order', 'family', 'genus', 'species']:
    df = pd.read_csv(
            f'../whole_data/vdWGS_microbiom/{tax}_abundance_table_vdWGS.unfiltered.tsv', 
            sep='\t',
            index_col = 'sampleid'
        )        
    if len(blacklist) == 0:
        blacklist = df[df.sum(axis=1) == 0].index
        
    df = df[~df.index.isin(blacklist)]
    label = f'{dataset_name}_{tax}'
    df = df.fillna(0)

    df_normalized = df.div(df.sum(axis=1), axis=0) 
    if dataset_name == "abundance":
        data_orig[label] = df_normalized[df_normalized.sum(axis=1) != 0]
    else:
        data_prop = df_normalized
        pseudo_count = 0.5 * data_prop[data_prop > 0].min().min()  # половина минимального ненулевого значения
        data_prop = data_prop + pseudo_count
        
        # CLR-преобразование
        clr_data = pd.DataFrame(
            clr(data_prop),
            index=data_prop.index,
            columns=data_prop.columns
        )
        data_orig[label] = clr_data
            
dataframe = pd.read_csv(f'../data/RA/pivot_g_our.csv', sep=',')
label = 'g_our'
data_orig[label] = dataframe.drop('Unnamed: 0', axis=1)
processed_root = 'data_processed'
pca_root = './results/pca'


