from pathlib import Path
import yaml
import pandas as pd
from skbio.stats.composition import multiplicative_replacement, clr
from scipy.spatial.distance import pdist, squareform

def load_config(config_path):
    """Загружает конфигурацию из YAML-файла."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config