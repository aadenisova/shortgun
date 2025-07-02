from skbio.stats.composition import clr
from scipy.spatial.distance import pdist, squareform

import pandas as pd

from pathlib import Path

from utils.prepare import prepare_for_clr
from utils.config import load_config

config = load_config('config.yaml')

data_dir = Path(config['data_dir'])
output_dir = Path(config['output_dir'])
output_dir.mkdir(parents=True, exist_ok=True)

tax_levels = config.get('tax_levels', ['order', 'family', 'genus', 'species'])

for tax in tax_levels:
    file_path = data_dir / f"{tax}_abundance_table_vdWGS.unfiltered.tsv"
    if not file_path.exists():
        print(f"[WARNING] Файл не найден: {file_path}")
        continue

    df = pd.read_csv(file_path, sep='\t', index_col='sampleid')
    
    df_clr_ready = prepare_for_clr(df)
    
    clr_data = clr(df_clr_ready.values)
    
    aitchison = pd.DataFrame(
        squareform(pdist(clr_data, metric='euclidean')),
        index=df_clr_ready.index,
        columns=df_clr_ready.index
    )
    
    # Сохранение результатов
    output_file = output_dir / f"{tax}_aitchison_distance.csv"
    aitchison.to_csv(output_file)
    print(f"[INFO] Сохранено: {output_file}")