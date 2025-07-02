import pandas as pd
from scipy.spatial.distance import braycurtis, pdist, squareform, jensenshannon

preprocessed_data = dict()
dataset_name = 'abundance'

PATH = "/home/jupyter-denisova2/microbiome/Human-Gut-Microbiome-Analysis"

metric_abbreviation = "BC" # 'jensenshannon', 'braycurtis'

abbreviation_to_name = {
    "BC": "braycurtis",
    "JS": "jensenshannon",
}

for tax in ['order', 'family', 'genus', 'species']:
    df = pd.read_csv(
            f'{PATH}/data_processed_pat0.7/{dataset_name}_{tax}.csv', 
            sep=',',
        )        

    label = f'{dataset_name}_{tax}'

    data = df.copy()
    index = data["sampleid"]

    data = data.drop(columns = ["sampleid"])
    data_norm = data.div(data.sum(axis=1) + 1e-12, axis=0)

    pair_table = pd.DataFrame(
        squareform(
            pdist(
                data.values, 
                metric=abbreviation_to_name[metric_abbreviation]
            )
        ),
        index=data_norm.index,
        columns=data_norm.index
    )

    pair_table.index = index
    pair_table.to_csv(f"pair_distances/{metric_abbreviation}_{tax}.csv", index = True)