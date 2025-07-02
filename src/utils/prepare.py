import pandas as pd
import numpy as np
from skbio.stats.composition import multiplicative_replacement

def prepare_for_clr(df,
                    min_sample_sum=1e-5,
                    min_taxa_presence=0.05,
                    min_median_abundance=1e-6,
                    max_single_taxon_frac=0.95,
                    normalize=True,
                    verbose=True):
    from skbio.stats.composition import multiplicative_replacement
    
    original_shape = df.shape
    df = df.copy()

    if normalize:
        df = df.div(df.sum(axis=1), axis=0)

    # Удаляем сэмплы с почти нулевой суммой
    df = df[df.sum(axis=1) > min_sample_sum]

    # Удаляем таксоны, которые встречаются слишком редко
    min_presence = int(df.shape[0] * min_taxa_presence)
    df = df.loc[:, (df > 0).sum(axis=0) >= min_presence]

    # Удаляем таксоны с низким медианным abundance
    df = df.loc[:, df.median(axis=0) > min_median_abundance]

    # Удаляем сэмплы, где один таксон доминирует
    max_frac = df.max(axis=1)
    df = df[max_frac < max_single_taxon_frac]

    # ❗ Удаляем строки, где все значения нули (иначе multip_replacement упадёт)
    df = df[(df > 0).any(axis=1)]

    if df.shape[0] == 0:
        raise ValueError("После фильтрации не осталось ни одного валидного сэмпла.")

    # Преобразование — замена нулей
    df_replaced = multiplicative_replacement(df.values)
    df_replaced = pd.DataFrame(df_replaced, index=df.index, columns=df.columns)

    if verbose:
        print(f"[INFO] Отфильтровано: {original_shape} → {df_replaced.shape}")
    
    return df_replaced
