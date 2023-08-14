import pandas as pd
import more_itertools as itx


def filter_participants(participant_path, **filters):
    df = pd.read_csv(participant_path, sep="\t")
    for col, values in filters.items():
        if col not in df:
            continue
        df = df[df[col].isin(map(df[col].dtype.type, itx.always_iterable(values)))]
    return df
