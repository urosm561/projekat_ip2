import pandas as pd
import numpy as np

from typing import List, Optional, Tuple
from dataclasses import dataclass

from sklearn.base import BaseEstimator, TransformerMixin

COLUMNS_TO_DROP = ["left_end", "right_end", "left_start", "right_start", "sequence_id", "length", "protein_name", "repeat_type"]

def preprocess_data(
        data: str | pd.DataFrame, 
        seed: int, 
        train_size: Optional[float] = None, 
        valid_size: Optional[float] = None, 
        test_size: Optional[float] = None,
        columns_to_drop: List[str] = COLUMNS_TO_DROP 
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Prepares the data by dropping `columns_to_drop`, and splitting
    it into train, validation and test sets according to specified sizes.

    Returns these datasets.
    """

    # Various checks for split sizes.
    eps = 1e-5

    if not train_size and not valid_size and not test_size:
        train_size = 0.7
        valid_size = 0.2
        test_size = 0.1

    elif not train_size and not valid_size:
        train_size = 7/9 * (1 - test_size)
        valid_size = 2/9 * (1 - test_size)

    elif not train_size and not test_size:
        train_size = 7/8 * (1 - valid_size)
        test_size = 1/8 * (1 - valid_size)

    elif not valid_size and not test_size:
        valid_size = 2/3 * (1 - train_size)
        test_size = 1/3 * (1 - train_size)

    elif not train_size:
        train_size = 1 - valid_size - test_size
        train_size = train_size if abs(train_size) > eps else 0

    elif not valid_size:
        valid_size = 1 - train_size - test_size
        valid_size = valid_size if abs(valid_size) > eps else 0

    elif not test_size:
        test_size = 1 - train_size - valid_size
        test_size = test_size if abs(test_size) > eps else 0

    else:
        if abs(train_size + valid_size + test_size - 1) > eps:
            print(train_size + valid_size + test_size - 1)
            raise ValueError("Split sizes must add up to 1.")

    # Reading data.
    if isinstance(data, str):
        data = pd.read_csv(data)
    elif not isinstance(data, pd.DataFrame):
        raise ValueError("data type is wrong!")
    
    # Dropping columns.
    data.drop(columns = columns_to_drop, inplace = True)

    # Cleaning up the repeats, we don't need both the left and right ones.
    if "left_repeat" in data.columns and "right_repeat" in data.columns:
        data.drop(columns = ["right_repeat"], inplace = True)
    
    # Make sure not to drop all the repeats and virus types as that is necessary for the models.
    if (("left_repeat" not in data.columns and "right_repeat" not in data.columns)
        or "virus_name" not in data.columns):
        print(data)
        raise ValueError("Removed important data from the dataframe.")

    # Since we now have only one repeat in the data we rename it.
    repeats_column_name = "left_repeat" if "left_repeat" in data.columns else "right_repeat"
    data["repeat"] = data[repeats_column_name]
    data.drop(columns = [repeats_column_name], inplace = True)

    # Reorder columns so that the repeat is at the beginning and virus type is at the end.
    middle_cols = [col for col in data.columns if col not in ["repeat", "virus_name"]]
    new_order = ["repeat"] + middle_cols + ["virus_name"]
    data = data[new_order]

    # Split the data into train, validation and test datasets.
    train_data = data.sample(frac = train_size, random_state = seed)
    valid_test_data = data.drop(train_data.index)

    valid_data = valid_test_data.sample(frac = valid_size / (valid_size + test_size), random_state = seed)
    test_data = valid_test_data.drop(valid_data.index)

    return train_data, valid_data, test_data

@dataclass
class Metrics:
    accuracy: float
    precision: float
    recall: float
    f1: float
    confusion_matrix: np.ndarray

AMINO_ACIDS = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
PAD = "_"
CATS_21 = AMINO_ACIDS + [PAD]

class PositionalEncoder(BaseEstimator, TransformerMixin):
    def __init__(self, col_name: str = "repeat", width: int = 7, pad_char: str = PAD):
        self.col_name = col_name
        self.width = width
        self.pad_char = pad_char

    def fit(self, X, y = None):
        return self
    
    def transform(self, X):
        s = X[self.col_name].astype(str).values
        out = []

        for seq in s:
            seq = seq[-self.width : ]
            seq = (self.pad_char * (self.width - len(seq))) + seq
            out.append(list(seq))

        cols = [f"p{i + 1}" for i in range(self.width)]
        return pd.DataFrame(out, columns = cols)

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = None
    valid_size = None
    test_size = None
    columns_to_drop = ["left_end", "right_end", "left_start", "right_start", "sequence_id", "length", "protein_name", "repeat_type"]
    train, valid, test = preprocess_data(data_path, seed, train_size, valid_size, test_size)
    print(train)
    print(valid)
    print(test)