import pandas as pd
import numpy as np

from typing import List, Optional, Tuple
from dataclasses import dataclass

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix

from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler

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
    train_valid_data, test_data = train_test_split(data, test_size = test_size, random_state = seed, stratify = data["virus_name"])

    if valid_size > 0:
        valid_ratio = valid_size / (valid_size + train_size)
        train_data, valid_data = train_test_split(train_valid_data, test_size = valid_ratio, random_state = seed, stratify = train_valid_data["virus_name"])
    else:
        valid_data = train_valid_data.iloc[0:0]
        train_data = train_valid_data

    return train_data, valid_data, test_data

def rebalance_data(data: pd.DataFrame, seed: int, strategy: str = "over"):
    X = data[["repeat"]]
    y = data["virus_name"]

    if strategy == "over":
        ros = RandomOverSampler(random_state = seed)
        Xr, yr = ros.fit_resample(X, y)
    elif strategy == "under":
        rus = RandomUnderSampler(random_state = seed)
        Xr, yr = rus.fit_resample(X, y)
    else:
        return data
    
    out = pd.concat([Xr, yr], axis = 1)
    out.columns = ["repeat", "virus_name"]

    return out

@dataclass
class Metrics:
    accuracy: float
    precision: float
    recall: float
    f1: float
    confusion_matrix: np.ndarray

    def __str__(self) -> str:
        cm = np.asarray(self.confusion_matrix, dtype = int)
        lines = [" ".join(f"{v:5d}" for v in row) for row in cm]
        cm_str = "\n".join(lines)

        return (
            "Metrics:\n\n"
            f"  accuracy    : {self.accuracy:.4f}\n"
            f"  precision   : {self.precision:.4f}\n"
            f"  recall      : {self.recall:.4f}\n"
            f"  f1          : {self.f1:.4f}\n\n"
            f"  confusion_matrix:\n{cm_str}"
        )
    
def calculate_metrics(y_test, y_pred) -> Metrics:
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average = "weighted")
    recall = recall_score(y_test, y_pred, average = "weighted")
    f1 = f1_score(y_test, y_pred, average = "weighted")
    cm = confusion_matrix(y_test, y_pred, labels = VIRUS_NAMES)

    return Metrics(accuracy, precision, recall, f1, cm)

AMINO_ACIDS = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
PAD = "_"
ALPHABET = AMINO_ACIDS + [PAD]
VIRUS_NAMES = ["BCOV", "BAT_SARS", "MERS", "SARS1"]

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
    # data_path = "proteins/data/virus_protein_repeats.csv"
    # seed = 561
    # train_size = 0.8
    # valid_size = 0
    # test_size = 0.2
    # columns_to_drop = ["left_end", "right_end", "left_start", "right_start", "sequence_id", "length", "protein_name", "repeat_type"]
    # train, valid, test = preprocess_data(data_path, seed, train_size, valid_size, test_size)
    # print(train["virus_name"].value_counts())
    # print(valid["virus_name"].value_counts())
    # print(test["virus_name"].value_counts())

    df = pd.DataFrame({"repeat": ["ANQ"]})
    pos = PositionalEncoder(width = 7)

    df = pos.transform(df)
    print(df)