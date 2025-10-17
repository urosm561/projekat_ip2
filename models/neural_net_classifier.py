import numpy as np
import pandas as pd
import torch.nn as nn
import seaborn as sns 
import matplotlib.pyplot as plt

import torch

from torch.utils.data import TensorDataset, DataLoader
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from typing import List, Tuple, Dict

from preprocess_data import preprocess_data, PositionalEncoder, CATS_21, Metrics, VIRUS_NAMES

class SimpleMLP(nn.Module):
    def __init__(self, in_dim: int, num_classes: int, p_drop: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 256),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(64, num_classes)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)

def class_weights_from_train(y: np.ndarray, VIRUS_NAMES: List[str]) -> torch.Tensor:
    counts = pd.Series(y).value_counts().reindex(VIRUS_NAMES).fillna(0).astype(int)
    inv = 1.0 / counts.replace(0, np.nan)
    inv = inv / inv.max()
    inv = inv.fillna(0.0)
    return torch.tensor(inv.values, dtype=torch.float32)

def build_features(train_df: pd.DataFrame, test_df: pd.DataFrame, max_aa_length: int = 7):
    pos = PositionalEncoder(width = max_aa_length)
    ohe = OneHotEncoder(
        categories = [CATS_21] * max_aa_length,
        handle_unknown = "error",
        sparse_output = False
    )
    prep = Pipeline([
        ("pos", pos),
        ("ohe", ColumnTransformer(
            transformers = [("onehot", ohe, [f"p{i + 1}" for i in range(max_aa_length)])],
            remainder = "drop"
        ))
    ])
    X_train = prep.fit_transform(train_df[["repeat"]])
    X_test = prep.transform(test_df[["repeat"]])
    return X_train.astype(np.float32), X_test.astype(np.float32)

def encode_labels(train_df: pd.DataFrame, test_df: pd.DataFrame, target_col: str) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str,int]]:
    class2idx = {c: i for i, c in enumerate(VIRUS_NAMES)}
    y_train = train_df[target_col].astype(str).map(class2idx).values
    y_test = test_df[target_col].astype(str).map(class2idx).values
    return y_train, y_test

def train_simple_nn_classifier(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    seed: int,
    max_aa_length: int = 7,
    batch_size: int = 1024,
    lr: float = 1e-3,
    epochs: int = 20,
    p_drop: float = 0.1,
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
) -> Tuple[nn.Module, Metrics]:
    torch.manual_seed(seed)
    np.random.seed(seed)

    X_train, X_test = build_features(train_df, test_df, max_aa_length=max_aa_length)
    y_train, y_test = encode_labels(train_df, test_df, target_col="virus_name")

    in_dim = X_train.shape[1]
    num_classes = len(VIRUS_NAMES)

    Xtr = torch.from_numpy(X_train)
    ytr = torch.from_numpy(y_train).long()
    Xte = torch.from_numpy(X_test)
    yte = torch.from_numpy(y_test).long()

    train_loader = DataLoader(TensorDataset(Xtr, ytr), batch_size=batch_size, shuffle=True)
    test_loader  = DataLoader(TensorDataset(Xte, yte), batch_size=batch_size, shuffle=False)

    model = SimpleMLP(in_dim=in_dim, num_classes=num_classes, p_drop=p_drop).to(device)
    weights = class_weights_from_train(train_df["virus_name"].astype(str).values, VIRUS_NAMES).to(device)
    criterion = nn.CrossEntropyLoss(weight = weights)
    opt = torch.optim.AdamW(model.parameters(), lr = lr)

    model.train()
    for _ in range(epochs):
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            loss = criterion(logits, yb)

            opt.zero_grad()
            loss.backward()
            opt.step()

    model.eval()
    preds, trues = [], []
    with torch.no_grad():
        for xb, yb in test_loader:
            xb = xb.to(device)
            logits = model(xb)
            pred = logits.argmax(dim=1).cpu().numpy()
            preds.append(pred)
            trues.append(yb.numpy())

    y_pred = np.concatenate(preds)
    y_true = np.concatenate(trues)

    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average = "weighted", zero_division=0)
    recall  = recall_score(y_true, y_pred, average = "weighted", zero_division=0)
    f1 = f1_score(y_true, y_pred, average = "weighted", zero_division=0)
    cm = confusion_matrix(y_true, y_pred, labels = list(range(num_classes)))

    return model, Metrics(accuracy, precision, recall, f1, cm)

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.8
    valid_size = 0
    test_size = 0.2
    columns_to_drop = ["left_end", "right_end", "left_start", "right_start", "sequence_id", "length", "protein_name", "repeat_type"]
    train, _, test = preprocess_data(data_path, seed, train_size, valid_size, test_size)

    model, metrics = train_simple_nn_classifier(train, test, seed)

    print(metrics)

    sns.heatmap(metrics.confusion_matrix, annot = True, fmt = "d", cmap = "inferno",
                xticklabels = VIRUS_NAMES, yticklabels = VIRUS_NAMES)
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.show()
