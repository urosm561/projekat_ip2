import numpy as np
import pandas as pd
import torch.nn as nn

import torch

from torch.utils.data import TensorDataset, DataLoader
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import ParameterGrid
from typing import List, Tuple, Dict, Optional

from preprocess_data import preprocess_data, PositionalEncoder, ALPHABET, Metrics, VIRUS_NAMES, calculate_metrics, rebalance_data, log_training, split_into_proteins

class SimpleMLP(nn.Module):
    def __init__(
            self, 
            in_dim: int, 
            num_classes: int, 
            n1: int = 256,
            n2: int = 128,
            n3: int = 64,
            p_drop: float = 0.1
    ):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, n1),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(n1, n2),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(n2, n3),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(n3, num_classes)
        )

    def forward(
            self, 
            x: torch.Tensor
    ) -> torch.Tensor:
        return self.net(x)

def class_weights_from_train(
        y: np.ndarray, 
        VIRUS_NAMES: List[str]
) -> torch.Tensor:
    """Calculate the class weights, because the classes are heavily imbalanced."""
    counts = pd.Series(y).value_counts().reindex(VIRUS_NAMES).fillna(0).astype(int)
    inv = 1.0 / counts.replace(0, np.nan)
    inv = inv / inv.max()
    inv = inv.fillna(0.0)
    return torch.tensor(inv.values, dtype=torch.float32)

def build_features(
        train_data: pd.DataFrame, 
        test_data: pd.DataFrame, 
        max_aa_length: int = 7
):
    """
    Prepare the data by positionally encoding amino acids up to the `max_aa_length`,
    and then one-hot encoding these columns. We get `max_aa_length` x 21 columns.
    """
    pos = PositionalEncoder(width = max_aa_length)
    ohe = OneHotEncoder(
        categories = [ALPHABET] * max_aa_length,
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

    X_train = prep.fit_transform(train_data[["repeat"]])

    if not test_data.empty:
        X_test = prep.transform(test_data[["repeat"]])
    else:
        X_test = np.empty((0, X_train.shape[1]), dtype = X_train.dtype)

    return X_train.astype(np.float32), X_test.astype(np.float32)

def encode_labels(
        train_data: pd.DataFrame, 
        test_data: pd.DataFrame, 
        target_col: str
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str,int]]:
    """
    Encode the labels in the simplest possible way.
    """
    class2idx = {c: i for i, c in enumerate(VIRUS_NAMES)}
    y_train = train_data[target_col].astype(str).map(class2idx).values
    y_test = test_data[target_col].astype(str).map(class2idx).values
    return y_train, y_test

def train_nn_classifier(
    train_data: pd.DataFrame,
    test_data: pd.DataFrame,
    seed: int,
    max_aa_length: int = 7,
    resample_strategy: Optional[str] = None,
    batch_size: int = 1024,
    lr: float = 1e-3,
    epochs: int = 20,
    p_drop: float = 0.1,
    n1: int = 256,
    n2: int = 128, 
    n3: int = 64,
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
    patience = 5
) -> Tuple[nn.Module, Metrics]:
    """
    Train the neural network and then evaluate it.

    Args:
        train_data (pd.DataFrame): Data to train the nn on.
        test_data (pd.DataFrame): Data to evaluate the nn on.
        seed (int): Control randomness.
        max_aa_length (int, optional): Maximum allowed number of amino acids in a repeat. Defaults to 7.
        batch_size (int, optional): Batch size for the DataLoaders. Defaults to 1024.
        lr (float, optional): Learning rate for the Adam optimizer. Defaults to 1e-3.
        epochs (int, optional): Number of epochs for training. Defaults to 20.
        p_drop (float, optional): Dropout parameter. Defaults to 0.1.
        device (str, optional): torch device.

    Returns:
        Tuple[nn.Module, Metrics]: Trained network and all the relevant metrics.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    train_data = rebalance_data(train_data, seed = seed, strategy = resample_strategy)

    X_train, X_test = build_features(train_data, test_data, max_aa_length=max_aa_length)
    y_train, y_test = encode_labels(train_data, test_data, target_col="virus_name")

    in_dim = X_train.shape[1]
    num_classes = len(VIRUS_NAMES)

    Xtr = torch.from_numpy(X_train)
    ytr = torch.from_numpy(y_train).long()
    Xte = torch.from_numpy(X_test)
    yte = torch.from_numpy(y_test).long()

    train_loader = DataLoader(TensorDataset(Xtr, ytr), batch_size=batch_size, shuffle=True)
    test_loader  = DataLoader(TensorDataset(Xte, yte), batch_size=batch_size, shuffle=False)

    model = SimpleMLP(in_dim=in_dim, num_classes=num_classes, p_drop=p_drop, n1=n1, n2=n2, n3=n3).to(device)
    weights = class_weights_from_train(train_data["virus_name"].astype(str).values, VIRUS_NAMES).to(device)
    criterion = nn.CrossEntropyLoss(weight = weights)
    opt = torch.optim.AdamW(model.parameters(), lr = lr)

    best_val_loss = float("inf")
    epochs_since_improvement = 0

    metrics = None

    model.train()
    for epoch in range(epochs):
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            loss = criterion(logits, yb)

            opt.zero_grad()
            loss.backward()
            opt.step()

        model.eval()
        preds, trues = [], []
        val_loss = 0
        val_samples = 0

        if len(test_loader) > 0:

            with torch.no_grad():
                for xb, yb in test_loader:
                    xb = xb.to(device)
                    logits = model(xb)
                    loss = criterion(logits, yb).sum()
                    val_loss += loss.item()

                    pred = logits.argmax(dim=1).cpu().numpy()
                    preds.append(pred)
                    trues.append(yb.numpy())
                    val_samples += yb.size(0)

            y_pred = np.concatenate(preds)
            y_test = np.concatenate(trues)

            y_pred = [VIRUS_NAMES[i] for i in y_pred.tolist()]
            y_test = [VIRUS_NAMES[i] for i in y_test.tolist()]

            metrics = calculate_metrics(y_test, y_pred)

            val_loss /= val_samples

            # print(f"{epoch}/{epochs} - loss: {val_loss}")

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                epochs_since_improvement = 0
            else:
                epochs_since_improvement += 1

                if epochs_since_improvement == patience:
                    break

    return model, metrics

def choose_best_nn(
        train_data, 
        valid_data, 
        test_data, 
        seed, 
        max_aa_length, 
        resample_strategy,    
        p_drop,
        n1,
        n2, 
        n3,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
        batch_size: int = 1024,
        lr: float = 1e-3,
        epochs: int = 20,
        patience = 5,
        results_folder = "results"
):
    grid = ParameterGrid({
        "max_aa_length": max_aa_length,
        "resample_strategy": resample_strategy,
        "p_drop": p_drop,
        "n1": n1,
        "n2": n2,
        "n3": n3
    })
    
    best_val_f1 = -1.0
    best_val_metrics = None
    best_val_params = None
    rows = []

    for params in grid:
        _, metrics = train_nn_classifier(train_data, valid_data, seed, **params, 
                                         device = device, batch_size = batch_size, 
                                         lr = lr, epochs = epochs, patience = patience)
        f1 = metrics.f1
        
        rows.append({
            **params,
            "accuracy": metrics.accuracy,
            "precision": metrics.precision,
            "recall": metrics.recall,
            "f1": metrics.f1
        })

        if f1 > best_val_f1:
            best_val_f1 = f1
            best_val_params = params
            best_val_metrics = metrics

    leaderboard = pd.DataFrame(rows).sort_values("f1", ascending = False).reset_index(drop = True)

    _, metrics = train_nn_classifier(pd.concat([train_data, valid_data], ignore_index = True), 
                                     test_data, seed, **best_val_params, 
                                     device = device, batch_size = batch_size, 
                                     lr = lr, epochs = epochs, patience = patience)
    
    model, _ = train_nn_classifier(pd.concat([train_data, valid_data, test_data], ignore_index = True), 
                                     pd.DataFrame(columns = ["repeat", "virus_name"]), seed, **best_val_params, 
                                     device = device, batch_size = batch_size, 
                                     lr = lr, epochs = epochs, patience = patience)
    
    log_training(results_folder, "nn", leaderboard, best_val_params, metrics, seed)

    return model, metrics, best_val_params

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.7
    valid_size = 0.2
    test_size = 0.1
    
    spike_data, nucleocapsid_data = split_into_proteins(data_path)
    
    train_spike, valid_spike, test_spike = preprocess_data(spike_data, seed, train_size, valid_size, test_size)
    train_nucleocapsid, valid_nucleocapsid, test_nucleocapsid = preprocess_data(nucleocapsid_data, seed, train_size, valid_size, test_size)

    model_nucleocapsid, metrics_nucleocapsid, params_nucleocapsid = choose_best_nn(
        train_nucleocapsid, valid_nucleocapsid, test_nucleocapsid, seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", None],
        p_drop = [0.1],
        n1 = [256, 128],
        n2 = [128, 64],
        n3 = [64, 32],
        results_folder = "results/nucleocapsid" 
    )

    model_spike, metrics_spike, params_spike = choose_best_nn(
        train_spike, valid_spike, test_spike, seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", None],
        p_drop = [0.1],
        n1 = [256, 128],
        n2 = [128, 64],
        n3 = [64, 32],
        results_folder = "results/spike" 
    )
