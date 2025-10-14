from models.data_split import prepare_data
from utils.embeddings import protvec_encode_repeat

from sklearn.preprocessing import LabelEncoder

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from torch.utils.data import Dataset, DataLoader

class RepeatDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype = torch.float32)
        self.y = torch.tensor(y, dtype = torch.long)

    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

SEED = 561

train, valid, test = prepare_data(data = "proteins/data/virus_protein_repeats.csv", seed = SEED)
X_train = np.array([protvec_encode_repeat(repeat) for repeat in train["repeat"]])

le = LabelEncoder()
y_train = le.fit_transform(train["virus_name"])

train_dataset = RepeatDataset(X_train, y_train)

class ProtVecClassifier(nn.Module):
    def __init__(self, input_dim = 100, hidden_dim = 128, num_classes = 4):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.dropout = nn.Dropout(0.1)
        self.fc2 = nn.Linear(hidden_dim, num_classes)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        return self.fc2(x)
    
def train_model(model, dataloader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    for X_batch, y_batch in dataloader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)

        optimizer.zero_grad()
        outputs = model(X_batch)
        loss = criterion(outputs, y_batch)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(dataloader)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

train_loader = DataLoader(train_dataset, batch_size = 64, shuffle = True)
model = ProtVecClassifier(input_dim = 100, hidden_dim = 128, num_classes = len(le.classes_))
model.to(device)

optimizer = optim.Adam(model.parameters(), lr = 1e-3)
criterion = nn.CrossEntropyLoss()

for epoch in range(1000):
    loss = train_model(model, train_loader, optimizer, criterion, device)
    print(f"Epoch {epoch + 1}/20 - Loss: {loss:.4f}")