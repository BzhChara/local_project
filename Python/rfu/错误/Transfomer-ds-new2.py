import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

# Device configuration
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Define feature extractors
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, vocab_size=21, embedding_dim=128, hidden_dim=256, num_layers=2):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True, dropout=0.2)
        self.layer_norm = nn.LayerNorm(hidden_dim)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, 128),
            nn.GELU(),
            nn.Dropout(0.3)
        )

    def forward(self, x):
        x = self.embedding(x)
        lstm_out, _ = self.lstm(x)
        return self.fc(self.layer_norm(lstm_out[:, -1, :]))

class GlycanFeatureExtractor(nn.Module):
    def __init__(self, input_dim=342, embedding_dim=128):
        super().__init__()
        self.fc = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(256, embedding_dim)
        )
        
        self.attention = nn.MultiheadAttention(embedding_dim, 4, dropout=0.2)
        self.fc_out = nn.Sequential(
            nn.LayerNorm(embedding_dim),
            nn.Linear(embedding_dim, 128),
            nn.GELU()
        )

    def forward(self, x):
        x = self.fc(x.float())
        x = x.unsqueeze(1)
        attn_out, _ = self.attention(x, x, x)
        return self.fc_out(attn_out.squeeze(1))

# Transformer model
class TransformerRegressionModel(nn.Module):
    def __init__(self, protein_extractor, glycan_extractor, d_model=256, nhead=8, num_layers=4):
        super().__init__()
        self.protein_extractor = protein_extractor
        self.glycan_extractor = glycan_extractor

        self.protein_proj = nn.Sequential(
            nn.Linear(128, d_model),
            nn.LayerNorm(d_model),
            nn.Dropout(0.2)
        )
        
        self.glycan_proj = nn.Sequential(
            nn.Linear(128, d_model),
            nn.LayerNorm(d_model),
            nn.Dropout(0.2)
        )
        
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=512, dropout=0.2, activation='gelu', batch_first=True)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers)
        
        self.regressor = nn.Sequential(
            nn.Linear(d_model, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(256, 1)
        )

        self._init_weights()

    def _init_weights(self):
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.kaiming_normal_(module.weight, nonlinearity='relu')
                nn.init.constant_(module.bias, 0)

    def forward(self, protein_input, glycan_input):
        protein_feat = self.protein_proj(self.protein_extractor(protein_input))
        glycan_feat = self.glycan_proj(self.glycan_extractor(glycan_input))

        features = torch.stack([protein_feat, glycan_feat], dim=1)

        transformer_out = self.transformer(features)

        aggregated = transformer_out.mean(dim=1) + transformer_out.amax(dim=1)
        return self.regressor(aggregated).squeeze()

# Dataset class
class CustomDataset(Dataset):
    def __init__(self, protein_data, glycan_data, targets):
        self.protein_data = torch.tensor(protein_data, dtype=torch.long)
        self.glycan_data = torch.tensor(glycan_data, dtype=torch.float32)
        self.targets = torch.tensor(targets, dtype=torch.float32)

    def __len__(self):
        return len(self.targets)

    def __getitem__(self, idx):
        return self.protein_data[idx], self.glycan_data[idx], self.targets[idx]

# Load and split data
combined_data = pd.read_csv("combined_encode.csv", header=None)
X = combined_data.iloc[:, :-1].values
y = combined_data.iloc[:, -1].values

# Proper three-way split
X_temp, X_test, y_temp, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.25, random_state=42)  # 0.25*0.8 = 0.2

# Split features
def split_features(X):
    return X[:, :2710], X[:, 2710:2710+342]

X_protein_train, X_glycan_train = split_features(X_train)
X_protein_val, X_glycan_val = split_features(X_val)
X_protein_test, X_glycan_test = split_features(X_test)

# Create datasets
train_dataset = CustomDataset(X_protein_train, X_glycan_train, y_train)
val_dataset = CustomDataset(X_protein_val, X_glycan_val, y_val)
test_dataset = CustomDataset(X_protein_test, X_glycan_test, y_test)

# DataLoaders
batch_size = 64
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size)
test_loader = DataLoader(test_dataset, batch_size=batch_size)

# Initialize model
protein_extractor = ProteinFeatureExtractor().to(device)
glycan_extractor = GlycanFeatureExtractor().to(device)
model = TransformerRegressionModel(protein_extractor, glycan_extractor).to(device)

optimizer = optim.AdamW(model.parameters(), lr=1e-4, weight_decay=0.01)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, 
    mode='min', 
    factor=0.5, 
    patience=3,
    min_lr=1e-6
)
criterion = nn.HuberLoss()

# Training loop with early stopping
def train_model(model, train_loader, val_loader, epochs=200):
    best_loss = float('inf')
    patience = 10
    no_improve = 0
    
    for epoch in range(epochs):
        # Training phase
        model.train()
        train_loss = 0
        for proteins, glycans, labels in train_loader:
            proteins = proteins.to(device).long()
            glycans = glycans.to(device)
            labels = labels.to(device)
            
            optimizer.zero_grad()
            outputs = model(proteins, glycans)
            loss = criterion(outputs, labels)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            train_loss += loss.item() * len(labels)
        
        # Validation phase
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for proteins, glycans, labels in val_loader:
                proteins = proteins.to(device).long()
                glycans = glycans.to(device)
                labels = labels.to(device)
                
                outputs = model(proteins, glycans)
                val_loss += criterion(outputs, labels).item() * len(labels)
        
        # Update learning rate
        scheduler.step(val_loss)
        
        # Calculate metrics
        train_loss /= len(train_dataset)
        val_loss /= len(val_dataset)
        print(f"Epoch {epoch+1:02d} | "f"Train Loss: {train_loss:.4f} | "f"Val Loss: {val_loss:.4f}", flush=True)
        
        # Early stopping check
        if val_loss < best_loss:
            best_loss = val_loss
            no_improve = 0
            torch.save(model.state_dict(), 'best_model.pth')
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"Early stopping at epoch {epoch+1}", flush=True)
                break

# Run training
train_model(model, train_loader, val_loader)

# Final evaluation
def evaluate_model(model, test_loader):
    model.load_state_dict(torch.load('best_model.pth', weights_only=True))
    model.eval()
    
    predictions, actuals = [], []
    with torch.no_grad():
        for proteins, glycans, labels in test_loader:
            proteins = proteins.to(device).long()
            glycans = glycans.to(device)
            
            preds = model(proteins, glycans).cpu().numpy()
            predictions.extend(preds)
            actuals.extend(labels.numpy())
    
    # Calculate metrics
    mae = mean_absolute_error(actuals, predictions)
    mse = mean_squared_error(actuals, predictions)
    r2 = r2_score(actuals, predictions)
    
    print(f"\nFinal Test Metrics:", flush=True)
    print(f"MAE: {mae:.4f}", flush=True)
    print(f"MSE: {mse:.4f}", flush=True)
    print(f"R²: {r2:.4f}", flush=True)
    
    # Residual plot
    plt.figure(figsize=(10, 6))
    residuals = np.array(actuals) - predictions
    plt.scatter(predictions, residuals, alpha=0.5)
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Residual Plot')
    plt.axhline(0, color='red', linestyle='--')
    plt.savefig('result2.png', bbox_inches='tight')
    plt.close()

# Run evaluation
evaluate_model(model, test_loader)