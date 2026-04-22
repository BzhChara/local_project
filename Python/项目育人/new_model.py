import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from evaluation_reg import *
import matplotlib.pyplot as plt
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, roc_auc_score, recall_score, roc_curve, precision_score, f1_score
from sklearn.metrics import matthews_corrcoef
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.utils.data.sampler import WeightedRandomSampler

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


data = pd.read_csv('dataset_finger_1280.csv')
X = data.iloc[:, 2:]
y = data.iloc[:, 1]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
class DiabetesDataset(Dataset):
    def __init__(self,X,y):
        self.len = X.shape[0]
        self.x_data = torch.from_numpy(X.values)
        self.y_data = torch.from_numpy(y.values)
    def __getitem__(self, index):
        return self.x_data[index],self.y_data[index]
    def __len__(self):
        return self.len

train_dataset = DiabetesDataset(X_train,y_train)
# positive = []
# negative = []
# for x_data, y_data in train_dataset:
#     if y_data == 1:
#         positive.append(y_data)
#     else:
#         negative.append(y_data)
# positive_count = len(positive)/(len(positive)+len(negative))
# negative_count = len(negative)/(len(positive)+len(negative))
# positive_weight = 1/positive_count
# negative_count = 1/negative_count
# weights_train = [positive_weight if y_data == negative_count else 1 for x_data, y_data in train_dataset]
# sampler_train = WeightedRandomSampler(weights_train,len(train_dataset), replacement=True)
# train_loader = DataLoader(train_dataset, batch_size=128, sampler=sampler_train)
train_loader = DataLoader(dataset=train_dataset,batch_size=128,shuffle=True,num_workers=2,drop_last=True)


test_dataset = DiabetesDataset(X_test,y_test)
test_loader = DataLoader(dataset=test_dataset,batch_size=128,shuffle=True,num_workers=2,drop_last=True)
# for i, data in enumerate(train_loader, 0):
#     inputs, lables = data


## methods:如果label为1，那么对应的该类别被取出来的概率是另外一个类别的2倍


# weights_test = [2 if y_data == 1 else 1 for x_data, y_data in test_dataset]
# sampler_test = WeightedRandomSampler(weights_test,len(test_dataset), replacement=True)
# test_loader = DataLoader(test_dataset, batch_size=128, sampler=sampler_test)

# X_train_numpy = X_train.values
# X_test_numpy = X_test.values
# y_train_numpy = y_train.values
# y_test_numpy = y_test.values
# #
# #
# #
# prot_X_train_numpy = X_train_numpy[:, 0:1280]
# glycan_X_train_numpy = X_train_numpy[:, 1280:1685]
# prot_X_test_numpy = X_test_numpy[:, 0:1280]
# glycan_X_test_numpy = X_test_numpy[:, 1280:1685]
# print(prot_X_train_numpy.shape)
# print(glycan_X_train_numpy.shape)
#
#
# prot_X_train_tensor = torch.tensor(prot_X_train_numpy, dtype=torch.float32).to(device)
# glycan_X_train_tensor = torch.tensor(glycan_X_train_numpy, dtype=torch.float32).to(device)
# y_train_tensor = torch.tensor(y_train_numpy, dtype=torch.float32).view(-1, 1).to(device)
# prot_X_test_tensor = torch.tensor(prot_X_test_numpy, dtype=torch.float32).to(device)
# glycan_X_test_tensor = torch.tensor(glycan_X_test_numpy, dtype=torch.float32).to(device)
# y_test_tensor = torch.tensor(y_test_numpy, dtype=torch.float32).view(-1, 1).to(device)


class GlycoproteinProphet(nn.Module):
    def __init__(self):
        super(GlycoproteinProphet, self).__init__()
        self.prot_fc1 = nn.Linear(1280, 64)
        self.prot_fc2 = nn.Linear(64, 32)
        self.prot_dropout1 = nn.Dropout(0.3)
        self.prot_dropout2 = nn.Dropout(0.2)
        self.bn_prot1 = nn.BatchNorm1d(64)
        self.bn_prot2 = nn.BatchNorm1d(32)
        self.activation_fn = nn.GELU()
        self.glycan_fc1 = nn.Linear(402, 128)
        self.glycan_lstm = nn.LSTM(402, 64, 2, batch_first=True)
        self.conv1 = nn.Conv1d(128, 64, 1)
        self.glycan_rnn = nn.RNN(64, 64, 2)
        self.glycan_f2 = nn.Linear(64, 32)
        self.bn_glycan1 = nn.BatchNorm1d(32)

        self.bn_fc1 = nn.Linear(64, 32)
        self.bn_fc2 = nn.Linear(32, 12)
        self.bn_fc3 = nn.Linear(12, 1)
        self.bn_relu = nn.ReLU()
        # Attention
        self.W_query = nn.Linear(403, 64)
        self.W_key = nn.Linear(403, 64)
        self.W_value = nn.Linear(403, 64)
        self.softmax = nn.Softmax(dim=1)
        self.attention_glycan_fc1 = nn.Linear(467, 32)

    def forward(self, prot_X_train_tensor, glycan_X_train_tensor):
        prot_X_train_tensor = prot_X_train_tensor.float()
        glycan_X_train_tensor = glycan_X_train_tensor.float()
        x = self.prot_fc1(prot_X_train_tensor)
        prot1 = self.bn_prot1(self.prot_dropout1(x))
        prot2 = self.bn_prot2(self.prot_dropout2(self.prot_fc2(prot1)))
        # LSTM
        # glycan1,_= self.glycan_lstm(self.glycan_fc1(glycan_X_train_tensor))
        # glycan_X_train_tensor = glycan_X_train_tensor.transpose(0,1)
        # CNN
        # glycan1 = self.conv1(glycan_X_train_tensor)
        # glycan1 = glycan1.transpose(0, 1)
        # glycan1,_=self.glycan_rnn(glycan1)
        # glycan1 = glycan1.transpose(0, 1)
        # glycan1 = self.conv2(glycan1)
        # glycan1 = glycan1.transpose(0, 1)
        # glycan1 = self.conv1(self.glycan_fc1(glycan_X_train_tensor).transpose(0,1)).transpose(0,1)
        # #glycan2 = self.pool(glycan1).transpose(0,1)
        # #glycan2 = glycan2.view(glycan2.size(0), -1)
        # glycan3 = self.bn_glycan1(self.activation_fn(self.glycan_f2(glycan1)))
        # LSTM+RNN
        # glycan1, _ = self.glycan_lstm(glycan_X_train_tensor)
        # # glycan2,_ = self.glycan_rnn(glycan1)
        # glycan3 = self.bn_glycan1(self.glycan_f2(glycan1))
        # Attention
        query = self.W_query(glycan_X_train_tensor)
        key = self.W_key(glycan_X_train_tensor)
        value = self.W_value(glycan_X_train_tensor)
        key = key.transpose(0, 1)
        scores = torch.matmul(query, key)
        #attention_weights = torch.exp(scores)/torch.sum(scores, dim=1, keepdim=True)
        attention_weights = self.softmax(scores)
        #attention_weights = F.normalize(scores, p=2, dim=1)
        weighted_values = torch.matmul(attention_weights, value)
        #print(weighted_values.shape)
        glycan_feature = torch.cat((glycan_X_train_tensor, weighted_values), dim=1)
        glycan = self.attention_glycan_fc1(glycan_feature)
        h_n = torch.cat((prot2, glycan), 1)
        #x = F.sigmoid(self.bn_fc3(self.bn_relu(self.bn_fc2(self.bn_relu(self.bn_fc1(h_n))))))
        x = self.bn_fc3(self.bn_relu(self.bn_fc2(self.bn_relu(self.bn_fc1(h_n)))))
        return x


GlycoproteinProphet = GlycoproteinProphet().to(device)
# 定义损失函数和优化
#criterion = nn.BCELoss().to(device)
criterion = nn.MSELoss().to(device)
optimizer = optim.Adam(GlycoproteinProphet.parameters(), lr=0.001)

# 学习率调度器
scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.5)

# 保存每个epoch的训练损失和验证损失


num_epochs = 450

train = []
for epoch in range(num_epochs):
    print(1)
    train_losses = []
    #train_pred = []
    #train_labels = []
    GlycoproteinProphet.train().to(device)
    for data in train_loader:
        optimizer.zero_grad()
        inputs,lables = data
        lables = lables.float().to(device)
        prot_X_train_tensor = inputs[:, 0:1280].to(device)
        #print(prot_X_train_tensor.shape)
        glycan_X_train_tensor = inputs[:, 1280:1685].to(device)
        #print(glycan_X_train_tensor.shape)
        outputs = GlycoproteinProphet(prot_X_train_tensor,glycan_X_train_tensor).squeeze().to(device)
        loss = criterion(outputs, lables).to(device)
        #print(epoch,loss.item())
        loss.backward()
        optimizer.step()
        #scheduler.step() # 更新学习率
        train_losses.append(loss.item())# 记录训练损失
        train.append(np.mean(train_losses))
    #scheduler.step()
    if (epoch + 1) % 10 == 0:
        print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {np.mean(train_losses):.4f}')

result_test = []
lables_test = []
GlycoproteinProphet.eval().to(device)
with torch.no_grad():
    for data in test_loader:
        inputs, lables = data
        lables = lables.float().to(device)
        prot_X_test_tensor = inputs[:, 0:1280].to(device)
        # print(prot_X_train_tensor.shape)
        glycan_X_test_tensor = inputs[:, 1280:1685].to(device)
        # print(glycan_X_train_tensor.shape)
        outputs = GlycoproteinProphet(prot_X_test_tensor, glycan_X_test_tensor).squeeze().to(device)
        result_test.append(outputs.tolist())
        lables_test.append(lables.tolist())
my_result = [item for sublist in result_test for item in sublist]
my_labels = [item for sublist in lables_test for item in sublist]
evaluate_regression2(np.array(my_labels), np.array(my_result))
y_test_df = pd.DataFrame({'y_test':my_labels})
#y_test_df.to_csv('y_test_320.csv')
y_pred_df = pd.DataFrame({'y_pred': my_result})
#y_pred_df.to_csv('y_pred_320.csv')
# 将两个dataframe合并为一个
y_df = pd.concat([y_test_df, y_pred_df], axis=1)
# 将dataframe导出为csv文件
y_df.to_csv('Attention_1280.csv', index=False)



