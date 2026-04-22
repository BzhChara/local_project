import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# 数据加载
data_path = "combined_encode_esm_no-z-score.csv"
df = pd.read_csv(data_path, header=None)

# 特征和标签分离
X = df.iloc[:, :1280 + 314].values  
y = df.iloc[:, -1].values   

# 对目标变量进行标准化
scaler = StandardScaler()
y_scaled = scaler.fit_transform(y.reshape(-1, 1)).flatten()

# 数据集划分
X_train, X_test, y_train, y_test = train_test_split(
    X, y_scaled, test_size=0.2, random_state=42
)

# 转换为LightGBM数据集格式
train_data = lgb.Dataset(X_train, label=y_train)
test_data = lgb.Dataset(X_test, label=y_test, reference=train_data)

# 设置模型参数
params = {
    'objective': 'regression',  # 回归任务
    'metric': 'rmse',          # 评估指标
    'boosting_type': 'gbdt',    # 基础提升树类型
    'learning_rate': 0.1,       # 学习率
    'num_leaves': 31,          # 叶子数量
    'max_depth': -1,           # 树深度无限制
    'min_data_in_leaf': 20,    # 叶子最小样本数
    'feature_fraction': 0.8,   # 特征采样比例
    'bagging_fraction': 0.8,   # 数据采样比例
    'bagging_freq': 5,         # 每5次迭代执行采样
    'verbose': 0               # 不输出训练过程
}

# 训练模型
model = lgb.train(
    params,
    train_data,
    num_boost_round=2000,
    valid_sets=[test_data],
    callbacks=[
        lgb.early_stopping(stopping_rounds=50),  # 50轮无改善则停止
        lgb.log_evaluation(50)
    ]
)

# 预测和评估
y_pred = model.predict(X_test)
y_pred_original = scaler.inverse_transform(y_pred.reshape(-1, 1))  # 还原到原始量纲

mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"MSE: {mse:.4f}")
print(f"RMSE: {rmse:.4f}")
print(f"MAE: {mae:.4f}")
print(f"R² Score: {r2:.4f}")