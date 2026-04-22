import h2o
from h2o.automl import H2OAutoML
from sklearn.preprocessing import QuantileTransformer
from sklearn.metrics import root_mean_squared_error, mean_squared_error, mean_absolute_error, r2_score
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 初始化H2O
h2o.init(max_mem_size='12G')

# 加载数据
train = h2o.import_file("single_protein_prediction_train2.csv", header=0)
test = h2o.import_file("single_protein_prediction_test2.csv", header=0)

# 定义目标变量列名
protein_cols = train.columns[:1280]        # 前1280列为蛋白质特征
oligo_cols = train.columns[1280:1280+342]  # 接下来的342列为寡糖特征
y_col = train.columns[-1]                  # 最后一列为目标变量RFU

# 使用分位数变换（映射到正态分布）
quantile_transformer = QuantileTransformer(output_distribution='normal', random_state=42)

# 仅对训练集拟合
train_rfu_quantile = quantile_transformer.fit_transform(train[y_col].as_data_frame(use_multi_thread=True).values)
test_rfu_quantile = quantile_transformer.transform(test[y_col].as_data_frame(use_multi_thread=True).values)

# 替换原有的RFU_scaled列
train_final = train[protein_cols + oligo_cols].cbind(
    h2o.H2OFrame(train_rfu_quantile, column_names=["RFU_scaled"])
)
test_final = test[protein_cols + oligo_cols].cbind(
    h2o.H2OFrame(test_rfu_quantile, column_names=["RFU_scaled"])
)

# 训练AutoML
aml = H2OAutoML(max_runtime_secs=1800*2, 
                nfolds=5, 
                seed=42, 
                stopping_metric="RMSE",
                sort_metric="RMSE",
                project_name='rfu_prediction')

aml.train(
    x=protein_cols + oligo_cols,  # 使用所有原始特征
    y="RFU_scaled",
    training_frame=train_final
)

# 查看模型性能
print(aml.leaderboard)

# 预测测试集
pred_scaled = aml.leader.predict(test_final)

# 获取标准化后的真实值
true_scaled = test_final["RFU_scaled"].as_data_frame(use_multi_thread=True).values.flatten()
pred_scaled = pred_scaled.as_data_frame(use_multi_thread=True).values.flatten()

# 计算指标（标准化空间）
print("\n标准化空间评估指标:")
print(f"RMSE: {root_mean_squared_error(true_scaled, pred_scaled):.4f}")
print(f"MSE: {mean_squared_error(true_scaled, pred_scaled):.4f}")
print(f"MAE: {mean_absolute_error(true_scaled, pred_scaled):.4f}")
print(f"R²: {r2_score(true_scaled, pred_scaled):.4f}")

# ==== 新增的保存结果代码 ==== 
# 生成序列号（从1开始）
indices = np.arange(1, len(pred_scaled) + 1)

# 创建包含序列号和预测值的DataFrame
results = pd.DataFrame({
    'index': indices,
    'predicted_RFU_scaled': pred_scaled
})

# 保存为CSV文件（空格分隔，无列名）
results.to_csv('single_protein_prediction_result.csv', 
               sep=' ', 
               header=False, 
               index=False)
# ========================

# 绘制标准化残差图
plt.figure(figsize=(10,6))
plt.scatter(pred_scaled, pred_scaled - true_scaled, alpha=0.5)
plt.axhline(y=0, color='r', linestyle='--')
plt.xlabel("Predicted RFU (Z-score)")
plt.ylabel("Residuals (Z-score)")
plt.title("Standardized Residual Plot")
plt.show()

# 释放资源
h2o.cluster().shutdown()