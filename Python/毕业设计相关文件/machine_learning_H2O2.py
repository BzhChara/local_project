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
train = h2o.import_file("generalize_train.csv", header=0)

test_bacteria = h2o.import_file("generalize_test_bacteria.csv", header=0)
test_LysM = h2o.import_file("generalize_test_LysM.csv", header=0)
test_toxin = h2o.import_file("generalize_test_toxin.csv", header=0)
test_virus = h2o.import_file("generalize_test_virus.csv", header=0)

# 定义目标变量列名
protein_cols = train.columns[:1280]        # 前1280列为蛋白质特征
oligo_cols = train.columns[1280:1280+342]  # 接下来的342列为寡糖特征
y_col = train.columns[-1]                  # 最后一列为目标变量RFU

# 初始化分位数转换器
quantile_transformer = QuantileTransformer(output_distribution='normal', random_state=42)

# 只在训练集上拟合转换器
train_rfu = train[y_col].as_data_frame(use_multi_thread=True).values
train_rfu_quantile = quantile_transformer.fit_transform(train_rfu)

# 转换训练集
train_final = train[protein_cols + oligo_cols].cbind(
    h2o.H2OFrame(train_rfu_quantile, column_names=["RFU_scaled"])
)

# 定义处理测试集的函数
def process_test_set(test_h2o):
    # 转换测试集的RFU值
    test_rfu = test_h2o[y_col].as_data_frame(use_multi_thread=True).values
    test_rfu_quantile = quantile_transformer.transform(test_rfu)
    
    # 创建处理后的测试集
    return test_h2o[protein_cols + oligo_cols].cbind(
        h2o.H2OFrame(test_rfu_quantile, column_names=["RFU_scaled"])
    )

# 处理所有测试集
test_bacteria_final = process_test_set(test_bacteria)
test_LysM_final = process_test_set(test_LysM)
test_toxin_final = process_test_set(test_toxin)
test_virus_final = process_test_set(test_virus)

# 训练AutoML
aml = H2OAutoML(max_runtime_secs=1800*2, 
                # max_models=5,
                nfolds=5, 
                seed=42, 
                stopping_metric="RMSE",
                sort_metric="RMSE",
                # include_algos=["GBM", "GLM"],
                # exclude_algos=["XGBoost"],
                project_name='rfu_prediction')

aml.train(
    x=protein_cols + oligo_cols,  # 使用所有原始特征
    y="RFU_scaled",
    training_frame=train_final
)

# 查看模型性能
print(aml.leaderboard)

# 评估时需要分别处理每个测试集
def evaluate_model(model, test_set, test_name):
    pred_scaled = model.predict(test_set)
    
    true_scaled = test_set["RFU_scaled"].as_data_frame(use_multi_thread=True).values.flatten()
    pred_scaled = pred_scaled.as_data_frame(use_multi_thread=True).values.flatten()
    
    # 计算指标（标准化空间）
    print(f"\n{test_name} 评估指标:")
    print(f"RMSE: {root_mean_squared_error(true_scaled, pred_scaled):.4f}")
    print(f"MSE: {mean_squared_error(true_scaled, pred_scaled):.4f}")
    print(f"MAE: {mean_absolute_error(true_scaled, pred_scaled):.4f}")
    print(f"R²: {r2_score(true_scaled, pred_scaled):.4f}")

# 对每个测试集进行评估
evaluate_model(aml.leader, test_bacteria_final, "Bacteria Test")
evaluate_model(aml.leader, test_LysM_final, "LysM Test") 
evaluate_model(aml.leader, test_toxin_final, "Toxin Test")
evaluate_model(aml.leader, test_virus_final, "Virus Test")

# 释放资源
h2o.cluster().shutdown()