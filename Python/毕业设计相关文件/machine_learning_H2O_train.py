import h2o
from h2o.automl import H2OAutoML
from sklearn.preprocessing import QuantileTransformer
import joblib
import os

# 初始化H2O
h2o.init(max_mem_size='16G')

# 加载数据
data = h2o.import_file("combined_encode_esm_no-z-score.csv", header=0)

# 定义目标变量列名
protein_cols = data.columns[:1280]
oligo_cols = data.columns[1280:1280+342]
y_col = data.columns[-1]

# 数据分割
train, test = data.split_frame(ratios=[0.8], seed=42)

# 使用分位数变换
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

"""# 训练AutoML
aml = H2OAutoML(max_runtime_secs=3600, 
                nfolds=5, 
                seed=42, 
                stopping_metric="RMSE",
                sort_metric="RMSE",
                keep_cross_validation_predictions=True,
                keep_cross_validation_models=True,
                project_name='rfu_prediction')

aml.train(
    x=protein_cols + oligo_cols,
    y="RFU_scaled",
    training_frame=train_final
)

# 保存最佳模型
best_model = aml.leader
model_path = h2o.save_model(best_model, path="./saved_models", force=True)
print(f"模型已保存至: {model_path}")
"""
# 保存测试集
test_path = "./data/test_data.csv"
h2o.export_file(test_final, test_path)
print(f"测试集已保存至: {test_path}")

# 保存背景帧（从训练集中抽取100个样本）
if train_final.nrow > 100:
    safe_ratio = min(100 / train_final.nrow, 0.99)
    splits = train_final.split_frame(ratios=[safe_ratio], seed=42)
    background_frame = splits[0]
    bg_path = "./data/background_frame.csv"
    h2o.export_file(background_frame, bg_path)
    print(f"背景帧已保存至: {bg_path}")

# 保存列名信息
import joblib
column_info = {
    "protein_cols": protein_cols,
    "oligo_cols": oligo_cols,
    "y_col": y_col
}
joblib.dump(column_info, "./data/column_info.pkl")
print("列信息已保存")

# 关闭H2O
h2o.shutdown(prompt=False)