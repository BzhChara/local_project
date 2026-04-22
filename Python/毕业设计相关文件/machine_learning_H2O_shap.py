import h2o
import shap
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import joblib

# 初始化H2O
h2o.init(max_mem_size='16G')

# 加载模型
model_path = "./saved_models/StackedEnsemble_BestOfFamily_2_AutoML_1_20250602_70214"  # 替换为实际路径
best_model = h2o.load_model(model_path)

# 加载测试数据
test_path = "./data/test_data.csv"
test_final = h2o.import_file(test_path)

# 加载背景帧
bg_path = "./data/background_frame.csv"
background_frame = h2o.import_file(bg_path)

# 加载列信息
column_info = joblib.load("./data/column_info.pkl")
protein_cols = column_info["protein_cols"]
oligo_cols = column_info["oligo_cols"]

# 选择测试样本
test_sample = test_final[:2000] 

try:
    # 尝试使用背景帧计算贡献值
    contribs = best_model.predict_contributions(test_sample, background_frame=background_frame)
except Exception as e:
    print(f"使用背景帧计算贡献值时出错: {e}")
    try:
        # 如果失败，尝试不使用背景帧
        contribs = best_model.predict_contributions(test_sample)
    except Exception as e2:
        print(f"无法计算SHAP贡献值: {e2}")
        print("跳过SHAP分析部分")
        h2o.shutdown(prompt=False)
        exit()

# 处理SHAP贡献值
shap_df = contribs.as_data_frame(use_multi_thread=True)
shap_values = shap_df.iloc[:, :-1].values
feature_names_used = shap_df.columns[:-1].tolist()

# 1. 绘制SHAP摘要图
shap_explainer = best_model.explain(test_sample, background_frame=background_frame)
shap_summary = shap_explainer.visualize("summary")
shap_summary.savefig('shap_summary.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. 个体特征贡献分析
# 计算平均绝对SHAP值
mean_abs_shap = np.abs(shap_values).mean(axis=0)

# 确定特征类型
protein_set = set(protein_cols)
oligo_set = set(oligo_cols)
feature_types = [
    'Protein' if f in protein_set else 
    'Oligosaccharide' if f in oligo_set else 
    'Unknown'
    for f in feature_names_used
]

# 创建特征重要性DataFrame
feature_importance = pd.DataFrame({
    'Feature': feature_names_used,
    'Mean|SHAP|': mean_abs_shap,
    'Type': feature_types
})

# 获取前20个最重要的特征
top_20_features = feature_importance.sort_values('Mean|SHAP|', ascending=False).head(20)

# 统计前20中蛋白质和寡糖特征的数量
protein_count = top_20_features[top_20_features['Type'] == 'Protein'].shape[0]
oligo_count = top_20_features[top_20_features['Type'] == 'Oligosaccharide'].shape[0]

# 打印统计结果
print(f"\n前20位关键驱动因子统计:")
print(f"蛋白质特征数量: {protein_count}")
print(f"寡糖特征数量: {oligo_count}")

# 可视化前20个特征
plt.figure(figsize=(12, 8))
sns.barplot(x='Mean|SHAP|', y='Feature', hue='Type', data=top_20_features,
            palette={'Protein': '#1f77b4', 'Oligosaccharide': '#2ca02c'})
plt.title('Top 20 Features by Mean Absolute SHAP Value', fontsize=16, fontweight='bold')
plt.xlabel('Mean(|SHAP Value|)', fontsize=14)
plt.ylabel('Feature', fontsize=14)
plt.xticks(fontsize=12)  
plt.yticks(fontsize=12)  
plt.tight_layout()
plt.savefig('top_20_features.png', dpi=300)
plt.close()

# 3. 全局特征集贡献分析
protein_mask = [f in protein_set for f in feature_names_used]
oligo_mask = [f in oligo_set for f in feature_names_used]

protein_shap_sum = np.abs(shap_values[:, protein_mask]).sum()
oligo_shap_sum = np.abs(shap_values[:, oligo_mask]).sum()
total_shap_sum = protein_shap_sum + oligo_shap_sum

protein_contribution = protein_shap_sum / total_shap_sum * 100
oligo_contribution = oligo_shap_sum / total_shap_sum * 100

print("\n全局特征集贡献分析:")
print(f"蛋白质特征总贡献: {protein_contribution:.2f}%")
print(f"寡糖特征总贡献: {oligo_contribution:.2f}%")

plt.figure(figsize=(12, 8))
wedges, text_labels, autotexts = plt.pie([protein_contribution, oligo_contribution], 
        labels=['Protein Features', 'Oligosaccharide Features'],
        colors=['#1f77b4', '#2ca02c'], autopct='%1.1f%%',
        startangle=90, explode=(0.05, 0), shadow=False,textprops={'fontsize': 14})

# 设置自动百分比文本的字体大小和样式
for autotext in autotexts:
    autotext.set_fontsize(12)         # 设置百分比文本字体大小
    autotext.set_fontweight('bold')   # 设置文本加粗

plt.title('Global Contribution to Prediction', fontsize=16, fontweight='bold')
plt.savefig('global_contribution.png', dpi=300)
plt.close()

# 释放资源
h2o.cluster().shutdown()