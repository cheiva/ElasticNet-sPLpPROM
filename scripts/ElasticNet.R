library(tidyverse)
# 采用Elastic Net
library(caret)
library(caret)
library(glmnet)
library(pROC)
# 混淆矩阵绘制
library(ggplot2)
library(caret)
library(dplyr)
library(tidyr)

### 数据预处理 ###
#注：原始数据受版权保护，不可外传，已替换为模拟数据
metabolomics_path <- "data/raw/代谢组学.xlsx"  
sample_info_path <- "data/raw/代谢组学样本说明.xlsx"

# 数据读取与整合
metabolomics <- t(readxl::read_excel(metabolomics_path))
colnames(metabolomics) <- metabolomics[1, ]
metabolomics <- as.data.frame(metabolomics[-1, ]) %>% 
  mutate(Sample = rownames(.))

sample_info <- readxl::read_excel(sample_info_path, sheet = 1, skip = 3) %>% 
  select(Sample = '样本名称*', 
         Group = '备注') %>%
  mutate(Group = as.character(Group))

combined <- metabolomics %>% 
  left_join(sample_info, by = "Sample")

# 创建分类标签
combined <- combined %>% 
  mutate(Class = case_when(
    grepl("^ctrl", Sample) ~ "normal",
    grepl("^pPROM", Sample) ~ "pPROM",
    grepl("^sPTB", Sample) ~ "sPL",
    TRUE ~ NA_character_
  ))  %>%
  filter(!is.na(Class)) 

# 划分数据标签
combined <- combined %>% 
  mutate(DataType = ifelse(grepl("^ABB", Group), "test", "train")) %>% 
  filter(DataType %in% c("train", "test"))

# 转换数据类型
combined[, 1:(ncol(combined)-4)] <- apply(combined[,  1:(ncol(combined)-4)], 2, as.numeric)
combined$Class <- factor(combined$Class, levels = c("normal", "sPL", "pPROM"))

### 划分训练/测试集 ### 
train_data <- combined %>% filter(DataType == "train") %>% select(-Sample, -Group, -DataType)
test_data <- combined %>% filter(DataType == "test") %>% select(-Sample, -Group, -DataType)

# 处理缺失值
train_medians <- apply(train_data[, -ncol(train_data)], 2, median, na.rm = TRUE)
train_data[, -ncol(train_data)] <- sapply(1:(ncol(train_data)-1), function(i) {
  ifelse(is.na(train_data[, i]), train_medians[i], train_data[, i])
})

test_data[, -ncol(test_data)] <- sapply(1:(ncol(test_data)-1), function(i) {
  ifelse(is.na(test_data[, i]), train_medians[i], test_data[, i])
})

# 将初步处理的(未经模型特化标准化)数据保存备用，以对比其他模型效果
write.csv(test_data, "data/processed/test_data.csv")
write.csv(train_data, "data/processed/train_data.csv")

### 数据标准化 ###

train_classes <- train_data$Class
test_classes <- test_data$Class

train_features <- as.matrix(train_data[, -ncol(train_data)])
test_features <- as.matrix(test_data[, -ncol(test_data)])

scaler <- preProcess(train_features, method = c("center", "scale"))
train_scaled <- predict(scaler, train_features)
test_scaled <- predict(scaler, test_features)

### 特征预筛选（尽可能降维）###
var_filter <- nearZeroVar(train_scaled, freqCut = 95/5)
if(length(var_filter) > 0) {
  train_scaled <- train_scaled[, -var_filter]
  test_scaled <- test_scaled[, -var_filter]
}

### Elastic Net模型训练 ###
# 设置分层交叉验证
set.seed(123)
folds <- createMultiFolds(
  y = train_classes,
  k = 5,               
  times = 5            # 5折交叉验证 5次重复（比默认下调，减少过拟合情况）
)

ctrl <- trainControl(
  method = "repeatedcv",
  index = folds,        # 使用预生成的分层抽样折
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  selectionFunction = "best"
)

elastic_grid <- expand.grid(
  alpha = seq(0, 1, length = 11),
  lambda = 10^seq(-3, 2, length=50)  # 扩展lambda范围至更强的正则化
)

# Elastic Net模型训练
elnet_model <- train(
  x = train_scaled,
  y = train_classes,
  method = "glmnet",
  family = "multinomial",
  tuneGrid = elastic_grid,
  trControl = ctrl,
  metric = "Kappa"
)
### 模型评估 ###
# 训练集性能
train_pred <- predict(elnet_model, train_scaled)
confusionMatrix(train_pred, train_classes)

# 测试集预测
test_pred <- predict(elnet_model, test_scaled)
test_prob <- predict(elnet_model, test_scaled, type = "prob")

# 测试集评估
test_results <- confusionMatrix(test_pred, test_classes)
print(test_results)

# 多类别ROC曲线
roc_list <- lapply(levels(train_classes), function(cls){
  roc(response = as.numeric(test_classes == cls),
      predictor = test_prob[[cls]])
})
names(roc_list) <- levels(train_classes)

### 结果可视化 ###
plot_confusion_matrix <- function(cm, title_suffix = "", palette = "Greens") {

  cm_df <- as.data.frame(cm$table)

  total <- sum(cm_df$Freq)
  cm_df <- cm_df %>%
    mutate(Percent = round(Freq/total * 100, 1),
           Label = paste0(Freq, "\n(", Percent, "%)"))

  accuracy <- round(cm$overall["Accuracy"], 3)
  
  # 热力图
  ggplot(cm_df, aes(Reference, Prediction, fill = Freq)) +
    geom_tile(alpha = 0.8, color = "white") +
    geom_text(aes(label = Label), color = "black", size = 4) +
    scale_fill_distiller(palette = palette, direction = 1) +
    labs(title = paste("Confusion Matrix", title_suffix),
         subtitle = paste("Overall Accuracy =", accuracy),
         x = "Actual Class",
         y = "Predicted Class") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    guides(fill = guide_colorbar(title = "Count"))
}
plot_confusion_matrix(test_results, "(Test Set)", "Oranges")

# 变量重要性
var_imp <- varImp(elnet_model)$importance
var_imp <- var_imp[order(-rowMeans(var_imp)), ]
significant_features <- head(var_imp, 19)  # 只取了前19个

# 特征重要性
ggplot(data.frame(Feature=rownames(significant_features), 
                  Importance=rowMeans(significant_features)),
       aes(x=reorder(Feature, Importance), y=Importance)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() + 
  labs(title="Top 19 Important Features", y="Importance", x="Metabolites")

# 绘制多类别ROC曲线
colors <- rainbow(length(roc_list))  
plot(roc_list[[1]], col = colors[1], main = "ROC Curves")
for(i in 2:length(roc_list)) {
  lines(roc_list[[i]], col = colors[i])
}
legend("bottomright", legend = names(roc_list),
       col = colors, lwd = 2, bty = "n") 

### 模型保存与应用 ###
saveRDS(elnet_model, "results/model/preterm_predict_model.rds")

top_features <- rownames(significant_features)
write.csv(data.frame(Metabolites = top_features), 
          "results/figures/important_features.csv")
