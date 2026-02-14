
setwd("~/bioinformatics/ML_mtb/")


dst_crypritc_HIGH <- read.csv(file = "./data/raw/cryptic_phy_HIGH.csv")
drugs <- c("EMB", "INH","PZA","RIF","SM", "AMK","CM",
           "ETO","KM","LFX","MFX","OFX","PAS")
drug_compare <- c("EMB","INH","RIF","AMK","ETO","KM","LFX","MFX")
# drug_compare <- drug_compare[1:2] #debug
path_predict_cryptic <- paste0("./result/cryptic_project/", dst_crypritc_HIGH$sraid,".csv")


#### 下面是用模型预测 Cryptic project
model <- paste0(drugs,".rds")
model_show <- lapply(model, function(x){
  ak <- readRDS(paste0("./result/models/",x))
  do.call('rbind',lapply(ak,function(x){
    return(x$metrics)
  }))   
})
names(model_show) <- model
model_show

# ak <- readRDS("../FS_Mtb_2025/model/09_LFX_C1.rds")


if(FALSE){
  start_time <- Sys.time()
  filename <- paste0(dst_crypritc_HIGH$sraid,".vcf.gz")
  result_pred <- mclapply(seq_along(filename), function(i_f){
    cat('Sample:',i_f,"\n")
    filename_tmp <- filename[i_f]
    path_tmp <- paste0("../FS_Mtb_cryptic/data/vcf_files/",filename_tmp)
    file_tmp <- read.vcfR(path_tmp)
    pred_drug <- lapply(seq_along(drug_compare), function(i_d){
      drug_tmp <- drug_compare[i_d] 
      cat("Drug:",drug_tmp, "\n")
      path_model_drug <- paste0("./result/models/",drug_tmp,".rds")
      model_drug <- readRDS(path_model_drug)
      result_modelS <- lapply(1:5, function(i_m){
        model_tmp <- model_drug[[i_m]]
        model_fit <- model_tmp$fit
        if(is.null(model_fit)){
          result_return <- data.frame(drug = drug_tmp,
                                      prob_R = NA,prob_S = NA,
                                      Model = paste0("Model_S", model_tmp$stable))
          return(result_return)
        }
        
        feature_model_drug <- rownames(model_fit$importance)
        no <- match( (feature_model_drug %>% sub("rs","",.)),file_tmp@fix[,2])
        feature_drug_tmp <-  as.character(file_tmp@gt[no,2])
        feature_drug_tmp <- ifelse(grepl("1/1",feature_drug_tmp),1,0)
        names(feature_drug_tmp) <- feature_model_drug
        pred_drug_tmp <- predict(model_fit,as.data.frame(t(feature_drug_tmp)),type = "prob")
        result_return <- data.frame(pred_drug_tmp)
        colnames(result_return) <- c("prob_R", "prob_S")
        result_return <- data.frame(drug = drug_tmp,result_return,
                                    Model = paste0("model_S",model_tmp$stable))
        return(result_return)
      })
      result_modelS <- do.call('rbind',result_modelS)
      return(result_modelS)
    })
    pred_drug <- do.call('rbind',pred_drug)
    true_label <- dst_crypritc_HIGH[i_f,match(drug_compare,colnames(dst_crypritc_HIGH))]
    true_label <- as.character(true_label)
    pred_drug$true_label <- rep(true_label,each = 5)
    pred_drug$sra_id <- sub(".vcf.gz","", filename[i_f])
    pred_drug$fold <- NA
    # "sra_id","drug","fold","true_label","prob_R","prob_S" "Model"
    pred_drug <- pred_drug[,c(6,1,7,5,2:4)]
    write.csv(pred_drug,
              paste0("./result/cryptic_project/",sub(".vcf.gz", ".csv", filename_tmp)),
              row.names = F)
    return(pred_drug)
  }, mc.cores = 6)
  
  end_time <- Sys.time()
  
  # 计算时间差
  time_taken <- end_time - start_time
  print(time_taken)
  # Time difference of 1.627693 hours
}


calculate_metrics <- function(data, true_col = "true_label", 
                              prob_pos_col = "prob_R", 
                              prob_neg_col = "prob_S",
                              pos_class = "R",
                              threshold = 0.5) {
  
  # 加载必要的包
  if(!require(pROC)) {
    stop("需要安装 pROC 包: install.packages('pROC')")
  }
  if(!require(PRROC)) {
    stop("需要安装 PRROC 包: install.packages('PRROC')")
  }
  
  # 获取真实标签和预测概率
  actual <- factor(data[[true_col]], levels = c("R", "S"))
  prob_R <- data[[prob_pos_col]]
  prob_S <- data[[prob_neg_col]]
  
  # 根据阈值生成预测标签
  predictions <- factor(ifelse(prob_R >= threshold, "R", "S"), 
                        levels = c("R", "S"))
  
  # 构建混淆矩阵
  confusion <- table(Predicted = predictions, Actual = actual)
  
  # 提取 TP, FN, FP, TN (R为正类)
  if(all(c("R", "S") %in% rownames(confusion)) && 
     all(c("R", "S") %in% colnames(confusion))) {
    TP <- as.numeric(confusion["R", "R"])  # 预测R,实际R
    FN <- as.numeric(confusion["S", "R"])  # 预测S,实际R
    FP <- as.numeric(confusion["R", "S"])  # 预测R,实际S
    TN <- as.numeric(confusion["S", "S"])  # 预测S,实际S
  } else {
    # 处理缺失类别的情况
    TP <- ifelse("R" %in% rownames(confusion) && "R" %in% colnames(confusion),
                 as.numeric(confusion["R", "R"]), 0)
    FN <- ifelse("S" %in% rownames(confusion) && "R" %in% colnames(confusion),
                 as.numeric(confusion["S", "R"]), 0)
    FP <- ifelse("R" %in% rownames(confusion) && "S" %in% colnames(confusion),
                 as.numeric(confusion["R", "S"]), 0)
    TN <- ifelse("S" %in% rownames(confusion) && "S" %in% colnames(confusion),
                 as.numeric(confusion["S", "S"]), 0)
  }
  
  # 处理可能的NA值
  TP <- ifelse(is.na(TP), 0, TP)
  FN <- ifelse(is.na(FN), 0, FN)
  FP <- ifelse(is.na(FP), 0, FP)
  TN <- ifelse(is.na(TN), 0, TN)
  
  # 计算基本指标(添加除零保护)
  Recall <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  Precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  F1_score <- ifelse((Precision + Recall) == 0, 0, 
                     2 * (Precision * Recall) / (Precision + Recall))
  
  # Matthews Correlation Coefficient (MCC)
  MCC_numerator <- (TP * TN) - (FP * FN)
  term1 <- as.numeric(TP + FP)
  term2 <- as.numeric(TP + FN)
  term3 <- as.numeric(TN + FP)
  term4 <- as.numeric(TN + FN)
  MCC_denominator <- sqrt(term1 * term2 * term3 * term4)
  MCC <- ifelse(MCC_denominator == 0, 0, MCC_numerator / MCC_denominator)
  
  # 计算 AUC (ROC曲线下面积)
  # 将真实标签转换为数值: R=1, S=0
  actual_numeric <- ifelse(actual == "R", 1, 0)
  
  tryCatch({
    roc_obj <- roc(actual_numeric, prob_R, 
                   levels = c(0, 1), 
                   direction = "<",
                   quiet = TRUE)
    AUC <- as.numeric(auc(roc_obj))
  }, error = function(e) {
    warning("AUC计算失败: ", e$message)
    AUC <<- NA
  })
  
  # 计算 AUPRC (Precision-Recall曲线下面积)
  tryCatch({
    # 使用PRROC包计算AUPRC
    # 正类(R)的样本索引
    pos_indices <- which(actual == "R")
    neg_indices <- which(actual == "S")
    
    if(length(pos_indices) > 0 && length(neg_indices) > 0) {
      pr_obj <- pr.curve(scores.class0 = prob_R[pos_indices],
                         scores.class1 = prob_R[neg_indices],
                         curve = TRUE)
      AUPRC <- pr_obj$auc.integral
    } else {
      AUPRC <- NA
    }
  }, error = function(e) {
    warning("AUPRC计算失败: ", e$message)
    AUPRC <<- NA
  })
  
  # 返回结果
  results <- list(
    confusion_matrix = confusion,
    metrics = data.frame(
      TN = TN,
      FP = FP,
      FN = FN,
      TP = TP,
      Recall = round(Recall, 4),
      Precision = round(Precision, 4),
      Accuracy = round(Accuracy, 4),
      F1_score = round(F1_score, 4),
      MCC = round(MCC, 4),
      AUC = round(AUC, 4),
      AUPRC = round(AUPRC, 4),
      row.names = NULL
    ),
    roc_data = if(exists("roc_obj")) roc_obj else NULL,
    pr_data = if(exists("pr_obj")) pr_obj else NULL
  )
  
  return(results)
}





result_cryptic <- do.call('rbind',lapply(path_predict_cryptic,read.csv))

result_cryptic$Model <- gsub("^Model", "model", result_cryptic$Model)

set.seed(123)
# cvresult_drug_cryptic <- lapply(seq_along(drug_compare), function(i_d){
#   drug_tmp <- drug_compare[i_d]
#   cat("Drug:",drug_tmp,"\n")
#   col_type <- which(colnames(dst_crypritc_HIGH) %in% drug_tmp)
#   trainset <- data.frame(sra_id = dst_crypritc_HIGH$sraid, type = dst_crypritc_HIGH[,col_type])
#   
#   trainset$type <-factor(trainset$type, levels = c("R","S"))
#   no_D <- which(is.na(trainset$type))
#   trainset_omit <- trainset[-no_D,]
#   row.names(trainset_omit) <- NULL
#   folds <- createFolds(trainset_omit$sra_id, k = 5, list = TRUE, returnTrain = FALSE)
#   result_predict_tmp_drugomit <- result_predict_cryptic[which(result_predict_cryptic$drug == drug_tmp),]
#   
#   # 2. 执行5折交叉验证
#   cv_results <- lapply(seq_along(folds), function(i_f){
#     cat("\n========== Fold", i_f, "==========\n")
#     test_idx <- folds[[i_f]]
#     test_data <- trainset_omit[test_idx, ]
#     result <- result_predict_tmp_drugomit[match(test_data$sra_id, 
#                                                 result_predict_tmp_drugomit$sra_id),]
#     result$fold <- i_f
#     return(result)
#   })
#   # 2. 合并所有fold的结果
#   cv_results_df <- do.call(rbind, cv_results)
# })




metrics_cryptic <- lapply(seq_along(drug_compare), function(i_d){
  drug_tmp <- drug_compare[i_d]
  cat('Drug:',drug_tmp,"\n")
  result_sub <- result_cryptic[which(result_cryptic$drug == drug_tmp),]
  result_sub_split <- split.data.frame(result_sub,result_sub$Model)
  result_models <- lapply(seq_along(result_sub_split), function(i_split){
    set_tmp <- result_sub_split[[i_split]]
    result_tmp <- calculate_metrics(
      data = set_tmp,
      true_col = "true_label",
      prob_pos_col = "prob_R",
      prob_neg_col = "prob_S",
      pos_class = "R",
      threshold = 0.5
    )
    result_re <- result_tmp$metrics
    result_re$Model <- set_tmp$Model[1]
    return(result_re)
  })
  result_models <- do.call('rbind',result_models)
  result_models$Drug <- drug_tmp
  return(result_models)
})
metrics_cryptic <- do.call('rbind',metrics_cryptic)
metrics_cryptic <- split.data.frame(metrics_cryptic,metrics_cryptic$Model)
write.csv(metrics_cryptic[[3]],"./output/metric_test.csv", row.names = F)


##### 画图 cryptic project ######

# 雷达图
plot_radar_sector <- function(data, metric_name, fill_color, title) {
  library(dplyr)
  library(ggplot2)
  
  # 准备数据
  plot_data <- data %>%
    select(Drug, value = all_of(metric_name)) %>%
    mutate(
      Drug = factor(Drug, levels = Drug),
      id = row_number()
    )
  
  # 为每个扇形创建数据
  n_drugs <- nrow(plot_data)
  sector_data <- data.frame()
  
  for(i in 1:n_drugs) {
    # 每个扇形的4个顶点
    sector <- data.frame(
      drug_id = i,
      drug_name = plot_data$Drug[i],
      x = c(i - 0.5, i - 0.5, i + 0.5, i + 0.5),
      y = c(0, plot_data$value[i], plot_data$value[i], 0)
    )
    sector_data <- rbind(sector_data, sector)
  }
  
  # 绘图
  ggplot() +
    # 背景网格圆圈
    geom_hline(yintercept = seq(0.2, 1, 0.2), 
               color = "grey85", linewidth = 0.3) +
    # 扇形填充
    geom_polygon(data = sector_data,
                 aes(x = x, y = y, group = drug_id),
                 fill = fill_color,
                 alpha = 0.5,
                 color = fill_color,
                 linewidth = 0.2) +
    # 数值标签
    geom_text(data = plot_data,
              aes(x = id, y = value, label = sprintf("%.2f", value)),
              size = 3,
              hjust = 0.5,
              vjust = -0.5,
              color = "red",
              fontface = "bold") +
    # 药物标签
    scale_x_continuous(
      breaks = 1:n_drugs,
      labels = plot_data$Drug
    ) +
    # 极坐标
    coord_polar() +
    # Y轴范围
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    # 标题
    labs(title = title, x = NULL, y = NULL) +
    # 主题
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold", color = "grey30"),
      axis.text.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

data_rada_cryptic <- read.csv("./output/metric_test.csv")

p1 <- plot_radar_sector(data_rada_cryptic, "Accuracy", "#5B9BD5", "(a) Acc")
p2 <- plot_radar_sector(data_rada_cryptic, "MCC", "#FFC000", "(b) MCC")  # 改为MCC，使用黄色
p3 <- plot_radar_sector(data_rada_cryptic, "AUC", "#ED7D31", "(c) AUC")
p4 <- plot_radar_sector(data_rada_cryptic, "AUPRC", "#70AD47", "(d) AUPRC")


library(cowplot)

# 雷达图组合
metric_combined_cryptic <- plot_grid(p1, p2, p3,p4, ncol = 4, align = "h")


data_roc_cryptic <- result_cryptic[which(result_cryptic$Model == "model_S3"),]

f_plot_multi_drug_roc <- function(data, drug_list, title = "ROC Curves") {
  
  library(ggplot2)
  library(pROC)
  library(dplyr)
  
  # 颜色方案
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
              "#A65628", "#F781BF", "#999999")
  
  # 初始化图形
  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray50", linewidth = 0.5)
  
  # 存储 AUC 值
  auc_values <- numeric(length(drug_list))
  
  # 为每个药物绘制 ROC 曲线
  for(i in seq_along(drug_list)) {
    drug_name <- drug_list[i]
    
    # 筛选当前药物的数据
    drug_data <- data %>% filter(drug == drug_name)
    
    # 计算 ROC
    roc_obj <- roc(drug_data$true_label, drug_data$prob_S,
                   levels = c("R", "S"), quiet = TRUE)
    
    auc_values[i] <- as.numeric(auc(roc_obj))
    
    # 准备绘图数据
    roc_df <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Drug = drug_name
    )
    
    # 添加曲线
    p <- p + geom_path(data = roc_df,
                       aes(x = FPR, y = TPR),
                       color = colors[i],
                       linewidth = 1,
                       alpha = 0.8)
    
    cat(drug_name, "- AUC:", sprintf("%.4f", auc_values[i]), "\n")
  }
  
  # 图例标签
  legend_labels <- paste0(drug_list, " (AUC = ", 
                          sprintf("%.3f", auc_values), ")")
  
  # 计算图例位置
  n_drugs <- length(drug_list)
  legend_y <- seq(0.35, 0.05, length.out = n_drugs)
  
  # 添加图例
  for(i in 1:n_drugs) {
    # 图例线段
    p <- p + annotate("segment",
                      x = 0.55, xend = 0.62,
                      y = legend_y[i], yend = legend_y[i],
                      color = colors[i],
                      linewidth = 1,
                      alpha = 0.8)
    
    # 图例文字
    p <- p + annotate("text",
                      x = 0.64, y = legend_y[i],
                      label = legend_labels[i],
                      hjust = 0, size = 3.5)
  }
  
  # 完善图形 - 修复 margin 问题 ⬇️
  p <- p +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "False Positive Rate",
         y = "True Positive Rate",
         title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)  # 明确使用 ggplot2::margin
    )
  
  return(p)
}

f_plot_multi_drug_pr <- function(data, drug_list, title = "PR Curves") {
  
  library(ggplot2)
  library(PRROC)
  library(dplyr)
  
  # 颜色方案
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
              "#A65628", "#F781BF", "#999999")
  
  # 存储所有数据
  all_data <- list()
  auprc_values <- numeric(length(drug_list))
  
  # 计算每个药物的 PR 曲线
  for(i in seq_along(drug_list)) {
    drug_name <- drug_list[i]
    drug_data <- data %>% filter(drug == drug_name)
    
    pos_idx <- which(drug_data$true_label == "R")
    neg_idx <- which(drug_data$true_label == "S")
    
    if(length(pos_idx) > 0 && length(neg_idx) > 0) {
      pr_obj <- pr.curve(scores.class0 = drug_data$prob_R[pos_idx],
                         scores.class1 = drug_data$prob_R[neg_idx],
                         curve = TRUE)
      
      auprc_values[i] <- pr_obj$auc.integral
      
      all_data[[i]] <- data.frame(
        Recall = pr_obj$curve[, 1],
        Precision = pr_obj$curve[, 2],
        Drug = paste0(drug_name, " (AUPRC = ", sprintf("%.3f", auprc_values[i]), ")")
      )
      
      cat(drug_name, "- AUPRC:", sprintf("%.4f", auprc_values[i]), "\n")
    }
  }
  
  # 合并所有数据
  plot_data <- do.call(rbind, all_data)
  plot_data$Drug <- factor(plot_data$Drug, levels = unique(plot_data$Drug))
  
  # 绘图
  p <- ggplot(plot_data, aes(x = Recall, y = Precision, color = Drug)) +
    geom_path(linewidth = 1, alpha = 0.8) +
    scale_color_manual(values = colors[1:length(drug_list)]) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Recall (Sensitivity)",
         y = "Precision (PPV)",
         title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      # 图例在底部,单列垂直排列
      legend.position = "bottom",
      legend.direction = "vertical",  # 必须添加这个
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      legend.box = "vertical",  # 添加这个
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2),
                                ncol = 1))  # 单列
  
  return(p)
}


library(dplyr)
library(ggplot2)
library(patchwork)
# 定义药物分组
first_line_drugs <- c("EMB", "INH", "RIF")
second_line_drugs <- c("AMK", "ETO", "KM", "LFX", "MFX")
cat("\n一线药物:\n")
p_roc_1st <- f_plot_multi_drug_roc(
  data_roc_cryptic, 
  first_line_drugs, 
  "(e) ROC Curves - First-line Drugs"
)

p_roc_2nd <- f_plot_multi_drug_roc(
  data_roc_cryptic, 
  second_line_drugs, 
  "(f) ROC Curves - Second-line Drugs"
)


roc_combined_cryptic <- plot_grid(p_roc_1st, p_roc_2nd, ncol = 2, align = "h")


cat("\n=== 绘制 PR 曲线 ===\n")



p_pr_1st <- f_plot_multi_drug_pr(
  data_roc_cryptic, 
  first_line_drugs, 
  "(g) PR Curves - First-line Drugs"
);p_pr_1st

p_pr_2nd <- f_plot_multi_drug_pr(
  data_roc_cryptic, 
  second_line_drugs, 
  "(h) PR Curves - Second-line Drugs"
)
pr_combined_cryptic <- plot_grid(p_pr_1st, p_pr_2nd, ncol = 2, align = "h")
pr_combined_cryptic

plot_combined_cryptic <- plot_grid(metric_combined_cryptic,roc_combined_cryptic,
                                   pr_combined_cryptic, ncol = 1, align = "h")

ggsave("./output/metric_test.pdf", plot_combined_cryptic, width = 16, height = 16)

library(patchwork)

# 设置布局高度，考虑 PR 图的图例
plot_combined_cryptic <- 
  metric_combined_cryptic / 
  roc_combined_cryptic / 
  pr_combined_cryptic +
  plot_layout(heights = c(1, 1, 1.2685))  # PR 图需要更多空间

ggsave("./output/metric_test.pdf", plot_combined_cryptic, 
       width = 16, height = 16)




########### 下面是用模型预测 tbportal



dst <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_DST_sample_clean.csv")
X_genomic <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_X_genomic_005.csv")
drug_list <- function(drug_name) {
  # 创建数据结构（通常这部分会放在函数外部）
  drug_data <- list(
    df = data.frame(
      Abbrev = c("EMB","INH","PZA","RIF","SM",
                 "AMK","CM","ETO","KM","LFX",
                 "MFX","OFX","PAS"),
      Drug = c("ethambutol", "isoniazid", "pyrazinamide",
               "rifampicin", "streptomycin", "amikacin",
               "capreomycin", "ethionamide", "kanamycin",
               "levofloxacin", "moxifloxacin", "ofloxacin",
               "p_aminosalicylic_acid"),
      stringsAsFactors = FALSE
    )
  )
  
  drug_data$lookup <- c(
    setNames(drug_data$df$Drug, drug_data$df$Abbrev),
    setNames(drug_data$df$Abbrev, drug_data$df$Drug)
  )
  
  # 直接使用变量
  result <- drug_data$lookup[drug_name]
  names(result) <- NULL
  return(result)
}
set.seed(123)
cvresult_drug <- lapply(seq_along(drugs), function(i_d){
  drug_tmp <- drugs[i_d]
  model_rds <- readRDS(paste0("./result/models/",
                              drug_tmp,".rds"))
  feature_fit <- model_rds[[3]]$fit
  feature_tmp <- rownames(feature_fit$importance)
  col_feature <- which( colnames(X_genomic) %in%  feature_tmp)
  col_type <- which(colnames(dst) %in% drug_list(drug_tmp))
  trainset <-data.frame(X_genomic, type = factor(dst[,col_type], levels = c("R","S")))
  no_D <- which(is.na(trainset$type))
  trainset <- trainset[-no_D, ]
  rownames(trainset) <- NULL
  folds <- createFolds(trainset$type, k = 5, list = TRUE, returnTrain = FALSE)
  
  # 2. 执行5折交叉验证
  cv_results <- lapply(seq_along(folds), function(i_f){
    cat("\n========== Fold", i_f, "==========\n")
    test_idx <- folds[[i_f]]
    train_idx <- setdiff(1:nrow(trainset), test_idx)
    train_data <- trainset[train_idx, ]
    test_data <- trainset[test_idx, ]
    # 预测概率
    pred_prob <- predict(feature_fit , newdata = test_data, 
                         type = "prob")
    # 保存结果(包含真实标签)
    result <- data.frame(
      sra_id = test_data$sra_id,
      drug = drug_tmp,
      fold = i_f,
      true_label = test_data$type,  # 真实标签
      prob_R = pred_prob[, "R"],
      prob_S = pred_prob[, "S"]
    )
    
    return(result)
  })
  # 2. 合并所有fold的结果
  cv_results_df <- do.call(rbind, cv_results)
})

cvresult_drug <- do.call('rbind',cvresult_drug)

metrics_tbportal <- lapply(seq_along(drugs), function(i_d){
  drug_tmp <- drugs[i_d]
  cat('Drug:',drug_tmp,"\n")
  result_sub <- cvresult_drug[which(cvresult_drug$drug == drug_tmp),]
  result_tmp <- calculate_metrics(
    data = result_sub,
    true_col = "true_label",
    prob_pos_col = "prob_R",
    prob_neg_col = "prob_S",
    pos_class = "R",
    threshold = 0.5
  )
  result_re <- result_tmp$metrics
  result_re$Drug <- drug_tmp
  return(result_re)
})
metrics_tbportal <- do.call("rbind",metrics_tbportal)
write.csv(metrics_tbportal,"./output/metric_train.csv", row.names = F)

#### 画图。tbportal的雷达图
data_rada_tbportal <- read.csv("./output/metric_train.csv")
p1 <- plot_radar_sector(data_rada_tbportal, "Accuracy", "#5B9BD5", "(a) Acc")
p2 <- plot_radar_sector(data_rada_tbportal, "MCC", "#FFC000", "(b) MCC")  # 改为MCC，使用黄色
p3 <- plot_radar_sector(data_rada_tbportal, "AUC", "#ED7D31", "(c) AUC")
p4 <- plot_radar_sector(data_rada_tbportal, "AUPRC", "#70AD47", "(d) AUPRC")


library(cowplot)

# 雷达图组合
metric_combined_tbportal <- plot_grid(p1, p2, p3,p4, ncol = 4, align = "h")





head(cvresult_drug)

f_plot_cv_roc <- function(cv_results_df, species_name = NULL) {
  
  library(ggplot2)
  library(pROC)
  library(dplyr)
  
  # 检查是否有 fold 信息
  has_folds <- "fold" %in% names(cv_results_df) && 
    !all(is.na(cv_results_df$fold))
  
  if(has_folds) {
    # 原有的交叉验证逻辑
    # 1. 计算每个fold的ROC
    roc_list <- list()
    auc_values <- numeric(5)
    
    for(i in 1:5){
      fold_data <- cv_results_df[cv_results_df$fold == i, ]
      
      if(nrow(fold_data) > 0) {
        roc_obj <- roc(fold_data$true_label, fold_data$prob_S,
                       levels = c("R", "S"), quiet = TRUE)
        
        roc_list[[i]] <- roc_obj
        auc_values[i] <- as.numeric(auc(roc_obj))
        
        cat("Fold", i, "- AUC:", sprintf("%.4f", auc_values[i]), 
            "- Points:", length(roc_obj$sensitivities), "\n")
      }
    }
    
    # 2. Mean ROC
    roc_mean <- roc(cv_results_df$true_label, cv_results_df$prob_S,
                    levels = c("R", "S"), quiet = TRUE)
    mean_auc <- as.numeric(auc(roc_mean))
    
    cat("\nMean AUC:", sprintf("%.4f", mean_auc), "\n\n")
    
    # 3. 绘图
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
    
    p <- ggplot() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "gray50", linewidth = 0.5)
    
    # 绘制各个fold
    for(i in 1:5) {
      if(!is.null(roc_list[[i]])) {
        roc_obj <- roc_list[[i]]
        fold_df <- data.frame(
          FPR = 1 - roc_obj$specificities,
          TPR = roc_obj$sensitivities
        )
        
        p <- p + geom_path(data = fold_df,
                           aes(x = FPR, y = TPR),
                           color = colors[i],
                           alpha = 0.6,
                           linewidth = 0.8)
      }
    }
    
    # Mean ROC
    mean_df <- data.frame(
      FPR = 1 - roc_mean$specificities,
      TPR = roc_mean$sensitivities
    )
    
    p <- p + geom_path(data = mean_df,
                       aes(x = FPR, y = TPR),
                       color = "black",
                       linewidth = 1.2)
    
    # 标签和主题
    p <- p +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "False Positive Rate",
           y = "True Positive Rate",
           title = ifelse(is.null(species_name), "ROC Curves", species_name)) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0, face = "bold", size = 14),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
        plot.margin = ggplot2::margin(t = 1, r = 5, b = 1, l = 1, unit = "pt")
      )
    
    # 添加图例
    fold_labels <- paste0("Fold ", 1:5, " (AUC = ", sprintf("%.3f", auc_values), ")")
    mean_label <- paste0("Mean ROC (AUC = ", sprintf("%.3f", mean_auc), ")")
    
    legend_y <- seq(0.25, -0.02, length.out = 6)
    
    for(i in 1:5) {
      if(!is.null(roc_list[[i]])) {
        p <- p + annotate("segment",
                          x = 0.4, xend = 0.45,
                          y = legend_y[i], yend = legend_y[i],
                          color = colors[i],
                          linewidth = 0.8,
                          alpha = 0.6)
        
        p <- p + annotate("text",
                          x = 0.47, y = legend_y[i],
                          label = fold_labels[i],
                          hjust = 0, size = 3.2)
      }
    }
    
    p <- p + annotate("segment",
                      x = 0.35, xend = 0.4,
                      y = legend_y[6], yend = legend_y[6],
                      color = "black",
                      linewidth = 1.2)
    
    p <- p + annotate("text",
                      x = 0.42, y = legend_y[6],
                      label = mean_label,
                      hjust = 0, size = 3.2, fontface = "bold")
    
  } else {
    # 没有 fold 信息，只画整体 ROC
    cat("No fold information found. Plotting overall ROC curve.\n")
    
    roc_obj <- roc(cv_results_df$true_label, cv_results_df$prob_S,
                   levels = c("R", "S"), quiet = TRUE)
    
    auc_value <- as.numeric(auc(roc_obj))
    
    cat("Overall AUC:", sprintf("%.4f", auc_value), "\n\n")
    
    roc_df <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities
    )
    
    p <- ggplot() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "gray50", linewidth = 0.5) +
      geom_path(data = roc_df,
                aes(x = FPR, y = TPR),
                color = "#E41A1C",
                linewidth = 1.2) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "False Positive Rate",
           y = "True Positive Rate",
           title = ifelse(is.null(species_name), "ROC Curve", species_name)) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0, face = "bold", size = 14),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "gray90")
      )
    
    # 添加 AUC 标签
    auc_label <- paste0("AUC = ", sprintf("%.3f", auc_value))
    
    p <- p + annotate("text",
                      x = 0.6, y = 0.2,
                      label = auc_label,
                      hjust = 0, size = 4, fontface = "bold")
  }
  
  return(p)
}


# 方法1: 使用原有的 lapply 方式
plot_cv_tbportal <- lapply(seq_along(unique(cvresult_drug$drug)), function(i_d){
  drug_tmp <- unique(cvresult_drug$drug)[i_d]
  cat("\n=== Processing", drug_tmp, "===\n")
  
  # 筛选当前药物的数据
  drug_data <- cvresult_drug %>% filter(drug == drug_tmp)
  
  # 绘制 ROC 曲线
  p_tmp <- f_plot_cv_roc(drug_data, drug_tmp)
  
  return(p_tmp)
})


# plot_cv_tbportal[[1]]+labs(title = NULL)
# ggsave("./output_005/metric_cv_train_emb.pdf", plot_cv_tbportal[[1]]+labs(title = NULL), width = 3, height = 3)


plot_cv_combined_tbportal <- grid.arrange(
  grobs = plot_cv_tbportal, 
  ncol = 5,
  padding = unit(0, "line")  # 改为0
)

title_plot <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "(e)", 
           size = 5, fontface = "bold", hjust = 0) +
  theme_void() +
  theme(plot.margin = ggplot2::margin(t = 1, l = 1, b = 0, r = 0))

plot_cv_with_title <- grid.arrange(
  title_plot,
  arrangeGrob(grobs = plot_cv_tbportal, ncol = 5, padding = unit(0, "line")),
  ncol = 1,
  heights = c(0.05, 1)  # 标题占很小的高度
)

plot_combined_tbportal <- plot_grid(
  metric_combined_tbportal,
  plot_cv_with_title,
  ncol = 1,
  rel_heights = c(1, 2.5)  # 控制高度比例
)
ggsave("./output/metric_cv_train.pdf", plot_combined_tbportal, width = 16, height = 16)





