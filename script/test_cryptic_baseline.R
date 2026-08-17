#filename: test_cryptic_baseline.R
library(parallel)
library(randomForest)
library(vcfR)

setwd("~/bioinformatics/iTBDRS/")


dst_crypritc_HIGH <- read.csv(file = "./data/raw/cryptic_phy_HIGH.csv")
drugs <- c("EMB", "INH","PZA","RIF","SM", "AMK","CM",
           "ETO","KM","LFX","MFX","OFX","PAS")
drug_compare <- c("EMB","INH","RIF","AMK","ETO","KM","LFX","MFX")
# drug_compare <- drug_compare[1:2] #debug
path_predict_cryptic <- paste0("./result/cryptic_project/", dst_crypritc_HIGH$sraid,".csv")


#### 下面是用模型预测 Cryptic project

model_names <- c("rf","lasso_rf","mrmr_rf")


if(FALSE){
  start_time <- Sys.time()
  ## 加载模型文件
  cat("Loading all models...\n")
  all_models <- lapply(seq_along(drug_compare), function(i_d){
    drug_tmp <- drug_compare[i_d]
    path_model_drug <- paste0("./result/models_base/", drug_tmp, "_", model_names, ".rds")
    cat("Drug:",drug_tmp,"\n")
    models_list <- lapply(path_model_drug, function(path){
      if(file.exists(path)){
        model <- readRDS(path)
        return(model)
      } else {
        NULL
      }
    })
    names(models_list) <- model_names
    return(models_list)
  })
  names(all_models) <- drug_compare
  cat("All models loaded!\n")
  
  filename <- paste0(dst_crypritc_HIGH$sraid,".vcf.gz")
  result_pred <- lapply(seq_along(filename), function(i_f){
    cat('Sample:',i_f,"\n")
    filename_tmp <- filename[i_f]
    path_tmp <- paste0("../FS_Mtb_cryptic/data/vcf_files/",filename_tmp)
    file_tmp <- read.vcfR(path_tmp)
    pred_drug <- lapply(seq_along(drug_compare), function(i_d){
      drug_tmp <- drug_compare[i_d] 
      cat("Drug:",drug_tmp, "\n")
      # path_model_drug <- paste0("./result/models_base/",drug_tmp,"_",model_names ,".rds")
      
      result_baseline <- lapply(seq_along(model_names), function(i_m){
        # model_fit <- readRDS(path_model_drug[i_m])
        model_fit <- all_models[[drug_tmp]][[model_names[i_m]]] 
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
        # pred_drug_tmp <- predict(model,as.data.frame(t(feature_drug_tmp)),type = "prob")
        result_return <- data.frame(pred_drug_tmp)
        colnames(result_return) <- c("prob_R", "prob_S")
        result_return <- data.frame(drug = drug_tmp,result_return,
                                    Model = model_names[i_m])
        return(result_return)
      })
      
      result_modelS <- do.call('rbind',result_baseline)
      return(result_modelS)
    })
    pred_drug <- do.call('rbind',pred_drug)
    true_label <- dst_crypritc_HIGH[i_f,match(drug_compare,colnames(dst_crypritc_HIGH))]
    true_label <- as.character(true_label)
    pred_drug$true_label <- rep(true_label,each = 3)
    pred_drug$sra_id <- sub(".vcf.gz","", filename[i_f])
    pred_drug$fold <- NA
    # "sra_id","drug","fold","true_label","prob_R","prob_S" "Model"
    pred_drug <- pred_drug[,c(6,1,7,5,2:4)]
    write.csv(pred_drug,
              paste0("./result/cryptic_project_baseline/",sub(".vcf.gz", ".csv", filename_tmp)),
              row.names = F)
    return(pred_drug)
  })
  
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


#debug
path_predict_cryptic <- paste0("./result/cryptic_project_baseline/",dir("./result/cryptic_project_baseline/"))

result_cryptic <- do.call('rbind',lapply(path_predict_cryptic,read.csv))

# result_cryptic$Model <- gsub("^Model", "model", result_cryptic$Model)

set.seed(123)

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
# metrics_cryptic <- split.data.frame(metrics_cryptic,metrics_cryptic$Model)
write.csv(metrics_cryptic,"./output/metric_test_base.csv", row.names = F)


length(dir("./result/cryptic_project_baseline/"))
# ak <- readRDS("./result/models_base/AMK_lasso_rf.rds")
# ak$importance
# 
# ka <- readRDS("./result/models_005/AMK.rds")
# ka[[4]]$fit$importance
##### 画图 cryptic project ######

# 雷达图
data_line_cryptic_base <- read.csv("./output/metric_test_base.csv")
data_line_cryptic_base$Model <- sub("_","+",data_line_cryptic_base$Model)
data_line_cryptic_my <- read.csv("./output/metric_test.csv")
data_line_cryptic_my$Model <- "Our method"

data_plot_line <- rbind(data_line_cryptic_base,data_line_cryptic_my)
plot_metrics_facet <- function(data, 
                               colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),
                               line_size = 0.6,
                               point_size = 2,
                               auto_ylim = TRUE,
                               ncol = 2) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  # 准备数据 - 转换为长格式
  plot_data <- data %>%
    select(Drug, Model, Accuracy, MCC, AUC, AUPRC) %>%
    pivot_longer(cols = c(Accuracy, MCC, AUC, AUPRC),
                 names_to = "Metric",
                 values_to = "value") %>%
    mutate(
      Drug = factor(Drug, levels = unique(data$Drug)),
      Model = factor(Model, levels = unique(data$Model)),
      # 设置分面标签和顺序
      Metric = factor(Metric, 
                      levels = c("Accuracy", "MCC", "AUC", "AUPRC"),
                      labels = c("(a) Accuracy", "(b) MCC", "(c) AUC", "(d) AUPRC"))
    )
  
  # 绘图
  p <- ggplot(plot_data, aes(x = Drug, y = value, 
                             color = Model, group = Model)) +
    # 折线
    geom_line(linewidth = line_size, alpha = 0.8) +
    # 数据点
    geom_point(size = point_size, alpha = 0.8) +
    # 颜色
    scale_color_manual(values = colors) +
    # 分面
    facet_wrap(~ Metric, scales = ifelse(auto_ylim, "free_y", "fixed"), 
               ncol = ncol) +
    # 标题
    labs(x = NULL, y = NULL) +
    # 主题
    theme_minimal() +
    theme(
      # 🔥 分面标签 - 去掉背景阴影,和雷达图一致
      strip.text = element_text(size = 14, face = "bold", hjust = 0.5),  # 居中,和雷达图标题一样
      strip.background = element_blank(),  # 去掉背景
      
      # 🔥 坐标轴文字 - 和雷达图一致
      axis.text.x = element_text(size = 10, face = "bold", 
                                 color = "grey30", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9, color = "grey30"),
      
      # 🔥 图例 - 和雷达图一致
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      
      # 网格线
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      # 间距
      panel.spacing = unit(1, "lines"),
      
      # 背景
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

metric_combined_cryptic <- plot_metrics_facet(data_plot_line)
metric_combined_cryptic
ggsave("./output/metric_test_baseline.pdf",metric_combined_cryptic
       , width = 5, height = 5)
write.csv(data_plot_line,"./output/metric_test_baseline.csv",row.names = F)

data_plot_line_drug <- split.data.frame(data_plot_line,data_plot_line$Drug)

lapply(drug_compare, function(x)rowSums(data_plot_line_drug[[x]][,1:4]))
