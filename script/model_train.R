#filename: model_train.R

setwd("~/bioinformatics/ML_mtb/")
options(expressions = 500000)
drugs <- c("EMB", "INH","PZA","RIF","SM", "AMK","CM",
           "ETO","KM","LFX","MFX","OFX","PAS")


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
# 计算分类性能指标的函数 (R为阳性类)
calculate_metrics <- function(rf_model, data, response_var = "type") {
  
  # 获取预测值
  predictions <- predict(rf_model, newdata = data, type = "response")
  
  # 获取真实值
  actual <- data[[response_var]]
  
  # 确保是因子类型且水平一致
  actual <- factor(actual, levels = levels(predictions))
  
  # 构建混淆矩阵
  confusion <- table(Predicted = predictions, Actual = actual)
  
  # 提取 TN, FP, FN, TP (强制转换为数值类型避免溢出)
  if(all(c("R", "S") %in% rownames(confusion))) {
    TP <- as.numeric(confusion["R", "R"])
    FN <- as.numeric(confusion["S", "R"])
    FP <- as.numeric(confusion["R", "S"])
    TN <- as.numeric(confusion["S", "S"])
  } else {
    stop("混淆矩阵中缺少 'R' 或 'S' 类别")
  }
  
  # 处理可能的NA值
  TP <- ifelse(is.na(TP), 0, TP)
  FN <- ifelse(is.na(FN), 0, FN)
  FP <- ifelse(is.na(FP), 0, FP)
  TN <- ifelse(is.na(TN), 0, TN)
  
  # 计算各项指标(添加除零保护)
  Recall <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  Precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  F1_score <- ifelse((Precision + Recall) == 0, 0, 
                     2 * (Precision * Recall) / (Precision + Recall))
  
  # Matthews Correlation Coefficient (MCC) - 防止整数溢出
  MCC_numerator <- (TP * TN) - (FP * FN)
  
  # 分步计算,确保使用数值类型
  term1 <- as.numeric(TP + FP)
  term2 <- as.numeric(TP + FN)
  term3 <- as.numeric(TN + FP)
  term4 <- as.numeric(TN + FN)
  
  MCC_denominator <- sqrt(term1 * term2 * term3 * term4)
  MCC <- ifelse(MCC_denominator == 0, 0, MCC_numerator / MCC_denominator)
  
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
      row.names = NULL
    )
  )
  
  return(results)
}
# 逐步递减特征选择,找到MCC最优的特征组合
backward_selection_by_mcc <- function(data, response_var = "type",min_features = 5,verbose = TRUE){
  # 初始化
  current_features <- setdiff(colnames(data), response_var)
  best_mcc <- -Inf
  best_features <- NULL
  best_model <- NULL
  best_performance <- NULL
  performance_history <- list()
  
  iteration <- 1
  
  while(length(current_features) > min_features) {
    
    if(verbose) {
      cat(sprintf("\n=== Iteration %d: Testing %d features ===\n", 
                  iteration, length(current_features)))
    }
    
    # 使用当前特征训练模型
    subset_data <- data[, c(current_features, response_var)]
    fit_current <- randomForest(type ~ ., data = subset_data, importance = TRUE)
    
    # 在原始完整数据上评估性能
    perf_current <- calculate_metrics(fit_current, data, response_var)
    current_mcc <- perf_current$metrics$MCC
    
    if(verbose) {
      cat(sprintf("Current MCC: %.4f\n", current_mcc))
      cat(sprintf("Accuracy: %.4f, F1: %.4f, Precision: %.4f, Recall: %.4f\n",
                  perf_current$metrics$Accuracy,
                  perf_current$metrics$F1_score,
                  perf_current$metrics$Precision,
                  perf_current$metrics$Recall))
    }
    
    # 保存历史记录
    performance_history[[iteration]] <- list(
      n_features = length(current_features),
      features = current_features,
      mcc = current_mcc,
      all_metrics = perf_current$metrics,
      confusion_matrix = perf_current$confusion_matrix
    )
    
    # 更新最佳结果(只看MCC)
    if(current_mcc > best_mcc) {
      best_mcc <- current_mcc
      best_features <- current_features
      best_model <- fit_current
      best_performance <- perf_current
      
      if(verbose) {
        cat(sprintf("*** New best MCC: %.4f with %d features ***\n", 
                    best_mcc, length(best_features)))
      }
    }
    
    # 获取特征重要性(MeanDecreaseAccuracy)
    importance_scores <- importance(fit_current, type = 1)
    
    # 移除重要性最低的特征
    least_important <- names(which.min(importance_scores[, 1]))
    current_features <- setdiff(current_features, least_important)
    
    if(verbose) {
      cat(sprintf("Removed feature: %s (importance: %.4f)\n", 
                  least_important, 
                  min(importance_scores[, 1])))
    }
    
    iteration <- iteration + 1
  }
  
  # 最终结果
  if(verbose) {
    cat("\n" , rep("=", 60), "\n", sep="")
    cat("FINAL BEST MODEL\n")
    cat(rep("=", 60), "\n", sep="")
    cat(sprintf("Best MCC: %.4f\n", best_mcc))
    cat(sprintf("Number of features: %d\n", length(best_features)))
    cat("\nBest features:\n")
    print(best_features)
    cat("\nFull performance metrics:\n")
    print(best_performance$metrics)
    cat("\nConfusion matrix:\n")
    print(best_performance$confusion_matrix)
  }
  
  # 返回结果
  return(list(
    best_features = best_features,
    best_mcc = best_mcc,
    best_model = best_model,
    best_performance = best_performance,
    performance_history = performance_history,
    n_best_features = length(best_features)
  ))
}


dst <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_DST_sample_clean.csv")
X_genomic <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_X_genomic_005.csv")


set.seed(123)

all(X_genomic$sra_id == dst$sra_id)



library(randomForest)
library(caret)
library(pROC)
library(parallel)
# x_train <- X_genomic[,-(1:2)];y_train <- dst$levofloxacin
# no_D <- which(is.na(y_train))
# x_train <- x_train[-no_D,];y_train <- y_train[-no_D]
# mic_lfx <- get_mic(x_train,y_train)
get_mic <- function(x_train, y_train) {
  library(minerva)
  # 确保数据为数值型矩阵
  x_train <- as.matrix(x_train)
  y_train <- as.integer(ifelse(y_train =="S", "0", "1"))
  
  # 计算每个SNP的MIC值
  name_list <- colnames(x_train)
  mic_list <- numeric(length(name_list))
  
  # 添加进度显示
  cat("计算MIC值中...\n")
  pb <- txtProgressBar(min = 0, max = length(name_list), style = 3)
  
  for(i in seq_along(name_list)) {
    # 检查该SNP是否有变异
    if(var(x_train[,i]) > 0) {
      tryCatch({
        mic_list[i] <- mine(x_train[,i], y_train)$MIC
      }, error = function(e) {
        mic_list[i] <- NA
        message("\nSNP ", name_list[i], " 计算出错")
      })
    } else {
      mic_list[i] <- 0  # 对于无变异的SNP，设置MIC为0
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 创建结果数据框并排序
  gene_mic <- data.frame(
    SNP = name_list,
    MIC = mic_list,
    stringsAsFactors = FALSE
  )
  
  # 移除NA值并排序
  gene_mic <- na.omit(gene_mic)
  gene_mic <- gene_mic[order(gene_mic$MIC, decreasing = TRUE),]
  
  # 添加排名
  gene_mic$Rank <- 1:nrow(gene_mic)
  
  return(gene_mic)
}

## debug 
# i_d <- 1;i_f <- 1;i_imp <- 1

lapply(seq_along(drugs), function(i_d){
  drug_tmp <- drugs[i_d]
  cat('Drug:',drug_tmp, "\n")
  y_v <- dst[,which(colnames(dst) %in% drug_list(drug_tmp))]
  no_D <- which(is.na(y_v))
  y_tmp <- factor(y_v[-no_D], levels = c("R","S"))
  x_tmp <- X_genomic[-no_D,]
  set_drug <- data.frame(x_tmp,type = y_tmp)
  rownames(set_drug) <- NULL
  
  
  library(caret)
  fold_indices <- createFolds(set_drug$type, k = 5, list = FALSE)
  set_drug$fold <- fold_indices
  indices_save <- set_drug[,c(1:2,ncol(set_drug)+(-1:0))]
  write.csv(indices_save, 
            paste0("./result/TBportal/indices_fold/",drug_tmp,".csv"),
            row.names = F)
  mcc_kn_fold <- lapply(1:5, function(i_f){
    cat("\n",'Fold:',i_f,"\n")
    set_fold <- set_drug[which(set_drug$fold == i_f),]
    set_input <- set_fold[,-c(1:2,ncol(set_fold))]
    rf_model <- randomForest(
     type ~ .,                    # 公式
     data = set_input,           # 训练数据
     ntree = 10,                 # 树的数量
     mtry = sqrt(ncol(set_input)-1),  # 每次分裂考虑的变量数
     importance = TRUE,           # 关键参数:计算变量重要性
     proximity = FALSE,
     oob.prox = FALSE)
    importance_scores <- importance(rf_model, type = 1)
    # varImpPlot(rf_model, type = 1)  # MDA可视化
    no_D <- which(importance_scores[,1] <=0)

    set_input_filter <- set_input[,-no_D,drop = F]
    # set_input_filter <- set_input ## 去掉MDA筛选模块
    mic_scores <- get_mic(set_input_filter[, -ncol(set_input_filter)],
                          set_input_filter[, ncol(set_input_filter)])
    feature_topK <- mic_scores$SNP[1:5] #设置成5个最重要的snp
    cat("\n",'pair mcc 计算中. \n')
    pair_mcc <- lapply(seq_along(feature_topK), function(i_imp){
      feature_tmp <- feature_topK[i_imp]
      no_imp <- which(colnames(set_input) %in% feature_tmp)
      pair_mcc_sub <- mclapply(seq(ncol(set_input)-1), function(i_tmp){
        # print(i_tmp)
        if(no_imp == i_tmp){
          return(0)
        }
        fit_tmp <- randomForest(type ~ ., data = set_input[,c(no_imp,i_tmp,ncol(set_input))])
        performance <- calculate_metrics(fit_tmp, set_input, response_var = "type")
        return(performance$metrics$MCC)
      },mc.cores = 6)
      pair_mcc_sub <- do.call('cbind',pair_mcc_sub)
      return(pair_mcc_sub)
    })
    pair_mcc <- do.call('rbind',pair_mcc)
    colnames(pair_mcc) <- colnames(set_input)[-ncol(set_input)]
    rownames(pair_mcc) <- feature_topK
    cat("\n",'pair mcc 计算完成. \n')
    cb <- expand.grid(1:5,1:5);colnames(cb) <- c("k","n")
    mcc_kn <- lapply(seq(nrow(cb)), function(i_cb){
     k_tmp <- cb[i_cb,1]
     n_tmp <- cb[i_cb,2]
     top_features_list <- apply(pair_mcc[1:k_tmp, , drop = FALSE], 1, function(row) {
       top_indices <- order(row, decreasing = TRUE)[1:n_tmp]
       colnames(pair_mcc)[top_indices]
     })
     
     feature_pool <- unique(c(rownames(pair_mcc)[1:k_tmp],
                              as.vector(top_features_list)))
     subset_tmp <- set_input[, match(c(feature_pool,"type"),
                                     colnames(set_input)),drop = F]
     # fit_pool <- randomForest(type ~ ., data = subset_tmp)
     # performance_pool <- calculate_metrics(fit_pool, set_input, response_var = "type")
     
     result_kn_tmp <- backward_selection_by_mcc(
       data = subset_tmp,
       response_var = "type",
       min_features = 1,
       verbose = TRUE
     )
     result_return <- data.frame(bestfeature = paste(result_kn_tmp$best_features,
                                                     collapse = "_"),
                                 bestMCC = result_kn_tmp$best_mcc)
     return(result_return)
     # return(list(result_return,
     #             result_kn_tmp))
    })
    mcc_kn <- do.call('rbind',mcc_kn)
    return(data.frame(cb, MCC = mcc_kn, fold = i_f))
  })
  mcc_kn_fold <- do.call('rbind',mcc_kn_fold)
  mcc_kn_fold$drug <- drug_tmp
  write.csv(mcc_kn_fold, paste0("./result/TBportal/k-n-search/",
                   drug_tmp,".csv"),row.names = F)
  result_mean <- aggregate(MCC.bestMCC ~ k + n, data = mcc_kn_fold, FUN = mean)
  no_best <- which.max(result_mean[,3])
  best_features <- mcc_kn_fold[mcc_kn_fold$k == result_mean[no_best,1] & 
                                 mcc_kn_fold$n == result_mean[no_best,2], 3]
  best_features <- unlist(lapply(best_features, strsplit,split = "_"))
  
  feature_stable <- sort(table(best_features), decreasing = T)
  
  performance_stable <- lapply(5:1, function(i_s){
    fea_tmp <- names(feature_stable)[which(feature_stable >= i_s)]
    
    # 如果没有特征，返回空结构
    if(length(fea_tmp) == 0) {
      message(paste("⚠️  stable =", i_s, ": 没有特征"))
      
      empty_result <- list(
        confusion_matrix = matrix(NA, nrow = 2, ncol = 2,
                                  dimnames = list(Predicted = c("R", "S"),
                                                  Actual = c("R", "S"))),
        metrics = data.frame(
          TN = NA_real_,
          FP = NA_real_,
          FN = NA_real_,
          TP = NA_real_,
          Recall = NA_real_,
          Precision = NA_real_,
          Accuracy = NA_real_,
          F1_score = NA_real_,
          MCC = NA_real_
        ),
        stable = i_s,
        fit = NULL
      )
      return(empty_result)
    }
    
    set_input_stable <- set_drug[,match(c(fea_tmp, "type"),colnames(set_drug))]
    fit_tmp <- randomForest(type ~., set_input_stable)
    performance_tmp <- calculate_metrics(fit_tmp, set_input_stable, response_var = "type")
    performance_tmp$stable <- i_s
    performance_tmp$fit <- fit_tmp
    return(performance_tmp)
  })
  
  saveRDS(performance_stable, paste0("./result/models/",
                                     drug_tmp,".rds"))
  
})




f_plot_kn_heatmap <- function(kn_data, drug_name = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  # 计算每个 k-n 组合的平均 MCC
  heatmap_data <- kn_data %>%
    group_by(k, n) %>%
    summarise(
      mean_MCC = mean(MCC.bestMCC, na.rm = TRUE),
      sd_MCC = sd(MCC.bestMCC, na.rm = TRUE),
      n_folds = n(),
      .groups = "drop"
    )
  
  # 打印统计信息
  cat("Drug:", ifelse(is.null(drug_name), "Unknown", drug_name), "\n")
  cat("k range:", min(heatmap_data$k), "-", max(heatmap_data$k), "\n")
  cat("n range:", min(heatmap_data$n), "-", max(heatmap_data$n), "\n")
  cat("MCC range:", sprintf("%.4f - %.4f", 
                            min(heatmap_data$mean_MCC), 
                            max(heatmap_data$mean_MCC)), "\n\n")
  
  # 绘制热图
  p <- ggplot(heatmap_data, aes(x = k, y = n, fill = mean_MCC)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", mean_MCC)), 
              color = "black", size = 3) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7", 
      high = "#B2182B",
      midpoint = median(heatmap_data$mean_MCC),
      name = "Mean MCC"
    ) +
    scale_x_continuous(breaks = 1:5, expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:5, expand = c(0, 0)) +
    labs(
      x = "k (Number of Top SNP)",
      y = "n (Number of SNP pairs)",
      title = ifelse(is.null(drug_name), 
                     "Mean MCC Heatmap", 
                     drug_name)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    )
  
  return(p)
}

library(ggplot2)
library(patchwork)
library(cowplot)
library(pROC)

plot_heat <- lapply(seq_along(drugs), function(i_d){
  drug_tmp <- drugs[i_d]
  cat("Drug:",drug_tmp, "\n")
  kn_raw <- read.csv(paste0("./result/TBportal/k-n-search/",drug_tmp,".csv"))
  plot_tmp <- f_plot_kn_heatmap(kn_raw,drug_tmp)
})

# plot_heat[[1]]+labs(title = NULL)
# ggsave("./output_005/heatmap_train_emb.pdf", 
#        plot_heat[[1]]+labs(title = NULL), width = 3, height = 3)


plot_heat_tbportal <- grid.arrange(
  grobs = plot_heat, 
  ncol = 5,
  padding = unit(0, "line")  # 改为0
)

title_plot <- ggplot() + 
  annotate("text", x = 0, y = 0, label = "", 
           size = 5, fontface = "bold", hjust = 0) +
  theme_void() +
  theme(plot.margin = ggplot2::margin(t = 1, l = 1, b = 0, r = 0))

plot_heat_with_title <- grid.arrange(
  title_plot,
  plot_heat_tbportal,
  ncol = 1,
  heights = c(0.05, 1)
)
ggsave("./output/heatmap.pdf", plot_heat_with_title, width = 16, height = 8)



