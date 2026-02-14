#filename: model_train_baseline.R
# library(xgboost)
library(randomForest)
library(parallel)
set.seed(123)
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
calculate_metrics <- function(model, data, response_var = "type", model_type = "rf",
                              positive_class = "R", threshold = 0.5) {
  # 统一取出特征与标签（兼容：若 data 无该列名，则默认最后一列为标签）
  if (is.data.frame(data)) {
    if (response_var %in% colnames(data)) {
      X <- data[, !colnames(data) %in% response_var, drop = FALSE]
      y <- data[[response_var]]
    } else {
      X <- data[, -ncol(data), drop = FALSE]
      y <- data[[ncol(data)]]
      warning(sprintf("未找到列名 '%s'，已使用最后一列作为标签", response_var))
    }
  } else {
    stop("data 需要是 data.frame")
  }
  
  # 统一标签为因子，并确保水平包含正负类
  if (!is.factor(y)) y <- factor(y)
  # 若水平不含我们设定的正负类，按现有水平顺序确定：正类设为第二个
  if (!all(c(positive_class) %in% levels(y))) {
    if (nlevels(y) != 2) stop("当前实现仅支持二分类且必须能推断正类")
    positive_class <- levels(y)[min(2, nlevels(y))]
  }
  negative_class <- setdiff(levels(y), positive_class)[1]
  y <- factor(y, levels = c(negative_class, positive_class))
  
  # 根据模型类型预测
  if (model_type %in% c("xgboost","xgb")) {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("需要 xgboost 包用于 XGBoost 预测")
    }
    dtest <- xgboost::xgb.DMatrix(data = as.matrix(X))
    pred_prob <- predict(model, dtest)
    # 约定：pred_prob 为 label=1 的概率，对应正类 "R"
    # 若你的训练时编码确认为 1=R, 0=S，则保持如下；否则请调整 positive_class
    preds <- ifelse(pred_prob >= threshold, positive_class, negative_class)
    preds <- factor(preds, levels = c(negative_class, positive_class))
    p_pos <- pred_prob
  } else if (model_type %in% c("rf","randomforest")) {
    # 用概率 + 阈值生成预测（比 type="response" 更灵活）
    probs <- predict(model, newdata = data, type = "prob")
    # 若缺少概率（例如回归或某些设置），回退到 response
    if (is.null(probs) || !all(c(negative_class, positive_class) %in% colnames(probs))) {
      preds <- predict(model, newdata = data, type = "response")
      preds <- factor(preds, levels = c(negative_class, positive_class))
      # 无概率时无法算 AUC/阈值优化，但仍可算混淆矩阵与 MCC
      p_pos <- NULL
    } else {
      p_pos <- probs[, positive_class]
      preds <- ifelse(p_pos >= threshold, positive_class, negative_class)
      preds <- factor(preds, levels = c(negative_class, positive_class))
    }
  } else {
    stop("不支持的模型类型。请使用 'rf', 'randomforest', 'xgboost' 或 'xgb'")
  }
  
  # 混淆矩阵（确保维度完整）
  confusion <- table(Actual = y, Predicted = preds)
  # 填补缺失水平，防止索引不到
  if (!all(c(negative_class, positive_class) %in% rownames(confusion))) {
    confusion <- addmargins(confusion[levels(y), levels(y), drop = FALSE], FUN = NULL)
    confusion <- confusion[1:2, 1:2, drop = FALSE]
  }
  
  # 提取 TP/TN/FP/FN，并转为 numeric 避免整型溢出
  TP <- as.numeric(confusion[positive_class, positive_class])
  FN <- as.numeric(confusion[positive_class, negative_class])
  FP <- as.numeric(confusion[negative_class, positive_class])
  TN <- as.numeric(confusion[negative_class, negative_class])
  
  # 指标
  acc <- (TP + TN) / (TP + TN + FP + FN)
  precision <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  recall <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_ else {
    2 * precision * recall / (precision + recall)
  }
  
  # 稳健 MCC
  a <- TP + FP; b <- TP + FN; c <- TN + FP; d <- TN + FN
  mcc <- if (a == 0 || b == 0 || c == 0 || d == 0) {
    NA_real_
  } else {
    num <- (TP * TN) - (FP * FN)
    den <- sqrt(as.numeric(a) * as.numeric(b) * as.numeric(c) * as.numeric(d))
    if (!is.finite(den) || den == 0) NA_real_ else num / den
  }
  
  # AUC（有概率时才计算）
  auc <- NA_real_
  if (!is.null(p_pos) && requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(response = y, predictor = p_pos,
                         levels = c(negative_class, positive_class), direction = "<")
    auc <- as.numeric(pROC::auc(roc_obj))
  }
  
  out <- list(
    confusion_matrix = confusion,
    metrics = data.frame(
      TN = TN, FP = FP, FN = FN, TP = TP,
      Recall = round(recall, 4),
      Precision = round(precision, 4),
      Accuracy = round(acc, 4),
      F1_score = round(f1, 4),
      MCC = round(mcc, 4),
      AUC = round(auc, 4),
      row.names = NULL
    ),
    positive_class = positive_class,
    threshold = threshold
  )
  if (!is.null(p_pos)) out$probs <- p_pos
  return(out)
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

# ==================== 完整替换 Lasso 特征选择函数 ====================
library(glmnet)
library(randomForest)

lasso_feature_selection <- function(data, response_var = "type", alpha = 1, 
                                    nfolds = 5, selection_method = "lambda.1se",
                                    use_mcc = FALSE) {  # 默认 FALSE 保持原有行为
  
  # 准备数据
  X <- as.matrix(data[, !colnames(data) %in% response_var])
  y <- data[[response_var]]
  
  # 转换响应变量为0/1
  if (is.factor(y)) {
    # 确保 R=1, S=0
    y_numeric <- as.numeric(y == "R")
  } else {
    y_numeric <- y
  }
  
  if (use_mcc) {
    cat("  [使用 MCC 作为评估指标]\n")
    
    # 🔧 内嵌的 MCC 计算函数(仅用于 Lasso)
    calculate_mcc_internal <- function(actual_numeric, pred_numeric) {
      # actual_numeric 和 pred_numeric 都是 0/1
      # R=1(正类), S=0(负类)
      
      TP <- sum(actual_numeric == 1 & pred_numeric == 1)
      TN <- sum(actual_numeric == 0 & pred_numeric == 0)
      FP <- sum(actual_numeric == 0 & pred_numeric == 1)
      FN <- sum(actual_numeric == 1 & pred_numeric == 0)
      
      # 处理NA
      TP <- ifelse(is.na(TP), 0, TP)
      FN <- ifelse(is.na(FN), 0, FN)
      FP <- ifelse(is.na(FP), 0, FP)
      TN <- ifelse(is.na(TN), 0, TN)
      
      # 计算MCC
      MCC_numerator <- (TP * TN) - (FP * FN)
      term1 <- as.numeric(TP + FP)
      term2 <- as.numeric(TP + FN)
      term3 <- as.numeric(TN + FP)
      term4 <- as.numeric(TN + FN)
      MCC_denominator <- sqrt(term1 * term2 * term3 * term4)
      
      MCC <- ifelse(MCC_denominator == 0, 0, MCC_numerator / MCC_denominator)
      
      return(MCC)
    }
    
    # 手动交叉验证计算 MCC
    set.seed(123)
    n <- nrow(X)
    folds <- cut(seq(1, n), breaks = nfolds, labels = FALSE)
    
    # Lambda 序列
    lambda_seq <- exp(seq(log(0.001), log(1), length.out = 100))
    mcc_scores <- numeric(length(lambda_seq))
    
    for (i in seq_along(lambda_seq)) {
      fold_mccs <- numeric(nfolds)
      
      for (fold in 1:nfolds) {
        # 分割数据
        test_idx <- which(folds == fold)
        train_idx <- which(folds != fold)
        
        # 训练模型
        fit <- glmnet(
          x = X[train_idx, ],
          y = y_numeric[train_idx],
          family = "binomial",
          alpha = alpha,
          lambda = lambda_seq[i]
        )
        
        # 预测
        pred_prob <- predict(fit, newx = X[test_idx, ], type = "response")
        pred_numeric <- ifelse(pred_prob > 0.5, 1, 0)
        
        # 🔧 使用内嵌函数计算 MCC
        fold_mccs[fold] <- calculate_mcc_internal(y_numeric[test_idx], pred_numeric)
      }
      
      mcc_scores[i] <- mean(fold_mccs, na.rm = TRUE)
    }
    
    # 选择最优 lambda
    if (selection_method == "lambda.min") {
      best_idx <- which.max(mcc_scores)
    } else {
      # lambda.1se
      best_mcc <- max(mcc_scores, na.rm = TRUE)
      mcc_se <- sd(mcc_scores, na.rm = TRUE) / sqrt(length(mcc_scores))
      threshold <- best_mcc - mcc_se
      candidates <- which(mcc_scores >= threshold)
      best_idx <- candidates[which.max(lambda_seq[candidates])]
    }
    
    best_lambda <- lambda_seq[best_idx]
    best_mcc <- mcc_scores[best_idx]
    
    cat("  - 最优 MCC:", round(best_mcc, 4), "\n")
    
    # 用最优 lambda 训练最终模型
    cv_lasso <- glmnet(
      x = X,
      y = y_numeric,
      family = "binomial",
      alpha = alpha,
      lambda = best_lambda
    )
    
    # 添加 MCC 信息到结果
    extra_info <- list(
      best_mcc = best_mcc,
      mcc_scores = mcc_scores,
      lambda_seq = lambda_seq
    )
    
  } else {
    # 原始方法(不使用 MCC)
    set.seed(123)
    cv_lasso <- cv.glmnet(
      x = X,
      y = y_numeric,
      family = "binomial",
      alpha = alpha,
      nfolds = nfolds,
      type.measure = "deviance"
    )
    
    if (selection_method == "lambda.min") {
      best_lambda <- cv_lasso$lambda.min
    } else {
      best_lambda <- cv_lasso$lambda.1se
    }
    
    extra_info <- list()
  }
  
  # 提取非零系数的特征
  lasso_coef <- coef(cv_lasso, s = best_lambda)
  selected_features <- rownames(lasso_coef)[which(lasso_coef != 0)]
  selected_features <- selected_features[selected_features != "(Intercept)"]
  
  # 返回结果
  result <- c(
    list(
      cv_model = cv_lasso,
      best_lambda = best_lambda,
      selected_features = selected_features,
      n_features = length(selected_features),
      all_coefficients = as.matrix(lasso_coef)
    ),
    extra_info
  )
  
  return(result)
}

# ==================== 完整替换训练函数 ====================
train_lasso_rf <- function(data, response_var = "type", 
                           lasso_alpha = 1, 
                           lasso_nfolds = 5,
                           lambda_selection = "lambda.1se",
                           use_mcc = FALSE,  # 新增参数
                           rf_ntree = 500,
                           rf_mtry = NULL) {
  
  cat("步骤 1: Lasso 特征选择...\n")
  
  # Lasso特征选择
  lasso_result <- lasso_feature_selection(
    data = data,
    response_var = response_var,
    alpha = lasso_alpha,
    nfolds = lasso_nfolds,
    selection_method = lambda_selection,
    use_mcc = use_mcc
  )
  
  cat("  - 原始特征数:", ncol(data) - 1, "\n")
  cat("  - Lasso选择特征数:", lasso_result$n_features, "\n")
  cat("  - 最优 lambda:", round(lasso_result$best_lambda, 6), "\n")
  
  if (lasso_result$n_features == 0) {
    stop("Lasso没有选择任何特征!请尝试调整lambda或使用lambda.min")
  }
  
  # 准备RF训练数据
  selected_data <- data[, c(lasso_result$selected_features, response_var)]
  
  cat("\n步骤 2: 使用选定特征训练随机森林...\n")
  
  # 设置mtry
  if (is.null(rf_mtry)) {
    rf_mtry <- floor(sqrt(ncol(selected_data) - 1))
  }
  
  cat("  - RF 参数: ntree =", rf_ntree, ", mtry =", rf_mtry, "\n")
  
  # 训练随机森林
  rf_model <- randomForest(
    as.formula(paste(response_var, "~ .")),
    data = selected_data,
    ntree = rf_ntree,
    mtry = rf_mtry,
    importance = TRUE,
    proximity = FALSE,
    oob.prox = FALSE
  )
  
  cat("  - OOB 错误率:", round(rf_model$err.rate[rf_ntree, "OOB"] * 100, 2), "%\n")
  
  # 🔧 使用你原有的 calculate_metrics 计算 RF 的 MCC
  rf_metrics <- calculate_metrics(
    model = rf_model,
    data = selected_data,
    response_var = response_var,
    model_type = "rf"
  )
  
  cat("  - OOB MCC:", rf_metrics$metrics$MCC, "\n")
  
  # 返回完整结果
  result <- list(
    lasso_result = lasso_result,
    rf_model = rf_model,
    rf_metrics = rf_metrics,
    selected_data = selected_data,
    selected_features = lasso_result$selected_features,
    n_features = lasso_result$n_features
  )
  
  class(result) <- c("lasso_rf", "list")
  
  return(result)
}
# ==================== mRMR+ RF 特征选择 ====================
# 加载必要的包
# ==================== mRMR + RF with iFS (MCC-based) ====================
library(mRMRe)
library(randomForest)

# 计算 OOB-MCC（二分类）
.oob_mcc <- function(rf_model) {
  # 需要 randomForest 分类模型（有 votes 和 y）
  if (is.null(rf_model$votes) || is.null(rf_model$y)) return(NA_real_)
  classes <- colnames(rf_model$votes)
  if (length(classes) != 2) return(NA_real_)  # 仅支持二分类
  pos <- classes[2]
  neg <- classes[1]
  
  # OOB 概率（或投票比例）
  scores <- rf_model$votes[, pos]
  pred <- ifelse(scores >= 0.5, pos, neg)
  true <- rf_model$y
  
  # 统一因子水平
  pred <- factor(pred, levels = c(neg, pos))
  true <- factor(true, levels = c(neg, pos))
  
  # 计数并转换为 numeric，避免整数溢出
  tp <- as.numeric(sum(pred == pos & true == pos))
  tn <- as.numeric(sum(pred == neg & true == neg))
  fp <- as.numeric(sum(pred == pos & true == neg))
  fn <- as.numeric(sum(pred == neg & true == pos))
  
  # 分母（四项乘积）的稳健计算
  a <- tp + fp
  b <- tp + fn
  c <- tn + fp
  d <- tn + fn
  
  # 任何项为 0 → 分母为 0，返回 NA 或约定值
  if (a == 0 || b == 0 || c == 0 || d == 0) {
    return(NA_real_)
  }
  
  denom <- sqrt(a * b * c * d)  # numeric 乘法，不会整型溢出
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  
  ((tp * tn) - (fp * fn)) / denom
}

# iFS：前缀 1..K 上用 OOB-MCC 评估
.incremental_rf_eval_mcc <- function(data, response_var, ordered_features,
                                     rf_ntree = 500, rf_mtry = NULL,
                                     verbose = TRUE, seed = NULL) {
  K <- length(ordered_features)
  if (K < 1) stop("ordered_features 为空，无法进行 iFS")
  
  # mtry 设定
  if (is.null(rf_mtry)) {
    mtry_fun <- function(p) floor(sqrt(p))
  } else if (length(rf_mtry) == 1 && is.numeric(rf_mtry)) {
    mtry_const <- rf_mtry
    mtry_fun <- function(p) mtry_const
  } else {
    stop("rf_mtry 需要为 NULL 或 单个数值")
  }
  
  results <- data.frame(k = integer(), mcc = numeric(), stringsAsFactors = FALSE)
  if (!is.null(seed)) set.seed(seed)
  
  for (k in seq_len(K)) {
    feats_k <- ordered_features[1:k]
    df_k <- data[, c(feats_k, response_var), drop = FALSE]
    rf_k <- randomForest(
      as.formula(paste(response_var, "~ .")),
      data = df_k,
      ntree = rf_ntree,
      mtry = mtry_fun(k),
      importance = FALSE,
      proximity = FALSE,
      oob.prox = FALSE
    )
    val <- .oob_mcc(rf_k)
    results <- rbind(results, data.frame(k = k, mcc = val))
    if (verbose && (k %% max(1, floor(K/10)) == 0)) {
      cat(sprintf("  - iFS 进度: k=%d/%d, OOB-MCC=%.4f\n", k, K, val))
    }
  }
  
  # 最优 k（最大 OOB-MCC）
  best_idx <- which.max(results$mcc)
  best_row <- results[best_idx, , drop = FALSE]
  list(curve = results, best_k = best_row$k, best_mcc = best_row$mcc)
}

mrmr_feature_selection <- function(data, response_var = "type", 
                                   n_features = 30) {
  cat("准备 mRMR 数据...\n")
  
  # 准备数据
  feature_data <- data[, !colnames(data) %in% response_var, drop = FALSE]
  target <- data[[response_var]]
  
  # 转换目标变量为数值（mRMRe 需要数值矩阵；目标放第1列）
  if (is.factor(target)) {
    target_numeric <- as.numeric(target) - 1
  } else {
    target_numeric <- target
  }
  
  cat("  - 检查特征数据类型...\n")
  for (col in colnames(feature_data)) {
    if (is.factor(feature_data[[col]])) {
      feature_data[[col]] <- as.numeric(feature_data[[col]])
      cat("    转换因子列:", col, "\n")
    } else if (is.character(feature_data[[col]])) {
      feature_data[[col]] <- as.numeric(as.factor(feature_data[[col]]))
      cat("    转换字符列:", col, "\n")
    } else if (!is.numeric(feature_data[[col]])) {
      feature_data[[col]] <- as.numeric(feature_data[[col]])
      cat("    转换其他类型列:", col, "\n")
    }
  }
  
  feature_data <- as.data.frame(lapply(feature_data, as.numeric))
  mrmr_data <- data.frame(target = target_numeric, feature_data)
  mrmr_data <- as.data.frame(lapply(mrmr_data, as.numeric))
  
  cat("  - 数据维度:", nrow(mrmr_data), "行 x", ncol(mrmr_data), "列\n")
  cat("  - 所有列的类型:\n")
  print(sapply(mrmr_data, class))
  
  # 转换为 mRMRe 数据对象
  mrmr_data_obj <- tryCatch({
    mRMR.data(data = mrmr_data)
  }, error = function(e) {
    cat("错误详情:", e$message, "\n")
    stop("无法创建 mRMR.data 对象")
  })
  
  cat("\n运行 mRMR 特征选择...\n")
  cat("  - 候选特征数:", ncol(feature_data), "\n")
  cat("  - 目标选择特征数:", n_features, "\n")
  
  mrmr_result <- mRMR.classic(
    data = mrmr_data_obj,
    target_indices = 1,
    feature_count = min(n_features, ncol(feature_data))
  )
  
  # mRMR 有序索引与名称（顺序非常重要，用于 iFS）
  selected_indices <- solutions(mrmr_result)[[1]]
  ordered_features <- colnames(feature_data)[selected_indices]
  
  cat("  - 实际选择特征数:", length(ordered_features), "\n")
  
  feature_scores <- scores(mrmr_result)[[1]]
  feature_importance <- data.frame(
    Feature = ordered_features,
    Index = selected_indices,
    Score = feature_scores,
    Rank = seq_along(ordered_features)
  )
  feature_importance <- feature_importance[order(-feature_importance$Score), ]
  
  list(
    mrmr_result = mrmr_result,
    ordered_features = ordered_features,  # 供 iFS 使用
    n_features = length(ordered_features),
    feature_importance = feature_importance,
    original_feature_names = colnames(feature_data)
  )
}

train_mrmr_rf <- function(data, response_var = "type",
                          mrmr_n_features = 30,
                          rf_ntree = 500,
                          rf_mtry = NULL,
                          use_ifs = TRUE,
                          seed = NULL) {
  
  cat("步骤 1: mRMR 特征选择...\n")
  mrmr_result <- mrmr_feature_selection(
    data = data,
    response_var = response_var,
    n_features = mrmr_n_features
  )
  
  ordered_features <- mrmr_result$ordered_features
  cat("\nmRMR 有序特征（前若干）示例:\n")
  print(utils::head(ordered_features))
  
  if (use_ifs) {
    cat("\n步骤 2: iFS（以 OOB-MCC 选择最优前缀）...\n")
    eval_res <- .incremental_rf_eval_mcc(
      data = data,
      response_var = response_var,
      ordered_features = ordered_features,
      rf_ntree = rf_ntree,
      rf_mtry = rf_mtry,
      verbose = TRUE,
      seed = seed
    )
    best_k <- eval_res$best_k
    best_mcc <- eval_res$best_mcc
    cat(sprintf("  - 最优 k = %d, OOB-MCC = %.4f\n", best_k, best_mcc))
    final_features <- ordered_features[1:best_k]
    curve <- eval_res$curve
  } else {
    cat("\n步骤 2: 不使用 iFS，直接取前 K 个特征训练...\n")
    final_features <- ordered_features
    curve <- data.frame(k = length(final_features), mcc = NA_real_)
  }
  
  # 最终训练
  selected_data <- data[, c(final_features, response_var), drop = FALSE]
  final_mtry <- if (is.null(rf_mtry)) floor(sqrt(ncol(selected_data) - 1)) else rf_mtry
  
  cat("\n步骤 3: 使用最优特征子集训练最终随机森林...\n")
  cat("  - 最终特征数:", length(final_features), "\n")
  cat("  - RF 参数: ntree =", rf_ntree, ", mtry =", final_mtry, "\n")
  
  if (!is.null(seed)) set.seed(seed)
  rf_model <- randomForest(
    as.formula(paste(response_var, "~ .")),
    data = selected_data,
    ntree = rf_ntree,
    mtry = final_mtry,
    importance = TRUE,
    proximity = FALSE,
    oob.prox = FALSE
  )
  
  final_oob_mcc <- .oob_mcc(rf_model)
  cat(sprintf("  - 最终模型 OOB-MCC: %.4f\n", final_oob_mcc))
  
  result <- list(
    mrmr_result = mrmr_result,
    rf_model = rf_model,
    selected_data = selected_data,
    ordered_features = ordered_features,
    final_features = final_features,
    final_k = length(final_features),
    feature_importance_mrmr = mrmr_result$feature_importance,
    ifs_curve = curve,          # k vs OOB-MCC
    ifs_metric = "OOB_MCC",
    final_oob_mcc = final_oob_mcc
  )
  class(result) <- c("mrmr_rf", "list")
  return(result)
}

##################################### 
############## 主程序

dst <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_DST_sample_clean.csv")
X_genomic <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_X_genomic_005.csv")
# X_genomic <- read.csv("~/bioinformatics/ML_mtb/data/processed/01_X_genomic_0001.csv")

## debug 
# i_d <- 13;i_f <- 1;i_imp <- 1
metrics_base_train <- lapply(seq_along(drugs), function(i_d){
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
  # write.csv(indices_save,
  #           paste0("./result/TBportal/indices_fold/",drug_tmp,".csv"),
  #           row.names = F)
  model_fold <- mclapply(1:5, function(i_f){
    
    tryCatch({
      cat("\n",'Fold:',i_f,"\n")
      set_fold <- set_drug[which(set_drug$fold == i_f),]
      set_input <- set_fold[,-c(1:2,ncol(set_fold))]
      cat("\n",'rf model with out FS.')
      rf_model <- randomForest(
        type ~ .,                    # 公式
        data = set_input,           # 训练数据
        ntree = 500,                 # 树的数量
        mtry = sqrt(ncol(set_input)-1),  # 每次分裂考虑的变量数
        importance = TRUE,           # 关键参数:计算变量重要性
        proximity = FALSE,
        oob.prox = FALSE)
      rf_result <- calculate_metrics(rf_model,set_input)
      cat("\n",'begin lasso_rf_model',i_f,"\n")
      lasso_rf_model <- train_lasso_rf(
        data = set_input,
        response_var = "type",
        lasso_alpha = 1,
        lasso_nfolds = 5,
        lambda_selection = "lambda.min",
        use_mcc = TRUE,  # ← 启用 MCC
        rf_ntree = 500,
        rf_mtry = NULL
      )
      lasso_rf_result <- calculate_metrics(lasso_rf_model$rf_model,set_input)
      cat("\n",'begin mrmr_rf_model',i_f,"\n")
      
      mrmr_rf_model <- train_mrmr_rf(
        data = set_input,
        response_var = "type",
        mrmr_n_features = 30,
        rf_ntree = 500,
        rf_mtry = NULL,
        use_ifs = TRUE,
        seed = 2024
      )
      
      mrmr_rf_results <- calculate_metrics(
        model = mrmr_rf_model$rf_model,
        data = set_input,
        response_var = "type",
        model_type = "rf"
      )
      
      re_df <- do.call('rbind',list(rf_result$metrics,
                                    lasso_rf_result$metrics, mrmr_rf_results$metrics))
      re_df$fold <- i_f
      re_df$method <- c("rf","lasso_rf","mrmr_rf")
      re_set <- list(metrics = re_df,
                     rf_model = rf_model,
                     # xgb_model = gb_model,
                     lasso_rf_model = lasso_rf_model$rf_model,
                     mrmr_rf_model = mrmr_rf_model$rf_model)
      return(re_set)
    },
    error = function(e){
      error_msg <- sprintf(
        "\n❌❌❌ Fold %d 出错 ❌❌❌\n错误: %s\n位置: %s\n",
        i_f, 
        e$message,
        deparse(e$call))
    })
    
  }, mc.cores = 5)
  # mcc_kn_fold <- do.call('rbind',lapply(model_fold,function(x){x$metrics}))
  # mcc_kn_fold$drug <- drug_tmp
  # write.csv(mcc_kn_fold, paste0("./result/TBportal/basemodel/",
  #                               drug_tmp,"_5fold.csv"),row.names = F)
  # result_mean <- aggregate(MCC ~ method, data = mcc_kn_fold, FUN = mean)
  
  model_names <- c("rf","lasso_rf","mrmr_rf")
  layout_metric <- mclapply(seq_along(model_names), function(i_mn){
    model_name_tmp <- model_names[i_mn]
    
    cat('Model is:', model_name_tmp, "\n")
    feature_fold <- lapply(1:5, function(i_fmtmp){
      print(i_fmtmp)
      fold_model_tmp <- model_fold[[i_fmtmp]]
      if(any(class(fold_model_tmp) == "character")){
        feature_tmp <- NA
        return(feature_tmp)
      }
      fold_model_tmp <- fold_model_tmp[[i_mn + 1]]
      feature_tmp <- rownames(fold_model_tmp$importance)
      return(feature_tmp)
    })
    feature_m <- table(unlist(feature_fold))
    if(any(feature_m >=3)){
      feature_m <- names(feature_m[which(feature_m >= 3)])
    }else{
      feature_m <- names(feature_m[which(feature_m >= 1)])
    }
    
    set_layout <-set_drug[, match(c(feature_m,'type'), colnames(set_drug))]
    
    fit <- randomForest(
      type ~ .,                    # 公式
      data = set_layout,           # 训练数据
      ntree = 500,                 # 树的数量
      mtry = sqrt(ncol(set_layout)-1),  # 每次分裂考虑的变量数
      importance = TRUE,           # 关键参数:计算变量重要性
      proximity = FALSE,
      oob.prox = FALSE,
      keep.inbag = FALSE)
    
    # fit$y <- NULL                      # 🔥 新增
    # fit$predicted <- NULL              # 🔥 新增
    # fit$votes <- NULL                  # 🔥 新增(如果不需要查看投票)
    
    fit_result <- calculate_metrics(fit,set_layout)
    saveRDS(fit, paste0("./result/models_base/", drug_tmp,
                       "_", model_name_tmp,".rds"))
    re_fit <- fit_result$metrics
    re_fit$method <- model_name_tmp
    return(re_fit)
  },mc.cores = 3)
  layout_metric <- do.call('rbind',layout_metric)
  layout_metric$drug  <- drug_tmp
  return(layout_metric)
})
metrics_base_train <- do.call('rbind',metrics_base_train)
write.csv(metrics_base_train, "./output/metric_train_base.csv", row.names = F)


