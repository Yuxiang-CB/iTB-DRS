#!/usr/bin/env Rscript
# Predict resistance for 13 drugs from a single VCF or a TXT list of paths; only reads trained models.
# Dependencies: install.packages(c("randomForest", "vcfR", "ps"))
# Examples (run from project root):
# Rscript scripts/deploy_predict.R data/raw/Example/ERR2184230.vcf.gz
# Rscript scripts/deploy_predict.R sample.txt --outdir predictions
# Rscript scripts/deploy_predict.R --vcf sample.vcf --output results/pred.csv
# Each line in the TXT is one VCF path; relative paths are relative to the TXT directory; blank lines and # comment lines are ignored.
# By default, create a predictions folder in the current working directory; --outdir can specify the save directory.
# Batch CSV has one row per sample; the last column prediction_time_sec is the VCF reading and prediction time for that VCF.
# Models are loaded once before processing samples; loading time is reported separately and is not included in per-sample time.
#
# Keep training features consistent with 01_vcf_read.R: positions of records with FILTER == PASS are encoded as 1,
# positions absent or failing FILTER are encoded as 0; REF/ALT and GT are not distinguished.
# Input should be a single-sample variant-only VCF and use the same reference genome coordinates as training data.
# Default stable = 3, corresponding to the 3rd model saved by 02_model_train_mafs.R,
# --stable 2 corresponds to model_rds[[4]]$fit used by 03_model_test.R.
# prob_R >= threshold is called R (resistant), otherwise S (susceptible); default threshold is 0.5.
# Probability interface: https://search.r-project.org/CRAN/refmans/randomForest/html/predict.randomForest.html

drug_names <- c(
  EMB = "ethambutol", INH = "isoniazid", PZA = "pyrazinamide",
  RIF = "rifampicin", SM = "streptomycin", AMK = "amikacin",
  CM = "capreomycin", ETO = "ethionamide", KM = "kanamycin",
  LFX = "levofloxacin", MFX = "moxifloxacin", OFX = "ofloxacin",
  PAS = "p_aminosalicylic_acid"
)

project_directory <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script_arg)) return(getwd())
  script_path <- sub("^--file=", "", script_arg[1])
  if (!file.exists(script_path)) return(getwd())
  script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
  if (basename(script_dir) %in% c("script", "scripts")) dirname(script_dir) else script_dir
}

print_usage <- function() {
  cat(paste0(
    "Purpose: predict resistance to 13 drugs from a single VCF or a TXT list of paths.\n\n",
    "Usage: Rscript scripts/deploy_predict.R <sample.vcf[.gz] | sample.txt> [options]\n",
    "       Rscript scripts/deploy_predict.R --input sample.txt --outdir my_predictions\n",
    "       Rscript scripts/deploy_predict.R --vcf sample.vcf --output results/pred.csv\n\n",
    "Options:\n",
    "  --input <file>    Input VCF or TXT; can also be supplied as a positional argument\n",
    "  --vcf <file>      Compatibility alias for --input\n",
    "  --models <dir>    Model directory; defaults to models under the project root\n",
    "  --outdir <dir>    Save directory; defaults to predictions under the current working directory (created automatically)\n",
    "  --output <file>   CSV path; when used with --outdir, specify only the CSV filename\n",
    "  --stable <1-5>    Feature stability; default 3 (the 3rd saved model)\n",
    "  --threshold <0-1> Resistance threshold for prob_R; default 0.5; values equal to the threshold are called R\n",
    "  -h, --help        Show help\n\n",
    "Dependencies: randomForest, vcfR; ps can be used for cross-platform RAM statistics.\n",
    "TXT: one VCF path per line; relative paths are relative to the TXT directory.\n",
    "     Blank lines and comment lines starting with # are ignored; duplicate paths to the same VCF are processed only once.\n",
    "Input: single-sample variant-only VCF; reference genome coordinates must match the training data.\n",
    "Encoding: only positions with PASS records are 1; all others are 0; no filtering by GT.\n",
    "Output: a single VCF keeps one row per drug; a TXT batch CSV has one row per sample.\n",
    "        Default filenames are pred_<VCF name>.csv or pred_<TXT name>_batch.csv.\n",
    "        The last column prediction_time_sec is the time in seconds for VCF reading and prediction of 13 drugs.\n",
    "        Each drug model is loaded only once; loading time is printed separately.\n",
    "R = resistant; S = susceptible. RAM units are MiB, reporting process RSS.\n",
    "On Linux, the system RSS peak is also reported; on other platforms, sampled peaks for each processing stage are reported.\n"
  ))
}

parse_arguments <- function(args, project_dir) {
  opts <- list(input = NULL, models = file.path(project_dir, "models"),
               outdir = NULL, output = NULL, stable = "3", threshold = "0.5")
  option_names <- c("--input", "--vcf", "--models", "--outdir", "--output", "--stable", "--threshold")
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (startsWith(arg, "--") && grepl("=", arg, fixed = TRUE)) {
      key <- sub("=.*$", "", arg)
      value <- sub("^[^=]*=", "", arg)
      if (!key %in% option_names || !nzchar(value)) stop("Invalid argument: ", arg)
    } else if (arg %in% option_names) {
      key <- arg
      i <- i + 1L
      if (i > length(args) || startsWith(args[i], "--")) stop("Missing value for argument: ", key)
      value <- args[i]
    } else if (startsWith(arg, "-")) {
      stop("Unknown argument: ", arg, "; use --help for help.")
    } else {
      key <- "--input"
      value <- arg
    }
    field <- sub("^--", "", key)
    if (field == "vcf") field <- "input"
    if (field == "input" && !is.null(opts$input)) stop("Only one VCF or TXT input file can be specified at a time.")
    opts[[field]] <- value
    i <- i + 1L
  }
  if (is.null(opts$input)) stop("Please specify a VCF or TXT file; use --help for help.")
  stable <- suppressWarnings(as.numeric(opts$stable))
  if (!is.finite(stable) || !stable %in% 1:5) stop("--stable must be an integer from 1 to 5.")
  opts$stable <- as.integer(stable)
  opts$threshold <- suppressWarnings(as.numeric(opts$threshold))
  if (!is.finite(opts$threshold) || opts$threshold < 0 || opts$threshold > 1) {
    stop("--threshold must be a numeric value between 0 and 1.")
  }
  opts$is_batch <- grepl("[.]txt$", opts$input, ignore.case = TRUE)
  if (!is.null(opts$output) && !is.null(opts$outdir)) {
    if (basename(opts$output) != opts$output || grepl("[/\\\\]", opts$output)) {
      stop("When both --outdir and --output are specified, --output should contain only the CSV filename.")
    }
    opts$output <- file.path(opts$outdir, opts$output)
  }
  if (is.null(opts$outdir)) opts$outdir <- file.path(getwd(), "predictions")
  if (is.null(opts$output)) {
    stem <- if (opts$is_batch) sub("[.]txt$", "", basename(opts$input), ignore.case = TRUE) else sample_id(opts$input)
    filename <- paste0("pred_", stem, if (opts$is_batch) "_batch" else "", ".csv")
    opts$output <- file.path(opts$outdir, filename)
  }
  opts
}

sample_id <- function(path) {
  sub("[.]vcf([.]gz)?$", "", basename(path), ignore.case = TRUE)
}

read_samples <- function(input, is_batch) {
  if (!file.exists(input) || dir.exists(input) || file.access(input, 4) != 0) {
    stop("Input file does not exist or is not readable: ", input)
  }
  input <- normalizePath(input, winslash = "/", mustWork = TRUE)
  if (is_batch) {
    lines <- trimws(sub("^\ufeff", "", readLines(input, warn = FALSE, encoding = "UTF-8")))
    keep <- nzchar(lines) & !startsWith(lines, "#")
    paths <- lines[keep]
    line_numbers <- which(keep)
    if (!length(paths)) stop("No VCF paths found in TXT: ", input)
    # Allow quotes around paths; spaces and # inside paths are preserved as-is.
    first <- substr(paths, 1L, 1L)
    last <- substr(paths, nchar(paths), nchar(paths))
    quoted <- nchar(paths) >= 2L & first == last & first %in% c("\"", "'")
    paths[quoted] <- substr(paths[quoted], 2L, nchar(paths[quoted]) - 1L)
    if (any(!nzchar(paths))) stop("Empty file path found in TXT; line numbers: ", paste(line_numbers[!nzchar(paths)], collapse = ", "))
    paths <- path.expand(paths)
    absolute <- startsWith(paths, "/") | startsWith(paths, "\\") | grepl("^[A-Za-z]:", paths)
    paths[!absolute] <- file.path(dirname(input), paths[!absolute])
  } else {
    paths <- input
    line_numbers <- 1L
  }
  missing <- !file.exists(paths) | dir.exists(paths) | file.access(paths, 4) != 0
  if (any(missing)) {
    detail <- if (is_batch) paste0("line ", line_numbers[missing], ": ", paths[missing]) else paths[missing]
    stop("VCF file does not exist or is not readable:\n", paste(detail, collapse = "\n"))
  }
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  if (anyDuplicated(paths)) {
    cat(sprintf("Ignored %d duplicate VCF paths.\n", sum(duplicated(paths))))
    paths <- unique(paths)
  }
  data.frame(sra_id = sample_id(paths), vcf_path = paths, stringsAsFactors = FALSE)
}

# Actual resident memory of the process (RSS), not model object.size() or R memory allocation.
process_memory <- function() {
  rss <- peak <- NA_real_
  if (file.exists("/proc/self/status")) {
    status <- readLines("/proc/self/status", warn = FALSE)
    read_kib <- function(key) {
      line <- grep(paste0("^", key, ":"), status, value = TRUE)
      if (!length(line)) return(NA_real_)
      as.numeric(sub("^[^:]+:[[:space:]]*([0-9]+).*", "\\1", line[1])) / 1024
    }
    rss <- read_kib("VmRSS")
    peak <- read_kib("VmHWM")
  } else if (requireNamespace("ps", quietly = TRUE)) {
    info <- tryCatch(ps::ps_memory_info(), error = function(e) NULL)
    if (!is.null(info)) rss <- as.numeric(info[["rss"]]) / 1024^2
  } else if (.Platform$OS.type == "unix") {
    value <- suppressWarnings(tryCatch(
      system2("ps", c("-o", "rss=", "-p", as.character(Sys.getpid())),
              stdout = TRUE, stderr = FALSE), error = function(e) character()
    ))
    if (length(value)) rss <- suppressWarnings(as.numeric(trimws(value[1]))) / 1024
  }
  c(rss = rss, peak = peak)
}

read_pass_features <- function(path) {
  vcf <- vcfR::read.vcfR(path, verbose = FALSE)
  if (ncol(vcf@gt) > 2L) stop("VCF contains multiple samples; please provide a single-sample VCF.")
  fix <- vcf@fix
  if (length(unique(fix[, "CHROM"])) > 1L) {
    stop("VCF contains multiple contigs; the current model only uses rs<POS> and cannot distinguish different contigs.")
  }
  if (anyNA(fix[, "POS"]) || any(!grepl("^[1-9][0-9]*$", fix[, "POS"]))) {
    stop("VCF contains invalid POS coordinates.")
  }
  # Exactly consistent with training: presence of a PASS position is 1; do not use GT substrings such as 1/1.
  positions <- fix[which(fix[, "FILTER"] == "PASS"), "POS"]
  features <- if (length(positions)) unique(paste0("rs", positions)) else character()
  list(features = features, n_records = nrow(fix))
}

select_model <- function(object, stable, drug) {
  # Compatible with directly saved randomForest objects and stability lists saved by the training script.
  if (inherits(object, "randomForest")) return(object)
  if (!is.list(object)) stop(drug, ": RDS is not a recognized model structure.")
  entries <- if ("fit" %in% names(object)) list(object) else object
  matches <- which(vapply(entries, function(entry) {
    is.list(entry) && length(entry$stable) == 1L && isTRUE(entry$stable == stable)
  }, logical(1)))
  if (length(matches) != 1L) stop(drug, ": cannot find a unique stable = ", stable, " model.")
  fit <- entries[[matches]]$fit
  if (is.null(fit)) stop(drug, ": no available model for stable = ", stable, "; please choose another stability.")
  fit
}

load_model <- function(path, drug, stable, record_memory) {
  object <- readRDS(path)
  record_memory()
  fit <- select_model(object, stable, drug)
  if (!inherits(fit, "randomForest") || !identical(fit$type, "classification") ||
      is.null(fit$forest) || !setequal(fit$classes, c("R", "S"))) {
    stop(drug, ": a randomForest classification model containing forest and R/S classes is required.")
  }
  features <- names(fit$forest$ncat)
  if (is.null(features)) features <- rownames(fit$importance)
  if (!length(features) || anyNA(features) || anyDuplicated(features) ||
      any(!grepl("^rs[0-9]+$", features)) ||
      length(features) != length(fit$forest$ncat) || any(fit$forest$ncat != 1L)) {
    stop(drug, ": model features must be numeric rs<POS> sites.")
  }
  # Numeric sites are provided by newdata; remove the terms reference to the training loop environment,
  # so the cache only holds the selected model, not the training matrix and other stability models in that environment.
  if (!is.null(fit$terms)) attr(fit$terms, ".Environment") <- baseenv()
  list(fit = fit, features = features,
       Model = if (inherits(object, "randomForest")) "randomForest" else paste0("model_S", stable))
}

load_models <- function(paths, stable, record_memory) {
  models <- vector("list", length(drug_names))
  names(models) <- names(drug_names)
  for (i in seq_along(drug_names)) {
    drug <- names(drug_names)[i]
    cat(sprintf("[%02d/%02d] Loading %s model\n", i, length(drug_names), drug))
    models[[i]] <- load_model(paths[i], drug, stable, record_memory)
    invisible(gc(verbose = FALSE))
  }
  models
}

predict_drug <- function(model, drug, pass_features, opts, record_memory) {
  fit <- model$fit
  features <- model$features
  values <- as.integer(features %in% pass_features)
  newdata <- as.data.frame(matrix(values, nrow = 1L,
                                 dimnames = list(opts$sample_id, features)))
  # Select probability by class name to avoid reversed results if R/S column order changes.
  probs <- stats::predict(fit, newdata = newdata, type = "prob")
  if (nrow(probs) != 1L || !all(c("R", "S") %in% colnames(probs)) ||
      any(!is.finite(probs)) || any(probs < 0 | probs > 1) ||
      abs(sum(probs[1, c("R", "S")]) - 1) > 1e-8) {
    stop(drug, ": model did not return valid R/S probabilities.")
  }
  record_memory()
  data.frame(
    sra_id = opts$sample_id, drug = drug, drug_name = unname(drug_names[drug]),
    Model = model$Model,
    prediction = if (probs[1, "R"] >= opts$threshold) "R" else "S",
    prob_R = as.numeric(probs[1, "R"]), prob_S = as.numeric(probs[1, "S"]),
    threshold = opts$threshold, n_features = length(features),
    n_mutated_features = sum(values), stringsAsFactors = FALSE, row.names = NULL
  )
}

predict_sample <- function(sample, models, opts, record_memory) {
  started <- proc.time()[["elapsed"]]
  opts$sample_id <- sample$sra_id
  vcf_info <- read_pass_features(sample$vcf_path)
  record_memory()
  rows <- lapply(names(drug_names), function(drug) {
    predict_drug(models[[drug]], drug, vcf_info$features, opts, record_memory)
  })
  predictions <- do.call(rbind, rows)
  if (opts$is_batch) {
    result <- data.frame(sra_id = sample$sra_id, vcf_path = sample$vcf_path,
                         stable = opts$stable, threshold = opts$threshold,
                         stringsAsFactors = FALSE)
    for (i in seq_along(drug_names)) {
      drug <- names(drug_names)[i]
      result[[drug]] <- predictions$prediction[i]
      result[[paste0(drug, "_prob_R")]] <- predictions$prob_R[i]
      result[[paste0(drug, "_prob_S")]] <- predictions$prob_S[i]
    }
  } else {
    result <- predictions
  }
  # Both output formats place the single-VCF elapsed time in the last column, in seconds.
  elapsed <- round(proc.time()[["elapsed"]] - started, 6L)
  result$prediction_time_sec <- elapsed
  cat(sprintf("  VCF records: %d; PASS positions: %d; sample time: %.3f seconds\n",
              vcf_info$n_records, length(vcf_info$features), elapsed))
  result
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("-h", "--help"))) {
    print_usage()
    return(invisible(0L))
  }
  started <- proc.time()
  sampled_peak <- NA_real_
  record_memory <- function() {
    memory <- process_memory()
    if (is.finite(memory["rss"])) {
      sampled_peak <<- max(c(sampled_peak, memory["rss"]), na.rm = TRUE)
    }
    memory
  }
  on.exit({
    memory <- record_memory()
    elapsed <- proc.time() - started
    cat(sprintf("\nTotal time (including loading, reading, prediction, and output): %.3f seconds\n", elapsed[["elapsed"]]))
    cat(sprintf("CPU time: %.3f seconds\n", elapsed[["user.self"]] + elapsed[["sys.self"]]))
    if (is.finite(memory["rss"])) {
      cat(sprintf("RAM at end (process RSS): %.2f MiB\n", memory[["rss"]]))
      if (is.finite(memory["peak"])) {
        cat(sprintf("Peak RAM (system RSS peak): %.2f MiB\n", memory[["peak"]]))
      } else {
        cat(sprintf("Sampled peak RAM (RSS at each stage, not exact system peak): %.2f MiB\n", sampled_peak))
      }
    } else {
      cat("RAM (RSS): unavailable; install the ps package to enable cross-platform process memory statistics.\n")
    }
  }, add = TRUE)
  record_memory()
  opts <- parse_arguments(args, project_directory())
  samples <- read_samples(opts$input, opts$is_batch)
  model_paths <- file.path(opts$models, paste0(names(drug_names), ".rds"))
  missing <- !file.exists(model_paths) | dir.exists(model_paths) | file.access(model_paths, 4) != 0
  if (any(missing)) stop("Missing or unreadable models: ", paste(model_paths[missing], collapse = ", "))
  if (dir.exists(opts$output)) stop("Output path is a directory; please specify a CSV file: ", opts$output)
  protected_paths <- c(opts$input, samples$vcf_path, model_paths)
  if (file.exists(opts$output) && normalizePath(opts$output) %in% normalizePath(protected_paths)) {
    stop("Output CSV cannot overwrite input TXT, VCF, or model files.")
  }
  packages <- c("randomForest", "vcfR")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "; please run install.packages(c(\"randomForest\", \"vcfR\", \"ps\")).")
  }
  cat("Input file: ", normalizePath(opts$input), "\n", sep = "")
  cat("Model directory: ", normalizePath(opts$models), "\n", sep = "")
  cat("Output CSV: ", opts$output, "\n", sep = "")
  cat(sprintf("Stability: %d; resistance threshold: %.4f\n", opts$stable, opts$threshold))
  loading_started <- proc.time()[["elapsed"]]
  models <- load_models(model_paths, opts$stable, record_memory)
  loading_elapsed <- proc.time()[["elapsed"]] - loading_started
  cat(sprintf("Model loading complete: 13 drugs; loading time: %.3f seconds (loaded only once).\n", loading_elapsed))
  cat(sprintf("Number of samples to predict: %d\n", nrow(samples)))
  results <- vector("list", nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, , drop = FALSE]
    cat(sprintf("[%d/%d] Predicting sample %s\n", i, nrow(samples), sample$sra_id))
    results[[i]] <- predict_sample(sample, models, opts, record_memory)
    invisible(gc(verbose = FALSE))
  }
  predictions <- do.call(rbind, results)
  output_dir <- dirname(opts$output)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", output_dir)
  }
  utils::write.csv(predictions, opts$output, row.names = FALSE)
  cat("\nPrediction results (R = resistant; S = susceptible):\n")
  if (opts$is_batch) {
    columns <- c("sra_id", names(drug_names), "prediction_time_sec")
    print(utils::head(predictions[, columns], 10L), row.names = FALSE)
    if (nrow(predictions) > 10L) cat("Only the first 10 samples are shown; see the CSV for complete results.\n")
  } else {
    print(predictions[, c("drug", "prediction", "prob_R", "prob_S")], row.names = FALSE)
  }
  cat("\nResults saved to: ", normalizePath(opts$output), "\n", sep = "")
  invisible(0L)
}

if (sys.nframe() == 0L) {
  status <- tryCatch(main(), error = function(e) {
    cat("Error: ", conditionMessage(e), "\n", file = stderr(), sep = "")
    1L
  })
  quit(save = "no", status = status)
}
