setwd("~/bioinformatics/iTBDRS/")

drugs <- c("EMB", "INH","PZA","RIF","SM", "AMK","CM",
           "ETO","KM","LFX","MFX","OFX","PAS")


library(tidyverse)



get_genomic_coordinates <- function(gff_file_path, variant_list) {
  
  message("正在读取 GFF 文件...")
  
  # --- 1. 读取 GFF ---
  gff_data <- read_tsv(gff_file_path, comment = "#", col_names = FALSE, show_col_types = FALSE) %>%
    select(seqid = X1, type = X3, start = X4, end = X5, strand = X7, attributes = X9)
  
  # --- 2. 提取所有基因特征信息 ---
  gff_features <- gff_data %>%
    mutate(
      gene_name = str_extract(attributes, "(?<=gene=)[^;]+"),
      locus_tag = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
      parent_gene = str_extract(attributes, "(?<=Parent=gene-)[^;]+"),
      name = str_extract(attributes, "(?<=Name=)[^;]+"),
      product = str_extract(attributes, "(?<=product=)[^;]+"),
      
      gene_name = coalesce(gene_name, locus_tag, parent_gene, name),
      locus_tag_backup = locus_tag
    ) %>%
    filter(!is.na(gene_name)) %>%
    select(type, gene_name, locus_tag_backup, product, start, end, strand)
  
  # --- 3. 为每个基因创建完整的特征表 ---
  gene_regions <- gff_features %>%
    group_by(gene_name) %>%
    reframe(
      locus_tag = first(locus_tag_backup),
      gene_start = min(start, na.rm = TRUE),
      gene_end = max(end, na.rm = TRUE),
      strand = first(strand),
      
      has_cds = any(type == "CDS"),
      
      cds_start = {
        cds_rows <- start[type == "CDS"]
        if (length(cds_rows) > 0) min(cds_rows) else NA_real_
      },
      cds_end = {
        cds_rows <- end[type == "CDS"]
        if (length(cds_rows) > 0) max(cds_rows) else NA_real_
      },
      
      has_rRNA = any(type == "rRNA"),
      has_tRNA = any(type == "tRNA"),
      
      product_name = first(na.omit(product))
    ) %>%
    filter(is.finite(gene_start) & is.finite(gene_end)) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  # --- 4. 解析变体 ---
  variants_df <- tibble(input_variant = variant_list) %>%
    filter(!is.na(input_variant), input_variant != "-", input_variant != "NULL") %>%
    mutate(
      is_position = str_detect(input_variant, "^rs\\d+$"),
      
      gene_name = if_else(
        !is_position,
        str_extract(input_variant, "^[A-Za-z0-9]+"),
        NA_character_
      ),
      
      query_position = if_else(
        is_position,
        as.numeric(str_extract(input_variant, "(?<=^rs)\\d+")),
        NA_real_
      ),
      
      var_type = case_when(
        is_position ~ "position",
        str_detect(input_variant, "_p\\.") ~ "protein",
        str_detect(input_variant, "_c\\.") ~ "dna",
        str_detect(input_variant, "_n\\.") ~ "rna",
        str_detect(input_variant, "_r\\.") ~ "rna",
        TRUE ~ "unknown"
      ),
      
      temp_suffix = str_remove(input_variant, "^.*_[cpnr]\\."),
      extracted_num = as.numeric(str_extract(temp_suffix, "-?\\d+"))
    )
  
  # --- 5. 处理位置标记 ---
  position_variants <- variants_df %>%
    filter(is_position)
  
  if (nrow(position_variants) > 0) {
    position_results <- position_variants %>%
      rowwise() %>%
      mutate(
        matched_info = list({
          matches <- gene_regions %>%
            filter(gene_start <= query_position & gene_end >= query_position)
          
          if (nrow(matches) > 0) {
            match <- matches[1, ]
            
            # ✅ 计算基因内相对位置
            position_in_gene <- case_when(
              match$strand == "+" ~ query_position - match$gene_start + 1,
              match$strand == "-" ~ match$gene_end - query_position + 1,
              TRUE ~ NA_real_
            )
            
            region_type <- case_when(
              !is.na(match$cds_start) & !is.na(match$cds_end) & 
                query_position >= match$cds_start & query_position <= match$cds_end ~ {
                  
                  rel_pos <- case_when(
                    match$strand == "+" ~ query_position - match$cds_start + 1,
                    match$strand == "-" ~ match$cds_end - query_position + 1,
                    TRUE ~ NA_real_
                  )
                  
                  cds_length <- match$cds_end - match$cds_start + 1
                  
                  case_when(
                    is.na(rel_pos) ~ "CDS_coding",
                    rel_pos <= 30 ~ "CDS_5prime",
                    rel_pos > cds_length - 30 ~ "CDS_3prime",
                    TRUE ~ "CDS_coding"
                  )
                },
              
              match$has_rRNA ~ "rRNA",
              match$has_tRNA ~ "tRNA",
              
              !is.na(match$cds_start) & !is.na(match$cds_end) ~ {
                case_when(
                  match$strand == "+" & query_position < match$cds_start ~ "5prime_UTR",
                  match$strand == "+" & query_position > match$cds_end ~ "3prime_UTR",
                  match$strand == "-" & query_position > match$cds_end ~ "5prime_UTR",
                  match$strand == "-" & query_position < match$cds_start ~ "3prime_UTR",
                  TRUE ~ "gene_region"
                )
              },
              
              TRUE ~ "gene_region"
            )
            
            list(
              gene_name = match$gene_name,
              locus_tag = match$locus_tag,
              product_name = match$product_name,
              strand = match$strand,
              region_type = region_type,
              position = position_in_gene  # ✅ 基因内相对位置
            )
          } else {
            list(
              gene_name = NA_character_,
              locus_tag = NA_character_,
              product_name = NA_character_,
              strand = NA_character_,
              region_type = "intergenic",
              position = NA_real_  # ✅ 基因间区没有相对位置
            )
          }
        })
      ) %>%
      ungroup() %>%
      mutate(
        gene_name = map_chr(matched_info, "gene_name"),
        locus_tag = map_chr(matched_info, "locus_tag"),
        product_name = map_chr(matched_info, "product_name"),
        strand = map_chr(matched_info, "strand"),
        region_type = map_chr(matched_info, "region_type"),
        position = map_dbl(matched_info, "position"),  # ✅ 提取相对位置
        position_start = query_position,
        position_end = query_position
      ) %>%
      select(input_variant, gene_name, locus_tag, product_name, region_type, strand, var_type, position, position_start, position_end)
  } else {
    position_results <- tibble()
  }
  
  # --- 6. 处理 HGVS 格式变体 ---
  hgvs_variants <- variants_df %>%
    filter(!is_position)
  
  if (nrow(hgvs_variants) > 0) {
    calc_df <- hgvs_variants %>%
      left_join(gene_regions, by = "gene_name") %>%
      left_join(
        gene_regions %>% select(locus_tag, gene_start_lt = gene_start, gene_end_lt = gene_end, 
                                strand_lt = strand, product_name_lt = product_name,
                                cds_start_lt = cds_start, cds_end_lt = cds_end),
        by = c("gene_name" = "locus_tag")
      ) %>%
      mutate(
        gene_start = coalesce(gene_start, gene_start_lt),
        gene_end = coalesce(gene_end, gene_end_lt),
        strand = coalesce(strand, strand_lt),
        product_name = coalesce(product_name, product_name_lt),
        cds_start = coalesce(cds_start, cds_start_lt),
        cds_end = coalesce(cds_end, cds_end_lt)
      ) %>%
      select(-ends_with("_lt")) %>%
      mutate(
        # ✅ position 就是 HGVS 中的相对位置(extracted_num)
        position = extracted_num,
        
        rel_start = case_when(
          var_type %in% c("dna", "rna") ~ extracted_num,
          var_type == "protein" ~ (extracted_num * 3) - 2,
          TRUE ~ NA_real_
        ),
        rel_end = case_when(
          var_type %in% c("dna", "rna") ~ extracted_num,
          var_type == "protein" ~ extracted_num * 3,
          TRUE ~ NA_real_
        )
      )
    
    calc_abs_pos_vectorized <- function(s, idx, g_start, g_end) {
      case_when(
        is.na(g_start) ~ NA_real_,
        s == "+" & idx > 0 ~ g_start + idx - 1,
        s == "+" & idx < 0 ~ g_start + idx,
        s == "-" & idx > 0 ~ g_end - idx + 1,
        s == "-" & idx < 0 ~ g_end - idx,
        TRUE ~ NA_real_
      )
    }
    
    hgvs_results <- calc_df %>%
      mutate(
        pos_1 = calc_abs_pos_vectorized(strand, rel_start, cds_start, cds_end),
        pos_2 = calc_abs_pos_vectorized(strand, rel_end, cds_start, cds_end),
        position_start = pmin(pos_1, pos_2, na.rm = TRUE),
        position_end = pmax(pos_1, pos_2, na.rm = TRUE),
        
        region_type = case_when(
          var_type == "protein" ~ "CDS_coding",
          var_type %in% c("dna", "rna") & !is.na(cds_start) & !is.na(cds_end) & !is.na(rel_start) ~ {
            cds_length <- cds_end - cds_start + 1
            case_when(
              rel_start <= 30 ~ "CDS_5prime",
              rel_start > cds_length - 30 ~ "CDS_3prime",
              TRUE ~ "CDS_coding"
            )
          },
          TRUE ~ "unknown"
        )
      ) %>%
      select(input_variant, gene_name, locus_tag, product_name, region_type, strand, var_type, position, position_start, position_end)
  } else {
    hgvs_results <- tibble()
  }
  
  # --- 7. 合并结果 ---
  result <- bind_rows(hgvs_results, position_results)
  
  return(result)
}



# --- 再次运行你的数据 ---
# my_variants <- mu_imp_pos$mutation[1:3]
path_gff <- "./data/note/GCF_000195955.2_ASM19595v2_genomic.gff"

INFO <- lapply(seq_along(drugs), function(i_d){
  drug_tmp <- drugs[i_d]
  cat('\nDrug:', drug_tmp, "\n")
  fit_rds <- readRDS(paste0("./result/models_005/",drug_tmp,".rds")) 
  fit <- fit_rds[[3]]$fit
  mutation <- rownames(fit$importance)
  info_mutation <- get_genomic_coordinates(path_gff, mutation)
  re_info <- data.frame(Drug = drug_tmp, info_mutation)
  return(re_info)
})
INFO <- do.call('rbind',INFO)

write.csv(INFO,"./output/mutation_INFO.csv", row.names = F)

