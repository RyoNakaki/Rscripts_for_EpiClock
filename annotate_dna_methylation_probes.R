# 必要なパッケージの読み込み
library(dplyr)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
library(vroom)

# getBMをリトライする関数
getBM_with_retry <- function(attributes, filters, values, mart, max_attempts = 10, delay = 5) {
  attempt <- 1
  while (attempt <= max_attempts) {
    tryCatch({
      # getBMの実行
      result <- getBM(attributes = attributes, filters = filters, values = values, mart = mart)
      cat("getBM成功 (試行回数:", attempt, ")\n")
      return(result)
    }, error = function(e) {
      cat("getBM失敗 (試行回数:", attempt, "):", conditionMessage(e), "\n")
      attempt <<- attempt + 1
      Sys.sleep(delay) # 再試行前に遅延を挟む
    })
  }
  stop("getBMが", max_attempts, "回失敗しました。処理を中断します。")
}

# パス情報の指定

# output_dir_path <- "/Users/nakaki/data/annotate_custom_array/20241227/"
# input_data_path <- "/Users/nakaki/data/annotate_custom_array/20241227/beta_matrix.txt"
# ref_custom_array_path <- "/Users/nakaki/data/annotate_custom_array/ref/MethylationEPICv2.0/EPIC-8v2-0_A2.csv"
# ref_ewas_path <- "/Users/nakaki/data/20241202_edit_ewas/EWAS_Atlas_traits.tsv"
# ref_ensembl_mart_path <- "/Users/nakaki/data/20241113_custom_array/refs/ensembl_mart.rds"

output_dir_path <- "/Users/nakaki/data/annotate_custom_array/20250217/"
input_data_path <- "/Users/nakaki/data/annotate_custom_array/20250217/beta_matrix.txt"
ref_custom_array_path <- "/Users/nakaki/data/20241113_custom_array/refs/RhelixaMethylEpiclockV1_SS_20042798X391986_A1_All.csv" # カスタムアレイ用manifestファイル
ref_ewas_path <- "/Users/nakaki/data/20241202_edit_ewas/EWAS_Atlas_traits.tsv"
ref_ensembl_mart_path <- "/Users/nakaki/data/20241113_custom_array/refs/ensembl_mart.rds"

# vroom を使ってデータを読み込む
beta_matrix <- vroom(file = input_data_path, delim = "\t", col_names = TRUE)

# データを data.frame に変換しつつ、1列目を行名として設定
beta_matrix <- as.data.frame(beta_matrix)  # tibble から data.frame に変換
rownames(beta_matrix) <- beta_matrix[[1]]  # 1列目を行名に設定
beta_matrix <- beta_matrix[, -1]  # 1列目をデータから削除

# 1行目の要素が "cg" で始まる列を除外
beta_matrix <- beta_matrix[, !grepl("^cg", beta_matrix[1, ])]

# 「Cheek」という文字列を含む列名の列を除く
beta_matrix <- beta_matrix[, !grepl("Cheek", colnames(beta_matrix))]

# NAを8割以上含む行を除く
beta_matrix <- beta_matrix[rowMeans(is.na(beta_matrix)) < 0.8, ]

# IDが「cg」から始まるプローブ以外は除外
beta_matrix <- beta_matrix[!grepl("^cgBK", rownames(beta_matrix)), ]
beta_matrix <- beta_matrix[!grepl("^ch", rownames(beta_matrix)), ]
beta_matrix <- beta_matrix[!grepl("^ctl", rownames(beta_matrix)), ]

# 行ごとにZ-scoreを計算 (NAを無視する)
z_matrix <- t(apply(beta_matrix, 1, function(x) {
  # 有効な（NAでない）値のみを対象に計算
  valid_values <- x[!is.na(x)]
  
  # 有効な値がない場合、または標準偏差が0の場合は、NAで埋める
  if (length(valid_values) <= 1 || sd(valid_values) == 0) {
    rep(NA, length(x))
  } else {
    # 有効な値がある場合はZ-scoreを計算する
    (x - mean(valid_values)) / sd(valid_values)
  }
}))
colnames(z_matrix) <- colnames(beta_matrix)

# 平均値からの差分を計算 (NAを無視する)
diff_matrix <- t(apply(beta_matrix, 1, function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  ifelse(is.na(x), NA, x - mean_val)
}))
colnames(diff_matrix) <- colnames(beta_matrix)

# カスタムアレイのアノテーションファイルを読み込む
ref_annt <- read.csv(ref_custom_array_path)
colnames(ref_annt)[1] <- "Probe_ID_2"
colnames(ref_annt)[2] <- "Probe_ID"

# "_TC11"で終わるIDを持つプローブを除外
ref_annt <- ref_annt[!grepl("_TC11$", ref_annt$Probe_ID_2), ]

# 同じプローブの行を除外（Probe_IDが同じでProbe_ID_2が異なる行あり)
ref_annt <- ref_annt %>% distinct(Probe_ID, .keep_all = TRUE)

# CHR列の値に "chr" がついていなければ付ける
ref_annt$CHR <- ifelse(grepl("^chr", ref_annt$CHR), ref_annt$CHR, paste0("chr", ref_annt$CHR))

# EWAS-ATLAS情報の取得
ref_ewas <- read.table(file = ref_ewas_path, sep = "\t", header = T, row.names = NULL, quote = "")
colnames(ref_ewas)[1] <- "Probe_ID"

# サンプルごとにデータを出力
for (sample in colnames(beta_matrix)) {
  cat("Sample: ", sample, "\n")
  
  # 各サンプル列をデータフレームに変換
  sample_data <- data.frame(
    Probe_ID = rownames(beta_matrix),
    Beta_value = beta_matrix[, sample],
    Z_score = z_matrix[, sample],
    Delta_beta = diff_matrix[, sample]
  )
  
  # 平均差分の絶対値が0.05以上のプローブのみを選別
  sample_data <- subset(sample_data, abs(Delta_beta) > 0.05)
  
  # 行数が1,000を超えている場合に処理を実行
  if (nrow(sample_data) > 1000) {
    # Z_scoreの絶対値でソートし、上位1,000行を選択
    sample_data <- sample_data[order(abs(sample_data$Z_score), decreasing = TRUE), ][1:1000, ]
  }
  
  # Z-scoreに従いソーティング（降順）
  sample_data <- sample_data %>% arrange(desc(Z_score))
  
  # プローブの位置情報を紐付け
  sample_data <- sample_data %>% left_join(ref_annt, by = "Probe_ID")
  
  # プローブ位置情報を取得（近傍遺伝子情報の付与にて使用）
  ref_pos <- sample_data[, c("CHR", "MAPINFO", "MAPINFO")]
  colnames(ref_pos)[1] <- "Chr"
  colnames(ref_pos)[2] <- "Start"
  colnames(ref_pos)[3] <- "End"
  sample_data$Position <- paste(sample_data$CHR, sample_data$MAPINFO, sample_data$MAPINFO, sep = "_")
  
  # 近傍遺伝子情報の付与
  ref_gene <- build_annotations(genome = "hg38", annotations = c("hg38_basicgenes"))
  ref_pos <- makeGRangesFromDataFrame(ref_pos, na.rm = TRUE)
  ref_gene <- annotate_regions(regions = ref_pos, annotations = ref_gene, ignore.strand = TRUE)
  ref_gene <- data.frame(ref_gene)
  ref_gene$Position <- paste(ref_gene$seqnames, ref_gene$start, ref_gene$end, sep = "_")
  ref_gene <- ref_gene %>% distinct(Position, .keep_all = TRUE)
  ref_gene <- ref_gene[, c("seqnames", "start", "annot.id", "annot.symbol", "Position")]
  
  sample_data <- sample_data %>% left_join(ref_gene, by = "Position")
  
  # 出力列の絞り込み
  sample_data <- sample_data[, c("Probe_ID", "Beta_value", "Z_score", "Delta_beta", "seqnames", "start", "annot.id", "annot.symbol")]
  colnames(sample_data)[colnames(sample_data) == "seqnames"] <- "Chr"
  colnames(sample_data)[colnames(sample_data) == "start"] <- "Position"
  colnames(sample_data)[colnames(sample_data) == "annot.id"] <- "Gene_position"
  colnames(sample_data)[colnames(sample_data) == "annot.symbol"] <- "Gene_symbol"
  
  # 遺伝子名がNAでないプローブのみを選別
  sample_data <- sample_data[!is.na(sample_data$Gene_symbol), ]
  
  # GO情報の読み込み
  ref_go <- readRDS(ref_ensembl_mart_path)
  gene_id <- sample_data$Gene_symbol
  
  gene_descriptions <- getBM_with_retry(
    attributes = c("hgnc_symbol", "description"),
    filters = "hgnc_symbol",
    values = gene_id,
    mart = ref_go
  )
  
  ref_go <- getBM_with_retry(attributes = c("hgnc_symbol", "go_id", "name_1006", "namespace_1003"),
    filters = "hgnc_symbol",
    values = gene_id,
    mart = ref_go
  )
  
  ref_go <- ref_go %>%
    group_by(hgnc_symbol) %>%
    summarise(
      Molecular_Function = paste(name_1006[namespace_1003 == "molecular_function"], collapse = " | "),
      Cellular_Component = paste(name_1006[namespace_1003 == "cellular_component"], collapse = " | "),
      Biological_Process = paste(name_1006[namespace_1003 == "biological_process"], collapse = " | ")
    )
  
  ref_go <- left_join(ref_go, gene_descriptions, by = "hgnc_symbol")
  
  # 'description'列を'hgnc_symbol'列の直後に移動
  ref_go <- ref_go %>%
    relocate(description, .after = hgnc_symbol)
  
  # プローブとGO情報の紐付け
  colnames(ref_go)[1] <- "Gene_symbol"
  sample_data <- sample_data %>% left_join(ref_go, by = "Gene_symbol")
  
  colnames(sample_data)[colnames(sample_data) == "description"] <- "Description"
  colnames(sample_data)[colnames(sample_data) == "Molecular_Function"] <- "GO_molecular_function"
  colnames(sample_data)[colnames(sample_data) == "Cellular_Component"] <- "GO_cellular_component"
  colnames(sample_data)[colnames(sample_data) == "Biological_Process"] <- "GO_biological_process"
  
  # プローブとEWAS-ATLAS情報の紐付け
  sample_data <- sample_data %>% left_join(ref_ewas, by = "Probe_ID")
  colnames(sample_data)[colnames(sample_data) == "traits"] <- "EWAS_traits"
  
  # 空欄をNAで埋める
  sample_data[sample_data == ""] <- NA
  
  # ファイルに出力（CSV形式で保存）
  file_name <- paste0(output_dir_path, sample, ".csv")
  write.csv(sample_data, file = file_name, row.names = FALSE)
}

