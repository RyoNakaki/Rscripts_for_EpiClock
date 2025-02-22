# 必要なライブラリの読み込み
library(dplyr)

# 結合対象のタブ区切りテキストファイルのパスリストを指定
file_list <- c(
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/raw/EPAS0001/beta_matrix_EPAS0001.txt",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/raw/EPAS0002/beta_matrix_EPAS0002.txt",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/raw/EPAS0003/beta_matrix_EPAS0003.txt",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/raw/EPAS0004/beta_matrix_EPAS0004.txt",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/raw/EPAS0005/beta_matrix_EPAS0005.txt",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/raw/EPAS0006/beta_matrix_EPAS0006.txt"
)

# 出力ファイルのパスを指定
output_file <- "/Users/nakaki/Analysis/epiclock/custom/concat_beta_matrix/beta_matrix.txt"

# データを格納するリストを作成
data_list <- list()

# 各ファイルを読み込み、行名の順番を揃えてリストに格納
for (file in file_list) {
  # ファイルを読み込み、行名を rownames に設定
  data <- read.delim(file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
  
  # 最初のファイルの行名を基準とする
  if (length(data_list) == 0) {
    base_row_names <- rownames(data)
  } else {
    # 基準の行名に合わせて順番を並び替え（不足行は NA で補完）
    data <- data[match(base_row_names, rownames(data)), ]
  }
  
  # データフレームとしてリストに追加
  data_list[[length(data_list) + 1]] <- data
}

# 列方向に結合
combined_data <- do.call(cbind, data_list)

# 行名を保持
rownames(combined_data) <- base_row_names

# 行名を新しい列 "Probe_ID" として追加
combined_data <- data.frame(Probe_ID = rownames(combined_data), combined_data, check.names = FALSE)

# タブ区切りテキストとして出力
write.table(combined_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("ファイル結合が完了しました。出力ファイル:", output_file, "\n")