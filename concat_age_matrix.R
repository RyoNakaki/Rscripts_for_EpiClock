# ライブラリの読み込み
if (!require("data.table")) {
  install.packages("data.table")
  library(data.table)
}

# 絶対パスで指定するCSVファイルのリスト
csv_files <- c(
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/EPAS0001/age_matrix_EPAS0001.csv",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/EPAS0002/age_matrix_EPAS0002.csv",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/EPAS0003/age_matrix_EPAS0003.csv",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/EPAS0004/age_matrix_EPAS0004.csv",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/EPAS0005/age_matrix_EPAS0005.csv",
  "/Users/nakaki/Library/CloudStorage/GoogleDrive-nakaki@rhelixa.com/マイドライブ/Data/epiclock/custom/EPAS0006/age_matrix_EPAS0006.csv"
)

# 出力ファイルのパスを指定
# output_file <- "/Users/nakaki/Analysis/epiclock/concat_age_matrix/custom/age_matrix.csv"
output_file <- "/Users/nakaki/Analysis/epiclock/concat_age_matrix/EPIC/age_matrix.csv"

# 読み込んだデータを格納するリスト
data_list <- list()

# 各ファイルを読み込み、列名に従って並べ替え
for (file in csv_files) {
  cat("読み込み中：", file, "\n")
  # データの読み込み
  temp_data <- fread(file)
  
  # 最初のファイルの列名を基準にする
  if (length(data_list) == 0) {
    base_col_names <- colnames(temp_data)
  } else {
    # 基準列名にない列は除外し、順番を揃える
    temp_data <- temp_data[, ..base_col_names, with = FALSE]
  }
  
  # === バリデーション: Sex__c 列のチェック ===
  # Sex__c 列が存在する場合
  if ("Sex__c" %in% colnames(temp_data)) {
    # Male または Female 以外の要素を検出
    invalid_values <- unique(temp_data$Sex__c[!(temp_data$Sex__c %in% c("Male", "Female"))])
    
    # 不正な値が存在する場合
    if (length(invalid_values) > 0) {
      stop(paste(
        "エラー: ファイル", file, "に 'Sex__c' 列に不正な値があります：",
        paste(invalid_values, collapse = ", ")
      ))
    }
  } else {
    # Sex__c 列が存在しない場合もエラー
    stop(paste("エラー: ファイル", file, "に 'Sex__c' 列が存在しません。"))
  }
  
  # データをリストに追加
  data_list[[length(data_list) + 1]] <- temp_data
}

# リスト内のデータを結合
final_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

# データの確認
cat("結合後のデータプレビュー：\n")
print(head(final_data))

# 結合したデータをCSVとして出力
fwrite(final_data, output_file)

# 終了メッセージ
cat("結合したCSVファイルが", output_file, "として保存されました。\n")