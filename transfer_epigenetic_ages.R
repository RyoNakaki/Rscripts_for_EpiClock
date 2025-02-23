# 必要なパッケージの読み込み
library(ggplot2)
library(glmnet)

########################################
# 関数を定義
########################################

# データフレームの特定列の値を正規分布に従うことを仮定し変換し、新しい列として追加
add_transformed_column <- function(df, column_name, new_column_name, mean_ref, sd_ref) {
  # 平均と標準偏差を計算して変数に格納
  mean <- mean(df[,column_name], na.rm = TRUE)  # 平均
  sd <- sd(df[,column_name], na.rm = TRUE)     # 標準偏差
  
  # 結果の表示
  # cat("平均:", mean, "\n")
  # cat("標準偏差:", sd, "\n")
  
  # mean_dbを平均、sd_dbを標準偏差とする正規分布に当てはめ変換して、データフレームに新しい列として追加
  df[[new_column_name]] <- (df[,column_name] - mean) / sd * sd_ref + mean_ref
  
  return(df)
}

# スキャッタープロットを作成
make_scatter_plot <- function(output_file_path, df, column_name_x, column_name_y, output_name_x, output_name_y, mean_ref, sd_ref) {
  if (is.null(column_name_x)) { # x軸の指定がない場合はサンプル番号を要素としたベクトルを代入
    x <- seq_len(nrow(df))
  } else {
    x <- df[[column_name_x]]
  }
  
  y <- df[[column_name_y]]
  
  p <- ggplot(df, aes(x = x, y = y, color = Highlight)) +
    geom_point(size = 2) +                      # 点を描画
    geom_abline(intercept = mean_ref, slope = 0,        # 平均線を追加
                linetype = "dashed", color = "blue", size = 1) +
    geom_abline(intercept = mean_ref + sd_ref, slope = 0,        # +1σ線を追加
                linetype = "dashed", color = "#999999", size = 1) +
    geom_abline(intercept = mean_ref - sd_ref, slope = 0,        # -1σ線を追加
                linetype = "dashed", color = "#999999", size = 1) +
    scale_color_manual(
      values = c("Highlight" = "#EF7f00", "Other" = "gray")  # ハイライトとそれ以外の色を指定
    ) +
    coord_cartesian(ylim = c(mean_ref - sd_ref * 5.0, mean_ref + sd_ref * 5.0)) +        # 縦軸の範囲を設定
    labs(
      x = output_name_x,
      y = output_name_y,
      color = "Sample Type"
    ) +
    theme_minimal()
  
  # ファイルに保存
  ggsave(output_file_path, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
}

# 数値変換とスキャッタープロット作成を続けて実施
transform_and_plot <- function(output_file_path, df, column_name_x, column_name_y, output_name_x, output_name_y, mean_ref, sd_ref) {
  # ①  数値変換
  df <- add_transformed_column(
    df <- df,
    column_name <- column_name_y,
    new_column_name <- gsub("__c", "_t__c", column_name_y),
    mean_ref <- mean_ref,
    sd_ref <- sd_ref
  )
  
  # ② スキャッタープロット作成
  make_scatter_plot(
    output_file_path <- output_file_path,
    df <- df,
    column_name_x <- column_name_x,
    column_name_y <- gsub("__c", "_t__c", column_name_y),
    output_name_x <- output_name_x,
    output_name_y <- output_name_y,
    mean_ref <- mean_ref,
    sd_ref <- sd_ref
  )
  
  return(df)
}

# ElasticNetにより近似直線の傾きと切片を取得
run_elasticnet <- function(df, column_name_x, column_name_y) {
  # xとyを行列形式に変換
  x <- as.matrix(cbind(df[[column_name_x]], rep(1, nrow(df)))) # 説明変数
  y <- df[[column_name_y]] # 目的変数
  
  # alphaの候補リスト
  alpha_values <- seq(0, 1, by = 0.1)
  
  # グリッドサーチ用の結果格納
  results <- data.frame(alpha = numeric(), lambda = numeric(), mse = numeric())
  
  # グリッドサーチ
  for (alpha in alpha_values) {
    # Elastic Netの交差検証
    cv_fit <- cv.glmnet(x, y, alpha = alpha)
    
    # 最適なλとそのときのMSEを記録
    best_lambda <- cv_fit$lambda.min
    best_mse <- min(cv_fit$cvm)
    
    results <- rbind(results, data.frame(alpha = alpha, lambda = best_lambda, mse = best_mse))
  }
  
  # 最適なalphaを決定
  best_result <- results[which.min(results$mse), ]
  best_alpha <- best_result$alpha
  best_lambda <- best_result$lambda
  
  cat("最適なalpha:", best_alpha, "\n")
  cat("最適なlambda:", best_lambda, "\n")
  
  # 最適なalphaでElastic Netモデルを作成
  best_fit <- glmnet(x, y, alpha = best_alpha, lambda = best_lambda)
  
  return(best_fit)
}

# 線形に相関する要素間のスキャッタープロットの作成
transform_and_plot_linear <- function(output_file_path, df, column_name_x, column_name_y, output_name_x, output_name_y, fit, norm) {
  x <- as.matrix(cbind(df[[column_name_x]], rep(1, nrow(df))))
  
  # 近似直線の計算
  y_pred <- predict(fit, newx = x)
  
  # 近似直線をデータフレームに格納
  df[["y_pred"]] <- y_pred[, 1]  # 予測値を新しい列として追加
  
  # 係数を取得
  coef <- coef(fit)
  intercept <- coef[1]  # 切片
  slope <- coef[2]      # 傾き
  
  # ggplot2でスキャッタープロットと近似直線を作成（変換前）
  p <- ggplot(df, aes(x = df[, column_name_x], y = df[, column_name_y], color = Highlight)) +
    geom_point(size = 2) +                       # 元のデータポイント
    geom_line(aes(y = y_pred), linetype = "dashed", color = "blue", size = 1) +        # 近似直線
    labs(
      x = output_name_x,
      y = output_name_y,
      color = "Group"
    ) +
    scale_color_manual(values = c("Highlight" = "#EF7f00", "Other" = "gray")) +
    annotate("text", x = min(df[, column_name_x]), y = max(df[, column_name_y]),
             label = paste0("Intercept: ", round(intercept, 3), "\nSlope: ", round(slope, 3)),
             hjust = 0, vjust = 1, size = 4, color = "black") +  # 切片と傾き情報を追加
    theme_minimal()
  
  # プロットを表示
  # print(p)
  
  # ファイルに保存
  ggsave(gsub(".png", "_pre.png", output_file_path), plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  
  # EpiClockAgeの標準からの差分を暦年齢に加算
  if(norm == 1) {
    diff <- df[[column_name_y]] - (slope * df[[column_name_x]] + intercept)
    
    # 平均と標準偏差を計算して変数に格納
    mean <- mean(diff, na.rm = TRUE)  # 平均
    sd <- sd(diff, na.rm = TRUE)     # 標準偏差
    
    # mean_dbを平均、sd_dbを標準偏差とする正規分布に当てはめ変換して、データフレームに新しい列として追加
    df[[gsub("__c", "_t__c", column_name_y)]] <- (diff - mean) / sd
    
    # ggplot2でスキャッタープロットと近似直線を作成 (変換後)
    p <- ggplot(df, aes(x = df[, column_name_x], y = df[, gsub("__c", "_t__c", column_name_y)], color = Highlight)) +
      geom_point(size = 2) +                      # 点を描画
      geom_abline(intercept = 0.0, slope = 0,        # 平均線を追加
                  linetype = "dashed", color = "blue", size = 1) +
      geom_abline(intercept = 1.0, slope = 0,        # +1σ線を追加
                  linetype = "dashed", color = "#999999", size = 1) +
      geom_abline(intercept = -1.0, slope = 0,        # -1σ線を追加
                  linetype = "dashed", color = "#999999", size = 1) +
      scale_color_manual(
        values = c("Highlight" = "#EF7f00", "Other" = "gray")  # ハイライトとそれ以外の色を指定
      ) +
      coord_cartesian(ylim = c(-5.0, 5.0)) +        # 縦軸の範囲を設定
      labs(
        x = output_name_x,
        y = output_name_y,
        color = "Sample Type"
      ) +
      theme_minimal()
  } else{
    df[[gsub("__c", "_t__c", column_name_y)]] <- df[[column_name_y]] - (slope * df[[column_name_x]] + intercept) + df[[column_name_x]]
    
    # ggplot2でスキャッタープロットと近似直線を作成 (変換後)
    p <- ggplot(df, aes(x = df[, column_name_x], y = df[, gsub("__c", "_t__c", column_name_y)], color = Highlight)) +
      geom_point(size = 2) +                       # 元のデータポイント
      geom_abline(intercept = 0, slope = 1,        # 切片0、傾き1の直線を追加
                  linetype = "dashed", color = "blue", size = 1) +
      labs(
        x = output_name_x,
        y = output_name_y,
        color = "Group"
      ) +
      scale_color_manual(values = c("Highlight" = "#EF7f00", "Other" = "gray")) +
      theme_minimal()
  }
  
  # プロットを表示
  # print(p)
  
  # ファイルに保存
  ggsave(output_file_path, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  
  return(df)
}


########################################
# データの入出力先の指定と読み込み
########################################

# パス情報の指定
job_id <- "20250223_raw"
output_dir_path <- "/Users/nakaki/Analysis/epiclock/transfer_epigenetic_ages/custom/output/"
input_data_path <- "/Users/nakaki/Analysis/epiclock/concat_age_matrix/custom/age_matrix.csv"
input_list_path <- "/Users/nakaki/Analysis/epiclock/transfer_epigenetic_ages/sample_llist.csv"


# インプットデータを読み込む
df <- read.csv(input_data_path)
sample_list <- read.csv(input_list_path)

# リストに該当する要素をラベリング
df[["Highlight"]] <- ifelse(df[["Name"]] %in% sample_list[["Name"]], "Highlight", "Other")


########################################
# 結果の変換と出力（男女差なし）
########################################

# Hoverathクロック結果の変換と出力
#df <- transform_and_plot_linear(
#  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_HorvathClock.png"),
#  df <- df,
#  column_name_x <- "ChronologicalAge__c",
#  column_name_y <- "DNAmHorvathClock__c",
#  output_name_x <- "Chronological Age",
#  output_name_y <- "Horvath Clock",
#  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmHorvathClock__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
#  norm <- 0
#)

# PC Hoverathクロック結果の変換と出力
#df <- transform_and_plot_linear(
#  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_PCHorvathClock.png"),
#  df <- df,
#  column_name_x <- "ChronologicalAge__c",
#  column_name_y <- "DNAmPCHorvathClock__c",
#  output_name_x <- "Chronological Age",
#  output_name_y <- "Horvath Clock",
#  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmPCHorvathClock__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
#  norm <- 0
#)

# EpiClockAge結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_EpiClockAge.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmEpiclockAge__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "EpiClockAge",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmEpiclockAge__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 0
  # norm <- 1
)

# FitAge結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_FitAge.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmFitAge__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "FitAge",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmFitAge__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# ADM結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmADM.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmADM__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "ADM",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmADM__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# B2M結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmB2M.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmB2M__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "B2M",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmB2M__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# CystatinC結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmCystatinC.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmCystatinC__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "CystatinC",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmCystatinC__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# CystatinC結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmGDF15.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmGDF15__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "GDF15",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmGDF15__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# TIMP1結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmTIMP1.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmTIMP1__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "TIMP1",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmTIMP1__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# FRS結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmFRS.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmFRS__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "FRS",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmFRS__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# PACKYRS結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmPACKYRS.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmPACKYRS__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "PACKYRS",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmPACKYRS__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# Hb1AC結果の変換と出力
df <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmlogA1C.png"),
  df <- df,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmlogA1C__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "logA1C",
  fit <- run_elasticnet(df, "ChronologicalAge__c", "DNAmlogA1C__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)


# クロックリストの作成
clock_list <- c(
  "DNAmDunedinPACE",
  "DNAmLeptin",
  "DNAmPAI1",
  "DNAmlogCRP",
  "DNAmHGF",
  "DNAmIL6",
  "DNAmLE_Tibia",
  "DNAmLE_Patella",
  "DNAmSmoking_Philibert",
  "DNAmAlcohol_McCartney"
)

# 平均値リストの作成
mean_list <- c(
  1.0, # 0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
)

# 標準偏差リストの作成
sd_list <- c(
  0.1, # 1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0
)

# 全サンプルを対象に計算
for (i in seq_along(clock_list)) {
  df <- transform_and_plot(
    output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_", clock_list[i], ".png"),
    df <- df,
    column_name_x <- "ChronologicalAge__c",
    column_name_y <- paste0(clock_list[i], "__c"),
    output_name_x <- "Chronological Age",
    output_name_y <- paste0(clock_list[i], " (transformed)"),
    mean_ref <- mean_list[i],
    sd_ref <- sd_list[i]
  )
}


########################################
# 結果の変換と出力（男女差あり）
########################################

# 男女でデータフレームの分割
df_male <- df[df[["Sex__c"]] == "Male", ]  # 男性に該当するサンプルを抽出
df_female <- df[df[["Sex__c"]] == "Female", ]  # 女性に該当するサンプルを抽出

# Gait結果の変換と出力（男性）
df_male <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmGait_Male.png"),
  df <- df_male,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmGait__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "Gait",
  fit <- run_elasticnet(df_male, "ChronologicalAge__c", "DNAmGait__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# Grip結果の変換と出力（男性）
df_male <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmGrip_Male.png"),
  df <- df_male,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmGrip__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "Grip",
  fit <- run_elasticnet(df_male, "ChronologicalAge__c", "DNAmGrip__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# VO2結果の変換と出力（男性）
df_male <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmVO2_Male.png"),
  df <- df_male,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmVO2__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "VO2",
  fit <- run_elasticnet(df_male, "ChronologicalAge__c", "DNAmVO2__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# FEV1結果の変換と出力（男性）
df_male <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmFEV1_Male.png"),
  df <- df_male,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmFEV1__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "FEV1",
  fit <- run_elasticnet(df_male, "ChronologicalAge__c", "DNAmFEV1__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# TL結果の変換と出力（男性）
df_male <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmTL_Male.png"),
  df <- df_male,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmTL__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "TL",
  fit <- run_elasticnet(df_male, "ChronologicalAge__c", "DNAmTL__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# Gait結果の変換と出力（女性）
df_female <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmGait_Female.png"),
  df <- df_female,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmGait__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "Gait",
  fit <- run_elasticnet(df_female, "ChronologicalAge__c", "DNAmGait__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# Grip結果の変換と出力（女性）v
df_female <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmGrip_Female.png"),
  df <- df_female,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmGrip__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "Grip",
  fit <- run_elasticnet(df_female, "ChronologicalAge__c", "DNAmGrip__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# VO2結果の変換と出力（女性）
df_female <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmVO2_Female.png"),
  df <- df_female,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmVO2__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "VO2",
  fit <- run_elasticnet(df_female, "ChronologicalAge__c", "DNAmVO2__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# FEV1結果の変換と出力（女性）
df_female <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmFEV1_Female.png"),
  df <- df_female,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmFEV1__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "FEV1",
  fit <- run_elasticnet(df_female, "ChronologicalAge__c", "DNAmFEV1__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)

# TL結果の変換と出力（女性）
df_female <- transform_and_plot_linear(
  output_file_path <- paste0(output_dir_path, job_id, "_scatter_plot_DNAmTL_Female.png"),
  df <- df_female,
  column_name_x <- "ChronologicalAge__c",
  column_name_y <- "DNAmTL__c",
  output_name_x <- "Chronological Age",
  output_name_y <- "TL",
  fit <- run_elasticnet(df_female, "ChronologicalAge__c", "DNAmTL__c"), # ElasticNetで学習をさせ、最適な切片と傾きを取得
  norm <- 1
)


########################################
# 変換結果の出力
########################################

df <- rbind(df_male, df_female)

# 中間処理に作成した列の除去
df <- df[, !names(df) %in% "Highlight"]
df <- df[, !names(df) %in% "y_pred"]

# ファイルに出力
write.csv(df, paste0(output_dir_path, job_id, "_transformed_result.csv"), row.names = FALSE)

