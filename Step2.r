library(terra)
library(readr)
library(dplyr)
library(stringr)
library(corrplot)
library(RColorBrewer)
library(sf)
library(caret)

# 1) 指定路径与读取栅格
bio_dir <- "E:/workspace/Amyelois.transitella/Amyelois.transitella.R/re-bio"
bio_files <- list.files(bio_dir, pattern = "\\.tif(f)?$", full.names = TRUE)
stopifnot(length(bio_files) > 0)
env <- rast(bio_files)
nm <- basename(bio_files)
nm <- gsub("\\.tif(f)?$", "", nm, ignore.case = TRUE)
nm_std <- ifelse(str_detect(tolower(nm), "bio\\D*\\d+"),
                 paste0("bio", str_extract(nm, "\\d+")),
                 ifelse(str_detect(tolower(nm), "elev|alt|elevation"), "elev", nm))
names(env) <- nm_std
env_sel <- env

# 2) 读取物种坐标并转为点（默认 WGS84，经纬度）
pts_csv <- "E:/workspace/Amyelois.transitella/Amyelois.transitella.R/AT_cleaned.csv"
raw <- readr::read_csv(pts_csv, show_col_types = FALSE)
find_col_name <- function(nms, patterns) {
  for (p in patterns) {
    m <- grep(p, nms, ignore.case = TRUE)
    if (length(m)) return(nms[m[1]])
  }
  return(NA_character_)
}
nms <- names(raw)
lon_col <- find_col_name(nms, c("^lon$", "longitude", "^x$", "经度"))
lat_col <- find_col_name(nms, c("^lat$", "latitude", "^y$", "纬度"))
if (is.na(lon_col) || is.na(lat_col)) {
  stop(paste0(
    "未能自动识别经纬度列名。\n",
    "当前列名：", paste(nms, collapse = ", "), "\n",
    "请手动把 lon_col/lat_col 改成正确的列名字符串，比如：\n",
    "lon_col <- 'Longitude'; lat_col <- 'Latitude'"
  ))
}
pts_df <- raw |>
  dplyr::transmute(lon = .data[[lon_col]],
                   lat = .data[[lat_col]]) |>
  tidyr::drop_na()

pts_sf  <- sf::st_as_sf(pts_df, coords = c("lon","lat"), crs = 4326)
if (!is.na(terra::crs(env_sel))) {
  pts_sf <- sf::st_transform(pts_sf, crs = terra::crs(env_sel))
}
pts_v <- terra::vect(pts_sf)


# 3) 按物种点提取环境变量并构建数据表
vals <- terra::extract(env_sel, pts_v, ID = FALSE)
vals <- na.omit(vals)
stopifnot(nrow(vals) > 2)  

# 4) 计算相关性矩阵（Pearson 与 Spearman）
cor_pearson  <- cor(as.data.frame(vals), use = "pairwise.complete.obs", method = "pearson")
cor_spearman <- cor(as.data.frame(vals), use = "pairwise.complete.obs", method = "spearman")


# 5) 可视化热图
library(corrplot)
library(RColorBrewer)
nv <- ncol(cor_pearson)
num_cex <- max(0.38, 0.95 - 0.03 * nv)   
lab_cex <- max(0.60, 1.05 - 0.02 * nv)   
pal_bwr <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

# ------- Pearson -------
corrplot(
  cor_pearson,
  method = "color",
  type   = "upper",
  order  = "original",
  diag   = TRUE,
  col    = pal_bwr(220),
  addCoef.col   = "#222222",
  number.cex    = num_cex,     
  number.digits = 1,           
  tl.col  = "black",
  tl.cex  = lab_cex,           
  tl.srt  = 45,
  tl.pos  = "lt",
  tl.offset = 0.6,             
  cl.pos  = "b",
  cl.length = 9,
  cl.cex  = 0.85,
  col.lim = c(-1, 1),
  title   = "Pearson Correlation",
  mar     = c(0, 0, 3, 0),
  cex.main = 1.6              
)

# ------- Spearman -------
corrplot(
  cor_spearman,
  method = "color",
  type   = "upper",
  order  = "original",
  diag   = TRUE,
  col    = pal_bwr(220),
  addCoef.col   = "#222222",
  number.cex    = num_cex,
  number.digits = 1,
  tl.col  = "black",
  tl.cex  = lab_cex,
  tl.srt  = 45,
  tl.pos  = "lt",
  tl.offset = 0.6,
  cl.pos  = "b",
  cl.length = 9,
  cl.cex  = 0.85,
  col.lim = c(-1, 1),
  title   = "Spearman Correlation",
  mar     = c(0, 0, 3, 0),
  cex.main = 1.6
)



# 6) 按阈值 |r| ≥ 0.8 进行变量筛除（默认按 Pearson）

# 选择pearson相关性用于筛选
cormat_for_filter <- cor_pearson   
to_drop_idx <- findCorrelation(abs(cormat_for_filter), cutoff = 0.8, names = FALSE, exact = TRUE)
vars_all    <- colnames(cormat_for_filter)
vars_drop   <- vars_all[to_drop_idx]
vars_keep   <- setdiff(vars_all, vars_drop)

vars_drop
vars_keep

# 7) 产出筛选结果、保存并画“筛后热图”
# 保存矩阵与变量清单
out_dir <- file.path(bio_dir, "_out"); dir.create(out_dir, showWarnings = FALSE)
write.csv(round(cor_pearson, 3),  file.path(out_dir, "cor_pearson.csv"),  row.names = TRUE)
write.csv(round(cor_spearman, 3), file.path(out_dir, "cor_spearman.csv"), row.names = TRUE)
writeLines(vars_keep, file.path(out_dir, "vars_keep_cut0.8.txt"))
writeLines(vars_drop, file.path(out_dir, "vars_drop_cut0.8.txt"))
# 筛后热图
cor_keep <- cormat_for_filter[vars_keep, vars_keep]
par(mfrow = c(1,1), mar = c(1.5, 1.5, 3, 1.5))
corrplot(
  cor_keep, method = "color", type = "upper", order = "original", diag = TRUE,
  col = pal_bwr(200), addCoef.col = "#222222", number.cex = 0.6, number.digits = 2,
  tl.col = "black", tl.cex = 0.95, tl.srt = 45, tl.pos = "lt",
  cl.pos = "b", cl.length = 11, cl.cex = 0.85, col.lim = c(-1,1),
  title = "Filtered (|r| < 0.80)", mar = c(0,0,2,0)
)

env_keep <- env_sel[[vars_keep]]
writeRaster(env_keep, filename = file.path(out_dir, "predictors_keep.tif"), overwrite = TRUE)

```
# SDMtunes筛选
## 8) 生成背景点并构建 SWD
install.packages("SDMtune")
install.packages("pROC")
library(SDMtune); library(pROC)

set.seed(42)

# 在筛后栅格 env_keep 的有效区域随机采样背景点
bg_v <- terra::spatSample(env_keep, size = 1000, method = "random", na.rm = TRUE, as.points = TRUE)
bg_sf <- sf::st_as_sf(bg_v)
bg_xy <- sf::st_coordinates(bg_sf) |> as.data.frame()
names(bg_xy) <- c("lon","lat")

# presence / background 坐标
p_xy <- pts_df[, c("lon","lat")] |> as.data.frame()
a_xy <- bg_xy[, c("lon","lat")]  |> as.data.frame()

# 构建 SDMtune 的 SWD
swd <- SDMtune::prepareSWD(
  species = "Amyelois.transitella",
  p = p_xy, a = a_xy,
  env = env_keep
)
swd

# 9) 数据拆分
set.seed(42)
datasets <- trainValTest(
  swd,
  val  = 0.2,
  test = 0.2,
  only_presence = TRUE,   
  seed = 42
)
train_swd <- datasets[[1]]
val_swd   <- datasets[[2]]
test_swd  <- datasets[[3]]


# 10) 训练基模型（Maxnet）
base_model <- train("Maxnet", data = train_swd)


# 11) 定义超参网格
h <- list(
  reg = seq(0.2, 5, 0.2),
  fc  = c("l", "lq", "lh", "lp", "lqp", "lqph")
)
set.seed(798)
opt_out <- optimizeModel(
  model  = base_model,
  hypers = h,
  metric = "auc",
  test   = val_swd,  
  pop    = 15,       
  gen    = 2,        
  seed   = 798
)

# 查看搜索结果与最优模型
opt_results <- opt_out@results
write.csv(opt_results, file.path(out_dir, "optimizeModel_results.csv"), row.names = FALSE)
best_model <- opt_out@models[[1]] 

# 12) 在独立测试集上评估
pred_test <- predict(best_model, data = test_swd, type = "cloglog")
labels    <- test_swd@pa
auc_val   <- pROC::auc(response = labels, predictor = pred_test)
cat("Test AUC (best_model):", round(as.numeric(auc_val), 4), "\n")
writeLines(paste("Test AUC (best_model):", round(as.numeric(auc_val), 4)),
           file.path(out_dir, "best_model_test_auc.txt"))

# 13) 变量重要性评估
permut_n <- 5
vi <- SDMtune::varImp(best_model, permut = permut_n)
vi <- vi[order(vi$Permutation_importance, decreasing = TRUE), ]
print(vi)
write.csv(vi, file.path(out_dir, "best_model_variable_importance.csv"), row.names = FALSE)

png(file.path(out_dir, "best_model_variable_importance.png"), width = 1600, height = 1000, res = 300)
par(mar = c(6,4,2,1))
barplot(height = vi$Permutation_importance,
        names.arg = vi$Variable,
        las = 2, cex.names = 0.9,
        ylab = "Permutation importance")
dev.off()

