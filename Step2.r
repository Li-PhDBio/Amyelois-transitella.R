library(terra)
library(rnaturalearth)
library(usdm)

#
bio_dir <- "~/workspace/Amyelois_transitella/AT/re-bio"  
bio_files <- file.path(bio_dir, sprintf("bio%d.tif", 1:19))
env_global <- rast(bio_files)
names(env_global) <- paste0("BIO", 1:19)

# ---- 2) Build M polygon (North America) and crop/mask ----
na_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  continent = "North America",
  returnclass = "sf"
)

na_v <- vect(na_sf)
na_v <- project(na_v, crs(env_global))

# 15N–50N
ext_lat <- ext(-170, -50, 15, 50)   # 经度范围覆盖北美即可
lat_box <- as.polygons(ext_lat, crs = crs(env_global))

# 
M_v <- intersect(na_v, lat_box)

# 
env_M <- crop(env_global, M_v, snap = "out")
env_M <- mask(env_M, M_v)

# ---- 可视化检查 ----
plot(env_M[[1]], main = "M region (North America, 15–50°N)")
plot(M_v, add = TRUE, border = "black", lwd = 2)

# ---- 2b) Quick visual check (M extent + one layer) ----
plot(env_M[["BIO1"]], main = "M region check (BIO1 masked by North America)")
plot(na_v, add = TRUE, lwd = 2)

# ---- 3) Sample background pixels inside M (same sample for corr & VIF) ----
set.seed(123)
n_bg <- 20000

# regular sampling is usually more stable than pure random
bg <- spatSample(env_M, size = n_bg, method = "regular",
                 na.rm = TRUE, as.points = FALSE, values = TRUE)

bg_df <- as.data.frame(bg)
bg_df <- bg_df[complete.cases(bg_df), ]


# ---- 4) Spearman correlation clustering (|r| >= corr_thr) ----
corr_thr <- 0.7

corr_mat <- suppressWarnings(cor(bg_df, method = "spearman", use = "pairwise.complete.obs"))

# cluster by distance = 1 - |r|
dist_mat <- as.dist(1 - abs(corr_mat))
hc <- hclust(dist_mat, method = "average")
clusters <- cutree(hc, h = 1 - corr_thr)

pick_rep <- function(vars, cm) {
  if (length(vars) == 1) return(vars)
  sub <- abs(cm[vars, vars, drop = FALSE])
  diag(sub) <- NA
  m <- colMeans(sub, na.rm = TRUE)     # avg |r| to others in the same cluster
  vars[which.min(m)]                   # keep the least redundant
}

cluster_list <- split(names(clusters), clusters)
selected_vars_corr <- unname(vapply(cluster_list, pick_rep, character(1), cm = corr_mat))

# ---- 4b) Visualization: correlation heatmap + dendrogram ----


# ---- 优化 Spearman 相关热图 ----
library(ggplot2)
library(reshape2)

# 转为长格式
corr_df <- reshape2::melt(corr_mat, varnames = c("Var1", "Var2"), value.name = "rho")

# 固定变量顺序（保持对角线）
corr_df$Var1_id <- as.numeric(corr_df$Var1)
corr_df$Var2_id <- as.numeric(corr_df$Var2)
p_corr <- ggplot(corr_df, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#2c7bb6",     # 蓝（负相关）
    mid = "white",
    high = "#d7191c",    # 红（正相关）
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  coord_fixed() +
  labs(
    title = "Spearman Correlation Matrix (BIO1–BIO19)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_corr)
# ---- 5) VIF step (VIF < vif_thr) on the SAME sampled dataframe ----
vif_thr <- 5

df_corr <- bg_df[, selected_vars_corr, drop = FALSE]
vif_res <- usdm::vifstep(df_corr, th = vif_thr) 
selected_vars_vif <- as.character(vif_res@results$Variables)

# 完整 VIF 结果
vif_table_full <- vif_res@results
vif_table_full$VIF <- as.numeric(vif_table_full$VIF)

# 按 VIF 从高到低排序
vif_table_full <- vif_table_full[order(vif_table_full$VIF, decreasing = TRUE), ]

print(vif_table_full)

# ---- 5b) Visualization: VIF barplot (final) ----
vif_table <- vif_res@results
vif_table$VIF <- as.numeric(vif_table$VIF)

# keep only selected (final)
vif_final <- vif_table[vif_table$Variables %in% selected_vars_vif, c("Variables", "VIF")]
vif_final <- vif_final[order(vif_final$VIF, decreasing = TRUE), ]

barplot(vif_final$VIF, names.arg = vif_final$Variables, las = 2,
        main = sprintf("Final variables VIF (th=%.1f)", vif_thr),
        ylab = "VIF")
abline(h = vif_thr, lty = 2)

# ---- 6) Final raster stack for biomod2 ----
env_M_final <- env_M[[selected_vars_vif]]

# ---- 7) Interpretability: print BIO meanings ----
bio_desc <- c(
  BIO1="Annual Mean Temperature",
  BIO2="Mean Diurnal Range (Mean of monthly (max temp - min temp))",
  BIO3="Isothermality (BIO2/BIO7) (* 100)",
  BIO4="Temperature Seasonality (standard deviation *100)",
  BIO5="Max Temperature of Warmest Month",
  BIO6="Min Temperature of Coldest Month",
  BIO7="Temperature Annual Range (BIO5-BIO6)",
  BIO8="Mean Temperature of Wettest Quarter",
  BIO9="Mean Temperature of Driest Quarter",
  BIO10="Mean Temperature of Warmest Quarter",
  BIO11="Mean Temperature of Coldest Quarter",
  BIO12="Annual Precipitation",
  BIO13="Precipitation of Wettest Month",
  BIO14="Precipitation of Driest Month",
  BIO15="Precipitation Seasonality (Coefficient of Variation)",
  BIO16="Precipitation of Wettest Quarter",
  BIO17="Precipitation of Driest Quarter",
  BIO18="Precipitation of Warmest Quarter",
  BIO19="Precipitation of Coldest Quarter"
)

cat("\n=== Correlation clustering selected (representatives) ===\n")
print(selected_vars_corr)
cat("\n=== VIF final selected ===\n")
print(selected_vars_vif)
cat("\n=== Final BIO meanings ===\n")
print(bio_desc[selected_vars_vif])

write.csv(
  vif_table_full,
  "VIF_full_results.csv",
  row.names = FALSE
)

write.csv(
  vif_table_final,
  "VIF_final_selected.csv",
  row.names = FALSE
