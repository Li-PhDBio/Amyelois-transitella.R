library(tidyverse)
library(terra)
library(sf)
library(CoordinateCleaner)
library(spThin)

base_dir <- "~/workspace/Amyelois_transitella/AT"
setwd(base_dir)

# --- 1. 读取数据 ---
raw_data <- read.csv("AT.csv", stringsAsFactors = FALSE)

# 假设列名为 gbifID, decimalLongitude, decimalLatitude, year
df <- raw_data %>%
  dplyr::select(species = organismName, 
                lon = decimalLongitude, 
                lat = decimalLatitude, 
                year = year, 
                source = institutionCode) %>%
  filter(!is.na(lon), !is.na(lat))

# --- 2. 基础清洗 (CoordinateCleaner) ---
# 移除零坐标、海域坐标（如果是陆生）、重复点
df_clean <- clean_coordinates(df, 
                              lon = "lon", lat = "lat",
                              species = "species",
                              tests = c("capitals", "centroids", "equal", "gbif", 
                                        "institutions", "zeros", "seas"), 
                              value = "clean")

# 移除 1970 年以前的记录
# 假设气候数据是 WorldClim 1970-2000
df_clean <- df_clean %>% filter(year >= 1970)

# 保存清洗后的数据
write.csv(df_clean, "AT_cleaned.csv", row.names = FALSE)
