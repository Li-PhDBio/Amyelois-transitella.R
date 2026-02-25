# Amyelois transitella 物种分布模型 (SDM)

> 基于 **biomod2** 的印度谷螟 (*Amyelois transitella*) 物种分布建模源代码

## 📋 项目简介

本项目使用 R 语言构建 *Amyelois transitella*（印度谷螟 / Navel Orangeworm）的物种分布模型（Species Distribution Model, SDM）。通过整合 GBIF 分布记录与 WorldClim 生物气候变量，利用 `biomod2` 包中的多种算法进行集成建模，预测该物种在当前及未来气候情景下的潜在分布范围。

## 🔬 分析流程

本项目按三步骤（Step 1 → Step 2 → Step 3）顺序运行：

### Step 1：物种分布数据清洗 (`Step1.r`)

- 读取 GBIF 来源的物种分布记录（`AT.csv`）
- 使用 `CoordinateCleaner` 包执行坐标清洗：
  - 移除首都坐标、质心坐标、零坐标、GBIF 总部坐标、机构坐标、海洋坐标等异常记录
- 过滤 1970 年以前的历史记录（与 WorldClim 气候基准期匹配）
- 输出清洗后的数据 `AT_cleaned.csv`

### Step 2：环境变量筛选 (`Step2.r`)

- 加载 19 个 WorldClim 生物气候变量（BIO1–BIO19）
- 构建研究区域 M 区（北美洲，15°N–50°N）并裁剪环境栅格
- 采样 20,000 个背景像素点用于变��筛选
- **Spearman 相关性聚类**：基于 |ρ| ≥ 0.7 阈值进行层次聚类，从每个聚类中选取冗余最低的代表变量
- **方差膨胀因子 (VIF) 筛选**：对相关性筛选后的变量进一步执行 VIF 逐步剔除（阈值 VIF < 5）
- 生成相关性热图与 VIF 柱状图等可视化结果
- 输出最终选定的环境变量栅格集

### Step 3：biomod2 建模与预测 (`Step3.r`)

- **数据准备**：使用 SRE 策略生成伪缺失点（5 组 × 1000 个）
- **交叉验证**：采用分层交叉验证（stratified, k=4）
- **单模型建模**：运行 10 种算法
  - CTA、ANN、FDA、GAM、GBM、MARS、MAXNET、RF、SRE、XGBOOST
  - 评估指标：TSS、AUC-ROC、Kappa
  - 变量重要性评估（3 次置换）
- **集成建模**：使用 6 种集成方法（EMmedian、EMmean、EMwmean、EMca、EMci、EMcv），TSS > 0.60 的模型入选
- **空间预测**：
  - 当前气候条件投影（全球尺度）
  - 未来气候情景投影（4 组）：
    | 情景 | 时间段 | SSP |
    |------|--------|-----|
    | Future2140126 | 2021–2040 | SSP1-2.6 |
    | Future4160126 | 2041–2060 | SSP1-2.6 |
    | Future2140585 | 2021–2040 | SSP5-8.5 |
    | Future4160585 | 2041–2060 | SSP5-8.5 |
- **报告输出**：生成 BIOMOD 摘要报告、ODMAP 协议文档和代码报告

## 📂 项目结构

```
Amyelois-transitella.R/
├── Step1.r                 # 物种分布数据清洗
├── Step2.r                 # 环境变量相关性分析与 VIF 筛选
├── Step3.r                 # biomod2 建模、集成与预测
├── LICENSE                 # MIT 许可证
└── README.md               # 项目说明文档
```

### 所需数据文件（未包含在仓库中）

```
AT/
├── AT.csv                  # GBIF 原始分布记录
├── re-bio/                 # WorldClim BIO1–BIO19 栅格文件（.tif）
├── bioclim_back/           # 当前气候背景变量（研究区裁剪后）
├── bioclim_world/          # 当前气候变量（全球尺度）
├── bioclim_126n/           # 未来气候 SSP1-2.6（2021–2040）
├── bioclim_126f/           # 未来气候 SSP1-2.6（2041–2060）
├── bioclim_585n/           # 未来气候 SSP5-8.5（2021–2040）
└── bioclim_585f/           # 未来气候 SSP5-8.5（2041–2060）
```

## 🛠️ 依赖环境

### R 包

| 包名 | 用途 |
|------|------|
| `tidyverse` | 数据处理 |
| `terra` | 栅格数据处理 |
| `sf` | 空间矢量数据 |
| `CoordinateCleaner` | GBIF 坐标清洗 |
| `spThin` | 空间稀疏化（已加载） |
| `rnaturalearth` | 获取北美洲边界 |
| `usdm` | VIF 方差膨胀因子分析 |
| `ggplot2` | 可视化 |
| `reshape2` | 数据重塑 |
| `biomod2` | 物种分布建模框架 |

### 安装依赖

```r
install.packages(c(
  "tidyverse", "terra", "sf", "CoordinateCleaner", "spThin",
  "rnaturalearth", "usdm", "ggplot2", "reshape2", "biomod2"
))
```

## 🚀 使用方法

1. 准备所需数据文件（GBIF 分布记录和 WorldClim 气候栅格数据）
2. 按顺序运行脚本：

```r
source("Step1.r")  # 数据清洗
source("Step2.r")  # 变量筛选
source("Step3.r")  # 建模与预测
```

> ⚠️ **注意**：请根据实际数据路径修改各脚本中的 `base_dir` 和文件路径。

## 📊 输出结果

- `AT_cleaned.csv` — 清洗后的分布记录
- `VIF_full_results.csv` — 完整 VIF 分析结果
- `VIF_final_selected.csv` — 最终选定变量的 VIF 值
- Spearman 相关性热图
- VIF 柱状图
- 模型评估图表（TSS / AUC / Kappa 箱线图）
- 变量重要性排序图
- 物种响应曲线
- 当前与未来分布预测地图
- BIOMOD 摘要报告 / ODMAP 协议

## 📄 许可证

本项目基于 [MIT License](LICENSE) 开源。

## 📖 参考

- [biomod2 官方文档](https://biomodhub.github.io/biomod2/)
- [WorldClim 气候数据](https://www.worldclim.org/)
- [GBIF 全球生物多样性信息网络](https://www.gbif.org/)
- Thuiller, W., et al. (2009). BIOMOD – a platform for ensemble forecasting of species distributions. *Ecography*, 32(3), 369-373.
