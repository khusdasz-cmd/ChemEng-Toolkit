# ChemEng-Toolkit

**化工数据分析 Python 工具箱** — 热重（TGA/DSC）、红外光谱（FT-IR）、热解动力学、精馏塔设计。

由化工专业本科生开发，用于本科科研和课程项目。

## 功能

- **TGA/DSC 处理** — 读取 TA Instruments 的原始 .xls 文件，绘制 TG/DTG 曲线，识别分解阶段，提取特征参数（T5%、T10%、Tmax、残炭率）
- **热解动力学** — Coats-Redfern 法，内置 15 种反应机理模型（F1-F3、R2-R3、D1-D4、A2-A4、P2-P4），支持分段活化能分析（Ea vs α）
- **FT-IR 光谱** — 读取 CSV 光谱数据，绘制垂直偏移堆叠图，自动标注特征吸收峰
- **出版级图表** — 基于 Matplotlib，统一配色方案和专业样式

## 快速开始

```bash
# 安装
pip install -e .

# 或用 conda
conda install -c conda-forge numpy pandas matplotlib scipy xlrd openpyxl
pip install -e .
```

### TGA 分析

```python
from chemeng_toolkit.thermal_analysis import TGAProcessor

# 加载数据
processor = TGAProcessor()
data = processor.load_from_xls("data/tga_sample/gt_100k.xls")

# 绘图
processor.plot_tg()
processor.plot_dtg()
processor.plot_dual_axis(">100k")

# 提取特征参数
summary = processor.summary_table()
```

### FT-IR 分析

```python
from chemeng_toolkit.thermal_analysis import FTIRProcessor

processor = FTIRProcessor(data_dir="data/ftir_sample")
processor.load_batch([
    ("sample.CSV", "Label", "#1f77b4"),
])

processor.plot_stacked(title="FT-IR Spectra")
```

### 动力学分析

```python
from chemeng_toolkit.thermal_analysis import coats_redfern, plot_ea_vs_alpha

results = coats_redfern(temp, weight)
print(f"最佳模型: {results[0]['code']}, E = {results[0]['E_kJ']:.2f} kJ/mol")
```

## 项目结构

```
ChemEng-Toolkit/
├── chemeng_toolkit/
│   ├── thermal_analysis/
│   │   ├── tga_processor.py      # TGA 数据加载与可视化
│   │   ├── kinetics.py           # Coats-Redfern 动力学分析
│   │   └── ftir_processor.py     # FT-IR 光谱处理
│   └── utils/
│       └── plot_helpers.py       # 共享样式和工具函数
├── examples/
│   ├── tga_example.ipynb         # TGA 分析示例教程
│   └── ftir_example.ipynb        # FT-IR 分析示例教程
└── data/
    ├── tga_sample/               # TGA 示例 .xls 文件
    └── ftir_sample/              # FT-IR 示例 .CSV 文件
```

## 背景

本工具包是作者在 **广东工业大学**（化学工程与工艺专业，预计 2026 年 6 月毕业）本科科研期间开发的，主要应用于 **木质素表征与生物质资源化** 研究。涉及课程知识包括：化工热力学、传递过程、反应工程、分离工程。

## 环境要求

- Python >= 3.9
- numpy、pandas、matplotlib、scipy
- xlrd（读取 .xls TGA 数据）
- openpyxl（读取 .xlsx 数据）

## 开源协议

MIT
