下面把 Krambeck & Betenbaugh (2005) 文中“模型直接使用或可反推得到”的参数整理成一张总表。  
为便于阅读，按**参数类别**分组，但列结构保持一致：

**列**：参数名称｜数据大小｜单位｜数据说明（物种 / 测量对象 / 假设依据）｜文献来源

---

## 1. 运输 / 反应器参数

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| Golgi 各室停留时间 τ | 5.56 | min | UB1997 中 kG=0.18 /min 的倒数；四个高尔基室相同 | Umana & Bailey 1997 |
| Golgi 各室停留时间 τ（对照模型） | 5 | min | Monica 等糖基化唾液酸化模型采用值 | Monica et al. 1997 |
| 高尔基总体积 | 10 | μm³ | 四个室均分，每室 2.5 μm³ | Umana & Bailey 1997 |
| 单室体积 vⱼ | 2.5 | μm³ | 用于由细胞 Vm 反推酶浓度 [Eᵢ] | 本文假设 |
| 总糖链浓度 cₜₒₜ（基础例） | 500 | μM | TPO 产量反推并考虑其他糖蛋白；拟合实验用 | Inoue et al. 1999；Kaszubska et al. 2000 |
| 总糖链浓度 cₜₒₜ（文中 base case） | 800 | μM | 蛋白产率 518 pmol/10⁶ cells/24h，每蛋白 1 个糖基 | 本文 |
| 总糖链浓度范围 | 150–2300 | μM | 糖蛋白产率 100–1500 pmol/10⁶ cells/24h 对应 | Umana & Bailey 1997 |
| 高产量情形 cₜₒₜ | 2000 | μM | 基础 500 μM 的 4 倍，模拟蛋白表达提升 | 本文 |

---

## 2. 酶在 Golgi 各室的分布分数

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| ManI / ManII / GnTI–GnTV 分布 | 0.15 / 0.45 / 0.30 / 0.10 | 分数 | CGC / MGC / TGC / TGN 四室 | Rabouille et al. 1995（HeLa） |
| GalT / SiaT 分布 | 0 / 0.05 / 0.20 / 0.75 | 分数 | 后高尔基室为主 | 本文沿用 UB1997 思路 |
| GnTE 分布 | 0 / 0.60 / 0.30 / 0.10 | 分数 | GnTE 需 GalT 先作用；避免第 4 室抑制唾液酸化 | 本文假设 |

---

## 4. 初始/假定酶浓度 [Eᵢ]（单室 2.5 μm³）

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| ManI | 1.04 | μM | 由 Vm=300 与 k_f 反推 | UB1997 + 本文 |
| ManII | 1.04 | μM | 同上 | UB1997 + 本文 |
| FucT | 1.00 | μM | 无 CHO Vm，假定典型值 | 本文假设 |
| GnTI | 3.03 | μM | Vm=450 | UB1997 |
| GnTII | 0.71 | μM | Vm=140 | UB1997 |
| GnTIII | 1.10 | μM | Vm=100（野生 CHO 无活性，但模型保留） | UB1997 |
| GnTIV | 0.36 | μM | Vm=10 | UB1997 |
| GnTV | 0.05 | μM | Vm=10 | UB1997 |
| GnTE | 1.00 | μM | 无 CHO Vm，假定 | 本文假设 |
| GalT | 0.44 | μM | Vm=580 | UB1997 |
| SiaT | 1.00 | μM | 无 CHO Vm，假定 | 本文假设 |

---

## 5. 拟合实验后的调整酶浓度（TPO / CHO）

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| ManI | 1.78 | μM | 拟合 Inoue TPO 数据 | 本文 |
| ManII | 1.32 | μM | 同上 | 本文 |
| FucT | 2.50 | μM | 使非岩藻糖化 <1% | 本文 |
| GnTI | 3.05 | μM | 基本不变 | 本文 |
| GnTII | 1.29 | μM | 上调 | 本文 |
| GnTIV | 3.62 | μM | 显著上调 | 本文 |
| GnTV | 0.40 | μM | 原假设过低 | 本文 |
| GnTE | 3.47 | μM | 上调以匹配 LacNAc 延伸 | 本文 |
| GalT | 0.66 | μM | 略上调 | 本文 |
| SiaT | 1.00 | μM | 保持不变 | 本文 |

---

## 7. 供体核苷酸在 Golgi 内浓度

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| 胞质 UDP-GlcNAc | 230 | μM | CHO 0.24 fmol/cell；胞质体 1050 μm³ | Baker et al. 2001 |
| 胞质 UDP-Gal | 95 | μM | 0.10 fmol/cell | Baker et al. 2001 |
| 胞质 CMP-NeuAc | 60 | μM | 0.06 fmol/cell | Baker et al. 2001 |
| Golgi UDP-GlcNAc（40×） | 9200 | μM | 主动运输 20–100×，取 40× | Perez & Hirschberg 1985；本文 |
| Golgi UDP-Gal（40×） | 3800 | μM | 同上 | 同上 |
| Golgi CMP-NeuAc（40×） | 2400 | μM | 同上 | 同上 |
| Golgi UDP-GlcNAc（100×） | 23000 | μM | 敏感性分析用 | 本文 |

---

## 3. 酶本征活性V_m与周转数 k_f

| 参数名称 | V_m 数据大小 | V_m 单位 | k_f 数据说明 | 文献来源 |
|---|---:|---|---|---|
| ManI 本征活性 | — | μmol/mg·min | 未给；MW 230 kDa | Tulsiani & Touster 1989（大鼠肝） |
| ManII 本征活性 | 6.75 | μmol/mg·min | MW 285 kDa → k_f=1924 min⁻¹ | Moremen et al. 1991（大鼠肝） |
| FucT 本征活性 | 4.37 | μmol/mg·min | 猪脑，底物 GlcNAc-GlcNAc-Asn-PABA；MW 58 → k_f=253 | Uozumi et al. 1996 |
| GnTI 本征活性 | 19.8 | μmol/mg·min | 大鼠肝，M5 底物；MW 50 → k_f=990 | Nishikawa et al. 1988 |
| GnTII 本征活性 | 27.5 | μmol/mg·min | 大鼠肾；MW 48 → k_f=1320 | Nishikawa et al. 1992 |
| GnTIII 本征活性 | 5.52 | μmol/mg·min | 大鼠肾；MW 110 → k_f=607 | Nishikawa et al. 1992 |
| GnTIV 本征活性 | 3.23 | μmol/mg·min | 牛肠；MW 58 → k_f=187 | Oguri et al. 1997 |
| GnTV 本征活性 | 18.80 | μmol/mg·min | 大鼠肾；MW 75 → k_f=1410 | Shoreibah et al. 1992 |
| GnTE 本征活性 | 0.27 | μmol/mg·min | 小牛血清；MW 90 → k_f=25 | Kawashima et al. 1993 |
| GalT 本征活性 | 19.8 | μmol/mg·min | 大鼠肝 GlcNAc；MW 440 → k_f=8712 | Bendiak et al. 1993 |
| SiaT 本征活性 | 10.3 | μmol/mg·min | 猪颌下腺 AFGP；MW 47 → k_f=484 | Sadler et al. 1979 |
| Och1 | 95*10^-3 | μmol/mg·min | MW 55.22 k_f=5.2 | Kitajima 2006 |
| Mnn9 | | | k_f=0.9 | Striebeck 2013 |
---


## 6. 受体/供体底物 Km 与 Kmd

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| ManI Km | 100 | μM | 大鼠肝 M8 | Tabas & Kornfeld 1979 |
| FucT Km | 25 | μM | 猪脑 | Uozumi et al. 1996 |
| FucT Kmd | 46 | μM | GDP-Fuc | Uozumi et al. 1996 |
| GnTI Km | 260 | μM | CHO M5 | Chaney & Stanley 1986 |
| GnTI Kmd | 170 | μM | UDP-GlcNAc | Chaney & Stanley 1986 |
| GnTII Km | 190 | μM | 大鼠肝 | Schachter et al. 1989 |
| GnTII Kmd | 960 | μM | UDP-GlcNAc | Schachter et al. 1989 |
| GnTIII Km | 190 / 230 | μM | 大鼠 / 鸡输卵管 | Taniguchi 1989；Schachter 1989 |
| GnTIII Kmd | 3100 / 1100 | μM | 同上 | 同上 |
| GnTIV Km | 3400 | μM | 大鼠 | Taniguchi 1989 |
| GnTIV Kmd | 8300 | μM | UDP-GlcNAc | Taniguchi 1989 |
| GnTV Km | 130 / 150 / 270 | μM | 人肺 / 重组人 / 大鼠 | Gu 1993；Sasai & Ikeda 2002；Taniguchi 1989 |
| GnTV Kmd | 3500 / 4000–4600 / 7600 | μM | 同上 | 同上 |
| GnTE Km | 700 / 3000 / 9000 / 600 | μM | Novikoff / 人 / 大鼠 | Schachter 1989；Piller & Cartron 1983；van den Eijnden 1988 |
| GnTE Kmd | 55 | μM | 人源 Galβ1-4GlcNAc | Piller & Cartron 1983 |
| GalT Km | 130 | μM | 大鼠肝 | Paquet et al. 1984 |
| GalT Kmd | 0（假设饱和） | μM | 文献缺，假定饱和 | 本文假设 |
| SiaT Km | 260 | μM | 大鼠肝 asialo α1-AGP | Williams et al. 1995 |
| SiaT Kmd | 57 | μM | CMP-NeuAc | Williams et al. 1995 |
| Och1 Km | 390 | μM |  | Romero 1988, Kitajima 2006(后者提供adjustment信息) |
| Och1 Kmd | 10 | μM | 病毒 | Xiang 2010 |
| Mnn9 Km | 6500 | μM |  | Striebeck 2013 |
| Mnn9 Kmd | 500 | μM |  | Striebeck 2013 |
---


## 8. 规则修正因子（adjustment factors）

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| ManII：Man=4 时 Km× | 0.50 | 倍数 | 底物规则调整 | UB1997 |
| GnTI：Br2=0 时 Km× | 5 | 倍数 | 分支状态影响 | 本文 |
| GnTII：Br2>1 或 Br1>1 时 Km× | 1.5 | 倍数 | 已有天线影响 | 本文 |
| GnTIII：Br1>0 时 Km× | 0.178 | 倍数 | GnTV 作用后 | Oguri et al. 1997 |
| GnTIV：Br3>0 时 Km× | 0.692 | 倍数 | 已有 Br3 | 本文 |
| GnTV：Br2=0 时 Km× | 4 | 倍数 | 无 Br2 | 本文 |
| GnTE：Br1/Br2 延伸快于 Br3/Br4 | 10× | 倍数 | 相对速率 | Schachter et al. 1989 |
| GnTE：Br1+Br2 同在时 | 4× | 倍数 | 协同延伸 | Schachter et al. 1989 |
| SiaT：Sia>1 后 Km× | 5 | 倍数 | 负电荷排斥，二唾液酸后变慢 | 本文（拟合加入） |
| GnTIV：GnTII 前作用 | 0.20 | 分数活性 | 仅有 Br4 时可作用 | Schachter et al. 1989 |
| GnTIV：Br1/Br2 半乳糖化后 | 0.66 | 分数活性 | 仍可变慢作用 | Schachter et al. 1989 |

---

## 9. 初始结构 / 入口条件

| 参数名称 | 数据大小 | 单位 | 数据说明 | 文献来源 |
|---|---:|---|---|---|
| 入口 M9 : M8 | 50 : 50 | 摩尔分数 | ER/CGN 处理后进入 CGC | Umana & Bailey 1997 |
| TPO 平均糖基数 | 5.4 | 个/蛋白 | 6 个位点大多占用 | Inoue et al. 1999 |
| 野生 CHO GnTIII | 0 | 活性 | 无内源 bisecting GlcNAc | Stanley & Campbell 1984 |

---

如果你愿意，我可以把上面内容再压缩成**“模型输入参数清单（可直接编程用）”**，只保留变量名、默认值、单位和是否可调三列。