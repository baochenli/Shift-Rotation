# Shift-Rotation
# Shift–Rotation wafer 绝对面形复现

来源：Su et al., *Absolute surface figure testing by shift-rotation method using Zernike polynomials*, Optics Letters 37(15), 3198–3200 (2012)，DOI 10.1364/OL.37.003198。原始论文位于上一级目录 `ol.37.003198.pdf`。

已用本机 MATLAB R2022a 实际运行；中文推导 `ShiftRotation_derivation.pdf` 由指定 TeX Live 2023 的 **pdfLaTeX** 编译。所有本项目输出都在本目录内，`simulation` 内的四个 XYZ 是模拟数据，不是实验数据。

## 快速运行

在 MATLAB 命令窗口执行：

```matlab
cd('F:\Chronicle\2026\论文复现\code');
test_ShiftRotation;       % 小网格正确性、输入、秩亏测试
demo_ShiftRotation(875);  % 四个模拟 XYZ -> G*x=R -> 面形 -> RMS
```

实验文件到来后，将四个文件命名为以下名称，放在本目录，核对文件表头、像素起点、y 轴方向和真实口径，然后运行 `run_experiment`：

1. `wafer_original.xyz`：原位置 W1
2. `wafer_xshift.xyz`：相对原位置沿 +x 平移 0.03 m 的 W2
3. `wafer_yshift.xyz`：相对原位置沿 +y 平移 0.03 m 的 W3
4. `wafer_rotation.xyz`：相对原位置绕中心顺时针旋转 1° 的 W4

**位置不是累积动作。四个数据必须保留同一固定探测器/参考面坐标。不能先将四个面形配准到 wafer 坐标后再直接套用本函数。** 假如仪器实际执行的是累积动作，需更换对应变换模型。

## 接口

```matlab
cfg.ShiftPixels = [-0.03 -0.03]/0.000362;
cfg.RotationDeg = 1;
cfg.GridSize = [875 875];
cfg.PixelPitch = 0.000362; % 已确认 m/pixel
cfg.IndexBase = 1;
cfg.YAxisSign = 1;
cfg.K = 81;
[wafer,out] = fun_ShiftRotation( ...
    'wafer_original.xyz','wafer_xshift.xyz', ...
    'wafer_yshift.xyz','wafer_rotation.xyz',cfg);
```

`wafer` 是 875×875 高度矩阵，单位 m，`wafer(y,x)`。原位置有效圆口径外为 NaN。需要纳米显示时用 `wafer*1e9`，不改变输入和保存结果的米单位。

物理平台移动和函数自变量移动相反：平台 +x/+y 各 0.03 m 对应 `ShiftPixels=[-82.8729281768,-82.8729281768]`。在测量侧观察、笛卡尔 +Y 向上的约定下，物理顺时针 -1° 对应自变量旋转 +1°。如果观察面相反或仪器图像已镜像，必须相应修改。代码不会将位移四舍五入到整数像素。

### 当前 XYZ 格式

三列分别为 `x_pixel y_pixel z_m`，无表头，默认坐标 1..875；每个像素一行，允许任意行序、缺失像素、NaN 高度。支持空格/制表符等 `readmatrix` 可识别的分隔符。

```text
1 1 NaN
2 1 NaN
...
438 438 1.23456789e-9
```

高度 0 是合法测量值。重复坐标、越界索引、额外列或非整数坐标会报错，避免静默误读。后续表头固定行数可设置 `HeaderLines`；复杂表头、额外列或物理坐标输入需调整 `sr_read_xyz.m`。假如仪器给出的是波前 OPD 而非表面高度，须先确认反射双程、波长和正负号等换算；当前按用户指定高度单位直接读取，不擅自除以 2。

### 常用选项

| 字段 | 默认值 | 含义 |
|---|---:|---|
| `ShiftPixels` | 必填 | x 与 y 两次动作的**自变量**有符号位移，像素 |
| `RotationDeg` | 必填 | 自变量旋转角，度 |
| `GridSize` | `[875 875]` | `[ny nx]` |
| `PixelPitch` | `0.000362` | m/pixel；物理位移应先除以此值 |
| `IndexBase` | `1` | 文件索引起点，可设为 0 |
| `CenterPixel` | `[438 438]` | **内部一基坐标**的中心，零基文件也遵守此约定 |
| `RadiusPixels` | `437` | 归一化圆口径半径，需按实际有效口径确认 |
| `YAxisSign` | `1` | 若图像行增加方向向下、物理 +Y 向上，则设 `-1` |
| `K` | `81` | 含前四个低阶模式在内的展开项数，支持 6..153 |
| `HeaderLines` | `0` | 跳过文本头部行数 |
| `IncludePiston` | `true` | 每组差分增加独立常数偏置；false 使用论文的 9 姿态项 |
| `FitStride` | `1` | 默认全像素拟合；增大可加快探索，残差仍按全像素计算 |
| `BlockSize` | `10000` | 分块 QR 的行块大小 |
| `ValidMasks` | 未设置 | 875×875×4 逻辑有效掩膜；可表示孔洞、遮挡或零填充无效点 |
| `Truth` | 未设置 | 独立真值矩阵或 XYZ 文件路径，单位 m、同一坐标系 |
| `ReturnSystem` | `false` | 小网格调试时返回 `out.G,out.R`；大系统受内存阈值限制 |

默认中心/圆直径来自像素网格假设，**不等于已经知道 wafer 的实际尺寸或旋转中心**。原文件有效区域及这两个几何参数应在实验前核对。真实 wafer 若偏心旋转，需标定中心或扩展仿射变换。

## 求解和低阶限制

论文式（2）印为 `tx1*Z1`，按其正文说明和式（3）—（9）修正为 `tx1*Z2`。模型使用有限差 `Z(Tq)-Z(q)`，不是微分，也不重采样测量高度。

默认解为 `[a5..aK, b_x(4), b_y(4), b_rot(4)]`，89 个未知数。四个姿态模式为活塞、x 倾斜、y 倾斜、power。论文只有后三个，默认额外活塞是为了容纳每次测量高度零点偏置，已在推导中明确标记为扩展。`IncludePiston=false` 使用论文的 9 姿态项，约束低阶后共86个未知量。

采用 `a1=a2=a3=a4=0` 的面形基准。**本方法无法恢复 wafer 的绝对活塞、整体倾斜和 power/曲率。** 输出是去这些不可辨识低阶后的、有限 Zernike 带宽内的绝对面形，不能表示为完整包含曲率的绝对高度。真实 power 需要其他标定方法补充。

Zernike 使用单位圆连续 RMS 归一化，明确编号见 `sr_zernike.m` 和 `out.modeTable`，不宣称与 Noll/Fringe 完全一致。81项包含完整 n≤11 及 n=12 的 m=0、m=2两项。论文未列出完整编号，后续与仪器系数对比必须按 `(n,m,cos/sin)` 映射。

逐组使用共同有效口径；不存在的像素不补零参加方程。分块增广 QR 保持原最小二乘问题，再对小矩阵列归一化并用 SVD 求解，不求正规方程逆矩阵。低阶约束后仍秩亏会报错，0°或180°旋转无法识别像散的测试已通过。

## 输出与 RMS

- `out.coefficients`：K个系数，单位 m，前四项为基准零值。
- `out.nuisance`：每列对应一组差分的姿态系数，m。
- `out.systemSize`、`out.rank`、`out.scaledConditionNumber`：拟合行列数、数值秩、列归一化后条件数。
- `out.surfaceRMS` / `surfacePV`：恢复面形本身去均值 RMS / PV，m。
- `out.residualRMS`：三组全像素差分方程残差 RMS，m。
- `out.allResidualRMS`：按像素数加权的总方程残差 RMS，m，**不是绝对精度**。
- `out.errorRMSRaw`：有独立真值时的直接逐点误差 RMS，m。
- `out.errorRMSAligned`：同一有效口径下去除差异中的活塞、倾斜和 power 后的误差 RMS，m。
- 没有 `Truth` 时误差 RMS 返回 NaN，绝不以方程残差代替。
- `out.referenceWithOriginalAlignment` 为 W1−wafer，仍含第一次姿态项，不能当作纯参考面。

三组差分因共享 W1 而含相关噪声；本次采用论文的普通最小二乘，未附加 GLS 白化或未经验证的置信区间。

## 已运行的结果

完整875×875、81项、随机种子20261001，单次测量加入0.1 nm独立高斯噪声：

| 指标 | 结果 |
|---|---:|
| 有效 G 尺寸 | 1,655,053 × 89 |
| 秩 | 89 |
| 列归一化后条件数 | 37.3141 |
| 真值对齐误差 RMS | **0.00505495 nm** |
| 差分残差 RMS | **0.14138460 nm** |
| 面形自身 RMS | 1.11495505 nm |

此外独立多项式、无噪声测试误差约1.50×10⁻¹⁴ nm；81项无噪声测试误差约6.27×10⁻¹⁵ nm，均属浮点数精度量级。测试覆盖不同掩膜、乱序/零基/表头、原论文9姿态项、基准对齐、分块求解与直接求解一致性以及退化旋转角检测。

这是模拟验证，不能复现或证明论文实验的0.13 nm重复性与0.2 nm方法间误差。81项高阶部分采用相同基合成，只验证模型内一致性；真实模型失配、支撑变形、坐标误差和非独立噪声仍需实验验证。

`simulation` 文件夹保存四个模拟输入、真值 MAT、`reconstruction.mat`、`wafer_reconstructed.xyz`、`validation.png` 和 `validation_summary.json`。`test_output` 保存小网格测试与 JSON 记录。`tmp` 保存编译/测试日志和排版检查图。

## 重编译推导 PDF

在 PowerShell 中运行：

```powershell
& 'F:\Chronicle\2026\论文复现\code\build_pdf.ps1'
```

脚本使用 `D:\texlive\2023\texlive\2023\bin\windows\pdflatex.exe` 连续编译两次生成目录和正文。已实际成功编译；Codex 内置编译器在本机报告平台标准目录不可用，因此以此本地 pdfLaTeX 输出为准。
