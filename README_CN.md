# CO–CO<sub>2</sub> 全维（9D）势能面

英文版：[README_EN.md](README_EN.md)

面向 CO + CO<sub>2</sub> 碰撞能量传递的 9 维全维势能面：用基本不变量神经网络（fundamental invariant neural network, FINN）拟合 65330 个 CCSD(T)-F12a/aug-cc-pVTZ（含 counterpoise 校正）ab initio 能量点，并配套给出 CO<sub>2</sub> 单体神经网络势、全部 ab initio 原始数据、第 5 次（最终）拟合产物与一个 Molpro 计算样例。

本仓库同时是论文的代码与数据存档：

> Yuquan Feng, Yihuang Jiang, Zhengyan Jiang, Jun Chen\*,
> **Theoretical study of energy transfer between CO and CO<sub>2</sub> based on a full-dimensional potential energy surface**,
> *Chinese Journal of Chemical Physics*, **2025**, 38(6): 875–883.
> DOI: [10.1063/1674-0068/cjcp2409121](https://doi.org/10.1063/1674-0068/cjcp2409121)

## Abstract

**英文摘要（原文，Chin. J. Chem. Phys. 2025, 38(6): 875–883）**

We present a comprehensive study of the CO2-CO interaction potential using a 9-dimensional (9D) potential energy surface (PES) constructed with the fundamental invariant-neural networks (FINN) method. The PES was generated from 65330 CCSD(T)-F12a/ aug-cc-pVTZ *ab initio* data points, with counterpoise correction, applied to eliminate basis set superposition error (BSSE). We performed quasi-classical trajectory simulations using this PES to investigate the rotational energy transfer dynamics. Our results reveal complex energy transfer mechanisms, with significant rotational excitation and relaxation dependent on collision energy and initial rotational states. The 9D PES provides a more accurate representation of the CO2-CO system, offering new insights into molecular dynamics and interaction mechanisms.

**中文摘要（原文）**

本文利用基本不变量-神经网络方法，构建了CO2与CO分子之间相互作用的9维势能面，并用于碰撞能量传递研究. 势能面基于65330个CCSD(T)-F12a/aug-cc-pVTZ高精度从头算数据，计算中采用了均衡校正以消除基组重叠误差. 以此势能面为基础，通过准经典轨线计算研究了CO2与CO分子间的转动能量转移动力学. 研究结果揭示了其中复杂的能量转移机制，表明产物的转动激发和弛豫性质较大地依赖于初始碰撞能. 这个新的全维势能面能精确描述CO2-CO体系，有望为分子动力学和相互作用机制提供新的理解.

**关键词（原文）**　CO, CO2, Potential energy surface, Reaction dynamics, Energy transfer

## 引用

代码、权重与数据的引用信息见上；`CITATION.cff` 收录了同样的元数据，GitHub 侧边栏可直接导出 BibTeX。

```bibtex
@article{Feng2025COCO2,
  author  = {Feng, Yuquan and Jiang, Yihuang and Jiang, Zhengyan and Chen, Jun},
  title   = {Theoretical study of energy transfer between CO and CO2 based on a full-dimensional potential energy surface},
  journal = {Chinese Journal of Chemical Physics},
  volume  = {38},
  number  = {6},
  pages   = {875--883},
  year    = {2025},
  doi     = {10.1063/1674-0068/cjcp2409121}
}
```

## 目录结构

```
.
├── src/                        Fortran 源码与网络权重
│   ├── pes_co2co_finn.F90      CO–CO2 9D PES，FINN 拟合式求值
│   ├── pes_co2_nn.F90          CO2 单体神经网络势
│   ├── energy1_dr.F            动力学接口：笛卡尔坐标/Å 输入，能量与梯度以 eV、eV/Å 输出
│   ├── pes_co2co_w05.dat       CO–CO2 网络权重（第 5 次拟合）
│   └── pes_co2_w36.dat         CO2 单体网络权重
└── data/
    ├── abinitio/               全部 ab initio 数据，9 组共 18 个文件
    │   ├── a-10000-energy.dat  a-10000-zcart
    │   ├── b-10000-*           c-10000-*
    │   ├── d-2500-*            e-3555-*
    │   ├── f-5840-*            g-12066-*
    │   └── h-12350-*           i-5172-*
    ├── nn-fit/                 第 5 次（最终）神经网络拟合产物
    │   ├── input               拟合控制参数
    │   ├── M05                 训练/验证掩码（65330 行，0/1）
    │   ├── W05.txt             拟合得到的权重，与 src/pes_co2co_w05.dat 逐字节相同
    │   ├── log                 Levenberg–Marquardt 训练日志（50 个 cycle × 1000 epoch）
    │   ├── total-65330.dat     训练集：13 个基本不变量 + 第 14 列能量（65330 行）
    │   └── energy_hist/        能量分布统计（eab.txt、statistics.sh 及其所用数据副本）
    └── examples/
        └── co2-co-315001.out   Molpro 输入/输出样例（CO2–CO，equilibrium geometry #315001）
```

两个 `.dat` 权重文件与 `.F90` 同目录，是因为源码通过 `#define pesdir "./"` 拼接运行时的权重文件路径（见 `pes_co2_nn.F90:278`、`pes_co2co_finn.F90:109`）；改这个宏即可把权重放到其他位置。

## 势能面构成

| 部分 | 拟合对象 | 网络结构 | 坐标 |
|---|---|---|---|
| CO<sub>2</sub> 单体 | UCCSD(T)-F12a / AVTZ（源码注释） | 3–15–15–1，tansig | 内坐标 (r₁, r₂, θ) |
| CO–CO<sub>2</sub> 9D 复合物 | 65330 个 CCSD(T)-F12a/aug-cc-pVTZ 点，含 counterpoise 校正 | 13–45–45–1，tansig | 13 个基本不变量 |

训练超参数见 `data/nn-fit/input`：`NIN=13, NHID=2, NL=45/45, RATIO=0.90, NLOOP=1000, NCYCLE=50, MAX_FAIL=10`。`data/nn-fit/W05.txt` 与 `src/pes_co2co_w05.dat` 内容相同，即仓库发布的 PES 权重来自第 05 次拟合。

## 单位与调用约定

源码内部能量以 **cm⁻¹** 表示（`cm=219474.5`，1 hartree = 219474.5 cm⁻¹）；`energy1_dr.F` 这一层为动力学接口做换算，对外给出 eV。

| 例程 | 输入 | 输出 |
|---|---|---|
| `pes_co2co(x, v, vco2, vco, vint, dvdr, dvco2, dvco)` | `x(3,5)` 笛卡尔坐标，Å，原子顺序 O–O–C–O–C | `v, vint, vco, vco2`：cm⁻¹；`dvdr`：cm⁻¹/Å |
| `pot_co2(n, r, v)` / `pot_co2_bba` / `pot_co2_cart` | 键长或笛卡尔坐标，Bohr | 能量：Hartree |
| `ENERGY1(Q, VOUT, DVDQ)` | `Q(3,NATOM)` 笛卡尔坐标，Å | `VOUT`：eV；`DVDQ`：eV/Å |

`ENERGY1` 内部把复合物相互作用能、CO 单体、CO₂ 单体三项相加后统一换算，因此可以直接嵌进分子动力学程序。权重文件末尾的 `rmse= …` 行与 `data/nn-fit/log` 中的单位标签写作 `mcm-1`，是拟合程序输出时的标注；能量数值本身与上表的 cm⁻¹ 约定一致。

## 数据说明

### `data/abinitio/`

9 组、合计 71483 个构型：

| 组 | 点数 | 组 | 点数 |
|---|---|---|---|
| a-10000 | 10000 | f-5840 | 5840 |
| b-10000 | 10000 | g-12066 | 12066 |
| c-10000 | 10000 | h-12350 | 12350 |
| d-2500 | 2500 | i-5172 | 5172 |
| e-3555 | 3555 | | |

- `*-energy.dat`：每行 `编号 + 4 列能量值`（列含义沿用原始项目的约定）。
- `*-zcart`：XYZ 格式，每帧 7 行 —— 原子数 `5`、编号，随后 5 行 `元素 x y z`（Å）。

### `data/nn-fit/`

- `total-65330.dat`：最终训练集，14 列 —— 前 13 列为基本不变量，第 14 列为相互作用能（cm⁻¹）。
- `M05`：第 5 次拟合的训练/验证划分掩码（`1` = 训练集，`0` = 验证集），行序与 `total-65330.dat` 一致。
- `log`：该次拟合的完整训练日志。
- `energy_hist/`：`eab.txt` 是按 100 cm⁻¹ 分箱的能量分布计数，`statistics.sh` 是生成它的 awk 脚本；`energy_hist/total-65330.dat` 是训练集副本，与上一级同名文件内容相同，脚本因此无需改路径即可运行。

### `data/examples/co2-co-315001.out`

一个完整的 Molpro 计算样例（CO<sub>2</sub>–CO，equilibrium geometry #315001），可对照 ab initio 输入格式与能量。

## 编译与运行

源码在神经网络前向计算中调用 BLAS（`dgemv`、`dgemm`、`ddot`），对应的编译环境是 **ifort + MKL**。以 `src/` 为工作目录：

```bash
# 编译目标文件
ifort -fpp -DNATOM=5 -O2 -c pes_co2co_finn.F90 energy1_dr.F pes_co2_nn.F90

# 与自己的驱动程序一起链接，BLAS 由 MKL 提供
ifort -mkl -o your_program your_program.o energy1_dr.o pes_co2co_finn.o pes_co2_nn.o
```

几点说明：

- `-fpp` 是必需的：`pes_co2co_finn.F90`、`pes_co2_nn.F90` 里有 `#define pesdir "./"`，`energy1_dr.F` 里用的是 `#     DEFINE NATOM 5` 这类旧式指令；`NATOM` 也可以由命令行 `-DNATOM=5` 给出，`NBOND` 由文件内那一行推得 10。
- `-mkl` 在链接阶段引入 MKL；改用系统 BLAS 时可替换为 `-lblas`。
- 可选优化：`-O2 -xHost`（`-xHost` 生成的代码只适配编译机的指令集）。
- macOS（ifort 经 Rosetta、MKL 在 `/opt/intel/mkl`）链接时还要加两个参数：`-Wl,-syslibroot,$(xcrun --show-sdk-path) -Wl,-rpath,$MKLROOT/lib`。
- 权重文件按 `pesdir`（默认 `./`）读取，可执行程序在 `src/` 目录下运行时即可直接读到位。

gfortran 对上面这类旧式预处理指令、以及 `pot_co2_cart` 中一处实参秩的检查严格得多，直接用 gfortran 无法编译，因此仓库以 ifort 的构建方式为准。

## 许可证

- `src/` 下的 Fortran 代码：MIT，见 [LICENSE](LICENSE)。
- `data/` 下的 ab initio 数据与网络权重：CC BY 4.0，署名信息即上文的文章引用。
