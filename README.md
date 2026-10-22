# CO–CO<sub>2</sub> 全维（9D）势能面

面向 CO + CO<sub>2</sub> 碰撞能量传递的 9 维全维势能面：用基本不变量神经网络（fundamental invariant neural network, FINN）拟合 65330 个 CCSD(T)-F12a/aug-cc-pVTZ（含 counterpoise 校正）ab initio 能量点，并配套给出 CO<sub>2</sub> 单体神经网络势、全部 ab initio 原始数据、第 5 次（最终）拟合产物与一个 Molpro 计算样例。

本仓库同时是论文的代码与数据存档：

> Yuquan Feng, Yihuang Jiang, Zhengyan Jiang, Jun Chen\*,
> **Theoretical study of energy transfer between CO and CO<sub>2</sub> based on a full-dimensional potential energy surface**,
> *Chinese Journal of Chemical Physics*, **2025**, 38(6): 875–883.
> DOI: [10.1063/1674-0068/cjcp2409121](https://doi.org/10.1063/1674-0068/cjcp2409121)

> **Abstract (EN)** A 9-dimensional potential energy surface for the CO<sub>2</sub>–CO van der Waals complex is constructed with the FINN method from 65 330 counterpoise-corrected CCSD(T)-F12a/aug-cc-pVTZ points, and used in quasi-classical trajectory simulations of rotational energy transfer. This repository hosts the Fortran evaluation routines, the fitted network weights, the full *ab initio* data set, the final (5th) fitting run, and a Molpro example.

---

## 引用

如果本仓库的代码、权重或数据对你有帮助，请引用上面的文章。仓库内附 `CITATION.cff`，GitHub 侧边栏可直接导出 BibTeX。

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
├── README.md
├── CITATION.cff
├── LICENSE
├── .gitignore
├── .gitattributes
├── src/                        Fortran 源码 + 运行时权重文件
│   ├── pes_co2co_finn.F90      CO–CO2 9D PES（FINN 拟合式求值）
│   ├── pes_co2_nn.F90          CO2 单体 NN 势
│   ├── energy1_dr.F            对外接口：Cartesian/Å 进，能量与梯度以 eV、eV/Å 出
│   ├── pes_co2co_w05.dat       CO–CO2 NN 权重（第 5 次拟合，被 pes_co2co_finn.F90 运行时读取）
│   └── pes_co2_w36.dat         CO2 单体 NN 权重（被 pes_co2_nn.F90 运行时读取）
└── data/
    ├── abinitio/               全部 ab initio 数据（9 组，共 18 个文件）
    │   ├── a-10000-energy.dat  a-10000-zcart
    │   ├── b-10000-*           c-10000-*
    │   ├── d-2500-*            e-3555-*
    │   ├── f-5840-*            g-12066-*
    │   ├── h-12350-*           i-5172-*
    ├── nn-fit/                 第 5 次（最终）NN 拟合产物
    │   ├── input               拟合控制参数
    │   ├── M05                 训练/验证掩码（65330 行，0/1）
    │   ├── W05.txt             拟合得到的权重（与 src/pes_co2co_w05.dat 逐字节相同）
    │   ├── log                 Levenberg–Marquardt 训练日志（50 个 cycle × 1000 epoch）
    │   ├── total-65330.dat     训练集：13 个基本不变量 + 第 14 列能量（65330 行）
    │   └── energy_hist/        能量分布统计（eab.txt、statistics.sh 及其所用的数据副本）
    └── examples/
        └── co2-co-315001.out   Molpro 输入/输出样例（CO2–CO eq geom #315001）
```

> `src/` 里的两个 `.dat` 权重文件必须与 `.F90` 放在同一目录：源码用 `#define pesdir "./"` 拼接运行时权重路径
> （`pes_co2_nn.F90:278`、`pes_co2co_finn.F90:109`）。若要放到别处，改这个宏或在编译时重定义。

## 势能面构成

| 部分 | 拟合对象 | 网络结构 | 说明 |
|---|---|---|---|
| CO<sub>2</sub> 单体 | UCCSD(T)-F12a / AVTZ（源码注释） | 3–15–15–1，tansig | 内坐标 (r₁, r₂, θ) |
| CO–CO<sub>2</sub> 9D 复合物 | 65330 个 CCSD(T)-F12a/aug-cc-pVTZ 点，含 counterpoise 校正 | 13–45–45–1，tansig | 13 个基本不变量，Levenberg–Marquardt 训练 |

- 训练超参数见 `data/nn-fit/input`：`NIN=13, NHID=2, NL=45/45, RATIO=0.90, NLOOP=1000, NCYCLE=50, MAX_FAIL=10`。
- `data/nn-fit/W05.txt` 与 `src/pes_co2co_w05.dat` 逐字节相同（md5 `af1b3eebc28a624b647cf2f73ae321c4`），
  也就是说本仓库发布的 PES 权重正是第 05 次拟合的结果。

## 单位与调用约定

源码内部一律用 **cm⁻¹**（`cm=219474.5`，1 hartree = 219474.5 cm⁻¹）；只有 `energy1_dr.F` 这一层为分子动力学接口换算成 eV。

| 例程 | 输入 | 输出 |
|---|---|---|
| `pes_co2co(x, v, vco2, vco, vint, dvdr, dvco2, dvco)` | `x(3,5)` 笛卡尔坐标，**Å**，原子顺序 **O–O–C–O–C** | `v, vint, vco, vco2`：cm⁻¹；`dvdr`：cm⁻¹/Å |
| `pot_co2(n, r, v)` / `pot_co2_bba` / `pot_co2_cart` | 键长或笛卡尔坐标，**Bohr** | 能量：Hartree |
| `ENERGY1(Q, VOUT, DVDQ)` | `Q(3,NATOM)` 笛卡尔坐标，**Å** | `VOUT`：eV；`DVDQ`：eV/Å |

`ENERGY1` 内部把复合物相互作用能、CO 单体、CO₂ 单体三部分相加后统一换算，因此可直接嵌进 MD 程序。

## 数据说明

### `data/abinitio/`

9 组、共 71483 个构型（原目录名 `all-data-69244` 沿用自原始工作目录）：

| 组 | 点数 | 组 | 点数 |
|---|---|---|---|
| a-10000 | 10000 | f-5840 | 5840 |
| b-10000 | 10000 | g-12066 | 12066 |
| c-10000 | 10000 | h-12350 | 12350 |
| d-2500 | 2500 | i-5172 | 5172 |
| e-3555 | 3555 | | |

- `*-energy.dat`：每行 `编号 + 4 列能量值`（列含义沿用原项目约定）。
- `*-zcart`：XYZ 格式，每帧 7 行 —— 原子数 `5`、编号、然后 5 行 `元素 x y z`（Å）。

### `data/nn-fit/`

- `total-65330.dat`：最终训练集，14 列 —— 前 13 列为基本不变量，第 14 列为相互作用能（cm⁻¹）。
- `M05`：第 5 次拟合所用的训练/验证划分掩码（`1` = 训练集，`0` = 验证集），行序与 `total-65330.dat` 一致。
- `log`：该次拟合完整日志（每 1000 epoch 一个 cycle，共 50 个 cycle）。
- `energy_hist/`：`eab.txt` 是按 100 cm⁻¹ 分箱的能量分布计数；`statistics.sh` 是生成它的 awk 脚本（在当前目录下直接 `sh statistics.sh` 即可重算）；`energy_hist/total-65330.dat` 是训练集副本，与上一级同名文件内容相同，保留是为了让脚本无需改路径即可运行。

### `data/examples/co2-co-315001.out`

一个完整的 Molpro 计算样例（CO<sub>2</sub>–CO，equilibrium geometry #315001），可用于对照 ab initio 输入格式与能量。

## 编译

代码是普通 Fortran（`.F90` 含 `#define`，需预处理）：

```bash
cd src
gfortran -cpp -O2 -DNATOM=5 -c pes_co2co_finn.F90
gfortran -cpp -O2 -DNATOM=5 -c energy1_dr.F
gfortran -cpp -O2 -c pes_co2_nn.F90
```

- `energy1_dr.F` 用的是 `#     DEFINE NATOM 5` 这种带缩进的旧式预处理指令，gfortran 会报
  `invalid preprocessing directive`；用命令行 `-DNATOM=5 -DNBOND=10` 代替（或先把这类行删掉）即可。
- 权重文件按运行时目录 `./` 读取，所以可执行程序要在 `src/` 下运行，或相应修改 `pesdir`。

## 已知问题

1. **`pes_co2_nn.F90` 第 40 行**：`pot_co2_cart` 里 `call pot_co2(1,r,v)` 把标量 `v` 传给了数组形参，
   gfortran ≥ 10 会报 `Rank mismatch in argument 'v'`。该 wrapper 在本仓库中没有调用者（`ENERGY1` 走的是
   `pes_co2co` 里返回的 `dvco2`），所以不影响 `energy1_dr.F` + `pes_co2co_finn.F90` 这条主链路；
   若要单独使用 `pot_co2_cart`，需要把它改成先取数组再取元素。
2. 权重文件末尾 `rmse= …` 行与 `data/nn-fit/log` 中的单位标签写作 **`mcm-1`**（原始拟合程序输出的是 `meV`，
   整理时只替换了标签、数值一位未动）；`src/pes_co2_w36.dat` 仍是原件写法 `meV`。
   无论在哪种写法下，**源码内部一律以 cm⁻¹ 为准**（`cm=219474.5`）。
3. `data/nn-fit/log` 约 2.9 MB / 4 万行，属于原始运行记录，保留以备考证。

## 仓库整理说明

本仓库由原始计算工作目录整理而来，只做目录归类与文件移动，数据内容未变。**唯一例外是单位标签**：
`src/pes_co2co_w05.dat`、`data/nn-fit/W05.txt`、`data/nn-fit/log` 三个文件里的 `meV` 已在本地统一替换为
`mcm-1`（仅改标签，数值未动），因此这三个文件的 md5 与 neckar 原件不同 —— 原件 md5 依次为
`9b8e167d…`、`9b8e167d…`、`1a643b71…`，需要回溯时可直接从原计算目录取回；其余 **29 个文件与原件逐字节一致**。

原路径与新路径对应关系：

| 原始路径（neckar） | 本仓库路径 |
|---|---|
| `/home/chenjun/pes_co2_co` | `src/` |
| `.../vdw-co-co2/3-nnfit/z-final-fengyq/new-nn/{M05,W05.txt,input,log,...}` | `data/nn-fit/` |
| `.../vdw-co-co2/2-abinitio/z-final-fengyq/all-data-69244` | `data/abinitio/` |
| `.../2-abinitio/mol_inp_31/1-10000-out/done-5001-6000/co2-co-315001.out` | `data/examples/co2-co-315001.out` |

原始目录里还有第 01–04、06–50 次拟合的中间权重（`M*`/`W*`），本仓库只保留最终的第 05 次。

## 大文件

`data/` 约占 55 MB（其中两个 16 MB 的 `total-65330.dat`）。直接 clone 即可，无需 LFS；
如果你 fork 后要频繁改动数据，建议对 `data/**` 启用 [Git LFS](https://git-lfs.com/)。

## 许可证

- `src/` 下的 Fortran 代码：MIT，见 [LICENSE](LICENSE)。
- `data/` 下的 ab initio 数据与网络权重：CC BY 4.0，使用请注明出处并引用上面的文章。

## 联系

陈军（Jun Chen），中国科学院福建物质结构研究所，<chenjun@fjirsm.ac.cn>
