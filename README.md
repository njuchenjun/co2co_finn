# CO–CO<sub>2</sub> Full-Dimensional (9D) Potential Energy Surface

中文版：[README_CN.md](README_CN.md)

A 9-dimensional potential energy surface for the CO + CO<sub>2</sub> collision, built for studying collisional energy transfer. The interaction energy is fitted with fundamental-invariant neural networks (FINN) to 65 330 CCSD(T)-F12a/aug-cc-pVTZ *ab initio* points (counterpoise-corrected). The repository also contains the neural-network potential of the CO<sub>2</sub> monomer, the complete *ab initio* data set, the fifth (final) fitting run, and one Molpro example.

This repository serves as the code and data archive of the paper:

> Yuquan Feng, Yihuang Jiang, Zhengyan Jiang, Jun Chen\*,
> **Theoretical study of energy transfer between CO and CO<sub>2</sub> based on a full-dimensional potential energy surface**,
> *Chinese Journal of Chemical Physics*, **2025**, 38(6): 875–883.
> DOI: [10.1063/1674-0068/cjcp2409121](https://doi.org/10.1063/1674-0068/cjcp2409121)

## Abstract

**Abstract (as published in Chin. J. Chem. Phys. 2025, 38(6): 875–883)**

We present a comprehensive study of the CO2-CO interaction potential using a 9-dimensional (9D) potential energy surface (PES) constructed with the fundamental invariant-neural networks (FINN) method. The PES was generated from 65330 CCSD(T)-F12a/ aug-cc-pVTZ *ab initio* data points, with counterpoise correction, applied to eliminate basis set superposition error (BSSE). We performed quasi-classical trajectory simulations using this PES to investigate the rotational energy transfer dynamics. Our results reveal complex energy transfer mechanisms, with significant rotational excitation and relaxation dependent on collision energy and initial rotational states. The 9D PES provides a more accurate representation of the CO2-CO system, offering new insights into molecular dynamics and interaction mechanisms.

**Keywords (as published)** CO, CO2, Potential energy surface, Reaction dynamics, Energy transfer

## Citation

Citation information for the code, the weights and the data is given above; `CITATION.cff` carries the same metadata, and the GitHub sidebar exports it as BibTeX.

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

## Repository layout

```
.
├── src/                        Fortran sources and network weights
│   ├── pes_co2co_finn.F90      CO–CO2 9D PES, evaluation of the FINN fit
│   ├── pes_co2_nn.F90          CO2 monomer neural-network potential
│   ├── energy1_dr.F            Dynamics interface: Cartesian coordinates/Å in, energy and gradients in eV, eV/Å out
│   ├── pes_co2co_w05.dat       CO–CO2 network weights (5th fit)
│   └── pes_co2_w36.dat         CO2 monomer network weights
└── data/
    ├── abinitio/               Complete ab initio data set, 9 batches, 18 files
    │   ├── a-10000-energy.dat  a-10000-zcart
    │   ├── b-10000-*           c-10000-*
    │   ├── d-2500-*            e-3555-*
    │   ├── f-5840-*            g-12066-*
    │   └── h-12350-*           i-5172-*
    ├── nn-fit/                 Fifth (final) neural-network fitting run
    │   ├── input               Fitting control parameters
    │   ├── M05                 Training/validation mask (65330 rows, 0/1)
    │   ├── W05.txt             Fitted weights, byte-identical to src/pes_co2co_w05.dat
    │   ├── log                 Levenberg–Marquardt training log (50 cycles × 1000 epochs)
    │   ├── total-65330.dat     Training set: 13 fundamental invariants + energy in column 14 (65330 rows)
    │   └── energy_hist/        Energy distribution statistics (eab.txt, statistics.sh and the data copy it reads)
    └── examples/
        └── co2-co-315001.out   Molpro input/output example (CO2–CO, equilibrium geometry #315001)
```

The two `.dat` weight files sit next to the `.F90` files because the sources build the runtime weight-file path from `#define pesdir "./"` (see `pes_co2_nn.F90:278` and `pes_co2co_finn.F90:109`); changing that macro moves the weights anywhere else.

## The potential energy surface

| Part | Fitted to | Network | Coordinates |
|---|---|---|---|
| CO<sub>2</sub> monomer | UCCSD(T)-F12a / AVTZ (source comment) | 3–15–15–1, tansig | internal (r₁, r₂, θ) |
| CO–CO<sub>2</sub> 9D complex | 65 330 CCSD(T)-F12a/aug-cc-pVTZ points, counterpoise-corrected | 13–45–45–1, tansig | 13 fundamental invariants |

Training hyperparameters are in `data/nn-fit/input`: `NIN=13, NHID=2, NL=45/45, RATIO=0.90, NLOOP=1000, NCYCLE=50, MAX_FAIL=10`. `data/nn-fit/W05.txt` and `src/pes_co2co_w05.dat` hold the same content, i.e. the PES weights released here come from the fifth fit.

## Units and calling conventions

Energies inside the sources are in **cm⁻¹** (`cm=219474.5`, 1 hartree = 219474.5 cm⁻¹); the `energy1_dr.F` layer converts to eV for the dynamics interface.

| Routine | Input | Output |
|---|---|---|
| `pes_co2co(x, v, vco2, vco, vint, dvdr, dvco2, dvco)` | `x(3,5)` Cartesian, Å, atom order O–O–C–O–C | `v, vint, vco, vco2`: cm⁻¹; `dvdr`: cm⁻¹/Å |
| `pot_co2(n, r, v)` / `pot_co2_bba` / `pot_co2_cart` | bond lengths or Cartesian, Bohr | energy: Hartree |
| `ENERGY1(Q, VOUT, DVDQ)` | `Q(3,NATOM)` Cartesian, Å | `VOUT`: eV; `DVDQ`: eV/Å |

`ENERGY1` adds the complex interaction energy, the CO monomer and the CO₂ monomer and converts once at the end, so it can be embedded directly in a molecular dynamics code. The unit label on the trailing `rmse= …` line of the weight files and in `data/nn-fit/log` reads `mcm-1`; it is the label written by the fitting program, and the numbers follow the cm⁻¹ convention of the table above.

## Data

### `data/abinitio/`

9 batches, 71 483 geometries in total:

| Batch | Points | Batch | Points |
|---|---|---|---|
| a-10000 | 10000 | f-5840 | 5840 |
| b-10000 | 10000 | g-12066 | 12066 |
| c-10000 | 10000 | h-12350 | 12350 |
| d-2500 | 2500 | i-5172 | 5172 |
| e-3555 | 3555 | | |

- `*-energy.dat`: each line is `index + 4 energy values` (column meanings follow the original project convention).
- `*-zcart`: XYZ format, 7 lines per frame — number of atoms `5`, the index, then 5 lines of `element x y z` (Å).

### `data/nn-fit/`

- `total-65330.dat`: the final training set, 14 columns — the first 13 are the fundamental invariants, column 14 is the interaction energy (cm⁻¹).
- `M05`: the training/validation split mask of the fifth fit (`1` = training set, `0` = validation set), row order identical to `total-65330.dat`.
- `log`: the complete training log of that fit.
- `energy_hist/`: `eab.txt` counts the energy distribution in 100 cm⁻¹ bins and `statistics.sh` is the awk script that produces it; `energy_hist/total-65330.dat` is a copy of the training set, identical to the file one level up, so the script runs without path edits.

### `data/examples/co2-co-315001.out`

A complete Molpro job (CO<sub>2</sub>–CO, equilibrium geometry #315001), useful for comparing the *ab initio* input format and energies.

## Build and run

The neural-network forward passes call BLAS (`dgemv`, `dgemm`, `ddot`), and the build environment is **ifort + MKL**. From `src/`:

```bash
# object files
ifort -fpp -DNATOM=5 -O2 -c pes_co2co_finn.F90 energy1_dr.F pes_co2_nn.F90

# link with your own driver; BLAS comes from MKL
ifort -mkl -o your_program your_program.o energy1_dr.o pes_co2co_finn.o pes_co2_nn.o
```

Some notes:

- `-fpp` is required: `pes_co2co_finn.F90` and `pes_co2_nn.F90` contain `#define pesdir "./"`, and `energy1_dr.F` uses the old-style directive `#     DEFINE NATOM 5`; `NATOM` can also be supplied on the command line as `-DNATOM=5`, while `NBOND` evaluates to 10 from the line inside the file.
- `-mkl` brings in MKL at link time; a system BLAS can be used instead via `-lblas`.
- Optional optimisation: `-O2 -xHost` (code built with `-xHost` only runs on the instruction set of the build machine).
- On macOS (ifort under Rosetta, MKL in `/opt/intel/mkl`) linking needs two extra flags: `-Wl,-syslibroot,$(xcrun --show-sdk-path) -Wl,-rpath,$MKLROOT/lib`.
- The weight files are read through `pesdir` (default `./`), so running the executable inside `src/` finds them directly.

gfortran is considerably stricter about the old-style preprocessing directives above and about one actual-argument rank check in `pot_co2_cart`, and fails to compile these sources, so the build recipe here follows ifort.

## License

- Fortran code under `src/`: MIT, see [LICENSE](LICENSE).
- *Ab initio* data and network weights under `data/`: CC BY 4.0, with attribution given by the article citation above.
