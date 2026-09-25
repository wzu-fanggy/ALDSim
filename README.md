# ALDSim

ALDSim is a lattice-based kinetic Monte Carlo (KMC) simulator for atomic layer
deposition (ALD) processes. It comprises two components:

- **SPPARKS-based engine (C++)** — built on the open-source
  [SPPARKS](http://www.cs.sandia.gov/~sjplimp/spparks.html) framework (Sandia
  National Laboratories). This is the engine used for the ZnS DEZ/H₂S ALD
  simulations in the associated manuscript, with the ZnS-specific model
  (`app_style ald/zns`), the input files for all simulated temperatures, and
  the modified SPPARKS routines required to build and run the simulations.
- **ALDSim Python desktop edition** (`aldsim-python/`) — the software
  registered with the China Copyright Protection Center, providing a PyQt5
  graphical interface with 2D/3D visualization of the same lattice KMC method.

## Software registration

The ALDSim software is registered with the China Copyright Protection Center
under the name **原子层沉积多尺度模拟软件 (ALDSim) V1.0**:

| Item | Value |
|---|---|
| Software name | 原子层沉积多尺度模拟软件 (ALDSim) V1.0 |
| Copyright holder | Wenzhou University (温州大学) |
| Registration no. | 2026SR0404643 |
| Registration date | March 9, 2026 |

## Contents

```
├── app_ald_zns.cpp / app_ald_zns.h     ZnS ALD application (app_style ald/zns)
├── diag_ald_zns.cpp / diag_ald_zns.h   ZnS diagnostics (species/mass counters)
├── app_ald.cpp / app_ald.h             Modified common ALD application
├── solve_linear.cpp / solve_linear.h   Modified linear solver (bounds checks)
├── aldsim-python/                      ALDSim Python desktop edition (PyQt5 GUI)
├── examples/ald-ZnS/                   Input files for the ZnS ALD simulations
│   └── ald-ZnS_300K … ald-ZnS_500K     One directory per temperature
├── MAKE/, Makefile*, STUBS/, build_serial.sh, build_mpi.sh
└── *.cpp / *.h                         SPPARKS core sources
```

## Build

Serial build (g++, no MPI required):

```bash
./build_serial.sh
```

MPI build:

```bash
./build_mpi.sh
```

## Run a ZnS ALD simulation

```bash
./spk_serial -in examples/ald-ZnS/ald-ZnS_300K/in_300.ald
```

Each temperature directory (`ald-ZnS_300K` … `ald-ZnS_500K` and
`ald-ZnS_500K_restart`) contains the input files used in the manuscript:

- `in_<T>.ald` — baseline 10-cycle DEZ/H₂S ALD run at temperature T;
- `in_<T>_calibrated.ald` — run with the calibrated kinetic parameter set;
- `in_<T>_optimized.ald` — run with the final optimized parameter set;
- `in_<T>_le2_inc.ald` — sensitivity run with the LE2 barrier raised by
  0.10 eV;
- `data.ald` — initial lattice configuration;
- `job_spk_serial` / `job_spk_mpi` — example SLURM job scripts.

## Citation

If you use ALDSim in your work, please cite the associated manuscript
(Kinetic Monte Carlo Study of Atomic Layer Deposition Process of Zinc Sulfide)
and the SPPARKS paper:

- S. J. Plimpton, C. Bystrom, and S. P. Webb, "SPPARKS: Stochastic Parallel
  PARticle Kinetic Simulator," available at
  http://www.cs.sandia.gov/~sjplimp/spparks.html.

## License

SPPARKS is distributed under the GNU General Public License; ALDSim inherits
this license. See `LICENSE` and the copyright headers in the source files.
