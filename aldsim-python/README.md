# ALDSim (Python desktop edition)

The Python desktop edition of ALDSim — the software registered with the
China Copyright Protection Center under the name **原子层沉积多尺度模拟软件
(ALDSim) V1.0** (registration no. 2026SR0404643, copyright holder: Wenzhou
University). It provides a PyQt5-based graphical interface with 2D/3D
visualization for lattice-based ALD kinetic Monte Carlo simulations.

```
├── main.py                  Application entry (PyQt5 GUI)
├── requirements.txt         Python dependencies
├── core/                    Simulation engine
│   ├── config.py            Simulation configuration
│   ├── cell.py              Lattice cell definitions
│   ├── layer.py             Lattice layer / species states
│   ├── react_arrival.py     Precursor arrival (gas impingement) events
│   ├── react_event.py       Surface reaction events
│   └── simulation_proc.py   KMC simulation procedure
└── ui/                      Graphical interface
    ├── main_window.py       Main window
    ├── view_2d.py / view_3d.py   2D / 3D lattice visualization
    └── dialog_data.py / dialog_params.py   Data & parameter dialogs
```

## Run

```bash
pip install -r requirements.txt
python main.py
```

Simulation results are stored in a local SQLite database (`data/ald.db`)
created at runtime; no data files are bundled in this repository.
