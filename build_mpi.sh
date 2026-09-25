#!/bin/bash
# Build SPPARKS with real MPI WITHOUT make.
# This cluster uses Intel oneAPI.  The MPI wrapper (mpiicpc/mpiicpx) here defaults
# to the classic 'icpc' compiler, which is NOT installed; the modern compiler is
# 'icpx'.  We point the wrapper at a compiler that exists via I_MPI_CXX.
# Run from the spparks-ALD-master root directory.
set -e
cd "$(dirname "$0")"

# --- Load Intel oneAPI if present (harmless if already loaded in your shell) ---
INTEL_SETVARS="/data/software/intel/oneapi/setvars.sh"
if [ -f "$INTEL_SETVARS" ]; then
  echo "== sourcing Intel oneAPI: $INTEL_SETVARS =="
  # shellcheck disable=SC1090
  source "$INTEL_SETVARS" >/dev/null 2>&1 || true
fi

# --- Pick an MPI C++ wrapper ---
if command -v mpiicpx >/dev/null 2>&1; then
  MPICXX=mpiicpx
elif command -v mpiicpc >/dev/null 2>&1; then
  MPICXX=mpiicpc
elif command -v mpicxx >/dev/null 2>&1; then
  MPICXX=mpicxx
else
  echo "No MPI C++ compiler wrapper found (tried mpiicpx/mpiicpc/mpicxx)."
  echo "On this cluster run:  source /data/software/intel/oneapi/setvars.sh"
  exit 1
fi

# --- Force the wrapper to use a C++ compiler that really exists ---
# Intel MPI wrappers read I_MPI_CXX to choose the underlying compiler.
if [ "$MPICXX" != "mpicxx" ]; then
  if command -v icpx >/dev/null 2>&1; then
    export I_MPI_CXX=icpx
    echo "== C++ compiler: icpx (modern Intel) =="
  elif command -v icpc >/dev/null 2>&1; then
    export I_MPI_CXX=icpc
    echo "== C++ compiler: icpc (classic Intel) =="
  elif command -v g++ >/dev/null 2>&1; then
    export I_MPI_CXX=g++
    echo "== C++ compiler: g++ (GNU) =="
  else
    echo "No C++ compiler found (tried icpx/icpc/g++)."
    exit 1
  fi
fi

echo "== MPI wrapper: $MPICXX =="
# Show what the wrapper would invoke (diagnostic only, does not compile).
$MPICXX -show 2>&1 | head -3 || true

# --- Compile all sources into Obj_mpi ---
echo "== compiling sources =="
mkdir -p Obj_mpi
cp *.cpp *.h Obj_mpi/
cd Obj_mpi
for f in *.cpp; do
  $MPICXX -O -std=c++11 -DSPPARKS_GZIP -c "$f"
done

# --- Link ---
echo "== linking spk_mpi =="
$MPICXX -O *.o -o ../spk_mpi
echo "OK: ../spk_mpi built"
