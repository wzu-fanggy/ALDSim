#!/bin/bash
# Build SPPARKS serial (g++, STUBS MPI) WITHOUT make
# Run from the spparks-ALD-master root directory.
set -e
cd "$(dirname "$0")"

# 1) build STUBS dummy MPI library (skip if libmpi.a already exists)
if [ ! -f STUBS/libmpi.a ]; then
  echo "== building STUBS/libmpi.a =="
  ( cd STUBS && g++ -O -c mpi.cpp && ar rs libmpi.a mpi.o )
fi

# 2) compile all sources into Obj_serial (mirrors make serial)
echo "== compiling sources =="
mkdir -p Obj_serial
cp *.cpp *.h Obj_serial/
cd Obj_serial
for f in *.cpp; do
  g++ -O -std=c++11 -DSPPARKS_GZIP -I../STUBS -c "$f"
done

# 3) link
echo "== linking spk_serial =="
g++ -O *.o ../STUBS/libmpi.a -o ../spk_serial
echo "OK: ../spk_serial built"
