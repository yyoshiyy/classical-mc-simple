# classical-mc-simple

Minimal, beginner-friendly classical spin Monte Carlo code (single process, no MPI).

The current implementation is in [`src/`](./src).

## Quick Start

1. Build:

```bash
cd src
cmake -S . -B build
cmake --build build -j
```

2. Run the provided 2D Ising sample (`L=16`):

```bash
cd samples/square_L16_Ising
../../build/MC_simple
```

3. Check results:

```bash
cat MC_simple_result.dat
```

## How To Use

Run `MC_simple` in a directory containing:

- `param.def`
- `lattice.def`
- `interaction.def`

CLI forms:

```bash
./MC_simple
./MC_simple my_param.def
./MC_simple my_param.def my_lattice.def my_interaction.def
```

Optional reproducibility seed:

```bash
MC_SIMPLE_SEED=12345 ./MC_simple
```

## What It Computes

Output file: `MC_simple_result.dat`

Columns:

1. `T`
2. `E_per_site`
3. `C_per_site`
4. `M2`
5. `acceptance`

Definitions:

- `E_per_site`: thermal average of energy per site
- `C_per_site`: `C = (<E^2> - <E>^2) / (N * T^2)`
- `M2`: `<(Mx^2 + My^2 + Mz^2) / N^2>`
- `acceptance`: Metropolis acceptance ratio

## Input Files

- `param.def`: simulation parameters (`Burn_in`, `Total_Step`, `Sample`, `num_temp`, `Ini_T`, `Delta_T`, `spin_dim`, `lambda`, `H`)
- `lattice.def`: lattice size (`L_x`, `L_y`, `L_z`) and orbital count (`orb_num`)
- `interaction.def`: bond list (`i_orb j_orb dx dy dz J`)

## Code Structure

Main source files are in `src/`:

- `main.c`: overall simulation flow (read input, run MC, average, write output)
- `memory.c`, `memory.h`: allocation and cleanup
- `input_parser.c`: parsing of `param.def`, `lattice.def`, `interaction.def`
- `lattice.c`: neighbor table construction, spin initialization, initial fields/energy
- `mc_update.c`: one Metropolis sweep per temperature slot
- `mc_def.h`: shared structs/prototypes
- `dSFMT.c`, `dSFMT.h`, `dSFMT-params.h`, `dSFMT-params19937.h`: random number generator
- `CMakeLists.txt`: build config

## Notes

- Single-process only (no MPI)
- Replica exchange and over-relaxation are intentionally omitted
- dSFMT parameter is fixed to `MEXP=19937` in this minimal version

## License

- This code is MIT licensed: see [`src/LICENSE`](./src/LICENSE)
- Included dSFMT files are BSD 3-Clause licensed: see [`src/LICENSE.txt`](./src/LICENSE.txt)
