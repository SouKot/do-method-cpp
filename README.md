# Dynamically Orthogonal (DO) Method — C++ Implementation

A parallel C++17 / Fortran framework for uncertainty quantification of
time-dependent PDEs using the **Dynamically Orthogonal (DO) field equations**.
Three model problems are included:

| Target | PDE | Parallelism |
|--------|-----|-------------|
| `BRGRDO` | 1-D viscous Burgers equation | Serial + MPI |
| `QGDO`  | Barotropic Quasi-Geostrophic vorticity equation | Serial only |
| `SWEDO` | Rotating Shallow-Water Equations (FVM) | Serial + MPI |

The QG problem uses the same parameters as in the
[JBIMAU/DO-methods/QG/matlab](https://github.com/JBIMAU/DO-methods) reference.

---

## Mathematical Background

The DO method represents a stochastic field as a low-rank decomposition

$$u(x,t,\omega) = \bar{u}(x,t) \;+\; \sum_{i=1}^{m} V_i(x,t)\; Y_i(t,\omega)$$

where

| Symbol | Meaning | Owner class |
|--------|---------|-------------|
| $\bar{u}$ | Deterministic mean field | `Mean` (Burger / QG / SWE) |
| $V \in \mathbb{R}^{N \times m}$ | Orthonormal stochastic basis | `BasisSolver` |
| $Y \in \mathbb{R}^{m \times S}$ | Stochastic coefficients (S realisations) | `Y_Stoch` |

At each time step the framework evolves three coupled sub-systems:

1. **Mean solve** — implicit θ-method + Newton iteration for $\bar{u}$.
2. **Coefficient update** — Newton iteration for $Y$ including noise injection
   ($\mathrm{d}W$).
3. **Basis update** — bordered linear solve for $V$ followed by
   M-orthonormalisation.

---

## Code Structure

```
do-method-cpp/
├── DO/                        # Model-independent stochastic solver core
│   ├── timedependent_do.cpp   # Driver for Burgers & QG (direct Newton)
│   ├── timedependent_swe.cpp  # Driver for SWE (NOX/LOCA adaptive stepper)
│   ├── BasisSolver.hpp/cpp    # Stochastic basis (V) solver
│   ├── StochSys.hpp/cpp       # Stochastic coefficient (Y) solver (Y_Stoch)
│   ├── LinearSolverWrapper.*  # Amesos / Amesos2 / Belos back-end abstraction
│   ├── StochasticState.hpp    # Shared-state struct coupling V and Y solvers
│   ├── NoiseGenerator.hpp     # Gaussian RNG (TRNG or Boost back-end)
│   ├── StochIO.hpp            # MatrixMarket I/O helpers
│   ├── DOTimeLoop.hpp         # Inline time-stepping helpers (runStochStep, …)
│   ├── DOUtils.hpp            # printnormMV and other small utilities
│   └── Interface.hpp          # Compile-time type alias: Mean → Burger/QG
│
├── Burger/                    # 1-D Burgers mean-field solver
│   ├── Mean.hpp/cpp           # Coordinator: Newton + θ-stepper
│   ├── PDEAssembler.hpp/cpp   # FD stencil operators, RHS & Jacobian assembly
│   └── ForcingProvider.hpp    # Time-dependent W(t) (header-only)
│
├── QG/                        # Quasi-Geostrophic mean-field solver
│   ├── Mean.hpp/cpp           # Coordinator: Newton (η²-step) + line-search
│   ├── PDEAssembler.hpp/cpp   # Wrapper around QG::QG (Fortran/C++ core)
│   ├── ForcingProvider.hpp    # Time-invariant W (header-only)
│   ├── QG.hpp/cpp             # QG spatial discretisation core
│   └── QG_Data.hpp/cpp        # QG grid data
│
├── swe/                       # Shallow-Water Equations (FVM, Fortran + C++)
│   ├── FVM_LocaInterface.*    # LOCA/NOX model evaluator
│   ├── ImplicitTimeStepper.*  # Adaptive implicit stepper
│   └── …                      # Fortran FVM kernels
│
├── general/                   # Shared infrastructure
│   ├── FVM_Domain.C           # Parallel domain decomposition
│   ├── NOX_Epetra_LinearSystem_Belos.C
│   ├── GaleriExt_Utils.C
│   └── Filestreams.C
│
├── includes/                  # Shared header files
├── paramsXML/                 # Example parameter files
│   ├── BurgerRun/             #   params.xml + StochasticParams.xml
│   ├── qgRun/
│   └── sweRun/
└── CMakeLists.txt
```

---

## Code Flow

### Burgers / QG driver (`timedependent_do.cpp`)

```
main()
 └─ runSimulation<Mean>(Comm, timer)       // Mean = Burger::Mean or QG::Mean
     │
     ├─ Load params.xml, StochasticParams.xml
     │
     ├─ Construct Mean model
     │   ├─ PDEAssembler   (FD operators / QG::QG wrapper)
     │   ├─ ForcingProvider (stochastic forcing W)
     │   └─ LinearSolverWrapper (Amesos/Amesos2/Belos)
     │
     ├─ Construct StochasticState (shared V, Y, E[∂/∂yy], bilinearTerm callback)
     │
     ├─ Construct BasisSolver (V solver)
     │   ├─ Read or randomise initial V
     │   └─ Build (M − dt·A), symbolic factorisation
     │
     ├─ Construct Y_Stoch (coefficient solver)
     │   ├─ Allocate Y, dW, moment matrices
     │   └─ Initialise noise generator (TRNG / Boost)
     │
     ╞══ TIME LOOP ═══════════════════════════════════════
     │  for t = t_start → t_end:
     │  │
     │  ├─ model->refreshForcing(t)          // update W(t)
     │  │
     │  ├─ runStochStep()                    // (DOTimeLoop.hpp)
     │  │   ├─ Y_Stoch::StochasticIterations()   // Newton for Y
     │  │   ├─ BasisSolver::computeBlocks(dt)     // rebuild (M − dt·A)
     │  │   ├─ BasisSolver::v_stoch_init(Vold)    // factorise & solve for V
     │  │   ├─ BasisSolver::TransferNorm()         // QR → transfer R into Y
     │  │   ├─ Y_Stoch::HBilinV()                 // bilinear expectations
     │  │   ├─ Y_Stoch::computeEyyTyT()           // E[yy^T]
     │  │   └─ Y_Stoch::computeEVyVy()            // E[VyVy] for mean RHS
     │  │
     │  ├─ model->NewtonSolver()             // advance ū one step
     │  │   ├─ ThetaStepper()                // θ-method residual
     │  │   ├─ assembleJacobian()            // J(u), ThetaJac = I − dt·θ·J
     │  │   ├─ LinSolve()                    // LinearSolverWrapper
     │  │   └─ RunBackTracking()             // (if needed)
     │  │
     │  ├─ Update detA (Jacobian copy for V-solver)
     │  ├─ model->setExpVyVy(…)              // pass E[VyVy] to mean
     │  └─ saveTimestepOutputs()             // periodic V, yT, mean snapshots
     │
     ╞══ POST-PROCESSING ════════════════════════════════
     └─ PostProcess(), write final solution, timing report
```

### SWE driver (`timedependent_swe.cpp`)

Follows the same structure but replaces the direct Newton mean solve with the
NOX/LOCA implicit time stepper (`ImplicitTimeStepper`), supports adaptive
step-size control, gradual wind-forcing startup, and uses the FVM domain
decomposition via `FVM::LocaInterface`.

### Compile-time model dispatch

The three executables share the **same** DO solver code.  The model PDE is
selected at compile time via preprocessor definitions set in `CMakeLists.txt`:

| Macro | `BRGRDO` | `QGDO` | `SWEDO` |
|-------|----------|--------|---------|
| `brgr` | 1 | 0 | 0 |
| `quasi_geo` | 0 | 1 | 0 |
| `need_locaInterface` | 0 | 0 | 1 |

`Interface.hpp` uses these macros to resolve the `Mean` type alias and
include the correct model headers.

---

## Class Diagram (Simplified)

```
┌─────────────────────────────────────────────────────────────┐
│                    timedependent_do.cpp                      │
│              runSimulation<MeanType>(Comm, timer)            │
└────────┬───────────────────┬───────────────────┬────────────┘
         │                   │                   │
         ▼                   ▼                   ▼
┌─────────────────┐  ┌──────────────┐   ┌──────────────────┐
│   Mean (Burger   │  │ BasisSolver  │   │    Y_Stoch       │
│    or QG)        │  │  (V solver)  │   │ (coeff. solver)  │
│                  │  │              │   │                  │
│ ┌──────────────┐ │  │ LinearSolver │   │  NoiseGenerator  │
│ │PDEAssembler  │ │  │   Wrapper    │   │                  │
│ ├──────────────┤ │  └──────┬───────┘   └────────┬─────────┘
│ │ForcingProvider│ │         │                    │
│ ├──────────────┤ │         └────────┬───────────┘
│ │LinearSolver  │ │                  │
│ │  Wrapper     │ │                  ▼
│ └──────────────┘ │       ┌──────────────────┐
└──────────────────┘       │ StochasticState   │
                           │  (shared V, Y,    │
                           │   E[∂/∂yy],       │
                           │   bilinearTerm)   │
                           └──────────────────┘
```

---

## Dependencies

| Library | Purpose | Required |
|---------|---------|----------|
| [Trilinos](https://trilinos.github.io/) (≥ 13.x) | Epetra, Amesos/Amesos2, Belos, NOX/LOCA, Anasazi, Ifpack, EpetraExt | Yes |
| [HYMLS](https://github.com/BIMAU/hymls) | Matrix I/O utilities, domain decomposition helpers | Yes |
| OpenMPI | Parallel communication | Yes |
| [Boost](https://www.boost.org/) *or* [TRNG](https://www.numbercrunch.de/trng/) | Gaussian random number generation | One of them |
| [Doxygen](https://www.doxygen.nl/) | API documentation generation | Optional |

Trilinos must be built with at least the following packages enabled:
`Epetra`, `EpetraExt`, `Amesos`, `Amesos2`, `Belos`, `Ifpack`, `Anasazi`,
`AztecOO`, `NOX`, `LOCA`, `ML`, `Galeri`, `Teuchos`, `Kokkos`.

---

## Build

Set the following environment variables:

```bash
export TRILINOS_HOME=/path/to/trilinos/install
export HYMLS_HOME=/path/to/hymls/install
# If using TRNG instead of Boost:
export TRNG_HOME=/path/to/trng/install
```

Then:

```bash
mkdir build && cd build
CXX=mpicxx cmake ..
```

Build individual targets:

```bash
make BRGRDO    # Stochastic Burgers
make QGDO      # Stochastic Quasi-Geostrophic
make SWEDO     # Stochastic Shallow-Water
make           # All three
```

### Generate API documentation

```bash
make docs      # Requires Doxygen; output in build/html/
```

---

## Running

Each target needs two XML parameter files in the working directory:

| File | Contents |
|------|----------|
| `params.xml` | Mean-solver parameters (grid size, viscosity, θ, solver package, time stepping) |
| `StochasticParams.xml` | DO solver parameters (basis dimension, stochastic iterations, forcing strength, sub-time-steps, Newton tolerances) |

The `paramsXML/` directory contains ready-to-use examples:

```bash
# Example: run the Burgers problem
mkdir run && cd run
cp ../paramsXML/BurgerRun/* .
ln -s ../build/BRGRDO brgrdo
# Solver XML files (e.g. Amesos_Klu.xml, KLU2.xml) must also be present
./brgrdo
```

For parallel runs (Burgers and SWE):

```bash
mpirun -np 4 ./brgrdo
```

### Key output files

| File pattern | Contents |
|--------------|----------|
| `MeanSol_<t>.text` | Mean-field snapshot |
| `mean_<t>.mm` | Mean-field (MatrixMarket) |
| `v_<t>.mm` | Stochastic basis V |
| `yT_<t>.mm` | Coefficient matrix Y^T |
| `tsEyy_<t>.mm` | Buffered E[yy^T] snapshots |
| `jac.mm`, `mass.mm` | Initial Jacobian and mass matrix |
| `timeProf.yml` | Teuchos timing profile (when enabled) |
| `time_dt_norm_Evyvy.txt` | Per-step residual norms |

---

## Parameter Reference

### `params.xml` — Model Parameters

| Parameter | Type | Burger default | QG default | Description |
|-----------|------|---------------|------------|-------------|
| `nx` | int | 128 | 200 | Grid points in x |
| `ny` | int | — | 200 | Grid points in y (QG only) |
| `mu` | double | 0.005 | — | Viscosity (Burgers only) |
| `Reynolds Number` | double | — | 40.0 | Reynolds number (QG only) |
| `theta` | double | 0.5 | 0.5 | θ-method parameter (0.5 = Crank–Nicolson) |
| `Solver Package` | string | `Amesos` | `Amesos` | `Amesos`, `Amesos2`, or `Belos` |
| `No. of vectors in stoch. forcing` | int | 1 | 2 | Columns of W |

### `StochasticParams.xml` — Stochastic Solver

| Parameter (under `StochModel`) | Type | Description |
|------|------|-------------|
| `Use Stochastic` | bool | Enable/disable stochastic components |
| `StochFrc Strength` | double | Scaling factor for noise forcing |
| `Max. Stoch subspace Dimension` | int | Number of DO basis modes m |
| `Stochastic Iterations` | int | Number of stochastic realisations S |
| `Number of Sub Time Steps` | int | Sub-stepping within dt |
| `Max Num of Iter` | int | Newton iteration limit (Y solver) |
| `Use BackTracking` | bool | Enable backtracking in Y Newton |
| `Tolerance of RHS` | double | Newton convergence tolerance |

---

## License

See repository for license information.
