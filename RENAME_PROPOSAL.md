# Rename Proposal — DO Method Codebase

> **How to use this document:**
> - Each table row has a checkbox column. Mark `[x]` to approve or add a comment.
> - Add your comments/alternatives in the **Your Comment** column.
> - Rows marked `⭐` are high-priority (clarity gain is large).
> - After review, return this file and I will implement the approved renames.

---

## 1. Class Renames

| # | Current | Proposed | File(s) | Rationale | Approve? | Your Comment |
|---|---------|----------|---------|-----------|----------|--------------|
| 1.1 ⭐ | `Y_Stoch` | `CoeffSolver` | StochSys.hpp/cpp, timedependent_do.cpp, timedependent_swe.cpp, DOTimeLoop.hpp | Aligns with `BasisSolver`; "Y_Stoch" is a math-variable name, not a class role | [x] | |
| 1.2 ⭐ | `Mean` (Burger) | `MeanSolver` | Burger/Mean.hpp/cpp, Interface.hpp, timedependent_do.cpp | Aligns with `BasisSolver` / `CoeffSolver`; `Mean` sounds like a data container | [x] | |
| 1.3 ⭐ | `Mean` (QG) | `MeanSolver` | QG/Mean.hpp/cpp, Interface.hpp, timedependent_do.cpp | Same as above; both live in different directories so no collision | [x] | |

> **Note:** If classes are renamed, the corresponding filenames (StochSys.hpp → CoeffSolver.hpp, Mean.hpp → MeanSolver.hpp) would also change. Mark below if you want filename renames too.

| # | Current File | Proposed File | Approve? | Your Comment |
|---|-------------|---------------|----------|--------------|
| 1.4 | StochSys.hpp / StochSys.cpp | CoeffSolver.hpp / CoeffSolver.cpp | [x] | |
| 1.5 | Burger/Mean.hpp / Mean.cpp | Burger/MeanSolver.hpp / MeanSolver.cpp | [x] | |
| 1.6 | QG/Mean.hpp / Mean.cpp | QG/MeanSolver.hpp / MeanSolver.cpp | [x] | |

---

## 2. Expectation-Value Variables (`Exp_` → `E` Prefix)

Convention: `E` + subscript variables, no underscores between E and subscript.
"D" reads as "divided by" (i.e. multiplied by inverse).

> **Doxygen comments**: Every renamed E-variable will get an inline Doxygen comment
> with its math meaning, e.g. `Teuchos::RCP<Epetra_MultiVector> Ezy; ///< $E[zy^T]$, cross-covariance (nW × m).`
> **Trailing `_`**: Private members get `_`, public members don't.

### 2.1 Y_Stoch (CoeffSolver)

| # | Current | Proposed | Type | Math Meaning | Approve? | Your Comment |
|---|---------|----------|------|-------------|----------|--------------|
| 2.1a ⭐ | `Exp_zy` | `Ezy` | `RCP<MultiVector>` | $E[zy^T]$ current step | [x] | We should comment somewhere which explain the math meaning for all the variables in this section and other (should it be doxygen or normal coment?).Also be careful if a variable is private or not. We nee to differentiate (using '_'?).  |
| 2.1b ⭐ | `Exp_zy_` | `EzyPrev` | `RCP<MultiVector>` | $E[zy^T]$ previous step | [x] | |
| 2.1c ⭐ | `Exp_yy_` | `Eyy` | `RCP<MultiVector>` | $E[yy^T]$ covariance (m×m) | [x] | |
| 2.1d | `Exp_VyVyy` | `EVyVyy` | `RCP<MultiVector>` | $E[\langle Vy, Vy \rangle y^T]$ | [x] | |
| 2.1e | `ExpzyDExpyy` | `EzyDEyy` | `RCP<MultiVector>` | $E[zy^T] \cdot (E[yy^T])^{-1}$ | [x] | |
| 2.1f | `ExpVyVyyDExpyy` | `EVyVyyDEyy` | `RCP<MultiVector>` | $E[VyVyy^T] \cdot (E[yy^T])^{-1}$ | [x] | |
| 2.1g ⭐ | `ExpDExpyy` | `EDEyy` | `RCP<MultiVector>` | Basis coupling: $(n_{sub} W E[zy^T] + \Delta t\, E[\langle Vy,Vy\rangle y^T]) / E[yy^T]$. Appears in $F_v = \Delta t\, AV + \texttt{EDEyy}$. | [x] | |
| 2.1h | `ExpYY` | `EYY` | `RCP<CrsMatrix>` | CrsMatrix form of $E[YY^T]$ | [x] | |
| 2.1i | `ExpVar` | — | `RCP<MultiVector>` | Variance workspace | — | **Moved to §7 (dead code)**: allocated but never used anywhere. Remove instead of rename. |
| 2.1j | `LocExpyyy` | `LocEyyy` | `RCP<MultiVector>` | Local $E[yyy^T]$ accumulator | [x] | |
| 2.1k | `GlobExpyyy` | `GlobEyyy` | `RCP<MultiVector>` | Global $E[yyy^T]$ | [x] | |
| 2.1l | `EyyTyT` | `EyyTyT` | `RCP<MultiVector>` | Already uses E-prefix — **keep** | [x] | |
| 2.1m | `exp_yy_` (private) | `EyyOld_` | `RCP<MultiVector>` | Previous-step $E[yy^T]$ — trailing `_` because private | [ ] |
| 2.1n | `rep_exp_vyvy` | `repEVyVy` | `RCP<MultiVector>` | Replicated $E[VyVy]$ for all samples | [[x]] | |
| 2.1o ⭐ | `expv4` | `EVyVy` | `RCP<Vector>` | $E[\langle Vy,Vy \rangle]$ diagonal vector | [[x]] | |
| 2.1p | `Vexpv4` | `VEVyVy` | `RCP<MultiVector>` | V times $E[\langle Vy,Vy \rangle]$ (if used) | [[x]] | |
| 2.1q | `map_expv4` | remove or `mapEVyVy` | `RCP<Map>` | Map for expv4 — may be unused | [x] | |

### 2.2 BasisSolver

| # | Current | Proposed | Type | Math Meaning | Approve? | Your Comment |
|---|---------|----------|------|-------------|----------|--------------|
| 2.2a | `Exp_zy` | `Ezy` | `RCP<MultiVector>` | $E[zy^T]$ | [x] | |
| 2.2b | `Exp_yy` | `Eyy` | `RCP<MultiVector>` | $E[yy^T]$ | [x] | |
| 2.2c | `Exp_yyy` | `Eyyy` | `RCP<MultiVector>` | $E[yyy^T]$ | [x] | |
| 2.2d | `exp_yy_inv` | `EyyInv` | `RCP<MultiVector>` | $(E[yy^T])^{-1}$ | [x] | |
| 2.2e | `expv3` | `EVyVyBasis` | `RCP<MultiVector>` | $E[\langle Vy,Vy \rangle]$ in basis context | [x] | |

### 2.3 StochasticState

| # | Current | Proposed | Type | Approve? | Your Comment |
|---|---------|----------|------|----------|--------------|
| 2.3a | `ExpDExpyy` | `EDEyy` | `RCP<MultiVector>` | [x] | |

### 2.4 Mean (Burger & QG) and PDEAssembler

| # | Current | Proposed | Type | Approve? | Your Comment |
|---|---------|----------|------|----------|--------------|
| 2.4a | `ExpVyVy_` | `EVyVy_` | `RCP<Vector>` | [x] | |
| 2.4b | `setExpVyVy()` | `setEVyVy()` | method | [x] | |
| 2.4c | `ExpVyVy` param in `assembleRHS()` | `EVyVy` | param | [x] | |

### 2.5 DOTimeLoop.hpp / timedependent_do.cpp

| # | Current | Proposed | Context | Approve? | Your Comment |
|---|---------|----------|---------|----------|--------------|
| 2.5a | `getExpDExpyy()` | `getEDEyy()` | Y_Stoch accessor | [x] | |
| 2.5b | `getEyy()` | `getEyy()` | Already short — **keep** | [x] | |
| 2.5c | `getEVyVy()` | `getEVyVy()` | Already E-prefix — **keep** | [x] | |

---

## 3. Cryptic / Single-Letter Variable Names

### 3.1 Y_Stoch (CoeffSolver)

| # | Current | Proposed | Type | What It Stores | Approve? | Your Comment |
|---|---------|----------|------|---------------|----------|--------------|
| 3.1a ⭐ | `Vnew` | `V` | `RCP<MultiVector>` | Stochastic basis (not "new"). **Public** member — no trailing `_`. | [x] | |
| 3.1b ⭐ | `B` | `W` | `RCP<MultiVector>` | Stochastic forcing basis — same as `W_` in BasisSolver | [x] | |
| 3.1c ⭐ | `udet` | `uMean` | `RCP<Vector>` | Mean-field solution vector | [x] | |
| 3.1d ⭐ | `ttt` | remove or `debugCounter` | `int` | Initialised to 1, purpose unclear | [x] | |
| 3.1e | `viv` | `bilinWork` | `RCP<MultiVector>` | Workspace for B(V_i, V_j) output | [x] | |
| 3.1f | `rszyy` | `sizeYY` | `int` | Stores m_*m_ | [x] | |
| 3.1g | `MyLDA` | `numSamples` | `int` | Local sample count = yTrans_->MyLength() | [x] | |
| 3.1h | `N_` | `nDOF` | `int` | Number of spatial DOFs = `uMean->GlobalLength()`. **Public**, no trailing `_`. Doxygen: `///< Total spatial degrees of freedom (grid points × variables per point).` | [x] | |
| 3.1i | `H` | `BilinTensor` | `RCP<MultiVector>` | Bilinear interaction tensor $B(V_i,V_j)$ projected to coefficient space (m × m²). **Not** a Hessian — it's the quadratic coupling between basis modes. | [x] |  |
| 3.1j | `Hn` | `BilinTensorN` | `RCP<MultiVector>` | Same tensor in physical space (N × m²). Doxygen will explain: `BilinTensor = V^T · BilinTensorN`, i.e. `BilinTensorN` is computed first via `Bilinear(V_i, V_j)`, then projected to get `BilinTensor`. | [x] | |
| 3.1k | `VHn` | `VtBilinTensorN` | `RCP<MultiVector>` | $V^T \cdot$ BilinTensorN | [x] | |
| 3.1l | `Rv` | `EyyInvView` | `RCP<MultiVector>` | $(E[yy^T])^{-1}$ as MultiVector — all three (`Rv`/`Ru`/`ru`) alias the **same memory** in different views. `View` suffix clarifies it's a view, not a copy. | [x] | |
| 3.1m | `Ru` | `EyyInvTeuchos` | `RCP<SerialDenseMatrix>` | $(E[yy^T])^{-1}$ as Teuchos dense matrix view | [x] | |
| 3.1n | `ru` | `EyyInvEpetra` | `SerialDenseMatrix` | $(E[yy^T])^{-1}$ as Epetra dense matrix view | [x] | |
| 3.1o | `eye` | `identity` | `RCP<MultiVector>` | m×m identity matrix | [x] | |
| 3.1p | `Utmp` | `workMxM` | `RCP<MultiVector>` | Temporary m×m workspace. `work` is standard in numerical LA (cf. LAPACK WORK arrays). | [x] | |
| 3.1q | `Vtemp` | `Vwork` | `RCP<MultiVector>` | Temporary copy of V. Same `work` convention. | [x] | |
| 3.1r | `YY` | `yyOuter` | `RCP<MultiVector>` | $y \otimes y$ outer-product workspace (m² × S). **S** = `numSamples` (the number of stochastic realisations, formerly `NumStochIter`). Each column stores the flattened m×m outer product $y_s y_s^T$ for one sample $s$, hence m² rows × S columns. | [x] | |
| 3.1s | `f_out` | `residual` | `RCP<MultiVector>` | Newton residual output | [x] | |
| 3.1t | `itrtr` | `debugWienerCount` | `int` | Counter for debug Wiener-process file reading | [x] | |

### 3.2 Newton Iterates in Y_Stoch

| # | Current | Proposed | What It Stores | Approve? | Your Comment |
|---|---------|----------|---------------|----------|--------------|
| 3.2a | `x0_` | `yOld` | Previous time-step Y values. **Public** — no trailing `_`. | [x] | |
| 3.2b | `x_` | `yCurr` | Current Newton iterate for Y. **Public** — no trailing `_`. | [x] | |
| 3.2c | `dx_` | `yDelta` | Newton correction ΔY. **Public** — no trailing `_`. | [x] | |
| 3.2d | `jac` | `jacView` | MultiVector **view** of the dense Jacobian (aliases `jacDense` memory). Avoids MV=Mass×V ambiguity. | [x] | |
| 3.2e | `Yjac` | `jacDense` | Dense Jacobian (SerialDenseMatrix) | [x] | |

### 3.3 BasisSolver

| # | Current | Proposed | Type | What It Stores | Approve? | Your Comment |
|---|---------|----------|------|---------------|----------|--------------|
| 3.3a | `V1` | `Vold` | `RCP<MultiVector>` | Previous-step basis | [x] | |
| 3.3b | `z` | `rhsWork` | `RCP<MultiVector>` | RHS workspace | [x] | |
| 3.3c | `Rvec` | `Rfactor` | `RCP<MultiVector>` | QR R-factor | [x] | |
| 3.3d | `eye` | `identity` | `RCP<MultiVector>` | Identity matrix | [x] | |
| 3.3e | `stochiter` | `stochIter` | `int` | Public iteration counter (add trailing `_`?) | [x] | |

---

## 4. Map Variable Names

| # | Current | Proposed | Class | What It Maps | Approve? | Your Comment |
|---|---------|----------|-------|-------------|----------|--------------|
| 4.1 | `map_x_` | `localMapY` | Y_Stoch | Local map of size m (Y coefficients) | [x] | |
| 4.2 | `map_z_` | `localMapZ` | Y_Stoch | Local map of size nW (noise) | [x] | |
| 4.3 | `map_mm` | `localMapYY` | Y_Stoch | Local map of size m*m | [x] | |
| 4.4 | `Tmap` | `stochMap` | Y_Stoch | Distributed map of size NumStochIter (samples) | [x] | |
| 4.5 | `map_expv4` | remove or `mapN` | Y_Stoch | Map for EVyVy — may be unused | [x] | |
| 4.6 | `map_` | `mapDOF` | BasisSolver | Distributed map of size nDOF (spatial DOFs). Matches the `nDOF` rename in CoeffSolver. Used for `identity`, `Eyy`, `Eyyy`, `EyyInv`. | [x] | |
| 4.7 | `y_map` | `noiseMap` | BasisSolver | Distributed map of size nW (number of forcing vectors). Despite the name, it's **not** a map for y — it's used for `Ezy` which has nW rows. Different from `localMapY` (size m). | [x] | |
| 4.8 | `locmap_` | `localMapM` | BasisSolver | Local map of size m. Same size as `localMapY` in CoeffSolver but `Epetra_LocalMap` vs distributed. Used for `Rfactor`. Could merge with a CoeffSolver map if both classes ever share, but keep separate for now. | [x] | |

---

## 5. Method Names (Inconsistent Casing / Unclear)

### 5.1 Y_Stoch Methods

| # | Current | Proposed | Rationale | Approve? | Your Comment |
|---|---------|----------|-----------|----------|--------------|
| 5.1a ⭐ | `y_rhs()` | `assembleRHS()` | snake_case → PascalCase; verb-first | [x] | |
| 5.1b ⭐ | `y_jac()` | `assembleJacobian()` | Same | [x] | |
| 5.1c | `NonLinRHS()` | `computeNonlinearRHS()` | Add verb prefix | [x] | |
| 5.1d | `LinCoeff()` | `computeLinearCoeff()` | Add verb prefix | [x] | |
| 5.1e | `HBilinV()` | `computeBilinTensor()` | Matches the variable rename `H`→`BilinTensor`. Computes the bilinear interaction tensor $B(V_i, V_j)$ and its projections. | [x] | |
| 5.1f | `jacBilinTerm()` | `computeBilinearJacobian()` | Verb-first, clearer | [x] | |
| 5.1g | `computeExpVVyVy()` | `computeEVyVy()` | Fix double-V typo, use E-prefix | [x] | |
| 5.1h | `computeExpVal()` | `computeExpectations()` | "Val" is vague | [x] | |
| 5.1i | `CreateLocMultiVec()` | `CreateLocalMultiVec()` | Spell out "Local" | [x] | |
| 5.1j | `CreateDistTransMultivec()` | `localToDistributed()` | Describes the operation: copies local-map Y/Z into distributed-map transposed layout. Paired with `CreateLocalMultiVec()` which goes the other direction. | [x] | |
| 5.1k | `PostProcess()` | `PostProcess()` | **Keep as-is.** It does more than second-moment: SVD of $YY^T$, writes deviations/variances to files, *and* computes the second moment. `PostProcess` is an accurate umbrella name. | [x] | |
| 5.1l | `RunBackTracking()` | `runBackTracking()` | Lowercase verb prefix for consistency | [x] | |
| 5.1m | `printTransNormMV()` | `printTransposedNorms()` | Drop the abbreviations: "Transposed" is clearer than "Trans", and the parameter type already tells you it's a MultiVector. | [x] | |

### 5.2 Y_Stoch Private Compute Methods

| # | Current | Proposed | Rationale | Approve? | Your Comment |
|---|---------|----------|-----------|----------|--------------|
| 5.2a | `computeEVyVy()` | `computeEVyVy()` | Already good with E-prefix | [x] | |
| 5.2b | `computeCrossVariance()` | `computeCrossVariance()` | Already clear | [x] | |
| 5.2c | `computeExpDExpyy()` | `computeEDEyy()` | E-prefix consistency | [x] | |
| 5.2d | `computeEVyVyy()` | `computeEVyVyy()` | Already uses E-prefix | [x] | |
| 5.2e | `computeEyyTyT()` | `computeEyyTyT()` | Already uses E-prefix | [x] | |

### 5.3 BasisSolver Methods

| # | Current | Proposed | Rationale | Approve? | Your Comment |
|---|---------|----------|-----------|----------|--------------|
| 5.3a | `getExp_yy()` | `getEyy()` | E-prefix consistency | [v] | |
| 5.3b | `get_y()` | `getY()` | Remove snake_case | [x] | |
| 5.3c | `set_y()` | `setY()` | Remove snake_case | [x] | |
| 5.3d | `get_frcStrenth()` | `getFrcStrength()` | Fix typo "Strenth" | [x] | |
| 5.3e | `set_frcStrength()` | `setFrcStrength()` | Already correct spelling | [x] | |

### 5.4 Mean Methods

| # | Current | Proposed | Rationale | Approve? | Your Comment |
|---|---------|----------|-----------|----------|--------------|
| 5.4a | `get_dim_W()` | `getDimW()` | Remove snake_case. `W` is used consistently: `W_` in BasisSolver, `W` in CoeffSolver (renamed from `B`), `W` as constructor param. It always means the stochastic forcing basis. | [x] | |
| 5.4b | `get_W()` | `getW()` | Same — `W` is the universal name for stochastic forcing in this codebase. | [x] | |

---

## 6. Constructor Parameter Names

| # | Current | Proposed | Class | Rationale | Approve? | Your Comment |
|---|---------|----------|-------|-----------|----------|--------------|
| 6.1 | `NumStochIter` | `numSamples` | Y_Stoch | It's the sample count S, not "iterations" | [x] | |
| 6.2 | `num_Subtime_Step` | `numSubSteps` | Y_Stoch | Mixed snake_case/CamelCase | [x] | |
| 6.3 | `Wb` | `W` | Y_Stoch | Matches BasisSolver convention. `W` = stochastic forcing basis everywhere. | [x] | |
| 6.4 | `uav` | `uMean` | Y_Stoch | Clearer than "u-average" | [x] | |
| 6.5 | `NormRHS` | `initialRHSNorm` | Y_Stoch | Distinguish from tolerance | [x] | |
| 6.6 | `numBackTrackingSteps` | `numBackTrackingSteps` | Y_Stoch | Fine but type is `double` — should it be `int`? | [x] | I agree |

---

## 7. Possibly Unused / Dead Members

These appear to be declared but never (or rarely) used. Confirm and remove?

| # | Name | Class | Type | Approve Removal? | Your Comment |
|---|------|-------|------|-----------------|--------------|
| 7.1 | `p_` | Y_Stoch | `RCP<Vector>` | [x] | ✅ Verified: declared in .hpp only, never assigned or read in any .cpp ||
| 7.2 | `W_graph_` | Y_Stoch | `RCP<CrsGraph>` | [x] | ✅ Verified: declared in .hpp only, never used |
| 7.3 | `stress_` | Y_Stoch | `RCP<CrsMatrix>` | [x] | ✅ Verified: declared in .hpp only, not used by SWE or any other target |
| 7.4 | `isInitialized_` | Y_Stoch | `bool` | [x] | ✅ Verified: declared only, never assigned or tested |
| 7.5 | `showGetInvalidArg_` | Y_Stoch | `bool` | [x] | ✅ Verified: declared only, never used |
| 7.6 | `numPrecRecomputes_` | Y_Stoch | `int` | [x] | ✅ Verified: declared only, never used |
| 7.7 | `sig_` | Y_Stoch | `double` | [x] | ✅ Verified: declared only, never assigned or read |
| 7.8 | `row_one` | Y_Stoch | `RCP<MultiVector>` | [x] | ✅ Verified: declared only, never allocated or used |
| 7.9 | `Vexpv4` | Y_Stoch | `RCP<MultiVector>` | [x] | ✅ Verified: declared only, never allocated or used |
| 7.10 | `map_expv4` | Y_Stoch | `RCP<Map>` | [x] | ✅ Verified: declared only, never allocated or used |
| 7.11 | `sol` | Y_Stoch | `RCP<SerialDenseMatrix>` | [x] | ✅ Verified: declared only, never allocated or used |
| 7.12 | `ttt` | Y_Stoch | `int` | [x] | ✅ Verified: initialised to 1 in init-list, never read anywhere |
| 7.13 | `ExpVar` | Y_Stoch | `RCP<MultiVector>` | [x] | ✅ Verified: allocated in constructor but never read or written to |

---

## 8. Naming Convention Summary (For Reference)

| Category | Convention | Example |
|----------|-----------|---------|
| Classes | `PascalCase`, role-based noun | `BasisSolver`, `CoeffSolver`, `MeanSolver` |
| Public methods | `camelCase` with verb prefix | `computeBlocks()`, `assembleRHS()` |
| Public members | `camelCase`, no trailing `_` | `Ezy`, `uMean`, `V` |
| Private members | `camelCase_` with trailing `_` | `sharedState_`, `iter_`, `domain_` |
| Expectation values | `E` + subscript, no underscores | `Eyy`, `Ezy`, `EzyDEyy` |
| Maps | descriptive + `Map` suffix | `localMapY`, `stochMap` |
| Booleans | `is`/`use`/`enable` prefix | `isConverged_`, `useNewton_` |
| Workspace temps | descriptive noun | `bilinWork`, `workMxM` |

---

**Please review and return with your approvals/comments. I will implement only the approved rows.**
