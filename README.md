# rheoSimpleInterFoam

`rheoSimpleInterFoam` is a custom OpenFOAM solver based on the volume of fluid (VoF) method (extended from `interFoam`). It incorporates elastic extra stress (`tauElastic`) computed via an algebraic second-order fluid model (Criminale-Ericksen-Filbey / CEF constitutive model).

## 1. Overview
* **Solver Name**: `rheoSimpleInterFoam`
* **Target Binary**: `$(FOAM_USER_APPBIN)/rheoSimpleInterFoam`
* **Governing Equations**: Adds the divergence of the elastic extra stress tensor `+ fvc::div(tauElastic)` to the right-hand side of the momentum equation (`UEqn.H`).
* **Constitutive Model**: Second-order fluid (CEF) algebraic approximation using Rivlin-Ericksen tensors ($A_1, A_2$) and normal stress coefficients ($\Psi_1, \Psi_2$), with phase-masking for two-phase interfaces.

---

## 2. Directory Structure
```text
my/
├── README.md                 # This compilation and user guide
└── rheoSimpleInterFoam/      # Custom solver directory
    ├── Make/
    │   ├── files             # Source files and target binary definitions
    │   └── options           # Include paths and linked libraries
    ├── calcTauElastic.H      # Algebraic calculation of elastic stress tensor
    ├── createFields.H        # Field allocations and transport property reads (Psi1, Psi2, tauElastic)
    ├── UEqn.H                # Momentum predictor equation with elastic stress divergence
    ├── rheoSimpleInterFoam.C # Main solver source file (entry point)
    └── *.H                   # Header files for VoF/MULES multi-phase algorithm
```

---

## 3. Prerequisites
1. **OpenFOAM Environment**:
   Ensure that your OpenFOAM environment variables are loaded in your current shell:
   ```bash
   # Example: OpenFOAM v2312
   source /usr/lib/openfoam/openfoam2312/etc/bashrc
   # Or from a local installation
   source ~/OpenFOAM/OpenFOAM-v2312/etc/bashrc
   ```

2. **Verify Environment Variables**:
   Confirm that the OpenFOAM project and user application binary directories are set:
   ```bash
   echo $WM_PROJECT_DIR      # OpenFOAM installation directory
   echo $FOAM_USER_APPBIN    # User binary output directory
   ```

3. **Ensure Output Directory Exists**:
   ```bash
   mkdir -p $FOAM_USER_APPBIN
   ```

---

## 4. Compilation Instructions

1. **Navigate to the Solver Directory**:
   ```bash
   cd rheoSimpleInterFoam
   ```

2. **Clean Previous Build Artifacts (Recommended)**:
   ```bash
   wclean
   ```

3. **Compile the Solver**:
   ```bash
   wmake
   ```
   To enable parallel compilation, use the `-j` flag:
   ```bash
   wmake -j
   ```

4. **Verify Installation**:
   Check that the binary was successfully installed in `$FOAM_USER_APPBIN`:
   ```bash
   which rheoSimpleInterFoam
   # Example output: $FOAM_USER_APPBIN/rheoSimpleInterFoam
   ```
   Run the help flag to verify executable execution:
   ```bash
   rheoSimpleInterFoam -help
   ```

---

## 5. Notes & Configuration

* **Main Source File (`rheoSimpleInterFoam.C`)**:
  * `Make/files` specifies `rheoSimpleInterFoam.C` as the primary compilation unit. Ensure this file is present in the solver root directory before running `wmake`.
* **Required Transport Properties**:
  * Case configurations require definition of first and second normal stress difference coefficients in `constant/transportProperties` (or within initial field dictionaries):
    * `Psi1`: First normal stress difference coefficient $[Pa \cdot s^2] = [kg/m]$
    * `Psi2`: Second normal stress difference coefficient $[Pa \cdot s^2] = [kg/m]$
* **OpenFOAM Version Compatibility**:
  * Include paths and library names in `Make/options` (such as `twoPhaseInter`, `immiscibleIncompressibleTwoPhaseMixture`) may vary depending on the specific OpenFOAM version and release track (OpenFOAM Foundation vs. OpenCFD / ESI). Adjust `Make/options` accordingly if header include errors occur during compilation.