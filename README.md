# **README**

## **Overview**

This repository contains the full analysis chain used to perform 1D and 2D signal extraction based on mass and decay-length (τ<sub>z</sub>) fits.
The workflow is implemented in ROOT and RooFit and is controlled through a master steering macro.

The repository is organized into:

* A **master script** to run the full workflow
* Helper modules for building PDFs, performing fits, and plotting
* Input configuration files defining binning and initial fit parameters
* Utilities for data I/O and plotting

---

# **1. Main Workflow**

## **`InputToResults.C`**

This is the *master file* coordinating the full analysis.
Run it using:

```bash
root -l -b -q 'InputToResults.C+(ispO, isMC, caseName, remakeDS, fitMass1D, fitTauz1D, fitNpTauz1D, fit2D, plotResults)'
```

### **Arguments**

| Argument      | Type     | Description                                                             |
| ------------- | -------- | ----------------------------------------------------------------------- |
| `ispO`        | `bool`   | Use p–O collision settings                                              |
| `isMC`        | `bool`   | Run on Monte Carlo sample                                               |
| `caseName`    | `string` | Tag for input/output files, useful for systematics                      |
| `remakeDS`    | `bool`   | Rebuild the RooDataSet from reduced tables                              |
| `fitMass1D`   | `bool`   | Perform 1D mass fit                                                     |
| `fitTauz1D`   | `bool`   | Perform 1D τ<sub>z</sub> fit                                            |
| `fitNpTauz1D` | `bool`   | Perform 1D τ<sub>z</sub> non-prompt signal fit after resolution fit     |
| `fit2D`       | `bool`   | Perform 2D mass–τ<sub>z</sub> fit                                       |
| `plotResults` | `bool`   | Produce and save summary plots                                          |

This script orchestrates:

* dataset building
* PDF construction
* performing fits (1D or 2D depending on flags)
* saving results
* producing final plots

---

# **2. Module Descriptions**

## **`SignalExtraction.C`**

Extracts fit results **for a single bin**.
This contains the core logic used by `InputToResults.C` when looping over bins.

---

## **`BuildPDF.C`**

Constructs the RooFit models used in the analysis and imports all needed parameters into the `RooWorkspace`.

Responsibilities include:

* defining PDFs for mass, τ<sub>z</sub> resolution, τ<sub>z</sub> background, and combined 2D fits
* loading initial values from input files
* assigning default values for uninitialized parameters
* ensuring model conventions (e.g. `_mass`, `_tauzRes`, `_tauzBkg`, `_tauzMass` suffixes)

---

## **`inputUtils.h`**

Contains utility functions related to **input and dataset creation**, including:

* reading input paths
* loading trees
* applying bin selections
* building RooDataSets

---

## **`plotUtils.h`**

Contains utilities dedicated to **plotting and visualization**, such as:

* applying default ROOT/ATLAS/LHC-style aesthetics
* beautifying RooPlot objects
* saving canvases with consistent formatting

---

# **3. The `inputFiles/` Directory**

The `inputFiles` folder contains configuration files that control:

* input paths for data or MC (e.g. `input_data_pO.txt`)
* bin edges
* initial parameter values for the fits

These files provide the initial state for *all* fits and must follow certain conventions.

---

# **4. Input File Conventions**

### **General format**

* Parameters must be separated by **`;`**
* Bin edges are specified using:

  ```
  pt; y; cent; chi2;
  ```

  with each entry using `min-max` format.

### **FitStat**

* The `fitStat` flag determines which fits should run.
* Set `fitStat = todo` to perform a new fit.
* If results already exist for a bin, new results are added without overwriting previous ones.

### **Model names**

* Any variable starting with `model` defines which PDF model to use.
* When introducing a new model, **add it to `BuildPDF.C`**, including:

  * model definition
  * default parameter values
  * handling inside the combined fits

### **Parameter format**

Parameters can be given in two ways:

#### **1. Fixed value**

```
[x]
```

→ parameter is fixed at value *x*

#### **2. Initial value with allowed range**

```
[x, xmin, xmax]
```

→ *x* is the starting value
→ *xmin–xmax* is the allowed range during the fit

### **Unused parameters**

* Input files may include more parameters than needed
* Only parameters relevant for the chosen model are used
* Default values (listed inside `BuildPDF.C`) are applied for any missing parameter

### **Required naming conventions**

Suffixes indicate the category of the parameter:

| Suffix      | Used for                            |
| ----------- | ----------------------------------- |
| `_mass`     | Mass PDFs (1D)                      |
| `_tauzRes`  | τ<sub>z</sub> resolution models     |
| `_tauzBkg`  | τ<sub>z</sub> background models     |
| `_tauzMass` | Combined 2D mass–τ<sub>z</sub> PDFs |

These suffixes are essential for correct parameter parsing and model construction.

---

# **5. Typical Workflow**

1. Edit binning/parameters in `inputFiles/...`
2. Run:

   ```bash
   root -l -b -q 'InputToResults.C+(true, false, "MyCase", true, true, true, true, true)'
   ```
3. Results, logs, and plots appear in the output directory tagged with `caseName`
4. Repeat with different parameter files for systematics

---

# **6. Contact / Contributions**

If you add new models or utilities:

* update `BuildPDF.C`
* describe default values
* follow suffix conventions

Contributions and improvements are welcome.

---