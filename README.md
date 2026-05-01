# AMULET: Acoustic Metastructure-Based Underwater Localization and Estimation Toolkit

**AMULET: Underwater Acoutic Metastructure for Direction-of-Arrival Estimation Underwater Using a Single Hydrophone**  
*Accepted to ACM/IEEE SenSys, 2026*  

---

## Overview

**AMULET** is an open-source research artifact accompanying our work on **acoustic metastructure-enabled Direction-of-Arrival (DoA) estimation** in underwater environments.  

This repository provides the **processed experimental dataset**, **reproducibility scripts**, and **3D metastructure designs** used to evaluate a compact, low-power approach to underwater localization that leverages **physics-driven spatial encoding** rather than large hydrophone arrays.

Our approach demonstrates that carefully designed acoustic metastructures can encode directional information directly into received signals, enabling accurate DoA estimation with minimal sensing hardware.

---

## Key Contributions

- **Metastructure-enabled DoA sensing** without large arrays  
- **Robust performance across environments** (tank, controlled indoor, and open water)  
- **Cross-environment generalization** (indoor calibration → outdoor deployment)  
- **Open dataset + full figure reproducibility pipeline**

---




### Models (`/models`)
- **STL files** for all metastructure designs used in experiments  
- Files are suitable for **direct 3D printing** (100% infill, no support)

---

### Code

#### `generate_all_results.m`
This is the **primary entry point** for reproducing results from the paper.

- Loads relevant subsets of the dataset  
- Generates all figures presented in the paper  
- Encapsulates the full evaluation pipeline  

---

### Data
- amulet_Data.mat contains all of the data within a labeled struct
- The included data is primarily pre-processed experimental signatures to drastically reduce file size, processing time, and code complexity
- Includes data from multiple environments:
  - Controlled tank experiments  
  - Indoor calibration setups  
  - Open-water (lake) deployments  

---

## Reproducing Results

### Requirements
- MATLAB (signal processing toolbox required)

### Steps

Clone the repository:
   ```bash
   git clone https://github.com/adbergey/AMULET.git
   cd AMULET
```

Or just download amulet_Data.mat and generate_all_results.m and place them in the same folder.
Then run the script to generate all of the figures from the paper.


### License

This project is licensed under the MIT License.


### Contact

For questions, issues, or collaboration inquiries, please contact adbergey@uw.edu

## Thanks for checking out our project!
