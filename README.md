# MRM_Processor

**MRM_Processor** is a graphical user interface (GUI)-based tool for automated processing of LC-MS/MS MRM data in targeted metabolomics.

This repository contains two main modules:

- **MRM Processor GUI**: used for sample data processing
- **MRM Curve**: used for calibration curve construction

---

## Functions

### MRM Processor GUI
MRM Processor GUI is designed for automated analysis of sample MRM data.  
The workflow includes:

- retention time (RT) correction
- chromatographic peak detection
- peak area integration
- analyte concentration calculation
- result export

### MRM Curve
MRM Curve is designed for calibration curve construction based on standard samples.  
The workflow includes:

- RT correction
- chromatographic peak detection
- peak area integration
- peak area ratio calculation
- calibration curve fitting
- output of regression parameters and calibration plots

---

## Input Files

Both modules require:

- **mzML files**
- an **Excel file containing MRM transition information**

Optional:

- an additional **Excel file containing manual integration ranges**

---

## Usage

### 1. MRM Curve: Calibration Curve Construction

#### Input files
Prepare:

- mzML files for calibration standards
- an Excel file containing MRM transition information
- optional manual integration file

#### Procedure
1. Click the first **Select** button to load the mzML files used for calibration curve construction.
2. The selected files will appear in the list on the left side of the button.
3. Enter the concentration of each calibration standard in the **Concentration** column.
4. Click the second **Select** button to load the Excel file containing the MRM transition information.
5. Optionally, click the third **Select** button to load the Excel file containing manual integration ranges.
6. Set the required thresholds under **Peak Picking Params**.
7. Set the required thresholds under **Curve Params**.
8. Click **Run** to start data processing.

#### Output
All output files will be saved in the folder containing the mzML files, including:

- a file containing all calibration curve data
- an integration range file for each analyte
- an individual calibration curve file for each analyte

---

### 2. MRM Processor GUI: Sample Data Analysis

#### Input files
Prepare:

- mzML files for samples
- an Excel file containing MRM transition information
- optional manual integration file

#### Procedure
1. Click the first **Select** button to load the mzML files to be analyzed.
2. The selected files will appear in the list on the left side of the button.
3. Click the second **Select** button to load the Excel file containing the MRM transition information.
4. Optionally, click the third **Select** button to load the Excel file containing manual integration ranges.
5. Set the required thresholds under **Set Params**.
6. Click **Run** to start data processing.

#### Output
All output files will be saved in the folder containing the mzML files, including:

- a file containing the overall analysis results for all samples
- an integration range file for each analyte

---

## Notes

- **MRM_Processor_Analyze** is intended for routine MRM sample data processing.
- **MRM_Processor_Curve** is intended for calibration curve construction and regression analysis.
- Manual integration can be performed by supplying an additional Excel file containing predefined integration ranges.

---

## Citation

If you use this tool in your work, please cite the corresponding publication.
