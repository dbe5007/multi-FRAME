# CANLab-Multivariate-Pipeline
This package is a data processing pipeline for performing multivariate imaging analyses using fMRI data. Currently the pacakge can can compute  Multivoxel Pattern Analysis (MVPA), Represetational Similarity Analysis (RSA), and Encoding-Retireval Similarity analysis (ERS). The package requires:

* [MATLAB R2017b](https://www.mathworks.com/products/matlab.html) (Working in >=2019a functionality)
* [SPM12](https://www.fil.ion.ucl.ac.uk/spm/)
* [CoSMoMVPA toolbox](http://www.cosmomvpa.org/) (Also available [on Github](https://github.com/CoSMoMVPA/CoSMoMVPA))

## Function Description

```createParams.m```
Creates a parameter MAT file with all relevant information for preprocessing model specification/estimation, and multivariate analyses. Also includes querying location and registration of regions of interest (ROIs) to be used for multivaraite analyses. Parameter file is named using processing choices made within the script to aid in file identification.

```estimateModel.m```
Estimates a single trial model for every trial of interest using SPM12. The resulting model assings a beta value for each trial to be used for multivariate analysis. Models are gunzipped to space space and will not impede classification analysis. 
Requires SPM

```preprocessData.m```
Performs 1st level preprocessing on functional data, including registration of functional images, slice scan time correction, and motion calculation. 

```runMVPAClassification.m```
Runs MVPA within singular ROIs or a Searchlight within mask provided by user. Currently employs SVM classifier for ROI level, and LDA for searchlight. Output is saved as tidyverse-formatted CSV file and as MAT files. Searchlight results are saved as weights NiFTi files. Summary CSV file is generated as well, concatnating all subject accuracies. 
Requires CoSMoMVPA

```runRSA.m```
Runs RSA within singular ROIs, Searchlight within mask provided by user, or ERS analyses. Output is saved as tidyverse-formatted CSV file and as MAT files. Summary CSV file is generated as well, concatnating all subject accuracies. 
Requires CoSMoMVPA

```specifyModel.m```
Creates single trial models for functional data across all runs. Requires behavioral file to be sourced for onsets and trial tag information

## Running the pipeline
THe package is currently equipped to either preprocess data in SPM12 or take data resulting from fMRIPrep. Below is order of scripts to be run in the pipeline assuming a new data project:


1. ```createParams.m```          - Dialog box querying directories, analysis choices, and preprocessing programs. All subsequent programs require a parameter file created from this script
2. ```preprocessData.m```        - Preprocess functional data via SPM12 or organize motion regressors from fMRIPrep.
3. ```specifyModel.m```          - Specifiy SPM model with information for trials of interest (i.e. trials to use for classification/representation).
4. ```estimateModel.m```         - Estimate model with every trials as its own individual regressor.
5. ```runMVPAClassification.m``` - Perform multivariate classification using previously derived model.


## Publications
Dennis, N. A., & Overman, A. A., Gerver, C. R., McGraw, K., Rowley, M. A., & Salerno, J. M. (2019). *Different types of associative encoding evoke differential processing in both younger and older adults: evidence from univariate and multivariate analyses.* Neuropsychologia, 135. https://doi.org/10.1016/j.neuropsychologia.2019.107240

Gerver, C. R., Overman, A. A., Babu, H. J., Hultman, C. E., & Dennis, N. A. (in press). *Examining the neural basis of congruent and incongruent configural contexts during associative retrieval.* Journal of Cognitive Neuroscience.

