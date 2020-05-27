# multi-FRAME
The multi-FRAME (**multi**variate-**F**mri **R**s**A** **M**vpa **E**rs) package is a data processing pipeline for performing multivariate imaging analyses using fMRI data. Currently the package can compute  Multivoxel Pattern Analysis (MVPA), Representational Similarity Analysis (RSA), and Encoding-Retrieval Similarity analysis (ERS). The package requires:

* [MATLAB R2017b](https://www.mathworks.com/products/matlab.html)
* [SPM12](https://www.fil.ion.ucl.ac.uk/spm/)
* [CoSMoMVPA toolbox](http://www.cosmomvpa.org/) (Also available [on Github](https://github.com/CoSMoMVPA/CoSMoMVPA))

## Main Package Functions

```createParams.m```
Creates a parameter MAT file with all relevant information for preprocessing model specification/estimation, and multivariate analyses. Also includes querying location and registration of regions of interest (ROIs) to be used for multivariate analyses. Parameter file is named using processing choices made within the script to aid in file identification.

```estimateModel.m```
Estimates a single trial model for every trial of interest using SPM12. The resulting model assigns a beta value for each trial to be used for multivariate analysis. Models are gzipped to save space and will not impede classification analysis. 
Requires SPM

```preprocessData.m```
Performs 1st level preprocessing on functional data, including registration of functional images, slice scan time correction, and motion calculation.

```maskRegistration```
Transforms and registers mask file (```.nii/.nii.gz```) to each subject single trial model for use in multivariate analyses. Files saved as ```nii.gz```.

```runMVPAClassification.m```
Runs MVPA within singular ROIs or a Searchlight within mask provided by user. Currently employs SVM classifier for ROI level, and LDA for searchlight. Output is saved as tidyverse-formatted CSV file and as MAT files. Searchlight results are saved as weights NiFTi files. Summary CSV file is generated as well, concatenating all subject accuracies. 
Requires CoSMoMVPA

```runRSA.m```
Runs RSA within singular ROIs, Searchlight within mask provided by user, or ERS analyses. Output is saved as tidyverse-formatted CSV file and as MAT files. Summary CSV file is generated as well, concatenating all subject accuracies. 
Requires CoSMoMVPA

```specifyModel.m```
Creates single trial models for functional data across all runs. Requires behavioral file to be sourced for onsets and trial tag information

## Directory Organization
The package is built around a specific hierarchal structure inspired by the [BIDS](https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy) organization. Below is a representation of the structure, with the '^' denoting parent directories that are created by the package:

```
projectTemplate
├── multivariate^
│   ├── masks^
│   |   └── params_fmriprep_localizer_MVPA_conditions_face_object_ROI^
│   |       ├── sub-y001
│   |       └── sub-y002
│   ├── models^
│   |   ├── SingleTrialModellocalizer^
│   |   |   ├── sub-y001
│   |   |   └── sub-y002
│   |   └── params_fmriprep_localizer_MVPA_conditions_face_object_ROI^
│   |   |   ├── sub-y001
│   |   |   └── sub-y002
│   ├── params_spm12_encoding_MVPA_conditions_face_object_ROI.mat
│   └── params_fmriprep_encoding_MVPA_conditions_face_object_ROI.mat
├── preprocessing
│   ├── fmriprep
│   │   ├── fmriprep
│   │   │   ├── logs
│   │   │   ├── sub-y001
│   │   │   └── sub-y002
│   │   ├── fmriprep_wf
│   │   ├── freesurfer
│   │   └── reportlets
│   └── spmPreprocessing^
│       ├── psfiles
│       ├── sub-y001
│       └── sub-y002
├── rawdata
│   ├── sub-y0001
│   │   ├── anat
│   │   │   └── sub-01_T1w.nii.gz
│   │   └── func
│   │       ├── sub-01_task-encoding_run-01_bold.nii.gz
│   │       ├── sub-01_task-encoding_run-01_events.tsv
│   │       ├── sub-01_task-encoding_run-02_bold.nii.gz
│   │       ├── sub-01_task-encoding_run-02_events.tsv
│   │       ├── sub-01_task-retrieval-01_bold.nii.gz
│   │       ├── sub-01_task-retrieval-01_events.tsv
│   │       ├── sub-01_task-retrieval-02_bold.nii.gz
│   │       ├── sub-01_task-retrieval-02_events.tsv
...
└── subjects.tsv

```

## Using the Package
Tee package is currently equipped to either preprocess data in SPM12 or take data resulting from fMRIPrep. Below is order of scripts to be run in the pipeline assuming a new data project:


1. ```createParams.m```          - Dialog box querying directories, analysis choices, and preprocessing programs. All subsequent programs require a parameter file created from this script
2. ```preprocessData.m```        - Preprocess functional data via SPM12 or organize motion regressors from fMRIPrep.
3. ```specifyModel.m```          - Specify SPM model with information for trials of interest (i.e. trials to use for classification/representation).
4. ```estimateModel.m```         - Estimate model with every trials as its own individual regressor.
5. ```maskRegistration```        - Register mask files to individual trial model data.
6. ```runMVPAClassification.m``` - Perform multivariate classification using previously derived model.


## Current Work Using the Package
### Publications
Dennis, N. A., & Overman, A. A., Gerver, C. R., McGraw, K., Rowley, M. A., & Salerno, J. M. (2019). *Different types of associative encoding evoke differential processing in both younger and older adults: evidence from univariate and multivariate analyses.* Neuropsychologia, 135. https://doi.org/10.1016/j.neuropsychologia.2019.107240

Gerver, C. R., Overman, A. A., Babu, H. J., Hultman, C. E., & Dennis, N. A. (in press). *Examining the neural basis of congruent and incongruent configural contexts during associative retrieval.* Journal of Cognitive Neuroscience.

### Presentations
Elbich, D., Adams, R.B., Kveraga, K., Dennis, N.A. (May, 2020). *Discriminability of Neural Patterns within the Magnocellular and Parvocellular Visual Pathways.* Poster presented at the Cognitive Neuroscience Society Conference, Boston, Massachusetts.

Chamberlain, J. D., Turney, I. C., Dennis, N. A., (May, 2020). *Neural Discriminability Increases in Older Adults Following Cognitive Training to Reduce False Memories.* Accepted to be presented at the annual Cognitive Aging Conference, Atlanta, GA (Postponed)

Chamberlain, J. D., Turney, I. C., Dennis, N. A., (May, 2020). *Encoding-Retrieval Similarity (ERS) of Perceptually Related Items and Their Relation to False Memories in Aging.* Presented at the Cognitive Neuroscience Society’s Virtual Annual Meeting, Boston, MA
 
Chamberlain, J. D., Hultman, C., Martinez, V., Carpenter, C., Overman, A., Dennis, N., (March, 2019). *Configuration Manipulation Impacts Neural Patterns in Medial Temporal Lobe in Associative Memory Retrieval.* Presented at the Cognitive Neuroscience Society’s Annual Meeting, San Francisco, CA

