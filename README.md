This repository contains the codebase used for analysis and figure generation for **"SKI complex loss renders 9p21.3-deleted or MSI-H cancers dependent on PELO"** (Borck et al., *Nature*, 2024).  

To reproduce the figures and/or analysis, follow these steps. 

1. Pull this repository, creating the folder ```pelo-manuscript```.
2. Set up a working environment with the following command, using python 3.9:
```
pip install -r pelo-manuscript/requirements.txt
```
4. (Optional) The processed data produced in this study is available on figshare at https://doi.org/0.6084/m9.figshare.27249741.v1[https://doi.org/0.6084/m9.figshare.27249741.v1], you can use "Download all" to get the dataset as a zipped folder.  Unzip it, rename the folder to ```data```, and place it in the ```pelo-manuscript``` directory.
5. Download the "Source Data" files from the manuscript, and set the ```SOURCE_DATA_FIG``` and ```SOURCE_DATA_EXT_FIG``` paths at the start of the notebook to the local paths to these files.

All notebooks in this repository begin with a codeblock which reads in the figshare dataset (if it is not already present), as well as files from other publications/resources, and then saves them locally to the ```data``` folder.  All outputs produced by the script are saved to the ```outputs``` folder and all figures are saved to the ```figures``` folder.

The file structure should look as follows:
```
pelo-manuscript/
| |_ data/ #either downloaded from figshare or fetched automatically
| |   |_ ... raw CSVs of data
| |_ outputs/ #created at runtime 
| |   |_ ... output files from analyses
| |_ figures #created at runtime
|     |_ ... figures as PNG and PDF
|_ Fig 1 + Ext Data Fig 1.ipynb
|_ Fig 2 + Ext Data Fig 2.ipynb
|_ Fig 3 + Ext Data Fig 3.ipynb
|_ Fig 4 + Ext Data Fig 4+5.ipynb
|_ Analysis for 9p21_3 Modifier Screen.ipynb #produces 9p21_3 logfold change/gene effect files
|_ Analysis for RNAseq.ipynb #produces deseq2/gsea_prerank files
|_ helper_funcs.py #functions and global variables used throughout notebooks
```

Note: Immunoblot panels will render empty as raw image files are not present.
