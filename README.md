# MRI_Phantom
Matlab code for segmenting and calculating T1rho or T2* from Philips MRI images of phantom

The MRI phantom was purchased from The Phantom Lab, Salem, New York, USA.  The phantom contains three pairs of vials with agarose concentrations of 2%, 3% and 4%.  The imaging was done with a 3T Philips machine.  The MRI series description must contain the spin lock times as this is not stored in the DICOM header.  The spin lock times must be between the characters "SL" (spin lock) and "ms" (milliseconds).  See dicom_lst.m and dicom_lst2.m for implementation.

Note that there are three main files that need to be run in sequences as the subsequent files depend on the outputs from the previous file.  Here's the order:

1. dicom_lst2
2. phantomr_plt2
3. phantomr_T1rho2

dicom_lst2 brings up a list of potential DICOM subdirectories (directories starting with s).  Select all of the DICOM subdirectories.  As long as there are not any subdirectory names starting with s (e.g. segmentations/), just click the "select all" button.

dicom_lst2 produces a table listing information about the series in the different subdirectories.  This information is also written to a MS-Excel spreadsheet, dicom_lst2.xlsx, and the variables are saved in a MAT file, dicom_lst2.mat.

phantomr_plt2 does the segmentation of the middle half of the slices.  It requires the function file circ_plt.m and dicom_lst2.mat.  The program processes all of the T1rho series (series with spin lock times in the Series Description).

phantomr_plt2 produces logical masks (vmsk - # pixel by 6 vials by # middle slices) that are saved in MAT files, phantomr?_plt2.mat.  ? is the series number.  Note the number of initial slices skipped is in the variable istrts.  The program also produces a PS file (phantomr?_plt1.ps) and PDF file (phantomr?_plt2.pdf) of some of the plots for each T1rho series.  I manually combine the files into a single PDF file.  (The PS file could be just converted to a second PDF file.)  Note that not all the slice plots are saved, so all of the slice plots should be checked before continuing (the program pauses between series).

phantomr_T1rho2 does the calculation of the T1rho values.  It requires the function file exp_fun1.m, dicom_lst2.mat and phantomr?_plt2.mat files.

phantomr_T1rho2 produces two MS-Excel spreadsheets, phantomr?_T1rho2.xlsx and phantomrs?_T1rho2.xlsx.  ? is the series number.  Vial results are written to phantomr?_T1rho2.xlsx and slice results are written to phantomrs?_T1rho2.xlsx.  The data and results are saved to MAT files, phantomr?_T1rho2.mat.

The phantom has vials with different concentrations of MnCl2 for getting different T2* values.  Similar to the T1rho processing, the files necessary for calculating T2* (for a single set of echo times in different series) are:

1. dicom_lst2
2. phantoms_plt2
3. phantoms_T2star2

phantoms_plt2 lets the user select the series that make up a complete set of T2* echo times.  phantoms_T2star2 processes all of the phantoms*_plt2.mat files in the current directory.

If a DICOMDIR file exists, dicom_lst.m, will read the DICOMDIR file and use the series information.

See the comments in the individual files for more information.
