 NCFLOWSPEED

This package includes custom MATLAB scripts for blood flow analysis for in vivo negative contrast imaging as described in:

Wu, Jung et al. Intravital fluorescence microscopy with negative contrast. 

The scripts were last revised and tested in MATLAB R2018a, using a 3.1GHz Intel Core i7 MacBookPro with 16Gb memory. This package also includes a sample image stack for running a demo and a copy of the DEMO output.

## Instructions for running Demo

1. Download package, place in MATLAB path. Add path.
2. Enter path .../MATLAB/NCFLOWSPEED
3. Run NCFLOWSPEED.m

## Description

NCFLOWSPEED takes a negative contrast, time-series image stack of a vascular network, and maps the blood flow speed. The script performs digital line scanning on multiple vessel segments throughout the vascular network, and the local flow speed along each segment is estimated and mapped. 

In addition to the final blood flow speed map, the script also saves the digital line scan image (OUTPUT/xDistyTime_LinFitSlopeQual) corresponding to each vessel segment (OUTPUT/HL_Mask), and the estimated slope of the line scan that is estimated by Radon Transform and proportional to the local flow speed.

![*NCFLOWSPEED Input*](https://github.com/juwellwwu/NCFLOWSPEED/blob/main/README%20Img/NCFLOWSPEED%20Input.gif?raw=true)

![*NCFLOWSPEED Output*](https://github.com/juwellwwu/NCFLOWSPEED/blob/main/README%20Img/NCFLOWSPEED%20Output.png?raw=true)