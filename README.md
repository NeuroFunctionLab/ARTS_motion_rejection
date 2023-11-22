# ARTS Motion Rejection

Automatic Rejection based on Tissue Signal (ARTS) is an algorithm for automatic detection and exclusion of motion-contaminated images to improve the precision of cerebral venous oxygenation quantification. This paper has been published on Magnetic Resonance Imaging: https://www.sciencedirect.com/science/article/pii/S0730725X23001935

## Table of Contents

- [Introduction](#introduction)
- [Getting Started](#getting-started)

## Introduction

Cerebral venous oxygenation (Yv) is a key parameter for the brain’s oxygen utilization and has been suggested to be a valuable biomarker in various brain diseases including hypoxic ischemic encephalopathy in neonates and Alzheimer's disease in older adults. T2-Relaxation-Under-Spin-Tagging (TRUST) MRI is a widely used technique to measure global Yv level and has been validated against gold-standard PET. However, subject motion during TRUST MRI scan can introduce considerable errors in Yv quantification, especially for noncompliant subjects. Therefore, we have developed an Automatic Rejection based on Tissue Signal (ARTS) algorithm to identify motion-contaminated TRUST images and exclude them from Yv quantification.

## Getting Started
### Step 1: open the GUI
Run "trustcode.m", the "trustcode" GUI shows up

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/cd556b18-8640-4d25-b105-a439c0f3a54e)

### Step 2: select TRUST data
Click on “Browse” button to choose the data to process. This GUI now supports three data formats: .REC, .DCM, .IMA. 
For .REC, directly choose the .REC file to process.
For .DCM and .IMA, select the first TRUST image of the data to process. Please make sure all the TRUST images to process are stored in the same folder, and no space is allowed in the path, e.g. please name the folder as “TRUST_sample” instead of “TRUST sample”.

Note: for the Philips Enhanced DICOM format, all TRUST images are stored in a single DICOM file. In this case, just select the DICOM file.

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/a10ecf1e-d17f-4876-b025-8a7dfb631bb3)

### Step 3: set the parameters
Enter the “Row number” (number of rows in the image), “Column number” (number of columns in the image) and “Dynamic number” (number of images), or click on “Default”.

Enter the “Hematocrit” or click on “Default”. If you click on “Default”, it will prompt a dialog window to let the user choose the gender of the subject. The default hematocrit for male, female and unknown are 0.42, 0.40, 0.41, respectively.

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/2cc814ba-c38e-4821-9e01-395d2ce69480)

Enter an array of “Effective TE” (in ms), or click on “Default”. The array of effective TEs must exactly match the actual order of effective TEs of the TRUST images to process. When you click on "Default", the GUI first asks you what is the inner loop of the sequence scan order. For data acquired with Siemens TRUST, select "eTE"; for data acquired with Philips TRUST, select "Rep"

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/009af1d2-e310-4fc1-aec1-8444c702d57f)

Then the GUI asks you to select the minimal eTE of the TRUST sequence. In most cases you should select "0.44"

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/43151c9f-9efc-4767-8877-e84538daf6b4)

Select "Neonate" or "Adult" based on your data. For threshold, we recommend 0.03 (see our paper for details: https://www.sciencedirect.com/science/article/pii/S0730725X23001935).

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/54cf4bec-9d40-4266-a5e6-1aa25f1424ae)

### Step 4: process the data
Click on “Process” button, in a few seconds, the GUI will first show the control images of each dynamic

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/4fe4c957-9101-474d-89e5-3b8ac51cafb7)


Then the GUI will show the difference image of each dynamic

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/de78af03-1827-4f25-a7a6-3bf74c65cd2a)


In the Matlab command window, the difference images that have been identified as motion-corrupted and subsequently excluded are listed. You can see that the ARTS algorithm correctly identifies Img5 as motion-corrupted

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/f0237d31-fbdf-4e7e-a774-9e364f716a15)


Finally, the GUI displays a figure of the averaged difference image (after excluding motion-corrupted images) to let the user draw the region of interest (ROI) for the superior sagittal sinus (SSS). Please make sure the ROI is large enough to cover the entire SSS, as shown below. After drawing the ROI, double click or right click and then select “Create Mask”.

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/b880ea45-d66d-4601-ae85-2633eefea624)

The results will be displayed in the “Results” panel (bottom row in the GUI). The results will also be saved as a .txt file in the same folder as the TRUST data

![image](https://github.com/NeuroFunctionLab/ARTS_motion_rejection/assets/142855046/b187eeb9-b757-485d-9df1-47542b1e7caf)

















