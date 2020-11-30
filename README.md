# matlab-portfolio
*Repository of some Matlab scripts I have written. These include project-specific scripts and more general tools that I use regularly.*

## FBP-pipeline.m
Script to process Florbetapir-PET scans from nii frames to SUVR. Adapted from a defunct Python script for Matlab with original features added.

## affine-transform.m
Custom version of an existing script to affine warp PET scans via their respective MRI.

## longitudinal-mri-pet.m
Original script that automates the previously manual task of longitudinal PET processing via SPM12 for clinical interpretation.

## warp-binarize-native-pons.m
Original custom script used for independent research using FDG-PET. Warps a template pons (reference region) ROI to patient-specific native space and binarizes the ROI with help of a white matter probability image to ensure proper fit.

## warp-via-dps-or-mri.m
Warps PET images to MNI-template space using either a patient-specific MRI or previously defined deformation parameters. Builds on previous scripts for added flexibility.
