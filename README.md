# BrCa Project

This repository contains the analysis code for the project.

## GenerateAnndata1.ipynb

This script creates the Anndata objects from QuPath output from cell segmentation. Furthermore, it QCs the samples, runs filtering based on total fluorescence level as well as DAPI levels and outputs h5ad files that can easily be read into an Anndata object.
