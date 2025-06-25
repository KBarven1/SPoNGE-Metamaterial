# Soft Porous Metamaterials using Inflation-Induced Buckling for Smart Actuation

This repository contains code files necessary to analyze SPoNGe metamaterials, including the alignment of adjacent cylindrical buckling modes within the structure's super cell.

The .inp files are run using Abaqus/CAE, which will produce the buckling modes using a series of static steps. The MATLAB code then acts to curve fit the plate displacements that are output from Abaqus to represent their displacements analytically. Lastly, the Mathematica code is used to predict the lobe alignment of adjacent pores using the curve fitted MATLAB plate displacements.

The .inp files presented here are representative of the process needed to run any geometry of the SPoNGe metamaterial.
