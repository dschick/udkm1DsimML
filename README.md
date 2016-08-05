# Introduction

The udkm1Dsim toolbox is a collection of classes and routines 
to simulate the structural dynamics and the according X-ray 
diffraction response in one-dimensional sample structures after 
ultrafast excitation. 
The toolbox provides the capabilities to define arbitrary layered 
structures on the atomic level including a rich database of 
element-specific physical properties. 
The excitation of ultrafast dynamics is represented by an 
N-temperature-model which is commonly applied for ultrafast 
optical excitations. 
Structural dynamics due to thermal stresses are calculated by 
a linear-chain model of masses and springs. 
The resulting X-ray diffraction response is computed by dynamical 
X-ray theory. 
The udkm1Dsim toolbox is highly modular and allows to introduce 
user-defined results at any step in the simulation procedure.

# Citation

Please cite the following article if you use the udkm1Dsim toolbox for your own publications:

D. Schick, A. Bojahr, M. Herzog, R. Shayduk, C. von Korff Schmising & M. Bargheer,
udkm1Dsim - A Simulation Toolkit for 1D Ultrafast Dynamics in Condensed Matter,
[Comput. Phys. Commun. 185, 651 (2014)](http://doi.org/10.1016/j.cpc.2013.10.009) [(preprint)](http://www.udkm.physik.uni-potsdam.de/medien/udkm1Dsim/udkm1DsimManuscriptPrePrint.pdf).

# Installation

Add the udkm1Dsim toolbox folder and all of its 
subfolder to your MATLAB searchpath by excecuting 
the command:

```
   addpath(genpath('path2udkm1Dsim'));
```
   
In order to use the udkm1Dsim documentation
before MATLAB 2013a (8.1) 
- open the MATLAB Preferences from the File menu. 
- click Help, and then select the All Products button.

since MATLAB 2013a (8.1)
- open the MATLAB Help and click on "Supplemental Software" on the lower 
  left of the help start page

Open the MATLAB Product Help and open the udkm1Dsim
toolbox documentation from the content listing.
Follow the Getting Started to get familar with the 
simulation workflow.

In order to run the example files change your current 
MATLAB directory to `path2udkm1Dsim/examples/`.
