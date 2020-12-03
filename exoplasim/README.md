EXOPLASIM
======

General Circulation Models Planet Simulator (PlaSim) and PUMA

This repository has all necessary components to run the models.

These models are research models used in Meteorology and Earth Sciences.
You need knowledge in Meteorology and skills in Linux, C and FORTRAN
for running the models and interpreting the results.

Please start reading the manuals located in:

puma/doc for the PUMA model
plasim/doc for the Planet Simulator

Then have a look at the README file for configuration and setup.

=======

This version of PlaSim has been modified to include additional modules
used for simulating climate on geological timescales (carbon-silicate 
weathering on land), different rotation states (i.e. tidally-locked), 
different atmosphere masses, and arbitrary input spectrum. 

The model has also been adapted to permit a rudimentary glacial model,
with an adaptive glacier mask and (mostly) unbounded land ice thickness
(included in surface geopotential height), but no ice flow, denudation,
or adaptive sea level.

The framework for the glacial model may be adapted in the future to enable
orogeny, such as might be expected from a coupled plate-tectonics model.

[![DOI](https://zenodo.org/badge/97154456.svg)](https://zenodo.org/badge/latestdoi/97154456)
