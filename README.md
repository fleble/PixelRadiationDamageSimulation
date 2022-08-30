# Radiation Damage Simulation 

## Overview
The code in data_processing.cpp is a `c++` implementation of the Hamburg model of radiation damage. It calculates the effective doping concentration based on an introduction of stable/decaying acceptors and the removal of donors due to irradiation. Moreover, it handles the annealing processes depending in the temperature. According to the effective doping concentration, the depletion voltage is calculated.
Besides the depletion voltage, the leakage current (inclusive alpha parameter) is calculated based on the Hamburg model. Again, leakage current increase due to irradiation and decrease due to annealing is handled. Two different options for the temperature averaging are available (`use_CMS_paper_alpha0`).


## Installation
To simplify the environment setting, we advise to install the code on `lxplus` or anywhere else where you have access to `cvmfs` (see Section `Setup environment`)
```
git clone git@github.com:fleble/PixelRadiationDamageSimulation.git
```

## Setup environment
After every new login, do:
```bash
source setenv.sh
```
This will only work if you have access to `cvmfs`.


## Usage

### Profile irradiation dose and temperature

In order to simulate all of these properties, a profile of irradiation dose and temperature is necessary. The profile has to have the following format:

`duration (seconds, integer) / temperature (K, float) / dose rate (neq/cm2s, integer)`

Different examples are available in the `temp_rad_data folder`. Besides that, a program to create new profiles from two given temperature and irradiation profiles is available in this folder (help included in the cpp file).

### Compiling simulation code

The program can be compiled with:
```bash
./compile.sh
```

### Executing simulation code

After compilation the program can be called with:
```bash
./data_processing temp_rad_data/my_profile.txt
```

By default, a variety of pdf files is generated and all corresponding root histograms are stored in a single root file. Furthermore, by default, the pdf files will be overwritten and the histograms will be added to the existing root file - this can be changed in the first lines of the code, a unique naming scheme based on the system clock is available!

Some variables might be changed in the first part of the code, apart from that, no adaptions should be necessary.

Attention:  
1. duration in the profile file HAS to be dividable by global timestep (parameter in the code, default=1) ->Integer result required for the ratio of duration and global time step

2. profile file may not have additional (even/especially empty) lines after the last line or in the beginning

3. the value for the global_layer_conversion has to be changed in order to fit to another layer (parameter transforms luminosity in neq) eg: B-Layer: 2.5e12, IBL 6e12

4. due to some computational reasons the final luminosity might not be correct. The fluence (in neq/cm2) however, should be exact.

5. to change from one detector type to another (e.g. for a different layer in the detector), change mainly 3 things: thickness (200 for IBL and 250 for PIXEL) and global_layer_conversion (global values) and Ndonor (different for PIXEL and IBL)
