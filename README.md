# Radiation Damage Simulation 

## Overview
The repo houses a c++ implementation of the Hamburg model of radiation damage ("simulation code") and a Python script to fit the depletion voltage data ("depletion voltage fit"). The simulation does not describe well the data (for non totally clear reason) and therefore empirical fits are employed to predict future depletion voltage evolution.

The simulation code calculates the effective doping concentration based on an introduction of stable/decaying acceptors and the removal of donors due to irradiation. Moreover, it handles the annealing processes depending in the temperature. According to the effective doping concentration, the depletion voltage is calculated.

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

In order to simulate all of these properties, a profile of irradiation dose and temperature is necessary. The profile is produced by the PixelMonitoring repo, see https://github.com/fleble/PixelMonitoring 

### Annealing and leakage current constants
Annealing constants are saved in the config file `config/annealing_constants.py` to avoid duplications.     
The leakage current constant are hard-coded in the simulation code. For now they are only used in the simulation code.

### Compiling and running simulation code

The simulation code can be compiled and run with:
```bash
./run_radiation_damage_simulation.sh
```
The sector name and time period must be adjusted to your needs.

#### Information about the simulation code
Histograms of the various quantities of interest (e.g. leakage current and depletion voltage!) are stored in a single root file.

Some variables might be changed in the first part of the code, apart from that, no adaptions should be necessary - TODO: Understand these variables and differentiate code and configuration.

Attention:  
1. duration in the profile file HAS to be dividable by global timestep (parameter in the code, default=1) ->Integer result required for the ratio of duration and global time step

2. profile file may not have additional (even/especially empty) lines after the last line or in the beginning

3. the value for the global_layer_conversion has to be changed in order to fit to another layer (parameter transforms luminosity in neq) eg: B-Layer: 2.5e12, IBL 6e12

4. due to some computational reasons the final luminosity might not be correct. The fluence (in neq/cm2) however, should be exact.

5. to change from one detector type to another (e.g. for a different layer in the detector), change mainly 3 things: thickness (200 for IBL and 250 for PIXEL) and global_layer_conversion (global values) and Ndonor (different for PIXEL and IBL)

### Running the depletion voltage fit

The depletion voltage code can be run with:
```bash
./run_depletion_voltage_fit.sh
```
The sector name and time period must be adjusted to your needs.

TODO: The depletion voltage fit code is incomplete and does not perform any fit yet...