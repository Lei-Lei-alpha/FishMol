# FishMol
A python molecular dynamics data analysis tool.

## Install
### Requirements
- python3: >=3.8

### Install
Fishmol is still under development, but you can install it locally:

1. Download the code from GitHub:
```bash
git clone https://github.com/Lei-Lei-alpha/fishmol.git
```
2. Change to the fishmol directory:
```bash
cd fishmol
```

3. Install fishmol
```bash
pip install -e ./
```
4. Run the following command in your terminal to see if fishmol is successfully installed:
```bash
fishmol
```

                      Welcome!                    ▄▄█▀                  FishMol
                                              ▄▄███▀                 version 0.0.1
          ○                                ▄▄█████▀                          ○
               ○                        ▄▄████████▄                         /
                                    ▄▄▄████████████▄                    ○--○
             ○                ▄▄▄████████████████████▄▄                     \ 
                        ▄▄▄██████▀███████████████████████▄▄                  ○--○           ▄
                    ▄▄████████████ ██████████████████████████▄▄             /           ▄▄█▀
             ○   ▄█████████████████ ████████████████████████████▄▄         ○         ▄███▀
              ▄██████████▀▀▀████████ ██████████████████████████████▄▄             ▄█████▀
             ▄████████▀   ○  ▀███████ █████████████████████████████████▄▄      ▄▄█████▀
             ■▄███████▄      ▄███████ █████████████████████████████████████▄▄▄███████▀
              ▀█████████▄▄▄█████████ ███████████████████████████████████████████████
                ▀██████████████████ ████████████████████████████████████████████████▄
                  ▀██████████████▀▄███████████████████████████████████▀▀▀▀▀████████████▄
                    ▀▀█████████▀▄███████████████████████████████▀▀           ▀▀████████▀
                       ▀▀█████▄████████████████████████▀▀▀▀██▀                 ▀▀████▀
                            ▀▀▀█████████████████▀▀▀▀        ▀■                    ▀█▀
                                 ▀▀▀███████▀▀
                                     ▀████
                                        ▀▀▄                Contact: Lei.Lei@durham.ac.uk


## Usage
Please refer to docs for usage.
  - [I/O of trajectory files in `xyz` format](https://lei-lei-alpha.github.io/fishmol/trajectory_IO.html)
    - calibration of the trajectory file by fixing the centre of mass
    - Simple view function
    - Simple manipulate of trajectory: filter, select etc
  - [Calculation of mean square displacement (MSD) and diffusion coefficient](https://lei-lei-alpha.github.io/fishmol/MSD_diff_coeff.html)
  - [the anisotropy of diffusion by calculating the projection of MSD and diffusion coefficient along specified directions](https://lei-lei-alpha.github.io/fishmol/diff_aniso.html)
  - [hydrogen bond recognition and lifetime analysis](https://lei-lei-alpha.github.io/fishmol/H_bond.html)
  - [the reorientation lifetime of vectors](https://lei-lei-alpha.github.io/fishmol/VRD.html)
  - the distribution function of a range of scalars (including the radial distribution function, angular distribution function, dihedral angle distribution function)
  - 2D combined distribution functions
  - Van Hove correlation function
  - dimer lifetime correlation function
