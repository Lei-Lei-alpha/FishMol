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
Please refer to examples for usage.
  - I/O of trajectory files in `xyz` format
  - calibration of the trajectory file by fixing the centre of mass
  - the calculation of mean square displacement (MSD) and diffusion coefficient
  - the anisotropy of diffusion by calculating the projection of MSD and diffusion coefficient along specified directions
  - hydrogen bond recognition and lifetime analysis
  - the reorientation lifetime of vectors
  - the distribution function of a range of scalars (including the radial distribution function, angular distribution function, dihedral angle distribution function)
  - 2D combined distribution functions
  - Van Hove correlation function
  - dimer lifetime correlation function
