```
============================================================================================
||          +++++++                                                                       ||
||       +++++++         +++++++           +++++++  +++++     +     +++++++  +    +       ||
||      ++++             +        +     +     +     +   +    + +    +        +   +        ||
||    ++++++             +         +   +      +     +++++   +++++   +        + +          ||
||    ++++++++           +           +        +     + +    +     +  +        +   +        ||
||   ++++++++++          +++++++     +        +     +  +   +     +  +++++++  +    +       ||
|| +++++ +++++       <-------------------------------------------------------------->     ||
|| ++++    ++++                      Cyclone TRACKing  framework                          ||
|| +++++  +++++                                                                           ||
||  ++++++++++                                                                            ||
||   +++++++++                                                                            ||
||    +++++++                     EPhysLab (Environmental Physics Laboratory)             ||
||       +++++                             Universidade de Vigo                           ||
||    +++++++                     contact: albenis.perez.alarcon@uvigo.es                 ||
|| +++++++                                                                                ||
============================================================================================
```

# CyTRACK: Cyclone Tracking framework
CyTRACK is an open-source, comprehensive and user-friendly Python toolbox for detecting and tracking cyclones

[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

 [![Current Version: ](https://img.shields.io/badge/Current_Version-1.0.3-blue)](https://anaconda.org/tramo-ephyslab/cytrack)

If you use CyTRACK, please cite it as follows:


Pérez-Alarcón, A.; Coll-Hidalgo, P.; Trigo, R.M.; Nieto, R.; Gimeno, L. (2024). CyTRACK: An open-source and user-friendly python toolbox for detecting and tracking cyclones. Environmental Modelling & Software, 176, 106027. https://doi.org/10.1016/j.envsoft.2024.106027.


## Version History
[![Version: ](https://img.shields.io/badge/Version-1.0.3-blue)](https://anaconda.org/tramo-ephyslab/cytrack) updates CyTRACK for using the new release of the CDS API to download data from the ERA5 reanalysis.

# What do I need to get and run CyTRACK?

## To run CyTRACK, you need

   * [![python](https://img.shields.io/badge/Python-3-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
   * [![Git: ](https://img.shields.io/badge/Git--blue)](https://git-scm.com/)

and

 * [![Anaconda 3](https://img.shields.io/badge/Anaconda-3-green.svg)](https://www.anaconda.com/) (or similar to manage python packages)

or

  *  [![python](https://img.shields.io/badge/Python-3-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org) and the required modules on a cluster

## The packages required to run CyTRACK are:
  
```
- netCDF4
- numpy 
- scipy 
- mpi4py
- time
- datetime
- functools
- math 
- sys
- os
- matplotlib
- imp
- xarray
- sklearn
- argparse
```
# Installation

### First Method

1 - Using conda
```
conda install -c tramo-ephyslab cytrack
```

### Second Method
1 - Clone CyTRACK repository

  ```
git clone https://github.com/apalarcon/CyTRACK.git
  ```
2 - Verify you have installed all packages requiered for CyTRACK. If you use an Anaconda environment, please be sure you have activated the environment

3 - Copy the cytrack directory to your Anaconda instalation
```
cp -r cytrack path_to_anaconda_installation/lib/python3.x/site-packages/
````
To knnow the exactly path of your Anaconda installation (```patho_to_anaconda_installation/.../site-packages/```), you can follow these instructions:

```
$python
>>  from distutils.sysconfig import get_python_lib
>> print(get_python_lib())
```

### Third Method

1 - Clone CyTRACK repository.

 ```
git clone https://github.com/apalarcon/CyTRACK.git
  ```

2 - Verify you have installed all packages requiered for CyTRACK. If you use an Anaconda environment, please be sure you have activated the environment.

3 - run install_CyTRACK.sh.

### NOTE
If you have a problem with the mpi4py library, try these steps:

* Remove the ```mpi4py``` library ```conda remove mpi4py```
* Install the ```openmpi``` library ```conda install conda-forge::openmpi ```
* Install again the ```mpi4py``` library ```conda install mpi4py ```

If the problem continue (the problem is frequently related with the ```libmpi.so.12 ``` or similar), you can also try
* Search the mising library on your system and link it to your Anaconda lib path.
  ```
  ln -s path_to_missing_library/libmpi.so.12 patho_to_anaconda_installation/lib/
  ```
or 

* Contact your system administrator


# CyTRACK namelist file configuration

* For details and help on CyTRACK input file, type
```
python run_CyTRACK.py -cth t
```
or configure your code as follows:
```
import cytrack
cytrack.help()
```
You can get CyTRACK input file template by running:
```
python run_CyTRACK.py -gt t
```
or
```
import cytrack
cytrack.get_cytrack_inputs_template()
```


# Input data
CyTRACK can read and directly process input data from ERA5 and WRF-ARW model by typying source="ERA5" or source= "WRF" in the input file. Note that ERA5 input data must be downloaded with longitudes in 0-360 format. For other case of ERA5 data or other sources, please set source="CUSTOM". 

CyTRACK is configured for tracking tropical cyclones (TC), extratropical cyclones (EC), Mediterranean Cyclones (MC), subtropical cyclones (SC) and tropical-like-cyclones (TLC)

<b>* Revise CyTRACK help for more details</b>

# CyTRACK outputs
It is a comma-delimited text format following a similar to the HURDAT2 dataset (<a href="https://doi.org/10.1175/MWR-D-12-00254.1" target="blank"> Landsea and Franklin (2013) </a>) suported by the U.S. National Hurricane Center (NHC)
```
CyAL0132017, 63,
 Date,   hour,   latc,   lonc,    Pc,       Vmax,    Size,     Proci,      ROCI,     Core,  VTL   VTU    B
20170916, 06,   11.25,  -48.50,  1010.11,  52.68,   393.692,   1010.81,   339.477,   SLCC,  30,  -52,  3.0,
20170916, 12,   11.50,  -50.50,  1009.89,  49.47,   317.813,   1013.35,   434.489,   SLCC,  47,  -13,  3.0,
20170916, 18,   11.50,  -51.75,  1006.07,  56.79,   382.220,   1008.11,   192.095,   UDCC,  51,    0,  1.0,
20170917, 00,   12.25,  -53.25,  1007.04,  55.48,   351.333,   1010.95,   282.036,   SDWC,  54,   14,  0.0,
20170917, 06,   12.75,  -54.75,  1004.82,  61.51,   329.537,   1010.82,   515.375,   SDWC,  50,    4,  0.0,
20170917, 12,   13.00,  -56.25,  1004.13,  60.83,   410.666,   1012.79,   513.898,   SDWC,  42,   59,  0.0,
20170917, 18,   13.25,  -57.00,  1002.24,  43.33,   386.218,   1010.15,   470.297,   SDWC,  47,   54, -5.0,
20170918, 00,   14.25,  -58.00,  1001.37,  48.35,   409.628,   1011.21,   458.576,   SDWC,  39,   67, -2.0,
```


# Running CyTRACK

To run CyTRACK, we recomend the using of run_CyTRACK.py script

* On Linux PC
  
1 - By using run_CyTRACK.py
```
python run_CyTRACK.py -pf input_file
```
* Using MPI
```
mpirun -n N_proc python run_CyTRACK.py  -pf input_file
```


* On a HPC with Linux:

1 - Create a bash script (run_CyTRACK.sh). This example is valid for FINESTARRAE III cluster on Galician Supercomputing Center.
  
  ```
#!/bin/bash -l

#SBATCH --mem=64GB
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 7-00:00:00

module --purge
module load cesga/2020
module load miniconda3/4.9.2
conda activate envname


srun -n $SLURM_NTASKS  --mpi=pmi2 python run_CyTRACK.py -pf input_file
  ``` 
2 - Submit run script

```
sbatch run_CyTRACK.sh
```

# Testing CyTRACK
Please, follow the intructions at

[![Git: Testing_CyTRACK ](https://img.shields.io/badge/GitHub-Testing_CyTRACK-blue)](https://github.com/apalarcon/CyTRACK/tree/main/testing_CyTRACK)


# Contact and Support

* This code is not bug-free. Please report any bugs through 'Issues': https://github.com/apalarcon/CyTRACK/issues

or

* Contact to Albenis Pérez Alarcón:
  
  apalarcon1991[a]gmail.com
  
  albenis.perez.alarcon[a]uvigo.es


# LICENSE

This software is published under the GPLv3 license. This means: 
1. Anyone can copy, modify and distribute this software. 
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. If you dare build your business solely from this code, you risk open-sourcing the whole code base.
6. If you modify it, you have to indicate changes made to the code.
7. Any modifications of this code base MUST be distributed with the same license, GPLv3.
8. This software is provided without warranty.
9. The software author or license can not be held liable for any damages inflicted by the software.

# References
* Landsea, C. W., & Franklin, J. L. (2013). Atlantic hurricane database uncertainty and presentation of a new database format. Monthly Weather Review, 141(10), 3576-3592. https://doi.org/10.1175/MWR-D-12-00254.1
