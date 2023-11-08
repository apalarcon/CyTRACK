============================================================================================
||         +++++++                                                                        ||
||       +++++++         +++++++           +++++++  +++++     +     +++++++  +    +       ||
||      ++++             +        +     +     +     +   +    + +    +        +   +        ||
||    ++++++             +         +   +      +     +++++   +++++   +        + +          ||
||    ++++++++           +           +        +     + +    +     +  +        +   +        ||
||   ++++++++++          +++++++     +        +     +  +   +     +  +++++++  +    +       ||
|| +++++ +++++       <-------------------------------------------------------------->     ||
|| ++++    ++++                            Cyclone Tracking                               ||
|| +++++  +++++                                                                           ||
||  ++++++++++                                                                            ||
||   +++++++++                                                                            ||
||    +++++++                     EPhysLab (Environmental Physics Laboratory)             ||
||       +++++                             Universidade de Vigo                           ||
||    +++++++                     contact: albenis.perez.alarcon@uvigo.es                 ||
|| +++++++                                                                                ||
============================================================================================


```
============================================================================================
||         +++++++                                                                        ||
||       +++++++         +++++++           +++++++  +++++     +     +++++++  +    +       ||
||      ++++             +        +     +     +     +   +    + +    +        +   +        ||
||    ++++++             +         +   +      +     +++++   +++++   +        + +          ||
||    ++++++++           +           +        +     + +    +     +  +        +   +        ||
||   ++++++++++          +++++++     +        +     +  +   +     +  +++++++  +    +       ||
|| +++++ +++++       <-------------------------------------------------------------->     ||
|| ++++    ++++                            Cyclone Tracking                               ||
|| +++++  +++++                                                                           ||
||  ++++++++++                                                                            ||
||   +++++++++                                                                            ||
||    +++++++                     EPhysLab (Environmental Physics Laboratory)             ||
||       +++++                             Universidade de Vigo                           ||
||    +++++++                     contact: albenis.perez.alarcon@uvigo.es                 ||
|| +++++++                                                                                ||
============================================================================================
```

# CyTRACK
An open-source, comprehensive and user-friendly Python toolbox for detecting and tracking cyclones

# What do I need to get and run CyTRACK?

<table>
<thead>
<tr>
<th>python</th>
<th>status</th>
</tr>
</thead>
<tbody>
<tr>
<td>Python3</td>
<td> tested</td>
</tr>
</tbody>
</table>

## To run CyTRACK, you need

   * python 3
   * git

and

  *  anaconda (or similar to manage python packages)

or

  *  python 3 and the required modules on a cluster

## The packages required to run CyTRACK are:
  
```
- netCDF4
- numpy 
- scipy 
- mpi4py
- numpy 
- time
- struct
- datetime
- functools
- pathlib 
- gzip
- shutil
- math 
- fnmatch
- sys
- os
- matplotlib
- imp
```
# Installation

### First Method
  
1 - Clone CyTRACK repository

  ```
https://github.com/apalarcon/CyTRACK.git
  ```
2 - Verify you have installed all packages requiered for CyTRACK. If you use an Anaconda environment, please be sure you have activate the environment

3 - Copy the cytrack directory to your Anaconda instalation
```
cp -r cytrack path_to_anaconda_installation/lib/python3.x/site-packages/
````

### Second Method

1 - Clone CyTRACK repository.

2 - Verify you have installed all packages requiered for CyTRACK. If you use an Anaconda environment, please be sure you have activate the environment.

3 - run install_CyTRACK.sh.

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
Please note that you can get CyTRACK input parameters template by running:
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
A text-format file similar to the HURDAT2 dataset suported by the U.S. National Hurricane Centre


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


# Contact and Support

* This code is not bug-free. Please report any bugs through 'Issues': https://github.com/apalarcon/CyTRACK/issues

or

* Contact to Albenis Pérez Alarcón:
  
  apalarcon1991[a]gmail.com
  
  albenis.perez.alarcon[a]uvigo.es


# LICENSE
Copyright 2023 Albenis Pérez-Alarcón, Patricia Coll-Hidalgo, Ricardo M. Trigo, Raquel Nieto and Luis Gimeno

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
