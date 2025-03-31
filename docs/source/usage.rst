Run CyTRACK
=================================

How do I run CyTRACK?
----------------------

To run CyTRACK, create a file with the following code (could be **run_cytrack.py**)

.. code-block:: python

  import numpy as np
   import cytrack 

   args = cytrack.read_args()
   if args.cytrack_help:
      cytrack.help()
   elif args.get_template:
      cytrack.get_cytrack_inputs_template()
   else:
      cytrack.get_cytrack_main(args.parameterfile)

Once this code is created, CyTRACK can be run as follows:

On a Linux computer
----------------------

.. code-block:: bash

   mpirun -n N_proc python run_cytrack.py  -pf input_file

   e.g:  mpirun -np 4 python run_cytrack.py  -pf test_case.cfg

On a HPC with Linux
----------------------

Create a bash script (**run_cytrack.sh**). This example is valid for FINESTARRAE III cluster at the Galician Supercomputing Center.

.. code-block:: bash

   #!/bin/bash -l

   #SBATCH --mem=64GB
   #SBATCH -N 1
   #SBATCH -n 40
   #SBATCH -t 7-00:00:00

   module --purge
   module load cesga/2020
   module load miniconda3/4.9.2
   conda activate envname


   srun -n $SLURM_NTASKS  --mpi=pmi2 python run_cytrack.py -pf input_file


.. code-block:: bash
   sbatch run_cytrack.sh



Input and help
----------------------

You can also get the input file template and help by using the **run_cytrack.py** created above.

For input file template:

.. code-block:: bash

   python run_cytrack.py -gt t
   

For help on the input file:

.. code-block:: bash

   python run_cytrack.py -cth t



.. note::

   CyTRACK could run under Windows if you have installed the Anaconda distribution. This case has not been tested yet.
   If you have any problem, please contact us.