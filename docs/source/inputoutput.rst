 Input and Output data
=======================


Input Data
----------

CyTRACK can read and directly process input data from ERA5 and WRF-ARW model by typying source="ERA5" or source= "WRF" in the input file. Note that ERA5 input data must be downloaded with longitudes in 0-360 format. For other case of ERA5 data or other sources, please set source="CUSTOM".

CyTRACK is configured for tracking tropical cyclones (TC), extratropical cyclones (EC), Mediterranean Cyclones (MC), subtropical cyclones (SC) and tropical-like-cyclones (TLC)

.. seealso::

    CyTRACK help for more details


CyTRACK output
--------------

It is a comma-delimited text format following a similar to the HURDAT2 dataset ( Landsea and Franklin (2013) ) suported by the U.S. National Hurricane Center (NHC)

.. code:: bash

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

.. code:: bash

    latc: degrees (-90 to 90)
    lonc: degrees (-180 to 180)
    Pc: hPa
    Vmax: km/h
    Size: km
    Proci: hPa
    ROCI: km


.. seealso::

    Landsea, C.W., Franklin, J.L. (2013). Atlantic hurricane database (HURDAT2) 2.0. Technical Report, National Oceanic and Atmospheric Administration, National Weather Service, National Hurricane Center, Miami, FL.
    https://doi.org/10.5067/HURDAT2-AL022013
