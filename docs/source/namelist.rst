Input configuration file 
========================

Common parameters
------------------
.. code-block:: bash

    #============================================================================================================
    #||                                                                                                        ||
    #||              +++++++                                                                                   ||
    #||            +++++++         +++++++           +++++++  +++++     +     +++++++  +    +                  ||
    #||           ++++             +        +     +     +     +   +    + +    +        +   +                   ||
    #||         ++++++             +         +   +      +     +++++   +++++   +        + +                     ||
    #||         ++++++++           +           +        +     + +    +     +  +        +   +                   ||
    #||        ++++++++++          +++++++     +        +     +  +   +     +  +++++++  +    +                  ||
    #||       +++++ +++++       <-------------------------------------------------------------->               ||
    #||       ++++    ++++                           Cyclone TRACKing framework                                ||
    #||       +++++  +++++                                                                                     ||
    #||        ++++++++++                                                                                      ||
    #||         +++++++++        CyTRACK is free under the terms of the GNU General Public license             ||
    #||          +++++++                     EPhysLab (Environmental Physics Laboratory)                       ||
    #||            +++++                             Universidade de Vigo                                      ||
    #||         +++++++                      contact: albenis.perez.alarcon@uvigo.es                           ||
    #||      +++++++                                                                                           ||
    #============================================================================================================

    #CYTRACK INPUT PARAMETERS.
    #For details use python python run_CyTrack -cth t
    #For run CyTrack use python run_CyTrack.py -pf cytrack_inputs 
    #You can use your own input file following the instructions below 
    #------------------------------------------------------------------------------------------------------------
    #Print info during CyTrack runs ["True" / "False"]. Default value ['True']
    verbose="True"

    #Cyclone Type ["TC"/"EC"/"MC/TLC"/"SC/"].
    cyclone_type="EC"

    #============================================================================================================
    #CyTRACK Source information
    #============================================================================================================

    #Source of data ['WRF' / 'ERA5']
    source="ERA5"

    #source of data for tracking cyclones
    path_data_source="path to data source"

    #Only for WRF: wrf_domain ['wrfout_d01','wrfout_d02']
    wrfprefix="wrfout_d01"

    #Only for ERA5. The name of era5 files must be like this era_file_prefix_yyyymmdd_hh.nc or era_file_prefix_yyyymmddhh.nc. 
    #CyTRACK will automatically download ERA5 upper files if they are not found.
    #Prefix in the name of era file.
    era_file_prefix="prefix_"

    #Format of the date in ERA5 file ['yyyymmdd_hh' / 'yyyymmddhh']
    era_date_file_name='yyyymmdd_hh'

    #approximate data resolution in km
    model_res=30

    #Search regions ['NA','SA','NP','SP','SI','NI','EP','WP','NH','SH','GL']
    search_region="NA"

    #Search limits in the region [lonmin,latmin,lonmax,latmax] 
    search_limits=[lonmin,latmin,lonmax,latmax]


    #-------------------------------------------------------------------------------------------------------------
    #Including upper parameters to evaluate the cyclone thermal structure
    #Only if source = ERA5
    #-------------------------------------------------------------------------------------------------------------
    #Checking for upper level parameters ['yes' / 'no'].
    checking_upper_levels_parameters="yes"

    #Get VTL and VTU from linear regression. ["yes" / 'no']. Only if checking_upper_levels_parameters='yes'
    vtl_vtu_lr='yes'

    #Distance form storm center to compute cyclone phase space parameters. Only if Checking for upper level parameters = 'yes' 
    max_dist=500

    #path to upper level files. If source=WRF, set path_data_source_upper like to path_data_source
    path_data_source_upper="path to ERA5 upper level files"

    #Prefix for upper levels files. Only if source = ERA5. 
    #The name of era5 upper files must be like this era_upperfile_prefix_yyyymmdd_hh.nc or era_upperfile_prefix_yyyymmddhh.nc 
    #The date format will be the same as era_date_file_name
    #CyTRACK will automatically download ERA5 upper files if they are not found
    era_upperfile_prefix="upper_ERA5"


    #############################################################################
    #if source=="CUSTOM", please, define the following parameters:
    #############################################################################

    #Custom prefix in the name of custom source file.
    custom_file_prefix="uvmslp_ERA5"

    #Custom format of the date in ERA5 file ['yyyymmdd_hh' / 'yyyymmddhh']
    custom_date_file_name='yyyymmdd_hh'

    #Custom prefix for upper levels files. Only if source = CUSTOM 
    custom_upperfile_prefix="upper_ERA5"

    #Custom MSLP variable name
    custom_mslp_variable="msl"

    #Custom u-wind variable name
    custom_uwind_variable="u10"

    #Custom v-wind variable name
    custom_vwind_variable="v10"

    #Custom latitude_var_name
    custom_latitude_var="latitude"

    #Custom longitude_var_name
    custom_longitude_var="longitude"

    #Custom geopotential high variable name
    custom_geopotential_var_name="z"

    #Custom upper levels variable name
    custom_upper_level_variable_name="level"

    #Custom terrain_high_filename
    custom_terrain_high_filename="path to terrain high netcdf file"

    #custom terrain high variable name
    custom_terrain_high_var_name="z"



    #============================================================================================================
    #CyTRACK date configuration
    #============================================================================================================

    #Start date parameters  [yyyy mm dd hh]
    begin_year="yyyy" 
    begin_month="mm"
    begin_day="dd"
    begin_hour="hh"

    #End date parameters [yyyy mm dd hh]
    end_year="yyyy"
    end_month="mm"
    end_day="dd"
    end_hour="hh"


    #input file time_step, integer desde 1 hasta 6
    dt_h=6


    #Type of calendar in cyclone tracking. 365d to remove February 29. 366d to include February 29
    calendar="365d"


    #============================================================================================================
    #CyTRACK Output file information
    #============================================================================================================

    #path to save CyTRACK outputs
    path_out="path to save CyTRACk outputs"

    #path to save temporal files nedeed for CyTRACK runs
    tmp_dir="path to save CyTRACK temporal files"

    #Remove tmp_dir  ['yes' / 'no']. Default remove_tmp_dir='yes'
    remove_tmp_dir="no"


Tracking tracking extratropical cyclones (ECs)
-------------------------------------------------------------------

.. code-block:: bash

    #============================================================================================================
    #CyTRACK DEFAUL VALUES EXTRATROPICAL CYCLONES (ECs)
    #============================================================================================================
    #Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space. Only necessary if checking_upper_levels_parameters="yes".
    #Default value=0. Set core_criteria_length=-99 to match the full trajectory.
    core_criteria_length=0

    #Lower thermal wind threshold (VTL). Only necessary if checking_upper_levels_parameters="yes".
    #VTL<VTL_threshold. Default VTL_threshold=0
    VTL_threshold=0

    #Upper thermal wind threshold (VTU). Only necessary if checking_upper_levels_parameters="yes".
    #VTU<VTU_threshold. Default VTU_threshold=0
    VTU_threshold=0

    #B parameter. Only necessary if checking_upper_levels_parameters="yes".
    #|B|<Bhart_threshold. Default Bhart_threshold=10
    Bhart_threshold=10

    #Minimum wind speed in m/s threshold to consider a low pressure grid point as EC centre. Default max_wind_speed_threshold=0
    max_wind_speed_threshold=0

    #Outer ninimum wind speed in m/s threshold to consider compute the EC outer radius. Default outer_wind_speed_threshold=0
    outer_wind_speed_threshold=0

    # Minimum distance between two critical centers in km. Default filter_center_threshold=1000
    filter_center_threshold=1000


    #Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=1000
    dist_threshold=1000

    # Critical outer radius in km to considerer a low pressure point as critical center. Default critical_outer_radius=100
    critical_outer_radius=100

    #external search radius in km. Default rout=2000 km
    rout=2000

    #resolution for radial legs in km. Default dr_res=100 km
    dr_res=100

    #resolution of angle steps for radial legs in degrees. Default d_ang=10
    d_ang=10

    #Terrain filter in m. Set terrain_filter=0 to not apply terrain filter. Default terrain_filter=1000
    terrain_filter=1000

    #EC maximum intensity threshold in m/s along the full trajectory. intensity_threshold=0
    intensity_threshold=0

    #Threshold for EC lifetime in hours. Default dt_lifetime=48
    dt_lifetime=48

    #Relative vorticity threshold in 1/s to filter critical TCs centres. vorticity_threshold=1.45e-5
    vorticity_threshold=1.45e-5

    #Maximum slp treshold in hPa to filter EC centres. Deafult min_slp_threshold=1015 
    min_slp_threshold=1015

    #Minimum distance traveled in km by the system. Default minimum_distance_travelled=1000
    minimum_distance_travelled=1000

    #Radial distance (in degrees) for  checking the MSLP increase, default great_circle_distance=6.5
    great_circle_distance=6.5

    #Change in MSLP (in Pa) over a distance of great-circle-distance from the candidate point, default dmslp_great_circle_distance=200
    dmslp_great_circle_distance=200

    #Radius (in km) for computing the maximum surface winds, default radius_for_msw=250
    radius_for_msw=250

    #Dates before the specific date and hour to compute the average mslp. Default prev_days=14
    prev_days=14

    #Mean sea level pressure anomaly threshold in hPa to consideded a grid point as candidate for system centre. Default mslp_anomaly_threshold=-2.5
    mslp_anomaly_threshold=-2.5


Tracking tropical cyclones (TCs)
-------------------------------------------------------------------

.. code-block:: bash

    #============================================================================================================
    #CyTRACK DEFAUL VALUES FOR TROPICAL CYCLONES (TCs)
    #============================================================================================================
    #Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space. Only necessary if checking_upper_levels_parameters="yes".
    #Default value=3. Set core_criteria_length=-99 to match the full trajectory.
    core_criteria_length=3

    #Lower thermal wind threshold (VTL). Only necessary if checking_upper_levels_parameters="yes".
    #VTL>VTL_threshold. Default VTL_threshold=0
    VTL_threshold=0

    #Upper thermal wind threshold (VTU). Only necessary if checking_upper_levels_parameters="yes".
    #VTU<VTU_threshold. Default VTU_threshold=0
    VTU_threshold=0

    #B parameter. Only necessary if checking_upper_levels_parameters="yes".
    #|B|<Bhart_threshold. Default Bhart_threshold=10
    Bhart_threshold=10

    #Minimum wind speed in m/s threshold to consider a low pressure grid point as TC centre
    max_wind_speed_threshold=8

    #Outer ninimum wind speed in m/s threshold to consider compute the TC outer radius
    outer_wind_speed_threshold=6

    # Minimum distance between two critical centers in km. Default filter_center_threshold=400
    filter_center_threshold=400

    #Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=650
    dist_threshold=650

    # Critical outer radius in km to considerer a low pressure point as critical center. Default critical_outer_radius=100
    critical_outer_radius=100

    #resolution for radial legs in km. Default dr_res=100
    dr_res=100

    #resolution of angle steps for radial legs in degrees. Default d_ang=10
    d_ang=5

    #external search radius in km. Default rout=1000 km
    rout=1000

    #Terrain filter in m. Set terrain_filter=0 to not apply terrain filter. Default terrain_filter=0
    terrain_filter=0

    #TC maximum intensity threshold in m/s along the full trajectory. intensity_threshold=10
    intensity_threshold=10

    #Threshold for EC lifetime in hours. Default dt_lifetime=36
    dt_lifetime=36

    #Relative vorticity threshold in 1/s to filter critical TCs centres. vorticity_threshold=1.45e-5
    vorticity_threshold=1.45e-5

    #Maximum slp treshold in hPa to filter TC centres. Deafult min_slp_threshold=1015
    min_slp_threshold=1015

    #Minimum distance traveled in km by the system. Default minimum_distance_travelled=0
    minimum_distance_travelled=0

    #Radial distance (in degrees) for  cheking the MSLP increase, default great_circle_distance=5.5
    great_circle_distance=5.5

    #Change in MSLP (in Pa) over a distance of great-circle-distance from the candidate point, default dmslp_great_circle_distance=200
    dmslp_great_circle_distance=200

    #Radius (in km) for computing the maximum surface winds, default radius_for_msw=100
    radius_for_msw=100

    #Dates before the specific date and hour to compute the average mslp. Default prev_days=14
    prev_days=14

    #Mean sea level pressure anomaly threshold in hPa to consideded a grid point as candidate for system centre. Default mslp_anomaly_threshold=-2
    mslp_anomaly_threshold=-2


Tracking Mediterranean cyclones (MCs)
-------------------------------------------------------------------

.. code-block:: bash

    #============================================================================================================
    #CyTRACK DEFAUL VALUES FOR  MEDITERRANEAN CYCLONES (MCs).
    #============================================================================================================
    #Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space. Only necessary if checking_upper_levels_parameters="yes".
    #Default value=0. Set core_criteria_length=-99 to match the full trajectory.
    core_criteria_length=0

    #Lower thermal wind threshold (VTL). Only necessary if checking_upper_levels_parameters="yes".
    #VTL>VTL_threshold. Default VTL_threshold=0
    VTL_threshold=0

    #Upper thermal wind threshold (VTU). Only necessary if checking_upper_levels_parameters="yes".
    #VTU<VTU_threshold. Default VTU_threshold=0
    VTU_threshold=0

    #B parameter. Only necessary if checking_upper_levels_parameters="yes".
    #|B|<Bhart_threshold. Default Bhart_threshold=10
    Bhart_threshold=10


    #Minimum wind speed in m/s threshold to consider a low pressure grid point as MC centre
    max_wind_speed_threshold=0

    #Outer ninimum wind speed in m/s threshold to consider compute the MC outer radius
    outer_wind_speed_threshold=0

    # Minimum distance between two critical centers in km. Default filter_center_threshold=300
    filter_center_threshold=300

    #Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=400
    dist_threshold=400

    # Critical outer radius in km to considerer a low pressure point as critical center. Default critical_outer_radius=50
    critical_outer_radius=50

    #resolution for radial legs in km. Default dr_res=100
    dr_res=100

    #resolution of angle steps for radial legs in degrees. Default d_ang=10
    d_ang=10

    #external search radius in km. Default rout=1000 km
    rout=800

    #Terrain filter in m. Set terrain_filter=0 to not apply terrain filter. Default terrain_filter=1000
    terrain_filter=0

    #MC maximum intensity threshold in m/s along the full trajectory. Set intensity_threshold = 0 to not apply this criterion.  Default intensity_threshold=0
    intensity_threshold=0

    #Threshold for MC lifetime in hours. Default dt_lifetime=24
    dt_lifetime=24

    #Relative vorticity threshold in 1/s to filter critical MCs centres. vorticity_threshold=1.45e-5
    vorticity_threshold=1.45e-5

    #Maximum slp treshold in hPa to filter MC centres. Deafult min_slp_threshold=1015
    min_slp_threshold=1015

    #Minimum distance traveled in km by the system. Default minimum_distance_travelled=0
    minimum_distance_travelled=0

    #Radial distance (in degrees) for  cheking the MSLP increase, default great_circle_distance=5.5
    great_circle_distance=3

    #Change in MSLP (in Pa) over a distance of great-circle-distance from the candidate point, default dmslp_great_circle_distance=200
    dmslp_great_circle_distance=200

    #Radius (in km) for computing the maximum surface winds, default radius_for_msw=100
    radius_for_msw=100

    #Dates before the specific date and hour to compute the average mslp. Default prev_days=14
    prev_days=14

    #Mean sea level pressure anomaly threshold in hPa to consideded a grid point as candidate for system centre. Default mslp_anomaly_threshold=-2.5
    mslp_anomaly_threshold=-2.5



Trackingtropical-like cyclones (TLCs)
-------------------------------------------------------------------

.. code-block:: bash

    #============================================================================================================
    #CyTRACK DEFAUL VALUES FOR  MEDITERRANEAN TROPICAL-LIKE CYCLONES (TLCs).
    #============================================================================================================
    #Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space. Only necessary if checking_upper_levels_parameters="yes".
    #Default value=1. Set core_criteria_length=-99 to match the full trajectory.
    core_criteria_length=1

    #Lower thermal wind threshold (VTL). Only necessary if checking_upper_levels_parameters="yes".
    #VTL>VTL_threshold. Default VTL_threshold=0
    VTL_threshold=0

    #Upper thermal wind threshold (VTU). Only necessary if checking_upper_levels_parameters="yes".
    #VTU<VTU_threshold. Default VTU_threshold=0
    VTU_threshold=0

    #B parameter. Only necessary if checking_upper_levels_parameters="yes".
    #|B|<Bhart_threshold. Default Bhart_threshold=10
    Bhart_threshold=10

    #Minimum wind speed in m/s threshold to consider a low pressure grid point as TLC centre
    max_wind_speed_threshold=8

    #Outer ninimum wind speed in m/s threshold to consider compute the TLC outer radius
    outer_wind_speed_threshold=2

    # Minimum distance between two critical centers in km. Default filter_center_threshold=400
    filter_center_threshold=300

    #Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=400
    dist_threshold=400

    # Critical outer radius in km to considerer a low pressure point as critical center. Default critical_outer_radius=100
    critical_outer_radius=50

    #resolution for radial legs in km. Default dr_res=100
    dr_res=100

    #resolution of angle steps for radial legs in degrees. Default d_ang=10
    d_ang=10

    #external search radius in km. Default rout=1000 km
    rout=1000

    #Terrain filter in m. Set terrain_filter=0 to not apply terrain filter. Default terrain_filter=0
    terrain_filter=0

    #MC maximum intensity threshold in m/s along the full trajectory. intensity_threshold=10
    intensity_threshold=10

    #Minimum distance traveled in km by the system. Default minimum_distance_travelled=0 km
    minimum_distance_travelled=0

    #Threshold for TLC lifetime in hours. Default dt_lifetime=24
    dt_lifetime=24

    #TRelative vorticity threshold in 1/s to filter critical MCs centres. vorticity_threshold=1.45e-5
    vorticity_threshold=1.45e-5

    #Maximum slp treshold in hPa to filter TLC centres. Deafult min_slp_threshold=1010 
    min_slp_threshold=1013

    #Radial distance (in degrees) for  cheking the MSLP increase, default great_circle_distance=3
    great_circle_distance=3

    #Change in MSLP (in Pa) over a distance of great-circle-distance from the candidate point, default dmslp_great_circle_distance=200
    dmslp_great_circle_distance=200

    #Radius (in km) for computing the maximum surface winds, default radius_for_msw=100
    radius_for_msw=100

    #Dates before the specific date and hour to compute the average mslp. Default prev_days=14
    prev_days=14

    #Mean sea level pressure anomaly threshold in hPa to consideded a grid point as candidate for system centre. Default mslp_anomaly_threshold=-2.5
    mslp_anomaly_threshold=-2.5


Tracking subtropical cyclones (SCs)
-------------------------------------------------------------------

.. code-block:: bash

    #============================================================================================================
    #CyTRACK DEFAUL VALUES FOR SUBTROPICAL CYCLONES (SCs)
    #============================================================================================================
    #Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space. Only necessary if checking_upper_levels_parameters="yes".
    #Default value=7. Set core_criteria_length=-99 to match the full trajectory.
    core_criteria_length=7


    #Lower thermal wind threshold (VTL). Only necessary if checking_upper_levels_parameters="yes".
    #VTL>VTL_threshold. Default VTL_threshold=-50
    VTL_threshold=-50

    #Upper thermal wind threshold (VTU). Only necessary if checking_upper_levels_parameters="yes".
    #VTU<VTU_threshold. Default VTU_threshold=-10
    VTU_threshold=-10

    #B parameter. Only necessary if checking_upper_levels_parameters="yes".
    #|B|<Bhart_threshold. Default Bhart_threshold=25
    Bhart_threshold=25


    #Minimum wind speed in m/s threshold to consider a low pressure grid point as TC centre
    max_wind_speed_threshold=0

    #Outer ninimum wind speed in m/s threshold to consider compute the TC outer radius
    outer_wind_speed_threshold=2

    # Minimum distance between two critical centers in km. Default filter_center_threshold=400
    filter_center_threshold=400

    #Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=400
    dist_threshold=400

    # Critical outer radius in km to considerer a low pressure point as critical center. Default critical_outer_radius=100
    critical_outer_radius=0

    #resolution for radial legs in km. Default dr_res=100
    dr_res=100

    #resolution of angle steps for radial legs in degrees. Default d_ang=5
    d_ang=5

    #external search radius in km. Default rout=1000 km
    rout=1000

    #Terrain filter in m. Set terrain_filter=0 to not apply terrain filter. Default terrain_filter=0
    terrain_filter=0

    #TC maximum intensity threshold in m/s along the full trajectory. intensity_threshold=0
    intensity_threshold=0

    #Threshold for EC lifetime in hours. Default dt_lifetime=36
    dt_lifetime=36

    #TRelative vorticity threshold in 1/s to filter critical TCs centres. vorticity_threshold=1.45e-5
    vorticity_threshold=1.5e-5

    #Maximum slp treshold in hPa to filter tc centres. Deafult min_slp_threshold=1010 
    min_slp_threshold=1015

    #Minimum distance traveled in km by the system. Default minimum_distance_travelled=0
    minimum_distance_travelled=0

    #Radial distance (in degrees) for  cheking the MSLP increase, default great_circle_distance=5.5
    great_circle_distance=5.5

    #Change in MSLP (in Pa) over a distance of great-circle-distance from the candidate point, default dmslp_great_circle_distance=200
    dmslp_great_circle_distance=200

    #Radius (in km) for computing the maximum surface winds, default radius_for_msw=100
    radius_for_msw=100

    #Dates before the specific date and hour to compute the average mslp. Default prev_days=14
    prev_days=14

    #Mean sea level pressure anomaly threshold in hPa to consideded a grid point as candidate for system centre. Default mslp_anomaly_threshold=-2.5
    mslp_anomaly_threshold=-2.5


Getting input file template
---------------------------

You can get the input file template with default values using the following command:

.. code-block:: python

    import cytrack
    cytrack.get_cytrack_inputs_template()
 

Getting CyTRACK help on the input file
--------------------------------

To get the help of CyTRACK, you can use the following command:

.. code-block:: python

    import cytrack
    cytrack.help() 