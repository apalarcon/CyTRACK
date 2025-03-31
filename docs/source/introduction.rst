
CyTRACK: Cyclone Tracking framework
=================================
It is a new open-source, comprehensive and user-friendly Python toolbox for detecting and tracking cyclones in model and reanalysis datasets.


.. image:: _static/LogoV1.png
   :alt: CyTRACK logo
   :align: center
   :width: 400px


CyTRACK detects and tracks cyclones also using the MSLP. Like most algorithmic Lagrangian trackers, the procedure for detecting cyclones is divided into two parts: (i) critical centres detection
and (ii) pairing centres in continuous time steps. These procedures are described below.


Detecting critical centres
----------------------------

Critical centres are detected as mean sea level pressure (MSLP) minima and must satisfy the following conditions:

- The MSLP value is lower than a specified MSLP threshold (**min_slp_threshold**).
- The MSLP anomaly computed as the difference between the MSLP at time t0 and the mean MSLP in the previous N days (**prev_days**) is lower than a critical value (`mslp_anomaly_threshold`).
- The MSLP increases (**dmslp_great_circle_distance**) over a specified distance (**great_circle_distance**) from the candidate point.
- The maximum wind speed (**max_wind_speed_threshold**) around the critical centre (**radius_for_msw**) is higher than a predefined value.
- The surface relative vorticity (**vorticity_threshold**) is greater than a critical threshold.
-  The mean radial distance to the last closed isobar is higher than a critical value (**critical_outer_radius**). The radial distance to the last closed isobar is computed following the procedure developed by Rudeva and Gulev (2007).
-  When several centres exist in a critical radial distance (**filter_­center_threshold**), the centre that has the lowest MSLP is retained.

Additionally, CyTRACK allows discarding centres positioned over terrain higher than a critical high (**terrain_filter**).

.. seealso::

   Rudeva, I., Gulev, S.K., 2007. Climatology of cyclone size characteristics and their changes during the cyclone life cycle. Mon. Weather Rev. 135, 2568–2587. https://
doi.org/10.1175/MWR3420.1.



Paring cyclone centres in continuous time steps
---------------------------------------------

- Storm centres are linked together if they reoccur in the next time step (**t=t0 + dt**) within a critical distance (**dist_threshold**) from the previous low-pressure centre detected at time t0.
- If there are multiple identified centres within the critical distance, then the point with the lowest MSLP is chosen as the cyclone centre at the second time step. 
- CyTRACK also allows a one-time step gap. If no centre is detected at time t0+dt, CyTRACK searches for a candidate point at time t0+2dt. If at least one point is found at t0+2dt within a radial distance of two times dist_threshold from the cyclone centre at time t0, then the cyclone centre at time t0+dt is computed as the average latitude and longitude at time t0 and t0+2dt. To account for the “natural evolution” of the cyclone trajectory, the angle between the lines formed by the centres at time t0 and t0+dt and t0+dt and t0+2dt must be less than 10°. If the last condition is satisfied, the algorithm continues searching for the next centre at time t0+3dt; otherwise, the track ends at time t0. 
- After that, CyTRACK evaluates the lifetime (**dt_lifetime**) of the cyclone, the minimum distance travelled (**minimum_distance_travelled**) from genesis to dissipation, the maximum intensity (**intensity_threshold**) in terms of the maximum wind speed along the track.
- If **checking_upper_levels_parameters** is set to “yes” in the configuration file, CyTRACK classifies the cyclone core based on the thermal wind and thermal asymmetry according to the CPS.



For a more detailed understanding of CyTRACK, Please refer to 

.. seealso::

   Pérez-Alarcón, A.; Coll-Hidalgo, P.; Trigo, R.M.; Nieto, R.; Gimeno, L. (2024). CyTRACK: An open-source and user-friendly python toolbox for detecting and tracking cyclones. Environmental Modelling & Software, 176, 106027. https://doi.org/10.1016/j.envsoft.2024.106027.
