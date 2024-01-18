import numpy as np
import matplotlib.pylab as plt
import sys
import os
from netCDF4 import Dataset,num2date,date2num
from datetime import datetime
from scipy.interpolate import griddata
from cytrack.cytrack_functions import *
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning)
from mpi4py import MPI
import time
import functools
print = functools.partial(print, flush=True)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

start_time = time.time()

def get_cytrack_main(pathfile=""):
	if rank==0:
		disclaimer()
		print("\n")
		print("Using parameters from: " + pathfile)
		print("============================================================================================================\n")
	content = imp.load_source("", pathfile) 
	verbose = check_paths(content, "verbose")
	cyclone_type = check_paths(content, "cyclone_type")
	source = check_paths(content, "source")
	wrfprefix = check_paths(content, "wrfprefix")
	era_file_prefix = check_paths(content, "era_file_prefix")
	era_upperfile_prefix = check_paths(content, "era_upperfile_prefix")
	era_date_file_name = check_paths(content, "era_date_file_name")
	use_mslp_anomaly=check_paths(content, "use_mslp_anomaly")
	model_res = check_paths(content, "model_res")
	path_data_source = check_paths(content, "path_data_source")
	checking_upper_levels_parameters = check_paths(content, "checking_upper_levels_parameters")
	vtl_vtu_lr = check_paths(content, "vtl_vtu_lr") 
	max_dist = check_paths(content, "max_dist")
	core_criteria_length = check_paths(content, "core_criteria_length")
	VTL_threshold=check_paths(content, "VTL_threshold")
	VTU_threshold=check_paths(content, "VTU_threshold")
	Bhart_threshold=check_paths(content, "Bhart_threshold")
	path_data_source_upper = check_paths(content, "path_data_source_upper")
	search_limits = check_paths(content, "search_limits")
	search_region = check_paths(content, "search_region")
	begin_year = check_paths(content, "begin_year")
	begin_month = check_paths(content, "begin_month")
	begin_day = check_paths(content, "begin_day")
	begin_hour = check_paths(content, "begin_hour")
	end_year = check_paths(content, "end_year")
	end_month = check_paths(content, "end_month")
	end_day = check_paths(content, "end_day")
	end_hour = check_paths(content, "end_hour")
	calendar =  check_paths(content, "calendar")
	dt_h = check_paths(content, "dt_h")
	path_out = check_paths(content, "path_out")
	rout = check_paths(content, "rout")
	dr_res = check_paths(content, "dr_res")
	d_ang = check_paths(content, "d_ang")
	tmpdir = check_paths(content, "tmp_dir")
	filter_center_threshold = check_paths(content, "filter_center_threshold")
	critical_outer_radius = check_paths(content, "critical_outer_radius")
	dist_threshold = check_paths(content, "dist_threshold")
	remove_tmp_dir = check_paths(content, "remove_tmp_dir")
	min_slp_threshold = check_paths(content, "min_slp_threshold")
	dt_lifetime = check_paths(content, "dt_lifetime")
	minimum_distance_travelled = check_paths(content, "minimum_distance_travelled")
	terrain_filter = check_paths(content, "terrain_filter")
	prev_days = check_paths(content, "prev_days")
	mslp_anomaly_threshold = check_paths(content, "mslp_anomaly_threshold")
	max_wind_speed_threshold = check_paths(content, "max_wind_speed_threshold")
	outer_wind_speed_threshold = check_paths(content, "outer_wind_speed_threshold")
	intensity_threshold = check_paths(content, "intensity_threshold")
	vorticity_threshold = check_paths(content, "vorticity_threshold")
	great_circle_distance=check_paths(content,"great_circle_distance")
	dmslp_great_circle_distance=check_paths(content,"dmslp_great_circle_distance")
	radius_for_msw=check_paths(content,"radius_for_msw")
	plotting_maps=check_paths(content,"plotting_maps")
	
	custom_file_prefix=check_paths(content,"custom_file_prefix")
	custom_date_file_name=check_paths(content,"custom_date_file_name")
	custom_upperfile_prefix = check_paths(content,"custom_upperfile_prefix")
	custom_mslp_variable = check_paths(content,"custom_mslp_variable")
	
	custom_uwind_variable = check_paths(content,"custom_uwind_variable")
	custom_vwind_variable = check_paths(content,"custom_vwind_variable")
	
	custom_latitude_var = check_paths(content,"custom_latitude_var")
	custom_longitude_var = check_paths(content,"custom_longitude_var")
	custom_geopotential_var_name = check_paths(content,"custom_geopotential_var_name")
	custom_upper_level_variable_name = check_paths(content,"custom_upper_level_variable_name")
	custom_terrain_high_filename = check_paths(content,"custom_terrain_high_filename")
	custom_terrain_high_var_name = check_paths(content,"custom_terrain_high_var_name")
	
	
	
	
	filter_center_threshold, critical_outer_radius, dist_threshold, verbose, path_out, tmpdir, source, dt_h, rout,dr_res,d_ang,remove_tmp_dir,format_date_file_name,min_slp_threshold,dt_lifetime,minimum_distance_travelled,terrain_filter,prev_days,mslp_anomaly_threshold,search_limits,search_region,max_wind_speed_threshold,outer_wind_speed_threshold,intensity_threshold,vorticity_threshold,checking_upper_levels_parameters,max_dist,great_circle_distance,dmslp_great_circle_distance,radius_for_msw,path_data_source_upper,path_data_source,vtl_vtu_lr,core_criteria_length,VTL_threshold,VTU_threshold,Bhart_threshold,plotting_maps,use_mslp_anomaly, calendar=check_default_parameters(cyclone_type=cyclone_type,
																																										verbose=verbose, 
																																										source=source,
																																										path_out=path_out,
																																										model_res=model_res,
																																										dt_h=dt_h,
																																										tmpdir=tmpdir,
																																										filter_center_threshold=filter_center_threshold,
																																										critical_outer_radius=critical_outer_radius,
																																										dist_threshold=dist_threshold,
																																										rout=rout,
																																										dr_res=dr_res,
																																										d_ang=d_ang,
																																										remove_tmp_dir=remove_tmp_dir,
																																										era_date_file_name=era_date_file_name,
																																										min_slp_threshold=min_slp_threshold,
																																										dt_lifetime=dt_lifetime,
																																										minimum_distance_travelled=minimum_distance_travelled,
																																										terrain_filter=terrain_filter,
																																										prev_days=prev_days,
																																										mslp_anomaly_threshold=mslp_anomaly_threshold,
																																										search_limits=search_limits,
																																										search_region=search_region,
																																										max_wind_speed_threshold=max_wind_speed_threshold,
																																										outer_wind_speed_threshold=outer_wind_speed_threshold,
																																										intensity_threshold=intensity_threshold,
																																										vorticity_threshold=vorticity_threshold,
																																										checking_upper_levels_parameters=checking_upper_levels_parameters,
																																										vtl_vtu_lr=vtl_vtu_lr,
																																										max_dist=max_dist,
																																										great_circle_distance=great_circle_distance,
																																										dmslp_great_circle_distance=dmslp_great_circle_distance,
																																										radius_for_msw=radius_for_msw,
																																										path_data_source_upper=path_data_source_upper,
																																										path_data_source=path_data_source,
																																										core_criteria_length=core_criteria_length,
																																										VTL_threshold=VTL_threshold,
																																										VTU_threshold=VTU_threshold,
																																										Bhart_threshold=Bhart_threshold,
																																										plotting_maps=plotting_maps,
																																										use_mslp_anomaly=use_mslp_anomaly,
																																										calendar=calendar
																																										)
	use_mslp_anomaly=str2boolean(use_mslp_anomaly)

	verbose=str2boolean(verbose)
	remove_tmp_dir=str2boolean(remove_tmp_dir)
	checking_upper_levels_parameters=str2boolean(checking_upper_levels_parameters)
	vtl_vtu_lr=str2boolean(vtl_vtu_lr)
	plotting_maps=str2boolean(plotting_maps)
	
	
	pathoutput=path_out+"/"+program_name()+"_output"
	tmpdir=tmpdir+"/"+program_name()+"_tmpdir_"+source.upper()+"_"+search_region.upper()
	
	
	if rank==0:
		if not os.path.exists(pathoutput):os.makedirs(pathoutput)
		if not os.path.exists(tmpdir):os.makedirs(tmpdir)
	
	
	

		
	dates,hours=get_dates_vectors(year_case_init=begin_year,
		month_case_init=begin_month,
		day_case_init=begin_day,
		hour_case_init=begin_hour,
		year_case_end=end_year,
		month_case_end=end_month,
		day_case_end=end_day,
		hour_case_end=end_hour,
		dt_h=dt_h,
		prev_days=prev_days,
		calendar=calendar
		)

	i_bg=get_i_bg(prev_days)
	
	#print(dates)
	#dates=dates[i_bg:]
	#hours=hours[i_bg:]
	
	#sys.exit()
	checking_optimum_nproc(rank=rank, n_proc=mpisize, length_files=len(dates[i_bg:]))
	if rank==0 and verbose:
		print("\nTracking " + cyclone_type.upper()+"s over " + search_region +" region using the " + source  +" meteorological fields")
		print(program_name() +" is running using "+ str(int(mpisize)) + " CPUs")

		print_dafault_values(cyclone_type=cyclone_type,
				verbose=verbose, 
				path_out=pathoutput,
				dt_h=dt_h,
				tmpdir=tmpdir,
				filter_center_threshold=filter_center_threshold,
				critical_outer_radius=critical_outer_radius,
				dist_threshold=dist_threshold,
				rout=rout,
				dr_res=dr_res,
				d_ang=d_ang,
				remove_tmp_dir=remove_tmp_dir,
				min_slp_threshold=min_slp_threshold,
				dt_lifetime=dt_lifetime,
				minimum_distance_travelled=minimum_distance_travelled,
				terrain_filter=terrain_filter,
				start_date=dates[i_bg],
				start_hour=hours[i_bg],
				end_date=dates[-1],
				end_hour=hours[-1],
				prev_days=prev_days,
				mslp_anomaly_threshold=mslp_anomaly_threshold,
				search_limits=search_limits,
				search_region=search_region,
				max_wind_speed_threshold=max_wind_speed_threshold,
				outer_wind_speed_threshold=outer_wind_speed_threshold,
				intensity_threshold=intensity_threshold,
				vorticity_threshold=vorticity_threshold,
				great_circle_distance=great_circle_distance,
				dmslp_great_circle_distance=dmslp_great_circle_distance,
				radius_for_msw=radius_for_msw,
				checking_upper_levels_parameters=checking_upper_levels_parameters,
				vtl_vtu_lr=vtl_vtu_lr,
				core_criteria_length=core_criteria_length,
				VTL_threshold=VTL_threshold,
				VTU_threshold=VTU_threshold,
				Bhart_threshold=Bhart_threshold,
				use_mslp_anomaly=use_mslp_anomaly
				)
		
	
		
		print("\nPROCESSING:")
		print("============================================================================================================\n")
	
	if source.upper()=="WRF":
		input_files=get_wrf_files(dates=dates,hours=hours,wrfprefix=wrfprefix)
		source_prefix=wrfprefix
		source_upperprefix=wrfprefix
	elif source.upper()=="ERA5":
		input_files=get_era5_files(dates=dates,hours=hours,era_file_prefix=era_file_prefix,era_date_file_name=era_date_file_name)
		source_prefix=era_file_prefix
	
	elif source.upper()=="CUSTOM":
		input_files=get_custom_files(dates=dates,hours=hours,custom_file_prefix=custom_file_prefix,custom_date_file_name=custom_date_file_name)
		source_prefix=custom_file_prefix
	
	if source.upper()=="ERA5" and checking_upper_levels_parameters==True:
		input_files_upper=get_era5_files(dates=dates,hours=hours,era_file_prefix=era_upperfile_prefix,era_date_file_name=era_date_file_name)
		source_upperprefix=era_upperfile_prefix
		input_files_upper=input_files_upper[i_bg::]
	elif source.upper()=="CUSTOM" and checking_upper_levels_parameters==True:
		input_files_upper=get_custom_files(dates=dates,hours=hours,custom_file_prefix=custom_upperfile_prefix,custom_date_file_name=custom_date_file_name)
		source_upperprefix=custom_upperfile_prefix
		input_files_upper=input_files_upper[i_bg::]
	elif source.upper()=="WRF" and checking_upper_levels_parameters==True:
		input_files_upper=input_files[i_bg::]
	else:
		source_upperprefix=""
	
	if rank==0 and verbose:
		
		print("\nChecking input files")
		for i in range(0,len(input_files)):
			checkfiles= checking_input_files(pathfile=path_data_source, source_file=input_files[i],source=source,date=dates[i],hour=hours[i],flev="sfc")
		if source.upper()=="ERA5" and checking_upper_levels_parameters==True:
			upperdates=dates[i_bg::]
			upperhours=hours[i_bg::]
			for i in range(0,len(input_files_upper)):
				checkfiles_upper= checking_input_files(pathfile=path_data_source_upper, source_file=input_files_upper[i],source=source,date=upperdates[i],hour=upperhours[i],flev="upper")
		else:
			checkfiles_upper=True

		if checkfiles:
			print("    ---> | Surface Files : PASSED")
		if checkfiles_upper==True and checking_upper_levels_parameters==True:
			print("    ---> | Upper   Files : PASSED")
			print("----------------------------------------------------------------------------------------")



	dates=dates[i_bg:]
	hours=hours[i_bg:]
	input_files=input_files[i_bg:]


	comm.barrier()
	count = len(dates) // mpisize
	remainder = len(dates) % mpisize
	if rank < remainder:
		start = rank * (count + 1)
		stop = start + count + 1
	else:
		start = rank * count + remainder
		stop = start + count

	nproc_dates=dates[start:stop]
	nproc_input_files=input_files[start:stop]
	nproc_hours=hours[start:stop]
	if checking_upper_levels_parameters:
		nproc_input_files_upper=input_files_upper[start:stop]
	else:
		nproc_input_files_upper=np.empty_like(nproc_input_files)
		nproc_input_files_upper[:]=""

	if rank==0 and verbose:

		print("\nProcessing input files")
		print("----------------------------------------------------------------------------------------")
	time.sleep(1)
	if cyclone_type.upper() in ("TC","MC","TLC","SC","EC"):
		tracker_cyclones(cyclone_type=cyclone_type,
			source=source,
			idir=path_data_source,
			sourcefiles=nproc_input_files, 
			pathoutput=pathoutput,
			rout=rout,
			verbose=verbose,
			dates=nproc_dates,
			hours=nproc_hours,
			model_res=model_res,
			search_limits=search_limits,
			dr_res=dr_res,
			d_ang=d_ang,
			filter_center_threshold=filter_center_threshold,
			critical_outer_radius=critical_outer_radius,
			tmpdir=tmpdir,rank=rank, 
			search_region=search_region,
			min_slp_threshold=min_slp_threshold,
			terrain_filter=terrain_filter,
			prev_days=prev_days,
			mslp_anomaly_threshold=mslp_anomaly_threshold,
			source_filename_prefix=source_prefix,
			source_file_date_format=era_date_file_name,
			max_wind_speed_threshold=max_wind_speed_threshold,
			outer_wind_speed_threshold=outer_wind_speed_threshold,
			vorticity_threshold=vorticity_threshold,
			great_circle_distance=great_circle_distance,
			dmslp_great_circle_distance=dmslp_great_circle_distance,
			radius_for_msw=radius_for_msw,
			sourcefilesupper=nproc_input_files_upper,
			checking_upper_levels_parameters=checking_upper_levels_parameters,
			idir_upper=path_data_source_upper,
			plotting_maps=plotting_maps,
			use_mslp_anomaly=use_mslp_anomaly,
			custom_mslp_variable=custom_mslp_variable,
			custom_latitude_var=custom_latitude_var,
			custom_longitude_var=custom_longitude_var,
			custom_uwind_variable=custom_uwind_variable,
			custom_vwind_variable=custom_vwind_variable,
			custom_terrain_high_filename=custom_terrain_high_filename,
			custom_terrain_high_var_name=custom_terrain_high_var_name,
			)

	#if cyclone_type.upper()=="MC":
		#tracker_MC(source=source,
			#idir=path_data_source,
			#sourcefiles=nproc_input_files,
			#pathoutput=pathoutput,
			#rout=rout,
			#verbose=verbose,
			#dates=nproc_dates,
			#hours=nproc_hours,
			#model_res=model_res,
			#search_limits=search_limits,
			#dr_res=dr_res,
			#d_ang=d_ang,
			#filter_center_threshold=filter_center_threshold,
			#critical_outer_radius=critical_outer_radius,
			#tmpdir=tmpdir,rank=rank,
			#search_region=search_region,
			#min_slp_threshold=min_slp_threshold,
			#terrain_filter=terrain_filter,
			#prev_days=prev_days,
			#mslp_anomaly_threshold=mslp_anomaly_threshold,
			#source_filename_prefix=source_prefix,
			#source_file_date_format=era_date_file_name,
			#max_wind_speed_threshold=max_wind_speed_threshold,
			#outer_wind_speed_threshold=outer_wind_speed_threshold,
			#vorticity_threshold=vorticity_threshold,
			#great_circle_distance=great_circle_distance,
			#dmslp_great_circle_distance=dmslp_great_circle_distance,
			#radius_for_msw=radius_for_msw)





	comm.barrier()
	if rank==0:
		if verbose:

			print("\n\nLinking " + cyclone_type.upper() +" critical centers to contruct the full trajectory")
		sys_id=paring_centers(cyclone_type=cyclone_type,
				dates=dates,
				hours=hours,
				dist_threshold=dist_threshold,
				tmpdir=tmpdir,
				pathoutput=pathoutput,
				search_region=search_region,
				dt_h=dt_h,
				source=source,
				dt_lifetime=dt_lifetime,
				minimum_distance_travelled=minimum_distance_travelled,
				intensity_threshold=intensity_threshold,
				checking_upper_levels_parameters=checking_upper_levels_parameters,
				idir_upper=path_data_source_upper,
				source_upperprefix=source_upperprefix,
				era_date_file_name=era_date_file_name,
				search_limits=search_limits,
				max_dist=max_dist,
				vtl_vtu_lr=vtl_vtu_lr,
				core_criteria_length=core_criteria_length,
				VTL_threshold=VTL_threshold,
				VTU_threshold=VTU_threshold,
				Bhart_threshold=Bhart_threshold,
				custom_geopotential_var_name=custom_geopotential_var_name,
				custom_upper_level_variable_name=custom_upper_level_variable_name,
				custom_varlat=custom_latitude_var,
				custom_varlon=custom_longitude_var,
				custom_date_file_name=custom_date_file_name
				) 
				
		if verbose:
			print("\n " + program_name() +" found " + str(int(sys_id)) +" " + cyclone_type.upper() + "s")
			print("---------------------------------------------------------------------------------------")

	time.sleep(1)
	if rank==0 and remove_tmp_dir:
		os.system("rm -r "+tmpdir)


	elapsed_time = time.time() - start_time
	if rank==0:
		if verbose:
			print("\nRun time: %.2f seconds." % np.round(elapsed_time, 2))
		ending_credits()

