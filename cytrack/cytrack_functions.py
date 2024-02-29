import numpy as np
import sys
import os
from netCDF4 import Dataset,num2date,date2num
from datetime import datetime, timedelta
from scipy.interpolate import griddata
from sklearn.linear_model import LinearRegression
from scipy import interpolate
from numpy.core.numeric import normalize_axis_index
import imp
import argparse
import xarray as xr
import time
import math
import scipy.ndimage as sp
import matplotlib.collections as collections
import matplotlib.pylab as plt
import warnings
import functools
from scipy.interpolate import interp1d
from scipy.ndimage import maximum_filter, minimum_filter
print = functools.partial(print, flush=True)
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning) 


def program_name():
	return "CyTRACK"

def program_fullname():
	return "Cyclone TRACKing"


def get_currentversion():
	pathpkg = os.path.dirname(__file__)
	version_file = pathpkg+"/VERSION"
	with open(version_file) as vfile:
		version = vfile.readlines()[0].strip()
	return(version)


def get_lsatupdate():
        lupathpkg = os.path.dirname(__file__)
        version_upd = lupathpkg+"/LAST_UPDATE"
        with open(version_upd) as ufile:
                uversion = ufile.readlines()[0].strip()
        return(uversion)



def disclaimer():
	print("\n============================================================================================================")
	print("||                                                                                                        ||")
	print("||              +++++++                                                                                   ||")
	print("||            +++++++         +++++++           +++++++  +++++     +     +++++++  +    +                  ||")
	print("||           ++++             +        +     +     +     +   +    + +    +        +   +                   ||")
	print("||         ++++++             +         +   +      +     +++++   +++++   +        + +                     ||")
	print("||         ++++++++           +           +        +     + +    +     +  +        +   +                   ||")
	print("||        ++++++++++          +++++++     +        +     +  +   +     +  +++++++  +    +                  ||")
	print("||       +++++ +++++       <-------------------------------------------------------------->               ||")
	print("||       ++++    ++++                        " +   program_fullname() +" Version " +str(get_currentversion()) +"                               ||")
	print("||       +++++  +++++                                                                                     ||")
	print("||        ++++++++++                                                                                      ||")
	print("||         +++++++++    " +  program_name() + " Version " +str(get_currentversion())+ " is free under the terms of the GNU General Public license   ||")            
	print("||          +++++++                     EPhysLab (Environmental Physics Laboratory)                       ||")   
	print("||            +++++                             Universidade de Vigo                                      ||")
	print("||         +++++++                    contact: albenis.perez.alarcon@uvigo.es                             ||")
	print("||      +++++++                                                                                           ||")
	print("||                                                                         ** Last Update: " +    get_lsatupdate() +" **  ||")
	print("============================================================================================================")
	
def read_args():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"--parameterfile",
		"-pf",
		help="name of parameters file. To use " + program_name() + " run: python run_"+program_name()+".py -pf <your_input_file>",
		metavar="",
		type=str,
		default="ginputs",
	)
	parser.add_argument(
		"--cytrack_help",
		"-cth",
		help="Help for " + program_name() + " input parameters. Run: python run_"+program_name()+".py -cth t",
		metavar="",
		type=str2boolean,
		default=False,
	)
	parser.add_argument(
		"--get_template",
		"-gt",
		help="Get a template for the " + program_name() + " input parameters. Run: python run_"+program_name()+".py -gt t or python run_"+program_name()+".py --get_template t",
		metavar="",
		type=str2boolean,
		default=False,
	)
	args = parser.parse_args()
	return args


def get_cytrack_inputs_template():
	disclaimer()
	wpath=os.getcwd()
	pathpkg = os.path.dirname(__file__)
	os.system("cp "+pathpkg+"/cytrack_inputs.cfg " + wpath+"/cytrack_inputs_template.cfg")
	print("\n============================================================================================================")
	print("cytrack_inputs.cfg template was successfully copied to the working directory: " + wpath)
	print("Bye :)")
	print("============================================================================================================")


def ending_credits():
	print("\n\n============================================================================================================")
	print(program_name() +" Version " +str(get_currentversion()) + " has successfully finished")
	print("Bye :)")
	print("============================================================================================================")


def get_limits_by_region(search_region="",cyclone_type=""):
	if cyclone_type.upper()=="EC":
		if search_region.upper()=="NA":
			search_limits=[-120,20,25,75]
		elif search_region.upper()=="SA":
			search_limits=[-80,-75,30,25]
		elif  search_region.upper()=="NP":
			search_limits=[100,20,240,75]
		elif search_region.upper()=="SP":
			search_limits=[130,-75,280,-25]
		elif search_region.upper()=="SI":
			search_limits=[30,-75,130,-25]
		elif search_region.upper()=="NH":
			search_limits=[0,20,360,75]
		elif search_region.upper()=="SH":
			search_limits=[0,-75,360,-25]
		elif search_region.upper()=="GL":
			search_limits=[0,-90,360,90]
	if cyclone_type.upper()=="TC":
		if search_region.upper()=="AL":
			search_limits=[-110,0,10,55]
		elif search_region.upper()=="EP":
			search_limits=[180,0,270,50]
		elif search_region.upper()=="WP":
			search_limits=[100,0,180,55]
		elif search_region.upper()=="SP":
			search_limits=[130,-50,270,0]
		elif search_region.upper()=="NI":
			search_limits=[30,0,100,40]
		elif search_region.upper()=="SI":
			search_limits=[30,-50,130,0]
		elif search_region.upper()=="SA":
			search_limits=[-60,-45,20,0]
		elif search_region.upper()=="NH":
			search_limits=[0,0,360,60]
		elif search_region.upper()=="SH":
			search_limits=[0,-60,360,0]
		elif search_region.upper()=="GL":
			search_limits=[0,-60,360,60]
	if cyclone_type.upper() in ("MC", "TLC"):
			search_limits=[-5,20,50,50]
			
	return search_limits
def check_default_parameters(cyclone_type="",
				verbose='true', 
				source="ERA5",
				path_out="./",
				model_res=20,
				dt_h=6,
				tmpdir="./",
				filter_center_threshold=1500,
				critical_outer_radius=100,
				dist_threshold=1000,
				rout=2000,
				dr_res=100,
				d_ang=10,
				remove_tmp_dir='yes',
				era_date_file_name="yyyymmdd_hh",
				min_slp_threshold=1015,
				dt_lifetime=24,
				minimum_distance_travelled=1000,
				terrain_filter=1000,
				prev_days=14,
				mslp_anomaly_threshold=-3,
				search_limits=[None,None,None,None],
				search_region="",
				max_wind_speed_threshold=10,
				outer_wind_speed_threshold=2.5,
				intensity_threshold=17.5,
				vorticity_threshold=1e-5,
				checking_upper_levels_parameters=False,
				vtl_vtu_lr=False,
				core_criteria_length=0,
				VTL_threshold=0,
				VTU_threshold=0,
				Bhart_threshold=10,
				max_dist=500,
				great_circle_distance=5.5,
				dmslp_great_circle_distance=200,
				radius_for_msw=100,
				plotting_maps=False,
				path_data_source_upper="./",
				path_data_source="./",
				use_mslp_anomaly=True,
				calendar="366d"):
	
	
	if cyclone_type=="":
		print_error_message("cyclone_type is not defined")
		
	if search_region=="":
		search_region="GL"
		
	
	if plotting_maps == "" or plotting_maps==None:
		plotting_maps=False
	
	if search_limits==None or search_limits=="":
		search_limits=get_limits_by_region(search_region=search_region,cyclone_type=cyclone_type)
		
	if checking_upper_levels_parameters==None or checking_upper_levels_parameters=="":
		checking_upper_levels_parameters=False
	
	if vtl_vtu_lr==None or vtl_vtu_lr == "":
		vtl_vtu_lr=False
	
	if use_mslp_anomaly=="" or use_mslp_anomaly == None:
		use_mslp_anomaly=True

	
	if verbose!="no":
		verbose="yes"
	if source=="":
		source="ERA5"
	if path_out=="":
		wpath = os.getcwd()
		os.chdir(wpath)
		path_out=wpath
	if tmpdir=="":
		wpath = os.getcwd()
		os.chdir(wpath)
		tmpdir=wpath
		
	if calendar=="" or calendar==None:
		calendar="366d"
	
	if source.upper()=="ERA5":
		if path_data_source=="":
			swpath = os.getcwd()
			os.chdir(swpath+"/ERA5_input")
			path_data_source=swpath
		if checking_upper_levels_parameters==True and 	path_data_source_upper=="":
			wpath = os.getcwd()
			os.chdir(swpath+"/ERA5_input_upper")
			path_data_source_upper=swpath
	
	if source.upper()=="WRF" and path_data_source=="":
		print_error_message("ERROR: Please define input directory for WRF data files")
	
	if source.upper()=="WRF" and  checking_upper_levels_parameters==True and path_data_source_upper=="":
		path_data_source_upper=path_data_source
	
	if cyclone_type.upper()=="TLC":
		checking_upper_levels_parameters=True
		vtl_vtu_lr=True

		
	if dt_h==None or dt_h=="":
		dt_h = 6
	if prev_days==None or  prev_days=="":
		prev_days=14
	elif prev_days<=0:
		prev_days=0
		use_mslp_anomaly=False
	
	use_mslp_anomaly=str2boolean(use_mslp_anomaly)
	if use_mslp_anomaly==False:
		prev_days=0
	
	
		
	if cyclone_type.upper()=="EC":
		if mslp_anomaly_threshold==None or mslp_anomaly_threshold=="":
			mslp_anomaly_threshold=-2.5

		if filter_center_threshold==None or filter_center_threshold=="":
			filter_center_threshold=1000
		if critical_outer_radius==None or critical_outer_radius=="":
			critical_outer_radius=100
		if dist_threshold==None or dist_threshold=="None":
			dist_threshold=1000
		if rout==None or rout=="":
			rout=2000
		if dr_res==None or dr_res=="":
			dr_res=100
		if d_ang==None or d_ang=="":
			d_ang=10
		if remove_tmp_dir!="no":
			remove_tmp_dir="yes"
		if era_date_file_name=="":
			era_date_file_name="yyyymmdd_hh"
		if min_slp_threshold==None or min_slp_threshold=="":
			min_slp_threshold=1015

		if dt_lifetime==None or dt_lifetime=="":
			dt_lifetime=48
		if minimum_distance_travelled==None or minimum_distance_travelled=="":
			minimum_distance_travelled=1000

		if terrain_filter==None or terrain_filter=="":
			terrain_filter=1000
		elif terrain_filter<0:
			print_error_message("terrain_filter must be equal or higher than zero")

		if intensity_threshold==None or intensity_threshold=="":
			intensity_threshold=0,

		if max_wind_speed_threshold==None or max_wind_speed_threshold=="":
				max_wind_speed_threshold=0
		if outer_wind_speed_threshold==None or outer_wind_speed_threshold=="":
				outer_wind_speed_threshold=0
		if intensity_threshold==None or intensity_threshold=="":
				intensity_threshold=0,
		if 	vorticity_threshold==None or vorticity_threshold=="":
				vorticity_threshold=1.45e-5
		
		if max_dist==None or max_dist=="":
			max_dist=500
		
		if great_circle_distance==None or great_circle_distance=="":
			great_circle_distance=6.5
		
		if dmslp_great_circle_distance==None or dmslp_great_circle_distance == "":
			dmslp_great_circle_distance=200
			
		if radius_for_msw==None or radius_for_msw=="":
			radius_for_msw=250

		if core_criteria_length==None or core_criteria_length=="":
			core_criteria_length=0
		elif core_criteria_length!=-99 and core_criteria_length<0:
			core_criteria_length=0

		if VTL_threshold=="" or VTL_threshold==None:
			VTL_threshold=0
		if VTU_threshold=="" or  VTU_threshold==None:
			VTU_threshold=0
		if Bhart_threshold=="" or Bhart_threshold==None:
			Bhart_threshold=10




	elif cyclone_type.upper()=="TC":

		if mslp_anomaly_threshold==None or mslp_anomaly_threshold=="":
			mslp_anomaly_threshold=-2.

		if filter_center_threshold==None or filter_center_threshold=="":
			filter_center_threshold=650
		if critical_outer_radius==None or critical_outer_radius=="":
			critical_outer_radius=50
		if dist_threshold==None or dist_threshold=="None":
			dist_threshold=450
		if rout==None or rout=="":
			rout=1000
		if dr_res==None or dr_res=="":
			dr_res=100
		if d_ang==None or d_ang=="":
			d_ang=10
		if remove_tmp_dir!="no":
			remove_tmp_dir="yes"
		if era_date_file_name=="":
			era_date_file_name="yyyymmdd_hh"
		if min_slp_threshold==None or min_slp_threshold=="":
			min_slp_threshold=1015

		if dt_lifetime==None or dt_lifetime=="":
			dt_lifetime=36
		if minimum_distance_travelled==None or minimum_distance_travelled=="":
			minimum_distance_travelled=0

		if terrain_filter==None or terrain_filter=="":
			terrain_filter=0
		elif terrain_filter<0:
			print_error_message("terrain_filter must be equal or higher than zero")
			
		if max_wind_speed_threshold==None or max_wind_speed_threshold=="":
				max_wind_speed_threshold=8
		if outer_wind_speed_threshold==None or outer_wind_speed_threshold=="":
				outer_wind_speed_threshold=6
		if intensity_threshold==None or intensity_threshold=="":
				intensity_threshold=10,
		if 	vorticity_threshold==None or vorticity_threshold=="":
				vorticity_threshold=1.45e-5
		
		if max_dist==None or max_dist=="":
			max_dist=500
		
		if great_circle_distance==None or great_circle_distance=="":
			great_circle_distance=5.5
		
		if dmslp_great_circle_distance==None or dmslp_great_circle_distance == "":
			dmslp_great_circle_distance=200
			
		if radius_for_msw==None or radius_for_msw=="":
			radius_for_msw=100
			
		if core_criteria_length==None or core_criteria_length=="":
			core_criteria_length=3
		elif core_criteria_length!=-99 and core_criteria_length<0:
			core_criteria_length=3
		
		
		if VTL_threshold=="" or VTL_threshold==None:
			VTL_threshold=0
		if VTU_threshold=="" or  VTU_threshold==None:
			VTU_threshold=0
		if Bhart_threshold=="" or Bhart_threshold==None:
			Bhart_threshold=10
		
			
	elif cyclone_type.upper() in ("MC","TLC"):

		if mslp_anomaly_threshold==None or mslp_anomaly_threshold=="":
			mslp_anomaly_threshold=-2.5

		if filter_center_threshold==None or filter_center_threshold=="":
			filter_center_threshold=200
		if critical_outer_radius==None or critical_outer_radius=="":
			critical_outer_radius=50
		if dist_threshold==None or dist_threshold=="None":
			dist_threshold=200
		if rout==None or rout=="":
			rout=800
		if dr_res==None or dr_res=="":
			dr_res=100
		if d_ang==None or d_ang=="":
			d_ang=10
		if remove_tmp_dir!="no":
			remove_tmp_dir="yes"
		if era_date_file_name=="":
			era_date_file_name="yyyymmdd_hh"
		if min_slp_threshold==None or min_slp_threshold=="":
			min_slp_threshold=1015

		if dt_lifetime==None or dt_lifetime=="":
			dt_lifetime=24
		if minimum_distance_travelled==None or minimum_distance_travelled=="":
			minimum_distance_travelled=200

		if terrain_filter==None or terrain_filter=="":
			terrain_filter=0
		elif terrain_filter<0:
			print_error_message("terrain_filter must be equal or higher than zero")
			
		if max_wind_speed_threshold==None or max_wind_speed_threshold=="":
				max_wind_speed_threshold=0
		if outer_wind_speed_threshold==None or outer_wind_speed_threshold=="":
				outer_wind_speed_threshold=0
		if intensity_threshold==None or intensity_threshold=="":
				intensity_threshold=0
		if 	vorticity_threshold==None or vorticity_threshold=="":
				vorticity_threshold=1.45e-5
		
		if max_dist==None or max_dist=="":
			max_dist=300
		
		if great_circle_distance==None or great_circle_distance=="":
			great_circle_distance=3
		
		if dmslp_great_circle_distance==None or dmslp_great_circle_distance == "":
			dmslp_great_circle_distance=200
			
		if radius_for_msw==None or radius_for_msw=="":
			radius_for_msw=100
			
		if cyclone_type.upper() in ("TLC"):
			if core_criteria_length==None or core_criteria_length=="":
				core_criteria_length=1
			elif core_criteria_length!=-99 and core_criteria_length<0:
				core_criteria_length=1
				
			if VTL_threshold=="" or VTL_threshold==None:
				VTL_threshold=0
			if VTU_threshold=="" or  VTU_threshold==None:
				VTU_threshold=0
			if Bhart_threshold=="" or Bhart_threshold==None:
				Bhart_threshold=10
				
				
		elif cyclone_type.upper() in ("MC"):
			if core_criteria_length==None or core_criteria_length=="":
				core_criteria_length=0
			elif core_criteria_length!=-99 and core_criteria_length<0:
				core_criteria_length=0
				
			if VTL_threshold=="" or VTL_threshold==None:
				VTL_threshold=0
			if VTU_threshold=="" or  VTU_threshold==None:
				VTU_threshold=0
			if Bhart_threshold=="" or Bhart_threshold==None:
				Bhart_threshold=10
			
	elif cyclone_type.upper()=="SC":

		if mslp_anomaly_threshold==None or mslp_anomaly_threshold=="":
			mslp_anomaly_threshold=-2.5

		if filter_center_threshold==None or filter_center_threshold=="":
			filter_center_threshold=400
		if critical_outer_radius==None or critical_outer_radius=="":
			critical_outer_radius=0
		if dist_threshold==None or dist_threshold=="None":
			dist_threshold=400
		if rout==None or rout=="":
			rout=1000
		if dr_res==None or dr_res=="":
			dr_res=100
		if d_ang==None or d_ang=="":
			d_ang=10
		if remove_tmp_dir!="no":
			remove_tmp_dir="yes"
		if era_date_file_name=="":
			era_date_file_name="yyyymmdd_hh"
		if min_slp_threshold==None or min_slp_threshold=="":
			min_slp_threshold=1015

		if dt_lifetime==None or dt_lifetime=="":
			dt_lifetime=36
		if minimum_distance_travelled==None or minimum_distance_travelled=="":
			minimum_distance_travelled=0

		if terrain_filter==None or terrain_filter=="":
			terrain_filter=0
		elif terrain_filter<0:
			print_error_message("terrain_filter must be equal or higher than zero")
			
		if max_wind_speed_threshold==None or max_wind_speed_threshold=="":
				max_wind_speed_threshold=0
		if outer_wind_speed_threshold==None or outer_wind_speed_threshold=="":
				outer_wind_speed_threshold=2
		if intensity_threshold==None or intensity_threshold=="":
				intensity_threshold=0,
		if 	vorticity_threshold==None or vorticity_threshold=="":
				vorticity_threshold=1.5e-5
		
		if max_dist==None or max_dist=="":
			max_dist=500
		
		if great_circle_distance==None or great_circle_distance=="":
			great_circle_distance=5.5
		
		if dmslp_great_circle_distance==None or dmslp_great_circle_distance == "":
			dmslp_great_circle_distance=200
			
		if radius_for_msw==None or radius_for_msw=="":
			radius_for_msw=100
			
		if core_criteria_length==None or core_criteria_length=="":
			core_criteria_length=7
		elif core_criteria_length!=-99 and core_criteria_length<0:
			core_criteria_length=7
		
		if VTL_threshold=="" or VTL_threshold==None:
			VTL_threshold=-50
		if VTU_threshold=="" or  VTU_threshold==None:
			VTU_threshold=-10
		if Bhart_threshold=="" or Bhart_threshold==None:
			Bhart_threshold=25
		
		
			
	return filter_center_threshold, critical_outer_radius, dist_threshold, verbose, path_out, tmpdir, source, dt_h,rout,dr_res,d_ang,remove_tmp_dir,era_date_file_name, min_slp_threshold,dt_lifetime, minimum_distance_travelled,terrain_filter,prev_days,mslp_anomaly_threshold,search_limits,search_region,max_wind_speed_threshold,outer_wind_speed_threshold,intensity_threshold,vorticity_threshold,checking_upper_levels_parameters,max_dist,great_circle_distance,dmslp_great_circle_distance,radius_for_msw,path_data_source_upper,path_data_source, vtl_vtu_lr, core_criteria_length,VTL_threshold,VTU_threshold,Bhart_threshold,plotting_maps,use_mslp_anomaly, calendar

	

def print_dafault_values(cyclone_type="",
			verbose="yes", 
			path_out="./",
			dt_h=6,
			tmpdir=".",
			filter_center_threshold=None,
			critical_outer_radius=None,
			dist_threshold=None,
			rout=None,
			dr_res=None,
			d_ang=None,
			remove_tmp_dir="",
			min_slp_threshold=None,
			dt_lifetime=None,
			minimum_distance_travelled=None,
			terrain_filter=None,
			start_date="",
			start_hour="",
			end_date="",
			end_hour="",
			prev_days=None,
			mslp_anomaly_threshold=None,
			search_limits=[None,None,None,None],
			search_region="",
			max_wind_speed_threshold=10,
			outer_wind_speed_threshold=2.5,
			intensity_threshold=17.5,
			vorticity_threshold=1.45e-5,
			great_circle_distance=5.5,
			dmslp_great_circle_distance=200,
			radius_for_msw=100,
			checking_upper_levels_parameters=False,
			vtl_vtu_lr=False,
			core_criteria_length=0,
			VTL_threshold=0,
			VTU_threshold=0,
			Bhart_threshold=10,
			use_mslp_anomaly=True
			):
		
	if  str2boolean(verbose):
		print("\nParameters for "+program_name() + " runs")
		print("============================================================================================================")
		print("+ begin_date: " + start_date + " at " + start_hour.zfill(2) + " UTC")
		print("+ end_date:  " + end_date + " at " + end_hour.zfill(2) + " UTC")
		print("+ dt_h: " + str(int(dt_h)) + " hours")
		print("+ cyclone_type: "+ cyclone_type.upper())
		print("+ search_region " + search_region)
		print("+ search_limits :",search_limits)
		print("+ use_mslp_anomaly: ", use_mslp_anomaly)
		print("+ filter_center_threshold: " + str(filter_center_threshold) + " km")
		print("+ critical_outer_radius: " + str(critical_outer_radius) + " km")
		print("+ dist_threshold: " + str(dist_threshold) + " km")
		print("+ rout: " + str(rout) + " km")
		print("+ dr_res: " + str(dr_res) + " km")
		print("+ d_ang: " + str(d_ang) + " degrees")
		print("+ min_slp_threshold: " + str(min_slp_threshold) + " hPa")
		print("+ dt_lifetime: " + str(dt_lifetime)+ " hours")
		print("+ minimum_distance_travelled: " + str(minimum_distance_travelled) + " km")
		print("+ terrain_filter: " + str(terrain_filter) + " m")
		print("+ prev_days: " + str(prev_days) + " days")
		print("+ mslp_anomaly_threshold: " + str(mslp_anomaly_threshold) + " hPa")
		print("+ max_wind_speed_threshold: " + str(max_wind_speed_threshold) + " m/s")
		print("+ vorticity_threshold: " + str(vorticity_threshold) + " 1/s")
		print("+ outer_wind_speed_threshold: " + str(outer_wind_speed_threshold) + " m/s")
		print("+ intensity_threshold: " + str(intensity_threshold) + " m/s")
		print("+ great_circle_distance: " + str(great_circle_distance) + " degrees")
		print("+ dmslp_great_circle_distance: " + str(dmslp_great_circle_distance) + " Pa")
		print("+ radius_for_msw: " + str(radius_for_msw) + " km")
		print("+ checking_upper_levels_parameters: "+ str(checking_upper_levels_parameters))
		print("+ vtl_vtu_lr: " + str(vtl_vtu_lr))
		print("+ core_criteria_length: " + str(core_criteria_length))
		print("+ VTL_threshold : " + str(VTL_threshold))
		print("+ VTU_threshold: "  + str(VTU_threshold))
		print("+ Bhart_threshold: " + str(Bhart_threshold))
		print("+ Path to save "+ program_name() +" output: ", path_out)
		print("+ Path to save "+ program_name() +" temporal files: ", tmpdir)



def help():
	disclaimer()
	print("")
	print(program_name() + " HELP FOR INPUT PARAMETERS")
	print("------------------------------------------------------------------------------------------------------------")
	print("+ verbose              = 'True'or 'False'                -> Print info during " + program_name()+ " runs. Default value ['True']\n")
	
	print("\nStart date configuration")
	print("..........................................")
	print("+ begin_year           = 'yyyy'                          -> Year of the date for starting run")
	print("+ begin_month          = 'mm'                            -> Month of the date for starting run")
	print("+ begin_day            = 'dd'                            -> Day of the date for starting run")
	print("+ begin_hours          = 'hh'                            -> Hour of the date for starting run\n")
	
	print("\nEnd date configuration")
	print("..........................................")
	print("+ end_year             = 'yyyy'                          -> Year of the date for ending run")
	print("+ end_month            = 'mm'                            -> Month of the date for ending run")
	print("+ end_day              = 'dd'                            -> Day of the date for ending run")
	print("+ end_hour             = 'hh'                            -> Hour of the date for ending run\n")
	
	
	print("\nTime step for "+program_name() + " runs in hours")
	print("..........................................")
	print("+ dt_h                 = h                               -> Time step in hours. Default dt_h=6\n")
	
	print("\n Type of calendar for "+program_name() + " runs")
	print("+ calendar             = '365d' / '366d'                 -> Calendar: 365d: Remove February 29 in leap years\n")
	
	
	print("\nSystem tracking")
	print("..........................................")
	print("+ cyclone_type         = 'EC'/'TC'/'MC'/'TLC'/'SC'       -> Type of low pressure system for tracking: EC: extratropical cyclones")
	print("** Note: In the current version of " + program_name()+ " ("+str(get_currentversion())+")" + " is only available the tracking for extratropical cyclones\n" )
	print("+ use_mslp_anomaly     = 'yes'/'no'                      -> Apply MSLP anomaly to filter critical cyclones center. Default use_mslp_anomaly='yes'")
	
	
	
	print("\nSource information")
	print("..........................................")
	print("+ source               = 'WRF' or 'ERA5' or 'CUSTOM'     -> source of the mslp data for tracking cyclones")
	print("+ path_data_source     = 'path'                          -> Path to mslp dataset")
	print("+ wrfprefix            = 'wrfout_d01'                    -> Only for WRF. Prefix in WRF files before the date")
	print("+ era_file_prefix      = 'era_prefix'                    -> Only for ERA5. Prefix in ERA5 files before the date")
	print("+ era_date_file_name   = 'yyyymmdd_hh' or 'yyyymmddhh'   -> Only for ERA5. Format of date in era file name")
	
	print("\nOnly required for CUSTOM input data")
	print("custom_file_prefix     = 'custom_prefix'                 -> Custom source file name prefix")
	print("custom_date_file_name  = 'yyyymmdd_hh' or 'yyyymmddhh'   -> Format of date in era file name")
	print("custom_mslp_variable   = 'mslp'                          -> Name of MSLP variable in custom input file")
	print("custom_uwind_variable  = 'u10'                           -> Name of u-wind variable in custom input file")
	print("custom_vwind_variable  = 'v10'                           -> Name of v-wind variable in custom input file")
	print("custom_latitude_var    = 'latitude'                      -> Name of latitude variable in custom input file")
	print("custom_longitude_var   = 'longitude'                     -> Name of latitude variable in custom input file")
	print("custom_longitude_var   = 'longitude'                     -> Name of latitude variable in custom input file")
	print("custom_terrain_high_filename   = 'filename'              -> Path and name of file containing terrain high information in CUSTOM source")
	print("custom_terrain_high_var_name   = 'high'                  -> Name of Topography variable in custom_terrain_high_filename")
	print("** Note: Latitude and longitude in custom_terrain_high_filename must match latitude and longitude in custom source input data files" )
	
	print("\nOther tracking parameters")
	print("+ model_res            =  15                             -> Source data horizontal resolution in km")
	print("+ search_limits        = '[lonmin,latmin,lonmax,latmax]' -> Search limits")
	print("+ search_region        = 'NA'                            -> Search Region for EC (TC)")
	print("**EC (TC): Search region code and name")
	print("  NA (AL): North Atlantic")
	print("  SA (SA): South Atlantic")
	print("  NP (  ): North Pacific")
	print("  SP (SP): South Pacific")
	print("  SI (SI): South Indian Ocean")
	print("     (NI): North Indian Ocean")
	print("     (EP): Central and East Pacific Ocean")
	print("     (WP): South Indian Ocean")
	print("  NH (NH): North Hemisphere")
	print("  SH (SH): North Hemisphere")
	print("  GL (GL): Global")
	print("** Note: In the current version of " + program_name()+ " ("+str(get_currentversion())+")" + " it has been only tested for some (NA, MS, SA) regions\n" )
	
	
	print("\nChecking cyclone thermal structure")
	print(" checking_upper_levels_parameters    = ['yes'/'no']      -> Checking for upper level parameters ")
	print(" vtl_vtu_lr                          = ['yes'/'no']      -> Get VTL and VTU from linear regression. Only if checking_upper_levels_parameters='yes'")
	print(" max_dist                            = 500               -> Radius to compute the cyclone phase parameters. For TLCs max_dist=200")
	print(" path_data_source_upper              = path              -> Path to upper level files. If source=WRF, set path_data_source_upper like to path_data_source")
	print(" era_upperfile_prefix                = name_prefix       -> Only for ERA5. Prefix in upper ERA5 files before the date")
	print(" custom_upperfile_prefix              = name_prefix      -> Only for CUSTOM source. Prefix in upper CUSTOM files before the date")
	print(" custom_geopotential_var_name         = 'z'              -> Only for CUSTOM source. Name of gepotential variable in upper file ")
	print(" custom_upper_level_variable_name     = 'level'          -> Only for CUSTOM source. Name of level variable in upper file ")



	print("\n"+program_name() + " Output file information")
	print("..........................................")
	print("+ pathoutput            = 'path'                         -> Path to save "+program_name()+" outputs")
	print("+ tmp_dir               = 'path'                         -> Path to save temporal files of "+program_name())
	print("+ remove_tmp_dir        = 'yes' or 'no'                  -> Remove temporal files of "+program_name()+ " after runs ended\n")
	
	
	print("\nDefault parameters of "+program_name() + " for ECs")
	print("..........................................")
	print("+ filter_center_threshold      = 1000                    -> Minimum distance between two critical centers in km. Default filter_center_threshold=1000")
	print("+ critical_outer_radius        = 100                     -> Critical outer radius to considerer a low pressure point as critical center in km. Default critical_outer_radius=0")
	print("+ dist_threshold               = 1000                    -> Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=1000")
	print("+ rout                         = 2000                    -> External search radius in km. Default rout=2000 km")
	print("+ dr_res                       = 100                     -> Resolution of radial legs in km. Default dr_res=100 km")
	print("+ d_ang                        = 10                      -> resolution of angle steps for radial legs in degrees. Default d_ang=10 degrees")
	print("+ min_slp_threshold            = 1015                    -> Maximum slp treshold in hPa to filter EC centres. Deafult min_slp_threshold=1015 hPa")
	print("+ dt_lifetime                  = 48                      -> Threshold for EC lifetime in hours. Default dt_lifetime=48 hours")
	print("+ minimum_distance_travelled   = 1000                    -> Minimum distance traveled by the system in km . Default minimum_distance_travelled=1000 km")
	print("+ terrain_filter               = 1000                    -> Terrain filter in m. Set terrain_filter=0 to do not apply terrain filter. Default terrain_filter=1000 m")
	print("+ prev_days                    = 14                      -> Days before current time for computing the average mslp field. Default prev_days=14 ")
	print("+ mslp_anomaly_threshold       = -2.5                    -> MSLP anomaly threshold in hPa to considered a grid point as critical EC center")
	print("+ intensity_threshold          = 10                      -> Minimum intensity (m/s) of the EC along its track to save the track. Default intensity_threshold=0 m/s\n")
	print("+ core_criteria_length         = 0                       -> Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space" )
	print("+ VTL_threshold                = 0                       -> Lower thermal wind threshold (VTL). VTL<VTL_threshold. Default VTL_threshold=0")
	print("+ VTU_threshold                = 0                       -> Upper thermal wind threshold (VTU). VTU<VTU_threshold. Default VTU_threshold=0")
	print("+ Bhart_threshold              = 10                      -> B parameter of cyclone phase space. |B|>Bhart_threshold. Default Bhart_threshold=10")
	print("+ max_wind_speed_threshold     = 0                       -> Maximum wind speed (m/s) inside critical_outer_radius to considered a low pressure grid point as critical EC center.")
	print("                                                              Default max_wind_speed_threshold=0 m/s")
	print("+ outer_wind_speed_threshold   = 0                       -> Azimuthal wind speed threshold (m/s) to estimate the EC outer radius. Default outer_wind_speed_threshold=0 m/s")
	print("+ great_circle_distance        = 6.5                     -> Radial distance (degree) to check MSLP increase threshold. Default great_circle_distance=6.5 degree ")
	print("+ dmslp_great_circle_distance  = 200                     -> Change in MSLP (Pa) from the centro to a radial distance of great_circle_distance. Default dmslp_great_circle_distance=200 Pa")
	print("+ radius_for_msw               = 250                     -> Radius (in km) to compute de maximum wind speed. Default radius_for_msw=0 km")
	print("+ intensity_threshold          = 0                       -> Minimum intensity (m/s) of the TC along its track to save the track. Default intensity_threshold=0 m/s")
	print("+ vorticity_threshold          = 1.45e-5                 -> Relative vorticity threshold (1/s) to filter critical TC centers. Default vorticity_threshold=1.45e-5 1/s")
	

	print("\nDefault parameters of "+program_name() + " for TCs")
	print("..........................................")
	print("+ filter_center_threshold      = 400                     -> Minimum distance between two critical centers in km. Default filter_center_threshold=400")
	print("+ critical_outer_radius        = 100                     -> Critical outer radius to considerer a low pressure point as critical center in km. Default critical_outer_radius=50")
	print("+ dist_threshold               = 650                     -> Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=1000")
	print("+ rout                         = 1000                    -> External search radius in km. Default rout=1000 km")
	print("+ dr_res                       = 100                     -> Resolution for radial legs in km. Default dr_res=100 km")
	print("+ d_ang                        = 10                      -> Resolution of angle steps for radial legs in degrees. Default d_ang=10 degrees")
	print("+ min_slp_threshold            = 1015                    -> Maximum slp treshold in hPa to filter TC centres. Deafult min_slp_threshold=1015 hPa")
	print("+ dt_lifetime                  = 36                      -> Threshold for TC lifetime in hours. Default dt_lifetime=36 hours")
	print("+ minimum_distance_travelled   = 0                       -> Minimum distance traveled by the system in km . Default minimum_distance_travelled=0")
	print("+ terrain_filter               = 0                       -> Terrain filter in m. Set terrain_filter=0 to do not apply terrain filter. Default terrain_filter=0 m")
	print("+ prev_days                    = 14                      -> Days before current time for computing the average mslp field. Default prev_days=14 days")
	print("+ mslp_anomaly_threshold       = -2                      -> MSLP anomaly threshold in hPa to considered a grid point as critical TC center")
	print("+ max_wind_speed_threshold     = 8                       -> Maximum wind speed (m/s) inside critical_outer_radius to considered a low pressure grid point as critical TC center.")
	print("                                                            Default max_wind_speed_threshold=8 m/s")
	print("+ outer_wind_speed_threshold   = 6                       -> Azimuthal wind speed threshold (m/s) to estimate the TC outer radius. Default outer_wind_speed_threshold=6m/s")
	print("+ intensity_threshold          = 10                      -> Minimum intensity (m/s) of the TC along its track to save the track. Default intensity_threshold=10 m/s")
	print("+ vorticity_threshold          = 1.45e-5                 -> Relative vorticity threshold (1/s) to filter critical TC centers. Default vorticity_threshold=1.45e-5 1/s")
	print("+ great_circle_distance        = 5.5                     -> Radial distance (degree) to check MSLP increase threshold. Default great_circle_distance=5.5 degree ")
	print("+ dmslp_great_circle_distance  = 200                     -> Change in MSLP (Pa) from the centro to a radial distance of great_circle_distance. Default dmslp_great_circle_distance=200 Pa")
	print("+ radius_for_msw               = 100                     -> Radius (in km) to compute de maximum wind speed. Default radius_for_msw=100 km")
	print("+ core_criteria_length         = 3                       -> Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space" )
	print("+ VTL_threshold                = 0                       -> Lower thermal wind threshold (VTL). VTL<VTL_threshold. Default VTL_threshold=0")
	print("+ VTU_threshold                = 0                       -> Upper thermal wind threshold (VTU). VTU<VTU_threshold. Default VTU_threshold=0")
	print("+ Bhart_threshold              = 10                      -> B parameter of cyclone phase space. |B|>Bhart_threshold. Default Bhart_threshold=10")
	
	
	
	print("\nDefault parameters of "+program_name() + " for SCs")
	print("..........................................")
	print("+ filter_center_threshold      = 400                     -> Minimum distance between two critical centers in km. Default filter_center_threshold=400")
	print("+ critical_outer_radius        = 0                       -> Critical outer radius to considerer a low pressure point as critical center in km. Default critical_outer_radius=0")
	print("+ dist_threshold               = 400                     -> Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=1000")
	print("+ rout                         = 1000                    -> External search radius in km. Default rout=1000 km")
	print("+ dr_res                       = 100                     -> Resolution for radial legs in km. Default dr_res=100 km")
	print("+ d_ang                        = 10                      -> resolution of angle steps for radial legs in degrees. Default d_ang=10 degrees")
	print("+ min_slp_threshold            = 1015                    -> Maximum slp treshold in hPa to filter SC centres. Deafult min_slp_threshold=1015 hPa")
	print("+ dt_lifetime                  = 36                      -> Threshold for SC lifetime in hours. Default dt_lifetime=24 hours")
	print("+ minimum_distance_travelled   = 0                       -> Minimum distance traveled by the system in km . Default minimum_distance_travelled=1000 km")
	print("+ terrain_filter               = 0                       -> Terrain filter in m. Set terrain_filter=0 to do not apply terrain filter. Default terrain_filter=0 m")
	print("+ prev_days                    = 14                      -> Days before current time for computing the average mslp field. Default prev_days=14 days")
	print("+ mslp_anomaly_threshold       = -2.5                    -> MSLP anomaly threshold in hPa to considered a grid point as critical SC center")
	print("+ max_wind_speed_threshold     = 0                       -> Maximum wind speed (m/s) inside critical_outer_radius to considered a low pressure grid point as critical SC center.")
	print("                                                            Default max_wind_speed_threshold=0 m/s")
	print("+ outer_wind_speed_threshold   = 2                       -> Azimuthal wind speed threshold (m/s) to estimate the SC outer radius. Default outer_wind_speed_threshold=2m/s")
	print("+ intensity_threshold          = 0                       -> Minimum intensity (m/s) of the SC along its track to save the track. Default intensity_threshold=0 m/s")
	print("+ vorticity_threshold          = 1.5e-5                  -> Relative vorticity threshold (1/s) to filter critical SC centers. Default vorticity_threshold=1.5e-5 1/s")
	print("+ great_circle_distance        = 5.5                     -> Radial distance (degree) to check MSLP increase threshold. Default great_circle_distance=5.5 degree ")
	print("+ dmslp_great_circle_distance  = 200                     -> Change in MSLP (Pa) from the centro to a radial distance of great_circle_distance. Default dmslp_great_circle_distance=0 Pa")
	print("+ radius_for_msw               = 100                    -> Radius (in km) to compute de maximum wind speed. Default radius_for_msw=100 km")
	print("+ core_criteria_length         = 7                       -> Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space" )
	print("+ VTL_threshold                = -50                     -> Lower thermal wind threshold (VTL). VTL>VTL_threshold. Default VTL_threshold=-50")
	print("+ VTU_threshold                = -10                     -> Upper thermal wind threshold (VTU). VTU<VTU_threshold. Default VTU_threshold=-10")
	print("+ Bhart_threshold              = 25                      -> B parameter of cyclone phase space. |B|>Bhart_threshold. Default Bhart_threshold=20")
	
	
	print("\nDefault parameters of "+program_name() + " for MCs")
	print("..........................................")
	print("+ filter_center_threshold      = 300                     -> Minimum distance between two critical centers in km. Default filter_center_threshold=300")
	print("+ critical_outer_radius        = 50                      -> Critical outer radius to considerer a low pressure point as critical center in km. Default critical_outer_radius=50")
	print("+ dist_threshold               = 400                     -> Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=400")
	print("+ rout                         = 800                     -> External search radius in km. Default rout=800 km")
	print("+ dr_res                       = 100                     -> Resolution for radial legs in km. Default dr_res=100 km")
	print("+ d_ang                        = 10                      -> Resolution of angle steps for radial legs in degrees. Default d_ang=10 degrees")
	print("+ min_slp_threshold            = 1015                    -> Maximum slp treshold in hPa to filter MC centres. Deafult min_slp_threshold=1015 hPa")
	print("+ dt_lifetime                  = 24                      -> Threshold for MC lifetime in hours. Default dt_lifetime=24 hours")
	print("+ minimum_distance_travelled   = 200                     -> Minimum distance traveled by the system in km . Default minimum_distance_travelled=200 km")
	print("+ terrain_filter               = 0                       -> Terrain filter in m. Set terrain_filter=0 to do not apply terrain filter. Default terrain_filter=0 m")
	print("+ prev_days                    = 14                      -> Days before current time for computing the average mslp field. Default prev_days=14 days")
	print("+ mslp_anomaly_threshold       = -2.5                    -> MSLP anomaly threshold in hPa to considered a grid point as critical MC center")
	print("+ max_wind_speed_threshold     = 0                       -> Maximum wind speed (m/s) inside critical_outer_radius to considered a low pressure grid point as critical MC center.")
	print("                                                            Default max_wind_speed_threshold=0 m/s")
	print("+ outer_wind_speed_threshold   = 2                       -> Azimuthal wind speed threshold (m/s) to estimate the MC outer radius. Default outer_wind_speed_threshold=2m/s")
	print("+ intensity_threshold          = 0                       -> Minimum intensity (m/s) of the MC along its track to save the track. Default intensity_threshold=0 m/s")
	print("+ vorticity_threshold          = 1.45e-5                 -> Relative vorticity threshold (1/s) to filter critical MC centers. Default vorticity_threshold=1.45e-5 1/s")
	print("+ great_circle_distance        = 3                       -> Radial distance (degree) to check MSLP increase threshold. Default great_circle_distance=3 degree ")
	print("+ dmslp_great_circle_distance  = 200 Pa                  -> Change in MSLP (Pa) from the centro to a radial distance of great_circle_distance. Default dmslp_great_circle_distance=200 Pa")
	print("+ radius_for_msw               = 100                     -> Radius (in km) to compute de maximum wind speed. Default radius_for_msw=100 km")
	print("+ core_criteria_length         = 0                       -> Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space" )
	print("+ VTL_threshold                = 0                       -> Lower thermal wind threshold (VTL). VTL<VTL_threshold. Default VTL_threshold=0")
	print("+ VTU_threshold                = 0                       -> Upper thermal wind threshold (VTU). VTU<VTU_threshold. Default VTU_threshold=0")
	print("+ Bhart_threshold              = 10                      -> B parameter of cyclone phase space. |B|>Bhart_threshold. Default Bhart_threshold=10")
	
	
	print("\nDefault parameters of "+program_name() + " for TLCs")
	print("..........................................")
	print("+ filter_center_threshold      = 300                     -> Minimum distance between two critical centers in km. Default filter_center_threshold=300")
	print("+ critical_outer_radius        = 50                      -> Critical outer radius to considerer a low pressure point as critical center in km. Default critical_outer_radius=50")
	print("+ dist_threshold               = 400                     -> Maximum distance between centres (in km) in continuos time steps. Default dist_threshold=400")
	print("+ rout                         = 800                     -> External search radius in km. Default rout=800 km")
	print("+ dr_res                       = 100                     -> Resolution for radial legs in km. Default dr_res=100 km")
	print("+ d_ang                        = 10                      -> Resolution of angle steps for radial legs in degrees. Default d_ang=10 degrees")
	print("+ min_slp_threshold            = 1015                    -> Maximum slp treshold in hPa to filter TLC centres. Deafult min_slp_threshold=1015 hPa")
	print("+ dt_lifetime                  = 24                      -> Threshold for TLC lifetime in hours. Default dt_lifetime=24 hours")
	print("+ minimum_distance_travelled   = 200                     -> Minimum distance traveled by the system in km . Default minimum_distance_travelled=1000 km")
	print("+ terrain_filter               = 0                       -> Terrain filter in m. Set terrain_filter=0 to do not apply terrain filter. Default terrain_filter=0 m")
	print("+ prev_days                    = 14                      -> Days before current time for computing the average mslp field. Default prev_days=14 days")
	print("+ mslp_anomaly_threshold       = -2.5                    -> MSLP anomaly threshold in hPa to considered a grid point as critical TLC center")
	print("+ max_wind_speed_threshold     = 8                       -> Maximum wind speed (m/s) inside critical_outer_radius to considered a low pressure grid point as critical TLC center.")
	print("                                                            Default max_wind_speed_threshold=8 m/s")
	print("+ outer_wind_speed_threshold   = 2                       -> Azimuthal wind speed threshold (m/s) to estimate the TLC outer radius. Default outer_wind_speed_threshold=2m/s")
	print("+ intensity_threshold          = 10                      -> Minimum intensity (m/s) of the TLC along its track to save the track. Default intensity_threshold=10 m/s")
	print("+ vorticity_threshold          = 1.45e-5                 -> Relative vorticity threshold (1/s) to filter critical TLC centers. Default vorticity_threshold=1.45e-5 1/s")
	print("+ great_circle_distance        = 3                       -> Radial distance (degree) to check MSLP increase threshold. Default great_circle_distance=3 degree ")
	print("+ dmslp_great_circle_distance  = 200                     -> Change in MSLP (Pa) from the centro to a radial distance of great_circle_distance. Default dmslp_great_circle_distance=200 Pa")
	print("+ radius_for_msw               = 100                     -> Radius (in km) to compute de maximum wind speed. Default radius_for_msw=100 km")
	print("+ core_criteria_length         = 1                       -> Minimum time (time steps) in which the detected cyclone satisfies the thermal structure determined by the cyclone phase space" )
	print("+ VTL_threshold                = 0                       -> Lower thermal wind threshold (VTL). VTL<VTL_threshold. Default VTL_threshold=0")
	print("+ VTU_threshold                = 0                       -> Upper thermal wind threshold (VTU). VTU<VTU_threshold. Default VTU_threshold=0")
	print("+ Bhart_threshold              = 10                      -> B parameter of cyclone phase space. |B|>Bhart_threshold. Default Bhart_threshold=10")
	
	print("\n============================================================================================================")
	print("Bye :)")
	print("============================================================================================================")



def str2boolean(arg):
    if isinstance(arg, bool):
        return arg
    if arg.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif arg.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def print_error_message(message):
	print("=============================================================================================================")
	print ("ERROR: "+message)
	print("=============================================================================================================")
	raise SystemExit("Bye :)")

def checking_optimum_nproc(rank, n_proc, length_files):
	if n_proc > length_files:
		if rank==0:
			print_error_message("You are using " + str(int(n_proc)) + " processors. It exceeds the number of files to be processed.\n Use "+str(int(length_files))+" procesors")
		raise SystemExit("Bye :)")

def check_paths(pfile, path):
	try: 
		fpath = getattr(pfile, path)
	except:
		fpath = ""
	return fpath

def get_dates(date_init, date_end,dt_h):
	dates=[]
	delta=date_end-date_init
	for i in range(0,int(delta.total_seconds()/3600) + dt_h,dt_h):
		dates = np.append(dates,str(date_init + timedelta(hours=i)))
	return dates






def get_dates_vectors(year_case_init="",
		month_case_init="",
		day_case_init="",
		hour_case_init="",
		year_case_end="", 
		month_case_end="",
		day_case_end="",
		hour_case_end="",
		dt_h=None,
		prev_days=14,
		previous_dates=False,
		calendar="366d"
		):
	
	
	prev_days=np.abs(prev_days)
	
	if (dt_h==None or dt_h==0 or dt_h>6) and previous_dates==False:
		print_error_message("The time step is incorrect or zero. dt_h must be equal or lower than 6.\n"+program_name()+" Exit")

	
	if (year_case_init=="" or month_case_init=="" or day_case_init=="" or hour_case_init==""):
		print_error_message("Start date is incorrect in the configuration file\n"+program_name()+" Exit")

	else:
		if int(month_case_init)==1 or int(month_case_init)==3 or int(month_case_init)==5 or int(month_case_init)==7 or int(month_case_init)==8 or int(month_case_init)==10 or int(month_case_init)==12:
			if int(day_case_init)>31:
				print_error_message("The begin_day: " + day_case_init +" is incorrect in the configuration file")
		elif int(month_case_init)==4 or int(month_case_init)==6 or int(month_case_init)==9 or int(month_case_init)==11:
			if int(day_case_init)>30:
				print_error_message("The begin_day: " + day_case_init +" is incorrect in the configuration file")
		elif int(year_case_init)%4==0 and int(month_case_init)==2 and int(day_case_init)>29:
			print_error_message("The begin_day: " + day_case_init +" is incorrect in the configuration file")
		elif int(year_case_init)%4!=0 and int(month_case_init)==2 and int(day_case_init)>28:
			print_error_message("The begin_day: " + day_case_init +" is incorrect in the configuration file")

	if year_case_end=="" or month_case_end=="" or day_case_end=="" or hour_case_end=="":
		print_error_message("The end date info is incorrect in the configuration file")
	else:
		if int(month_case_end)==1 or int(month_case_end)==3 or int(month_case_end)==5 or int(month_case_end)==7 or int(month_case_end)==8 or int(month_case_end)==10 or int(month_case_end)==12:
			if int(day_case_end)>31:
				print_error_message("The end_date: " + day_case_end +" is incorrect in the configuration file")
		elif int(month_case_end)==4 or int(month_case_end)==6 or int(month_case_end)==9 or int(month_case_end)==11:
			if int(day_case_end)>30:
				print_error_message("The end_date: " + day_case_end +" is incorrect in the configuration file")

		elif int(year_case_end)%4==0 and int(month_case_end)==2 and int(day_case_end)>29:
			print_error_message("The end_date: " + day_case_end +" is incorrect in the configuration file")

		elif int(year_case_end)%4!=0 and int(month_case_end)==2 and int(day_case_end)>28:
			print_error_message("The end_date: " + day_case_end +" is incorrect in the configuration file")

	if int(month_case_init)<1 or int(month_case_init)>12 or int(month_case_end)<1 or int(month_case_end)>12:
		print_error_message("The month of the case dates is incorrect in the configuration file")

	if int(year_case_init)>int(year_case_end):
		print_error_message("Start year: "+ year_case_init +" is higher than end year: " + year_case_end)

	elif int(year_case_init)==int(year_case_end) and int(month_case_init)>int(month_case_end):
		print_error_message("It is the same year, but begin_month: "+ month_case_init +" is higher than month end: " + month_case_end)

	elif int(year_case_init)==int(year_case_end) and int(month_case_init)==int(month_case_end) and int(day_case_init)>int(day_case_end):
		print_error_message("It is the same year and same month, but begin_day: "+ day_case_init +" is higher than day end: " + day_case_end)

	elif int(year_case_init)==int(year_case_end) and int(month_case_init)==int(month_case_end) and int(day_case_init)== int(day_case_end) and int(hour_case_init)>int(hour_case_end):
		print_error_message("It is the same date, but begin_hour: "+ hour_case_init +" is higher than end_hour: " + hour_case_end)

	month_case_init=month_case_init.zfill(2)
	month_case_end=month_case_end.zfill(2)
	day_case_init=day_case_init.zfill(2)
	day_case_end=day_case_end.zfill(2)
	hour_case_init=hour_case_init.zfill(2)
	hour_case_end=hour_case_end.zfill(2)


	formatted_time = datetime.strptime(year_case_init+"-"+month_case_init+"-"+day_case_init+" "+hour_case_init+":00:00", "%Y-%m-%d %H:%M:%S")

	if calendar=="365d" and int(year_case_init)%4==0 and month_case_init=="03" and  int(day_case_init) in np.arange(1, prev_days+1,1):
		begin_time=str(formatted_time+timedelta(hours=int((prev_days+1)*(-24))))
		
	else:	
		begin_time=str(formatted_time+timedelta(hours=int(prev_days*(-24))))


	prev_year_init=begin_time.split(" ")[0].split("-")[0]
	prev_month_init=begin_time.split(" ")[0].split("-")[1]
	prev_day_init=begin_time.split(" ")[0].split("-")[2]
	prev_hour_init=begin_time.split(" ")[1].split(":")[0]
	
	prev_year_init=prev_year_init.zfill(4)
	prev_month_init=prev_month_init.zfill(2)
	prev_day_init=prev_day_init.zfill(2)
	prev_hour_init=prev_hour_init.zfill(2)


	if year_case_init+month_case_init+day_case_init+hour_case_init==year_case_end+month_case_end+day_case_end+hour_case_end and previous_dates==False:
		print_error_message(" Start date and End Date are the same \n"+program_name()+" Exit")
	else:
		date_init=datetime(int(prev_year_init),int(prev_month_init),int(prev_day_init),int(prev_hour_init),0,0)
		date_end=datetime(int(year_case_end),int(month_case_end),int(day_case_end),int(hour_case_end),0,0)
		dates_=get_dates(date_init, date_end,dt_h)
		dates=[]
		hours=[]
		for date_ in dates_:
			yyyymmdd=date_.split(" ")[0].split("-")[0]+date_.split(" ")[0].split("-")[1]+date_.split(" ")[0].split("-")[2]
			hh=date_.split(" ")[1].split(":")[0]
			
			if calendar=="365d" and  int( yyyymmdd[0:4])%4==0 and   yyyymmdd[4:8]=="0229":
				log_message="nothing to do"
			else:
			
				dates=np.append(dates,yyyymmdd)
				hours=np.append(hours,hh)

	return dates,hours


def ProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '+', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = "\r")
    if iteration == total: 
        print()


def get_era5_files(dates=[""],hours=[""],era_file_prefix="",era_date_file_name="yyyymmdd_hh"):
	erafiles=[]
	for i in range(0,len(dates)):
		year=dates[i][0:4].zfill(4)
		month=dates[i][4:6].zfill(2)
		day=dates[i][6:8].zfill(2)
		
		hour=hours[i].zfill(2)
		if era_date_file_name=="yyyymmdd_hh":
			nfame=year+month+day+"_"+hour
		elif era_date_file_name=="yyyymmdd_hh":
			nfame=year+month+day+hour
		else:
			print_error_message("Unrecognized ERA file era_date_file_name attribute in CyTRACK input paramters file\nera_date_file_name must be yyyymmdd_hh or yyyymmddhh\nRun << python run_CyTRACK.py - cth t >> for help")
		
		
		fname=era_file_prefix+"_"+nfame+".nc"
		erafiles=np.append(erafiles,fname)
	return erafiles



def get_custom_files(dates=[""],hours=[""],custom_file_prefix="",custom_date_file_name="yyyymmdd_hh"):
	customfiles=[]
	for i in range(0,len(dates)):
		year=dates[i][0:4].zfill(4)
		month=dates[i][4:6].zfill(2)
		day=dates[i][6:8].zfill(2)
		
		hour=hours[i].zfill(2)
		if custom_date_file_name=="yyyymmdd_hh":
			nfame=year+month+day+"_"+hour
		elif custom__date_file_name=="yyyymmdd_hh":
			nfame=year+month+day+hour
		else:
			print_error_message("Unrecognized CUSTOM file custom_date_file_name attribute in CyTRACK input paramters file\ncustom_date_file_name must be yyyymmdd_hh or yyyymmddhh\nRun << python run_CyTRACK.py - cth t >> for help")
		
		
		fname=custom_file_prefix+"_"+nfame+".nc"
		customfiles=np.append(customfiles,fname)
	return customfiles



def download_era5(erafile,year,month,day,hour):
	import cdsapi
	print(year, month, day, hour)
	c = cdsapi.Client()
	c.retrieve(
    		'reanalysis-era5-single-levels',
		{
        	'product_type': 'reanalysis',
		'variable':['10m_u_component_of_wind', '10m_v_component_of_wind', 'mean_sea_level_pressure',],
		'year': year,
		'month': month,
		'day': day,
		'time':hour+":00",
		'format':"netcdf",
		},
		erafile)


def download_era5_upper(erafile_upper,year,month,day,hour):
	import cdsapi
	print(year, month, day, hour)
	c = cdsapi.Client()
	c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
		'product_type': 'reanalysis',
		'variable':['geopotential'],
		'pressure_level': ['200', '250','300', '350','400', '450','500', '550', '600', '650','700', '750', '800', '850', '900',],
		'year': year,
		'month': month,
		'day': day,
		'time':hour+":00",
		'format':"netcdf",
		},
		erafile_upper)




def get_wrf_files(dates=[""],hours=[""], wrfprefix="d01"):
	wrfiles=[]
	for i in range(0,len(dates)):
		fname=wrfprefix+"_"+dates[i][0:4]+"-"+dates[i][4:6]+"-"+dates[i][6:8]+"_"+hours[i]+":00:00"
		wrfiles=np.append(wrfiles,fname)
	return wrfiles

def checking_input_files(pathfile, source_file, source,date,hour,flev='sfc'):
	if os.path.exists(pathfile+"/"+source_file):
		pass
	elif source.upper()=="WRF":
		print_error_message( pathfile+"/"+source_file+ " is not in the directory")
	
	elif source.upper()=="ERA5":
		print("--------------------------------------------------------------------------------------")
		print ("WARNING: ERA5 file: ", pathfile+"/"+source_file, " not found")
		print ("Trying to download it from " + source.upper() + " reanalysis")
		print("--------------------------------------------------------------------------------------")
		path=pathfile
		year=date[0:4].zfill(4)
		month=date[4:6].zfill(2)
		day=date[6:8].zfill(2)
		if flev=="sfc":
			download_era5(erafile=pathfile+"/"+source_file,year=year,month=month,day=day,hour=hour)
		elif flev=="upper":
			download_era5_upper(erafile_upper=pathfile+"/"+source_file,year=year,month=month,day=day,hour=hour)
	elif source.upper()=="CUSTOM":
		print_error_message( pathfile+"/"+source_file+ " is not in the directory")
	return True

def get_wrf_hgt(idir="",wrffile=""):
	ncwrffile=Dataset(idir+"/"+wrffile)
	hgt=ncwrffile.variables["HGT"][:]
	lat=ncwrffile.variables["XLAT"][0,:]
	lon=ncwrffile.variables["XLONG"][0,:]
	if len(hgt.shape)>2:
		hgt=hgt[0,:]
	return lat,lon,hgt


def get_i_bg(prev_days):
	if prev_days>0:
		i_bg=int(prev_days)*4
	else:
		i_bg=0
	return i_bg			


def get_mslp_anomaly(idir="./",
			source="",
			search_limits=[None,None,None,None],
			search_region="",
			prev_days=14,
			date="",
			hour="",
			source_filename_prefix="",
			source_file_date_format="",
			custom_mslp_variable="",
			custom_latitude_var="",
			custom_longitude_var=""):


	dates,hours=get_dates_vectors(year_case_init=date[0:4],
					month_case_init=date[4:6].zfill(2),
					day_case_init=date[6:8].zfill(2),
					hour_case_init=hour,
					year_case_end=date[0:4],
					month_case_end=date[4:6].zfill(2),
					day_case_end=date[6:8].zfill(2),
					hour_case_end=hour,
					dt_h=24,
					prev_days=prev_days,
					previous_dates=True
					 )

	dates=dates[:-1]
	hours=hours[:-1]

	if source.upper()=="ERA5":
		previous_files=get_era5_files(dates=dates,hours=hours,era_file_prefix=source_filename_prefix,era_date_file_name=source_file_date_format)

		mean_slp=0
		for index in range(0,len(previous_files)):
			varlist=get_era5_2dvar(idir=idir,erafile=previous_files[index],svariables=["msl"],search_limits=search_limits,search_region=search_region)
			eramslp=varlist[2,:]/100
			mean_slp=mean_slp+eramslp
			
	if source.upper()=="CUSTOM":
		previous_files=get_custom_files(dates=dates,hours=hours,custom_file_prefix=source_filename_prefix,custom_date_file_name=source_file_date_format)
		mean_slp=0
		for index in range(0,len(previous_files)):
			varlist=get_custom_2dvar(idir=idir,customfile=previous_files[index],svariables=[custom_mslp_variable],search_limits=search_limits,search_region=search_region, custom_latitude_var=custom_latitude_var, custom_longitude_var=custom_longitude_var)
			sourcemslp=varlist[2,:]
			
			if len(str(int(sourcemslp.max())))>=5:
				sourcemslp=sourcemslp/100
		
			mean_slp=mean_slp+sourcemslp
			
			
	if source.upper()=="WRF":
		previous_files=get_wrf_files(dates=dates,hours=hours,wrfprefix=source_filename_prefix)
		mean_slp=0
		for index in range(0,len(previous_files)):
			#wrflat,wrflon,wrfmslp=get_wrf_mslp(idir=idir,wrffile=previous_files[index],variables=["PSFC","T2","PHB","PH"])
			wrflat,wrflon,wrfmslp=get_wrf_mslp_new(idir=idir,wrffile=previous_files[index],variables=["PB","P","PHB","PH","T","QVAPOR"])
			
			
			
			mean_slp=mean_slp+wrfmslp


	mean_slp=mean_slp/len(previous_files)

	return mean_slp


def tracker_EC(source="",
		idir="./",
		sourcefiles="",
		pathoutput="./",
		rout=2000,
		verbose=True,
		dates=[""],
		hours=[""],
		model_res=20,
		search_limits=[None,None,None,None],
		dr_res=100,
		d_ang=10,
		filter_center_threshold=1500,
		critical_outer_radius=200,
		tmpdir="."+program_name()+"_tmpdir/",
		rank=0,
		search_region="",
		min_slp_threshold=1015,
		terrain_filter=1000,
		prev_days=14,
		mslp_anomaly_threshold=-3,
		source_filename_prefix="",
		source_file_date_format=""):

	if source.upper()=="ERA5":
		hgt_file=get_era5_hgtfile()
		varhgt=get_era5_2dvar(idir="",erafile=hgt_file,svariables=["z"],search_limits=search_limits,search_region=search_region)
		hgt_field=varhgt[2,:]/9.80665

		varu="u10"
		varv="v10"
		varlat="latitude"
		varlon="longitude"
		dx,dy=compute_dx_dy(lons=varhgt[1,:],lats=varhgt[0,:])
	elif source.upper()=="WRF":
		hgtlat,hgtlon,hgt_field=get_wrf_hgt(idir=idir,wrffile=sourcefiles[0])

		varu="U10"
		varv="V10"
		varlat="XLAT"
		varlon="XLONG"

		dx,dy=compute_dx_dy(lons=hgtlon,lats=hgtlat)

	if terrain_filter>0:
		hgt_field=hgt_field
	elif terrain_filter==0:
		hgt_field[:,:]=-1


	for index in range(0,len(sourcefiles)):
		fdate=dates[index]+hours[index]
		if verbose:
			print("   ---> | Processing "+source.upper()+ ": "+idir+"/"+sourcefiles[index])
		if source.upper()=="ERA5":
			varlist=get_era5_2dvar(idir=idir,erafile=sourcefiles[index],svariables=["msl"],search_limits=search_limits,search_region=search_region)
			sourcelat=varlist[0,:]
			sourcelon=varlist[1,:]
			sourcemslp=varlist[2,:]/100
		elif source.upper()=="WRF":
			#sourcelat,sourcelon,sourcemslp=get_wrf_mslp(idir=idir,wrffile=sourcefiles[index],variables=["PSFC","T2","PHB","PH"])
			sourcelat,sourcelon,sourcemslp=get_wrf_mslp_new(idir=idir,wrffile=sourcefiles[index],variables=["PB","P","PHB","PH","T","QVAPOR"])
			



			
		avg_mslp=get_mslp_anomaly(idir=idir,
						source=source,
						search_limits=search_limits,
						search_region=search_region,
						prev_days=prev_days,
						date=dates[index],
						hour=hours[index],
						source_filename_prefix=source_filename_prefix,
						source_file_date_format=source_file_date_format)

		mslp_anomaly=sourcemslp-avg_mslp

		clats,clons,croci,cclosedp,cpmin=get_low_EC_centers(lats=sourcelat,
						lons=sourcelon,
						mslp=sourcemslp,
						mslp_anomaly=mslp_anomaly,
						mslp_anomaly_threshold=mslp_anomaly_threshold,
						model_res=model_res,
						search_limits=search_limits,
						rout=rout,
						dr_res=dr_res,
						d_ang=d_ang,
						critical_outer_radius=critical_outer_radius,
						min_slp_threshold=min_slp_threshold,
						terrain_filter=terrain_filter,
						hgt_field=hgt_field)
		
		
		
		
		flats, flons,froci,fpmin, fclosedp,fmws,fouter_r=filter_centers(lats=np.array(clats),
								lons=np.array(clons),
								roci=np.array(croci),
								pmin=np.array(cpmin),
								closedp=np.array(cclosedp),
								filter_center_threshold=filter_center_threshold)


		fwind_speed=get_wind_speed(latsc=flats,
							lonsc=flons,
							radius=froci,
							wfile=sourcefiles[index],
							idir=idir,
							varu=varu,
							varv=varv,
							varlat=varlat,
							varlon=varlon,
							source=source.upper(),
							search_limits=search_limits,
							search_region=search_region,
							r_uv=False)
		
		centers_data=np.empty((len(flats),8))
		centers_data[:,0]=flats
		centers_data[:,1]=flons
		centers_data[:,2]=fpmin
		centers_data[:,3]=fwind_speed
		centers_data[:,4]=fclosedp
		centers_data[:,5]= froci
		centers_data[:,6]= 0
		centers_data[:,7]=1
		
		np.savetxt(tmpdir+"/critical_centers_"+fdate+".dat",centers_data)
		

def get_era5_2dvar(idir="./",
		erafile="",
		svariables=[""],
		search_limits=[None,None,None,None],
		search_region=""):
	
	
	ncera=Dataset(idir+"/"+erafile)
	eralat=ncera.variables["latitude"][:]
	eralon=ncera.variables["longitude"][:]




	for var in svariables:
		vars()[var]=ncera.variables[var][:]
		if len(vars()[var].shape)>2:
			vars()[var]=vars()[var][0,:]
	
	if search_region.upper() in ("NA","SA","AL","MS"):
		for i in range(0,len(svariables)):
			vars()[svariables[i]],nlon=convert_era5_matrix(vars()[svariables[i]],eralon,search_lon=180)

	else:
		nlon=np.copy(eralon)



	if search_region.upper() in ("NA","SA","AL","MS"):
		for i in range(0,len(nlon)):
			if nlon[i]>=180:
				nlon[i]=nlon[i]-360
	
	if len(nlon.shape):
		nlon,eralat=np.meshgrid(nlon,eralat)

	for i in range(0,len(svariables)):
		nera_lat,nera_lon, vars()["n"+svariables[i]]= era_subregion(lat=eralat,lon=nlon,var=vars()[svariables[i]],search_limits=search_limits)

	
	varlist=np.empty((len(svariables)+2,nera_lat.shape[0],nera_lat.shape[1]))
	varlist[0,:]=nera_lat
	varlist[1,:]=nera_lon
	ncera.close()	
	for i in range(0,len(svariables)):
		varlist[i+2,:]=vars()["n"+svariables[i]]
	return varlist 



def get_custom_2dvar(idir="./",
		customfile="",
		svariables=[""],
		search_limits=[None,None,None,None],
		search_region="",
		custom_latitude_var="latitude", 
		custom_longitude_var="longitude"):
	
	
	nc=Dataset(idir+"/"+customfile)
	customlat=nc.variables[custom_latitude_var][:]
	customlon=nc.variables[custom_longitude_var][:]


	for var in svariables:
		vars()[var]=nc.variables[var][:]
		if len(vars()[var].shape)>2:
			vars()[var]=vars()[var][0,:]
	
	
	#if len(customlon.shape)>1:
		#checklon=customlon[0,:]
	#else:
		#checklon=np.copy(customlon)
	
	#if search_region.upper() in ("NA","SA","AL","MS"):
		#for i in range(0,len(svariables)):
			#vars()[svariables[i]],nlon=convert_era5_matrix(vars()[svariables[i]],checklon,search_lon=180)

	#else:
		#nlon=np.copy(customlon)
	
	
	if len(customlon.shape):
		nlon,customlat=np.meshgrid(customlon,customlat)
		
		
	#if search_region.upper() in ("NA","SA","AL","MS"):
	for i in range(0,nlon.shape[0]):
		for j in range(0, nlon.shape[1]):
			if nlon[i,j]>=180:
				nlon[i,j]=nlon[i,j]-360	

	#for i in range(0,len(svariables)):
		#ncust_lat,ncust_lon, vars()["n"+svariables[i]]= era_subregion(lat=customlat,lon=nlon,var=vars()[svariables[i]],search_limits=search_limits)

	
	varlist=np.empty((len(svariables)+2,customlat.shape[0],customlat.shape[1]))
	varlist[0,:]=customlat
	varlist[1,:]=nlon
	nc.close()	
	for i in range(0,len(svariables)):
		varlist[i+2,:]=vars()[svariables[i]]
	return varlist 




def get_cumstom_hgt_data(filename, varname):
	cnc=Dataset(filename)
	
	hgt=cnc.variables[varname][:]
	
	if len(hgt.shape)>2:
		hgt=hgt[0,:]


	return hgt




def get_era5_3dvar(idir="./",
		erafile="",
		svariable="",
		varlevel="level",
		search_limits=[None,None,None,None],
		search_region="",
		fdate="",
		full=False
		):
	
	ncera=Dataset(idir+"/"+erafile)
	check_vars=[svariable, varlevel]
	variablenot=False
	for var in check_vars:
		if not var in ncera.variables.keys():
			print("		-WARNING: Variable " + svariable + " in " + idir+"/"+erafile+ " is missing")
			variablenot=True
	
	if  variablenot:
		print("		--------------------------------------------------------------------------------------")
		print ("		Trying to download it from EAR5 reanalysis")
		print("		--------------------------------------------------------------------------------------")
		year=fdate[0:4].zfill(4)
		month=fdate[4:6].zfill(2)
		day=fdate[6:8].zfill(2)
		hour=fdate[8:10].zfill(2)
		download_era5_upper(erafile_upper=idir+"/"+erafile,year=year,month=month,day=day,hour=hour)

	ncerau=Dataset(idir+"/"+erafile)
	eralat_=ncerau.variables["latitude"][:]
	eralon=ncerau.variables["longitude"][:]
	levels=ncerau.variables[varlevel][:]
	evar=ncerau.variables[svariable][:]/9.80665
	ncera.close()
	


	
	lonn,latt=np.meshgrid(eralon,eralat_)
	
	
	if len(evar.shape)>3:
		evar=evar[0,:]
	nevar=np.empty_like(evar)
	nevar[:,:,:]=0
	if search_region.upper() in ("NA","SA","AL","MS"):
		for i in range(0,evar.shape[0]):
			auxvar,nlon=convert_era5_matrix(evar[i,:],eralon,search_lon=180)
			nevar[i,:]=auxvar
	else:
		nlon=np.copy(eralon)

	if search_region.upper() in ("NA","SA","AL","MS"):
		for i in range(0,len(nlon)):
			if nlon[i]>=180:
				nlon[i]=nlon[i]-360

	
	if len(nlon.shape):
		nlon,eralat=np.meshgrid(nlon,eralat_)

	nera_lat,nera_lon, sevar= era_subregion(lat=eralat,lon=nlon,var=nevar[0,:],search_limits=search_limits)

	sub_evar=np.empty((len(levels),sevar.shape[0],sevar.shape[1]))
	sub_evar[:,:,:]=0

	for index in range(0,nevar.shape[0]):
		nera_lat,nera_lon, svar_aux = era_subregion(lat=eralat,lon=nlon,var=nevar[index,:],search_limits=search_limits)
		sub_evar[index,:]=svar_aux
	
	if full:
		return evar, latt, lonn, levels
	else:


		return sub_evar, nera_lat,nera_lon, levels



def get_custom_3dvar(idir="./",
		customfile="",
		svariable="",
		varlevel="level",
		search_limits=[None,None,None,None],
		search_region="",
		fdate="",
		full=False,
		varlat="", 
		varlon=""
		):
	
	nc=Dataset(idir+"/"+customfile)
	check_vars=[svariable, varlevel]
	variablenot=False
	for var in check_vars:
		if not var in nc.variables.keys():
			print("		-WARNING: Variable " + svariable + " in " + idir+"/"+customfile+ " is missing")
			variablenot=True
	
	ncerau=Dataset(idir+"/"+customfile)
	eralat_=ncerau.variables[varlat][:]
	eralon=ncerau.variables[varlon][:]
	levels=ncerau.variables[varlevel][:]
	
	evar=ncerau.variables[svariable][:]/9.80665
	nc.close()
	

	if len(eralat_.shape)<2:
	
		lonn,latt=np.meshgrid(eralon,eralat_)
	
	
	#if search_region.upper() in ("NA","SA","AL","MS"):
	for i in range(0,lonn.shape[0]):
		for j in range(0,lonn.shape[1]):
			if lonn[i,j]>=180:
				lonn[i,j]=lonn[i,j]-360
	
	
	
	
	if len(evar.shape)>3:
		evar=evar[0,:]
	#nevar=np.empty_like(evar)
	#nevar[:,:,:]=0
	#if search_region.upper() in ("NA","SA","AL","MS"):
		#for i in range(0,evar.shape[0]):
			#auxvar,nlon=convert_era5_matrix(evar[i,:],lonn[0,:],search_lon=180)
			#nevar[i,:]=auxvar
	#else:
		#nlon=np.copy(lonn[0,:])

	#if search_region.upper() in ("NA","SA","AL","MS"):
		#for i in range(0,len(nlon)):
			#if nlon[i]>=180:
				#nlon[i]=nlon[i]-360

	
	#if len(nlon.shape):
		#nlon,eralat=np.meshgrid(nlon,eralat_)

	#nera_lat,nera_lon, sevar= era_subregion(lat=eralat,lon=nlon,var=nevar[0,:],search_limits=search_limits)

	#sub_evar=np.empty((len(levels),evar.shape[0], evar.shape[1]))
	#sub_evar[:,:,:]=0

	#for index in range(0,evar.shape[0]):
		#nera_lat,nera_lon, svar_aux = era_subregion(lat=eralat,lon=nlon,var=nevar[index,:],search_limits=search_limits)
		#sub_evar[index,:]=svar_aux
	
	#if full:
	return evar, latt, lonn, levels
	#else:


		#return sub_evar, nera_lat,nera_lon, levels





def get_era5_hgtfile():
	pathpkg = os.path.dirname(__file__)
	hgt_file= pathpkg+"/ERA5_terrain_high.nc"
	return hgt_file


def convert_era5_matrix(matrix,lon,search_lon):
		
	index=np.where(lon==search_lon)
	index=int(index[0])
	aux_matrix=np.empty_like(matrix)
	aux_matrix[:,:index]=matrix[:,index:]
	aux_matrix[:,index:]=matrix[:,2:index+2]
	lon_aux=np.empty_like(lon)
	lon_aux[:index]=lon[index:]
	lon_aux[index:]=lon[:index]
	return aux_matrix,lon_aux


def era_subregion(lat=np.array(None),lon=np.array(None),var="",search_limits=[None,None,None,None]):
	latmin=search_limits[1]
	latmax=search_limits[3]
	lonmin=search_limits[0]
	lonmax=search_limits[2]
	
	
	l1=(lat>=latmin)&(lat<=latmax)
	l2=(lon>=lonmin)&(lon<=lonmax)

	m=l1&l2
	r,c=0,0

	for l in range(m.shape[0]):
		if m[l,:].sum():
			c=m[l,:].sum()
			break
	for l in range(m.shape[1]):
		if m[:,l].sum():
			r=m[:,l].sum()
			break

	nlon=lon[l1&l2].reshape(r,c)
	nlat=lat[l1&l2].reshape(r,c)
	nvar=var[l1&l2].reshape(r,c)

	return nlat,nlon,nvar


		

def get_wind_speed(latsc=[None],lonsc=[None],radius=[None],wfile="",idir="./",varu="U",varv="V",varlat="lat",varlon="lon",source="WRF",search_limits=[None,None,None,None],search_region="",r_uv=False):
	

	if source.upper()=="WRF":
		ncwfile=Dataset(idir+"/"+wfile)
		u=ncwfile.variables[varu][:]
		v=ncwfile.variables[varv][:]
		wlat=ncwfile.variables[varlat][:]
		wlon=ncwfile.variables[varlon][:]

		if len(u.shape)>2:
			u=u[0,:]
			v=v[0,:]
		if len(wlat.shape)>2:
			wlat=wlat[0,:]
			wlon=wlon[0,:]
	elif source.upper()=="ERA5":
		varlist=get_era5_2dvar(idir=idir,erafile=wfile,svariables=[varu,varv],search_limits=search_limits,search_region=search_region)
		wlat=varlist[0,:]
		wlon=varlist[1,:]
		u=varlist[2,:]
		v=varlist[3,:]
	
	elif source.upper()=="CUSTOM":
	
		varlist = get_custom_2dvar(idir=idir,customfile=wfile,svariables=[varu, varv],search_limits=search_limits,search_region=search_region, custom_latitude_var=varlat, custom_longitude_var=varlon)
		wlat=varlist[0,:]
		wlon=varlist[1,:]
		u=varlist[2,:]
		v=varlist[3,:]
		
	if r_uv==False:
		ff=np.sqrt(u**2+v**2)
		ffcentres=[]
		for i in range(0,len(latsc)):

			#dist=compute_grid_distance(lats=wlat,lons=wlon,latc=latsc[i],lonc=lonsc[i])
			dist=haversine(lon1=lonsc[i], lat1=latsc[i], lon2=wlon, lat2=wlat)
			if dist.min()>radius[i]:
				mws=-9999
			else:
				mws_=ff[dist<=radius[i]]
				mws=mws_.max()
			ffcentres=np.append(ffcentres,mws)
	
	if r_uv:
		return u,v
	else:
		return ffcentres



def tracker_cyclones(cyclone_type="",
		source="",
		idir="./",
		sourcefiles="",
		pathoutput="./",
		rout=2000,
		verbose=True,
		dates=[""],
		hours=[""],
		model_res=20,
		search_limits=[None,None,None,None],
		dr_res=100,
		d_ang=10,
		filter_center_threshold=800,
		critical_outer_radius=50,
		tmpdir="."+program_name()+"_tmpdir/",
		rank=0,
		search_region="",
		min_slp_threshold=1015,
		terrain_filter=1000,
		prev_days=14,
		mslp_anomaly_threshold=-3,
		source_filename_prefix="",
		source_file_date_format="",
		max_wind_speed_threshold=10,
		outer_wind_speed_threshold=2.5,
		vorticity_threshold=1.45e-5,
		great_circle_distance=5.5,
		dmslp_great_circle_distance=200,
		radius_for_msw=100,
		sourcefilesupper=[None],
		checking_upper_levels_parameters=False,
		idir_upper="./",
		plotting_maps=False,
		use_mslp_anomaly=True,
		custom_mslp_variable="",
		custom_latitude_var="",
		custom_longitude_var="",
		custom_uwind_variable="",
		custom_vwind_variable="",
		custom_terrain_high_filename="",
		custom_terrain_high_var_name=""):

	if source.upper()=="ERA5":
		hgt_file=get_era5_hgtfile()
		varhgt=get_era5_2dvar(idir="",erafile=hgt_file,svariables=["z"],search_limits=search_limits,search_region=search_region)
		hgt_field=varhgt[2,:]/9.80665

		varu="u10"
		varv="v10"
		varlat="latitude"
		varlon="longitude"

		dx,dy=compute_dx_dy(lons=varhgt[1,:],lats=varhgt[0,:])


	elif source.upper()=="WRF":
		hgtlat,hgtlon,hgt_field=get_wrf_hgt(idir=idir,wrffile=sourcefiles[0])

		varu="U10"
		varv="V10"
		varlat="XLAT"
		varlon="XLONG"
		dx,dy=compute_dx_dy(lons=hgtlon,lats=hgtlat)
		
	elif source.upper()=="CUSTOM":
		hgt_field =get_cumstom_hgt_data(custom_terrain_high_filename, custom_terrain_high_var_name)
		varu=custom_uwind_variable
		varv=custom_vwind_variable
		varlat=custom_latitude_var
		varlon=custom_longitude_var
		
		
		hgt_field=hgt_field/9.80665

	if terrain_filter>0:
		hgt_field=hgt_field
	elif terrain_filter==0:
		hgt_field[:,:]=-1

	for index in range(0,len(sourcefiles)):

		fdate=dates[index]+hours[index]
		if verbose:
			print("   ---> | Processing "+source.upper()+ ": "+idir+"/"+sourcefiles[index])
		if source.upper()=="ERA5":
			varlist=get_era5_2dvar(idir=idir,erafile=sourcefiles[index],svariables=["msl"],search_limits=search_limits,search_region=search_region)
			sourcelat=varlist[0,:]
			sourcelon=varlist[1,:]
			sourcemslp=varlist[2,:]/100


		elif source.upper()=="CUSTOM":
			varlist=get_custom_2dvar(idir=idir,customfile=sourcefiles[index],svariables=[custom_mslp_variable],search_limits=search_limits,search_region=search_region, custom_latitude_var=custom_latitude_var, custom_longitude_var=custom_longitude_var)
			sourcelat=varlist[0,:]
			sourcelon=varlist[1,:]
			sourcemslp=varlist[2,:]
			
			if len(str(int(sourcemslp.max())))>=5:
				sourcemslp=sourcemslp/100
			
			dx,dy=compute_dx_dy(lons=sourcelon,lats=sourcelat)
			

		elif source.upper()=="WRF":
			sourcelat,sourcelon,sourcemslp=get_wrf_mslp_new(idir=idir,wrffile=sourcefiles[index],variables=["PB","P","PHB","PH","T","QVAPOR"])


		if use_mslp_anomaly:	
			avg_mslp=get_mslp_anomaly(idir=idir,
							source=source,
							search_limits=search_limits,
							search_region=search_region,
							prev_days=prev_days,
							date=dates[index],
							hour=hours[index],
							source_filename_prefix=source_filename_prefix,
							source_file_date_format=source_file_date_format,
							custom_mslp_variable=custom_mslp_variable,
							custom_latitude_var=custom_latitude_var,
							custom_longitude_var=custom_longitude_var)

			mslp_anomaly=sourcemslp-avg_mslp
		else:
			mslp_anomaly=np.empty_like(sourcemslp)
			mslp_anomaly[:]=-9999
		
		if cyclone_type.upper() in ("TC","MC","TLC","SC","EC"):
			clats,clons,couter_r,cpmin,cmws,cclosedp,croci=get_low_centers(lats=sourcelat,
							lons=sourcelon,
							mslp=sourcemslp,
							mslp_anomaly=mslp_anomaly,
							mslp_anomaly_threshold=mslp_anomaly_threshold,
							use_mslp_anomaly=use_mslp_anomaly,
							model_res=model_res,
							search_limits=search_limits,
							rout=rout,
							dr_res=dr_res,
							d_ang=d_ang,
							critical_outer_radius=critical_outer_radius,
							min_slp_threshold=min_slp_threshold,
							terrain_filter=terrain_filter,
							hgt_field=hgt_field,
							sourcefile=sourcefiles[index],
							idir=idir,
							varu=varu,
							varv=varv,
							varlat=varlat,
							varlon=varlon,
							source=source.upper(),
							search_region=search_region,
							max_wind_speed_threshold=max_wind_speed_threshold,
							outer_wind_speed_threshold=outer_wind_speed_threshold,
							dx=dx,
							dy=dy,
							vorticity_threshold=vorticity_threshold,
							filter_center_threshold=filter_center_threshold,
							great_circle_distance=great_circle_distance,
							dmslp_great_circle_distance=dmslp_great_circle_distance,
							radius_for_msw=radius_for_msw)
		#elif cyclone_type.upper() in ("EC"):
			#clats,clons,couter_r,cpmin,cmws,cclosedp,croci=get_low_EC_centers(lats=sourcelat,
							#lons=sourcelon,
							#mslp=sourcemslp,
							#mslp_anomaly=mslp_anomaly,
							#mslp_anomaly_threshold=mslp_anomaly_threshold,
							#model_res=model_res,
							#search_limits=search_limits,
							#rout=rout,
							#dr_res=dr_res,
							#d_ang=d_ang,
							#critical_outer_radius=critical_outer_radius,
							#min_slp_threshold=min_slp_threshold,
							#terrain_filter=terrain_filter,
							#hgt_field=hgt_field,
							#sourcefile=sourcefiles[index],
							#idir=idir,
							#varu=varu,
							#varv=varv,
							#varlat=varlat,
							#varlon=varlon,
							#source=source.upper(),
							#search_region=search_region,
							#max_wind_speed_threshold=max_wind_speed_threshold,
							#outer_wind_speed_threshold=outer_wind_speed_threshold,
							#dx=dx,
							#dy=dy,
							#vorticity_threshold=vorticity_threshold,
							#filter_center_threshold=filter_center_threshold,
							#great_circle_distance=great_circle_distance,
							#dmslp_great_circle_distance=dmslp_great_circle_distance,
							#radius_for_msw=radius_for_msw)


		flats, flons,froci,fpmin,fclosedp,fmws,fouter_r=filter_centers(lats=np.array(clats),
								lons=np.array(clons),
								roci=np.array(croci),
								pmin=np.array(cpmin),
								closedp=np.array(cclosedp),
								ff=np.array(cmws),
								outer_r=np.array(couter_r),
								filter_center_threshold=filter_center_threshold)





		centers_data=np.empty((len(flats),8))
		centers_data[:,0]=flats
		centers_data[:,1]=flons
		centers_data[:,2]=fpmin
		centers_data[:,3]=fmws
		centers_data[:,4]=fclosedp
		centers_data[:,5]=froci
		centers_data[:,6]=fouter_r
		centers_data[:,7]=1	
		

		if plotting_maps:
			plotting_mslp_anomaly(sourcelons=sourcelon,
							sourcelats=sourcelat, 
							mslp=sourcemslp, 
							mslp_anomaly=mslp_anomaly, 
							search_limits=search_limits, 
							cenlons=flons, 
							cenlats=flats,
							search_region=search_region,
							fname=tmpdir+"/"+fdate+".png")
		
		
		np.savetxt(tmpdir+"/critical_centers_"+fdate+".dat",centers_data)
		


def compute_dx_dy(lons=np.array(None),lats=np.array(None)):
	dy=np.empty_like(lons)
	dy=dy[:-1,:]
	dy[:]=0
	dx=np.empty_like(lats)
	dx=dx[:,:-1]
	dx[:]=0

	for i in range(0,lats.shape[0]):
		for j in range(0,lats.shape[1]):
			if lons[i,j]>=180:
				lons[i,j]=lons[i,j]-360
			



	for i in range(0,lats.shape[0]-1):
			for j in range(0,lats.shape[1]-1):
						
				dist_dx= geo_distance(lat1=lats[i,j],
										lon1=lons[i,j],
										lat2=lats[i,j+1],
										lon2=lons[i,j+1])
				dist_dy=geo_distance(lat1=lats[i,j],
										lon1=lons[i,j],
										lat2=lats[i+1,j],
										lon2=lons[i+1,j])
			
				
				dx[i,j]=dist_dx
				dy[i,j]=dist_dy
					
	return dx*1000,dy*1000


def create_map(search_limits=[None, None, None, None], 
			 search_region="",):
	from cartopy import config
	from cartopy.util import add_cyclic_point
	import cartopy.feature as cfeature
	import cartopy.crs as ccrs
	from cartopy.mpl.geoaxes import GeoAxes
	from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
	import matplotlib.ticker as mticker
	if search_region in ("NA", "SA", "SI", "NI", "MS","AL"):
		center=0
	else:
		center=180
	
	
	
	min_lon,min_lat=(search_limits[0],search_limits[1])
	max_lon,max_lat=(search_limits[2],search_limits[3])
	paso_h=25
	

	
	crs = ccrs.PlateCarree()
	mapa=plt.subplot(1,1,1,projection=ccrs.PlateCarree(center) )
	mapa.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=1)
	mapa.add_feature(cfeature.STATES, linewidth=0.25)
	mapa.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
	
	
	gl = mapa.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='black', alpha=1, linestyle='--')
	lons=np.arange(math.ceil(min_lon),math.ceil(max_lon),paso_h)
	
	gl_lon_info=[]
	for clons in lons:
		if clons<180:
			gl_lon_info=np.append(gl_lon_info,clons)
		else:
			gl_lon_info=np.append(gl_lon_info,clons-360)

	gl_loc=[True,False,False,True]
	gl.ylabels_left = gl_loc[0]
	gl.ylabels_right = gl_loc[1]
	gl.xlabels_top = gl_loc[2]
	gl.xlabels_bottom = gl_loc[3]

	lons=np.arange(math.floor(min_lon-paso_h),math.ceil(max_lon+paso_h),paso_h)
	gl.xlocator = mticker.FixedLocator(lons)
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER
	gl.xlabel_style = {'size': 25, 'color': 'black'}
	gl.ylabel_style = {'size': 25,'color': 'black'}


	return mapa,crs


def plotting_mslp_anomaly(sourcelons=np.array([None]),
						  sourcelats=np.array([None]), 
						  mslp=np.array([None]), 
						  mslp_anomaly=np.array([None]), 
						  search_limits=[None, None, None, None], 
						  cenlons=np.array([None]), 
						  cenlats=np.array([None]),
						  search_region="",
						  fname="./"):
	
	
		
	plt.figure(figsize=(18,12))
	
	mapa,crs=create_map(search_limits=search_limits, search_region=search_region)
	
	lvs=np.arange(int(mslp.min()), int(mslp.max())+3,3)
	cs=mapa.contour(sourcelons,sourcelats,mslp,lvs, transform=crs)
	plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=17) 
	lvs=np.arange(-5,5,0.2)
	cf=mapa.contourf(sourcelons,sourcelats,mslp_anomaly,lvs, cmap=plt.cm.PiYG_r, extend="both", transform=crs)
	cb=plt.colorbar(cf, orientation="horizontal",pad=0.05,shrink=0.95)
	cb.set_label(label=" 14-day MSLP anomalies (hPa)", size=25,labelpad=15)
	cb.ax.tick_params(labelsize=25)
	mapa.scatter(cenlons,cenlats,color="b",s=90, marker="o",transform=crs)
	
	plt.savefig(fname, bbox_inches="tight", dpi=600)
	plt.close()



def first_derivative(f, axis=None, x=None, delta=None):
    
	n, axis, delta = _process_deriv_args(f, axis, x, delta)
	take = make_take(n, axis)

	slice0 = take(slice(None, -2))
	slice1 = take(slice(1, -1))
	slice2 = take(slice(2, None))
	delta_slice0 = take(slice(None, -1))
	delta_slice1 = take(slice(1, None))

	combined_delta = delta[delta_slice0] + delta[delta_slice1]
	delta_diff = delta[delta_slice1] - delta[delta_slice0]
	center = (- delta[delta_slice1] / (combined_delta * delta[delta_slice0]) * f[slice0]
				+ delta_diff / (delta[delta_slice0] * delta[delta_slice1]) * f[slice1]
				+ delta[delta_slice0] / (combined_delta * delta[delta_slice1]) * f[slice2])

	slice0 = take(slice(None, 1))
	slice1 = take(slice(1, 2))
	slice2 = take(slice(2, 3))
	delta_slice0 = take(slice(None, 1))
	delta_slice1 = take(slice(1, 2))

	combined_delta = delta[delta_slice0] + delta[delta_slice1]
	big_delta = combined_delta + delta[delta_slice0]
	left = (- big_delta / (combined_delta * delta[delta_slice0]) * f[slice0]
			+ combined_delta / (delta[delta_slice0] * delta[delta_slice1]) * f[slice1]
			- delta[delta_slice0] / (combined_delta * delta[delta_slice1]) * f[slice2])


	slice0 = take(slice(-3, -2))
	slice1 = take(slice(-2, -1))
	slice2 = take(slice(-1, None))
	delta_slice0 = take(slice(-2, -1))
	delta_slice1 = take(slice(-1, None))

	combined_delta = delta[delta_slice0] + delta[delta_slice1]
	big_delta = combined_delta + delta[delta_slice1]
	right = (delta[delta_slice1] / (combined_delta * delta[delta_slice0]) * f[slice0]
				- combined_delta / (delta[delta_slice0] * delta[delta_slice1]) * f[slice1]
				+ big_delta / (combined_delta * delta[delta_slice1]) * f[slice2])

	return np.concatenate((left, center, right), axis=axis)


def _process_deriv_args(f, axis, x, delta):
   
	n = f.ndim
	axis = normalize_axis_index(axis if axis is not None else 0, n)

	if f.shape[axis] < 3:
		raise ValueError('f must have at least 3 point along the desired axis.')

	if delta is not None:
		if x is not None:
			raise ValueError('Cannot specify both "x" and "delta".')

		delta = np.atleast_1d(delta)
		if delta.size == 1:
			diff_size = list(f.shape)
			diff_size[axis] -= 1
			delta = np.broadcast_to(delta, diff_size, subok=True)
		else:
			delta = _broadcast_to_axis(delta, axis, n)
	elif x is not None:
		x = _broadcast_to_axis(x, axis, n)
		delta = np.diff(x, axis=axis)
	else:
		raise ValueError('Must specify either "x" or "delta" for value positions.')

	return n, axis, delta


def calc_vorticity(u, v, dx=None, dy=None, x_dim=-1, y_dim=-2):
    
    dudy = first_derivative(u, delta=dy, axis=y_dim)
    dvdx = first_derivative(v, delta=dx, axis=x_dim)
    return dvdx - dudy

def _broadcast_to_axis(arr, axis, ndim):
    if arr.ndim == 1 and arr.ndim < ndim:
        new_shape = [1] * ndim
        new_shape[axis] = arr.size
        arr = arr.reshape(*new_shape)
    return arr

def make_take(ndims, slice_dim):
    def take(indexer):
        return tuple(indexer if slice_dim % ndims == i else slice(None) 
                     for i in range(ndims))
    return take


def compute_grid_distance(lats=np.array(None),lons=np.array(None),latc=None,lonc=None):
	dist=np.empty_like(lats)
	
	if lonc>180:
		lonc=lonc-360
	
	for i in range(0,lats.shape[0]):
		for j in range(0,lats.shape[1]):
			
			if lons[i,j]>180:
				lons[i,j]=lons[i,j]-360
			
			
			dist[i,j]=geo_distance(lat1=latc,
							lon1=lonc,
							lat2=lats[i,j],
							lon2=lons[i,j])
			

	return dist



def filter_centers(lats=np.array(None),
			lons=np.array(None),
			roci=np.array(None),
			pmin=np.array(None),
			closedp=np.array(None),
			ff=np.array([None]),
			outer_r=np.array([None]),
			filter_center_threshold=1000):

	
	
	if len(lats)>0:
		if ff[0]==None:
			ff=np.empty_like(lons)
			ff[:]=0
		if outer_r[0]==None:
			outer_r=np.empty_like(lons)
			outer_r[:]=0



		flats=[]
		flons=[]
		froci=[]
		fpmin=[]
		fclosedp=[]
		fmws=[]
		fouter_r=[]

		for i in range(0,len(lats)):
			
			dist=[]
			if np.isfinite(pmin[i]) and np.isfinite(lats[i]) and np.isfinite(lons[i]):
				if lons[i]>180:
					lons[i]=lons[i]-360
				coord_indices=[]
				npmin=[]
				check_roci=roci[i]
				dist_check=False
				for j in range(0,len(lats)):
					if lons[j]>180:
						lons[j]=lons[j]-360

					
					if np.isfinite(pmin[j]) and np.isfinite(lats[j]) and np.isfinite(lons[j]) :
						pdist= geo_distance(lat1=lats[i],
							lon1=lons[i],
							lat2=lats[j],
							lon2=lons[j])
						dist=np.append(dist,pdist)
						coord_indices=np.append(coord_indices,j)
						npmin=np.append(npmin,pmin[j])
						check_roci=np.append(check_roci,roci[j])
						dist_check=True

			dist=np.array(dist)
			if len(dist)>0 and dist_check==True:

				indices=np.where(dist<filter_center_threshold)
				if len(indices[0])==0 or len(indices[0])==1:
					flats=np.append(flats, lats[i])
					flons=np.append(flons,lons[i])
					froci=np.append(froci, roci[i])
					fpmin=np.append(fpmin, pmin[i])
					fclosedp=np.append(fclosedp,closedp[i])
					fmws=np.append(fmws,ff[i])
					fouter_r=np.append(fouter_r,outer_r[i])
					

				elif len(indices[0])>1:
					indices=indices[0]
					check_pmin=pmin[i]
					check_indices=i
					for m in range(0,len(indices)):
						check_pmin=np.append(check_pmin,npmin[int(indices[m])])
						check_indices=np.append(check_indices,coord_indices[int(indices[m])])
					check_pmin=np.array(check_pmin)
					i_npmin=np.where(check_pmin==check_pmin.min())
					pindex=int(i_npmin[0][0])
					flats=np.append(flats, lats[int(check_indices[pindex])])
					flons=np.append(flons, lons[int(check_indices[pindex])])
					froci=np.append(froci, roci[int(check_indices[pindex])])
					fpmin=np.append(fpmin, pmin[int(check_indices[pindex])])
					fclosedp=np.append(fclosedp,closedp[int(check_indices[pindex])])
					fmws=np.append(fmws,ff[int(check_indices[pindex])])
					fouter_r=np.append(fouter_r,outer_r[int(check_indices[pindex])])
					
					pmin[int(check_indices[pindex])]=np.nan
					lats[int(check_indices[pindex])]=np.nan
					lons[int(check_indices[pindex])]=np.nan
					for k in range(0,len(indices)):
						pmin[int(coord_indices[int(indices[k])])] = np.nan
						lats[int(coord_indices[int(indices[k])])] = np.nan
						lons[int(coord_indices[int(indices[k])])] = np.nan
			pmin[i]=np.nan
			lats[i]=np.nan
			lons[i]=np.nan
			dist=[]
	else:
			flats=[0]
			flons=[0]
			froci=[0]
			fpmin=[0]
			fclosedp=[0]
			fmws=[0]
			fouter_r=[0]
			fctype=[0]
	return flats, flons,froci,fpmin, fclosedp,fmws,fouter_r



def get_wrf_mslp(idir="./",wrffile="",variables=[""]):
	ncwrffile=Dataset(idir+"/"+wrffile)

	dimNx=ncwrffile.dimensions['west_east'].size
	dimNy=ncwrffile.dimensions['south_north'].size
	dimNz=ncwrffile.dimensions['bottom_top'].size

	wrflat=ncwrffile.variables["XLAT"][0,:]
	wrflon=ncwrffile.variables["XLONG"][0,:]
	for var in variables:
		vars()[var]=ncwrffile.variables[var][0,:]

	ncwrffile.close()
	

	wrfmslp=compute_mslp(psfc=vars()['PSFC'],T2=vars()['T2'],phb=vars()['PHB'],ph=vars()['PH'])

	return wrflat,wrflon,wrfmslp



def get_wrf_mslp_new(idir="./",wrffile="",variables=[""]):
	ncwrffile=Dataset(idir+"/"+wrffile)

	dimNx=ncwrffile.dimensions['west_east'].size
	dimNy=ncwrffile.dimensions['south_north'].size
	dimNz=ncwrffile.dimensions['bottom_top'].size

	wrflat=ncwrffile.variables["XLAT"][0,:]
	wrflon=ncwrffile.variables["XLONG"][0,:]

	TBASE=ncwrffile.variables["T00"][0]

	for var in variables:
		vars()[var]=ncwrffile.variables[var][0,:]

	ncwrffile.close()

	wrfmslp=slp(PB=vars()['PB'], P=vars()['P'], PHB=vars()['PHB'], PH=vars()['PH'], T=vars()['T'], QVAPOR=vars()['QVAPOR'],TBASE=TBASE)

	return wrflat,wrflon,wrfmslp


def slp(PB=np.array(None), P=np.array(None), PHB=np.array(None), PH=np.array(None), T=np.array(None), QVAPOR=np.array(None), TBASE=300.0):
	ptot = P + PB
	t = (TBASE+T) * (ptot/100000.) ** (287.04/1004.)
	# compute geopotential
	ph = (PHB + PH)/9.81
	nz,ny,nx = ph.shape
	# Unstagger geopotential
	ph[1:nz-1,:,:] = 0.5 * (ph[1:nz-1,:,:] + ph[2:nz,:,:])

	# These constants are from WRF-POST
	TC = 273.16+17.5
	PCONST = 10000.

	# Find lowest level that is at least PCONST above the surface
	klo = np.argmax(ptot < ptot[0]-PCONST, axis=0) - 1

	# Get the temperature, pressure and moisture at this level
	# and the level above
	klo[(klo-1 < 0)] = 0
	khi = klo+1
	nz,ny,nx = t.shape
	Tlo = t.reshape(nz,ny*nx)[klo.flatten(),range(ny*nx)].reshape(ny,nx)
	Plo = ptot.reshape(nz,ny*nx)[klo.flatten(),range(ny*nx)].reshape(ny,nx)
	Qlo = QVAPOR.reshape(nz,ny*nx)[klo.flatten(),range(ny*nx)].reshape(ny,nx)
	Thi = t.reshape(nz,ny*nx)[khi.flatten(),range(ny*nx)].reshape(ny,nx)
	Phi = ptot.reshape(nz,ny*nx)[khi.flatten(),range(ny*nx)].reshape(ny,nx)
	Qhi = QVAPOR.reshape(nz,ny*nx)[khi.flatten(),range(ny*nx)].reshape(ny,nx)
	nz,ny,nx = ph.shape
	Zhi = ph.reshape(nz,ny*nx)[khi.flatten(),range(ny*nx)].reshape(ny,nx)
	Zlo = ph.reshape(nz,ny*nx)[klo.flatten(),range(ny*nx)].reshape(ny,nx)
	# Virtual temperature correction
	Tlo = Tlo * (1.0+0.608 * Qlo)
	Thi = Thi * (1.0+0.608 * Qhi)

	p_at_pconst = ptot[0] - PCONST
	t_at_pconst = Thi-(Thi-Tlo)*np.log(p_at_pconst/Phi)*np.log(Plo/Phi)
	z_at_pconst = Zhi-(Zhi-Zlo)*np.log(p_at_pconst/Phi)*np.log(Plo/Phi)

	t_surf = t_at_pconst * (ptot[0]/p_at_pconst) ** (0.0065*287.04/9.81)
	t_sea_level = t_at_pconst + 0.0065*z_at_pconst

	# They call this the "Ridiculous MM5 test" for too hot temperatures
	l1 = t_sea_level < TC

	l2 = t_surf <= TC
	locs = np.bitwise_and(l2,~l1)
	t_sea_level = TC - 0.005*(t_surf-TC)**2
	t_sea_level[locs] = TC

	# Now final computation
	z_half_lowest = ph[0]
	mslp = ptot[0] * np.exp((2.*9.81*z_half_lowest)/\
						(287.04*(t_sea_level+t_surf)))
	#print np.max(out/100.,axis=None), np.min(out/100.,axis=None)
	# Smooth the output field and return in hPa
	#[gaussian_filter(out/100., sigma=3)]
	mslp=mslp/100
	return  mslp


def interp_levels(x, y, levels):
    shp = y.shape

    x = np.reshape(x, [shp[0],-1]).T
    y = np.reshape(y, [shp[0],-1]).T
    
    new_shape = (len(y), len(levels))

    values = np.empty( new_shape )
    
    for i, (xo, yo) in enumerate(zip(x, y)):
        interp_mod = interp1d(xo, yo, fill_value='extrapolate', axis=0)
        values[i] = interp_mod(levels)
    var_levels=values.T.reshape( new_shape[1:] + shp[1:])
    return var_levels



def interpolate_to_non_staggered_sigma_levels(non_staggered, staggered, matrix):
	
	new_matrix=np.empty((len(non_staggered), matrix.shape[1], matrix.shape[2]))
	new_matrix[:,:,:]=0
	for i in range(0, matrix.shape[1]):
		for j in range(0,matrix.shape[2]):
			
			idx = np.argsort(staggered)
			
			vals=matrix[:,i,j]
			
			new_matrix[:,i,j]=np.interp(non_staggered, staggered[idx], vals[idx])
		
	return new_matrix




def get_wrf_uppervar(idir="./",
					 wrffile="",
					 variables=[""],
					 plevels=np.array([None])):
	ncwrffile=Dataset(idir+"/"+wrffile)
	dimNx=ncwrffile.dimensions['west_east'].size
	dimNy=ncwrffile.dimensions['south_north'].size
	dimNz=ncwrffile.dimensions['bottom_top'].size

	wrflat=ncwrffile.variables["XLAT"][0,:]
	wrflon=ncwrffile.variables["XLONG"][0,:]
	znstag=ncwrffile.variables["ZNU"][0,:]
	zstag=ncwrffile.variables["ZNW"][0,:]
	
	Pp=ncwrffile.variables['P'][0,:]
	Pb=ncwrffile.variables['PB'][0,:] 
	
	P=Pp + Pb
	P=P/100
	
	if "PH" in variables:
		pgh=ncwrffile.variables[variables[0]][0,:]
		pghb=ncwrffile.variables["PHB"][0,:]
		aux_var=pghb+pgh
		var=interpolate_to_non_staggered_sigma_levels(znstag, zstag, aux_var)
		var=var/9.80665
	else:
		var=ncwrffile.variables[variable][0,:]
	
	ncwrffile.close()
	
	if plevels[0]==None:
		plevels=np.arange(1000,150,-50)
	
	
	nvar=interp_levels(P, var, plevels)
	
	return wrflat,wrflon,nvar,plevels





def compute_mslp(psfc=np.array(None),T2=np.array(None),phb=np.array(None),ph=np.array(None)):
	Rd=287.5
	gamma=0.00065
	g=9.8 
	hgp=phb+ph
	hgp=hgp/9.805

	psfc=psfc/100
	alpha=gamma*Rd/g
	
	exp1=hgp[0,:]/(Rd*T2)
	exp2=1-alpha*hgp[0,:]/(2*Rd*T2) + (1/3.)*(alpha*hgp[0,:]/(2*Rd*T2))**2

	mslp=psfc*np.exp(exp1*exp2)
	
	return mslp



def calc_relvort(u=np.array([None]), v=np.array([None]), lon=np.array([None]), lat=np.array([None])):

	pi = math.pi;
	pid = pi/180.;
	R_earth = 6371200.;
	lon2d, lat2d = lon, lat
	a_d2r = R_earth*pid 

	# Calculate hor distance in m
	dx = np.zeros([lat.shape[0], lat.shape[1]])
	dx[:,    0] = a_d2r*(lon2d[:,  1] - lon2d[:,   0])*np.cos(lat2d[:, 0]*pid)
	dx[:,   -1] = a_d2r*(lon2d[:, -1] - lon2d[:,  -2])*np.cos(lat2d[:,-1]*pid)
	dx[:, 1:-1] = a_d2r*(lon2d[:, 2:] - lon2d[:, :-2])*np.cos(lat2d[:, 1:-1]*pid)/2.

	# Calculate hor distance in m
	dy = np.zeros([lat.shape[0], lon.shape[1]])
	dy[   0, :] = a_d2r*(lat2d[ 1, :] - lat2d[  0, :])
	dy[  -1, :] = a_d2r*(lat2d[-1, :] - lat2d[ -2, :])
	dy[1:-1, :] = a_d2r*(lat2d[2:, :] - lat2d[:-2, :])/2.

	# Calculate derivatives of v in x direction
	dvdx = np.zeros([v.shape[0], v.shape[1]])
	dvdx[ :,    0] = (v[ :,  1] - v[ :,   0])/dx[ :,  0]
	dvdx[:,   -1] = (v[ :, -1] - v[ :,  -2])/dx[ :, -1]
	dvdx[ :, 1:-1] = (v[ :, 2:] - v[ :, :-2])/(2.*dx[:, 1:-1])

	# Calculate derivatives of u in y direction
	dudy = np.zeros([u.shape[0], u.shape[1]])

	dudy[  0, :] = (u[1, :]*np.cos(lat2d[1, :]*pid) \
					-u[ 0, :]*np.cos(lat2d[0, :]*pid))/(dy[0, :]*np.cos(lat2d[0, :]*pid))
	dudy[  -1, :] = (u[ -1, :]*np.cos(lat2d[-1, :]*pid) \
					-u[ -2, :]*np.cos(lat2d[-2, :]*pid))/(dy[-1, :]*np.cos(lat2d[-1, :]*pid))
	dudy [1:-1, :] = (u[ 2:, :]*np.cos(lat2d[2:, :]*pid) \
					-u[ :-2, :]*np.cos(lat2d[:-2, :]*pid))/(2.*dy[1:-1, :]*np.cos(lat2d[1:-1, :]*pid))
	return dvdx - dudy



def get_maxmin_mslp_points(lats=np.array([None]), lons=np.array([None]),  mslp=np.array([None]), mslp_anomaly=np.array([None]), hgt_field=np.array([None]), rel_vort=np.array([None]), extrema="min", nsize=25):
	data_ext = minimum_filter(mslp, nsize, mode='nearest')
	mxy, mxx = np.where(data_ext == mslp)


	nlats=[]
	nlons=[]
	nmslp=[]
	nmslp_anoms=[]
	nterrain=[]
	nrelvort=[]
	for i in range(len(mxy)):
		nlats=np.append(nlats, lats[mxy[i], mxx[i]])
		nlons=np.append(nlons, lons[mxy[i], mxx[i]])
		nmslp=np.append(nmslp, mslp[mxy[i], mxx[i]])
		nmslp_anoms=np.append(nmslp_anoms, mslp_anomaly[mxy[i], mxx[i]])
		nterrain=np.append(nterrain, hgt_field[mxy[i], mxx[i]])
		nrelvort=np.append(nrelvort, rel_vort[mxy[i], mxx[i]])
		
	return nlats, nlons, nmslp, nmslp_anoms, nterrain, nrelvort


def get_low_centers(lats=np.array([None]),
		lons=np.array([None]), 
		mslp=np.array([None]),
		mslp_anomaly=np.array([None]),
		mslp_anomaly_threshold=None,
		use_mslp_anomaly=True,
		model_res=None,
		search_limits=[None,None,None,None],
		rout=2000,
		dr_res=100,
		d_ang=10,
		critical_outer_radius=200,
		min_slp_threshold=1015,
		terrain_filter=1000,
		hgt_field=np.array([None]),
		sourcefile="",
		idir='./',
		varu="",
		varv="",
		varlat="",
		varlon="",
		source="",
		search_region="",
		max_wind_speed_threshold=10,
		outer_wind_speed_threshold=2.5,
		dx=20,
		dy=20,
		vorticity_threshold=1e-5,
		filter_center_threshold=800,
		great_circle_distance=5.5,
		dmslp_great_circle_distance=200,
		radius_for_msw=100):


	center_lats=[]
	center_lons=[]
	outer_r=[]
	outer_p=[]
	mcp=[]
	mws=[]
	croci=[]
	ctype=[]
	npi=int(111/model_res)
	great_circle_distance=great_circle_distance*111

	u,v=get_wind_speed(latsc=[0],
					lonsc=[0],
					radius=[critical_outer_radius],
					wfile=sourcefile,
					idir=idir,
					varu=varu,
					varv=varv,
					varlat=varlat,
					varlon=varlon,
					source=source.upper(),
					search_limits=search_limits,
					search_region=search_region,
					r_uv=True)
	


	rel_vort=calc_relvort(u=u, v=v, lon=lons, lat=lats)
	
	
	if search_region in ("SA","SH","SI","SP"):
		rel_vort=rel_vort*(-1)
	
	
	minlats, minlons, minmslp, minmslp_anoms, minterrain, minrelvort =get_maxmin_mslp_points(lats=lats, lons=lons,  mslp=mslp, mslp_anomaly=mslp_anomaly, hgt_field=hgt_field, rel_vort=rel_vort, extrema="min", nsize=25)
	
	bulk_amslp=np.copy(minmslp)
	
	bulk_amslp[minmslp>min_slp_threshold]=np.nan
	if use_mslp_anomaly:
		bulk_amslp[minmslp_anoms>mslp_anomaly_threshold]=np.nan
	if terrain_filter>0:
		bulk_amslp[minterrain>terrain_filter]=np.nan
	bulk_amslp[minrelvort<vorticity_threshold]=np.nan
	bulk_amslp[minlats>search_limits[3]]=np.nan
	bulk_amslp[minlats<search_limits[1]]=np.nan
	bulk_amslp[minlons>search_limits[2]]=np.nan
	bulk_amslp[minlons<search_limits[0]]=np.nan
	
	
	
	nnmslp=minmslp[np.isfinite(bulk_amslp)]
	nnlats=minlats[np.isfinite(bulk_amslp)]
	nnlons=minlons[np.isfinite(bulk_amslp)]
	
	
	
	#bulk_amslp=np.copy(mslp_anomaly)
	
	#bulk_amslp[bulk_amslp>mslp_anomaly_threshold]=np.nan
	#bulk_amslp[mslp>min_slp_threshold]=np.nan
	#if terrain_filter>0:
	#	bulk_amslp[hgt_field>terrain_filter]=np.nan
	#bulk_amslp[rel_vort<vorticity_threshold]=np.nan
	#bulk_amslp[lats>search_limits[3]]=np.nan
	#bulk_amslp[lats<search_limits[1]]=np.nan
	#bulk_amslp[lons>search_limits[2]]=np.nan
	#bulk_amslp[lons<search_limits[0]]=np.nan
	


	#nnmslp=mslp[np.isfinite(bulk_amslp)]
	#nnlats=lats[np.isfinite(bulk_amslp)]
	#nnlons=lons[np.isfinite(bulk_amslp)]
	
	
	nnroci=np.empty_like(nnlats)
	nnroci[:]=0
	
	nnclosedp=np.empty_like(nnlats)
	nnclosedp[:]=0

	nnmws=np.empty_like(nnlats)
	nnmws[:]=0

	nncouter_r=np.empty_like(nnlats)
	nncouter_r[:]=0

	nlats, nlons,nroci,nmslp,nclosedp,nmws,nouter_r=filter_centers(lats=nnlats,
									lons=nnlons,
									roci=nnroci,
									pmin=nnmslp,
									closedp=nnclosedp,
									ff=nnmws,
									outer_r=nncouter_r,
									filter_center_threshold=filter_center_threshold)

	for i in range(0, len(nmslp)):
		dmslp,nlatc, nlonc, npmin=compute_dmslp(latc=nlats[i],
							lonc=nlons[i],
							lats=lats,
							lons=lons,
							mslp=mslp,
							pmin=nmslp[i],
							dang=d_ang,
							dradius=model_res,
							search_radius=great_circle_distance+100,
							model_res=model_res,
							great_circle_distance=great_circle_distance,
							filter_center_threshold=.45*filter_center_threshold,
							search_limits=search_limits)

		
		if dmslp*100 - npmin*100 >= dmslp_great_circle_distance:
			
			fwind_speed=get_wind_speed(latsc=[nlatc],
					lonsc=[nlonc],
					radius=[radius_for_msw],
					wfile=sourcefile,
					idir=idir,
					varu=varu,
					varv=varv,
					varlat=varlat,
					varlon=varlon,
					source=source.upper(),
					search_limits=search_limits,
					search_region=search_region,
					r_uv=False)


			if fwind_speed[0] >= max_wind_speed_threshold:
				
				outer_size=compute_TC_size(latc=nlatc,
						lonc=nlonc,
						lats=lats,
						lons=lons,
						u=u,
						v=v,
						dradius=dr_res,
						dang=d_ang,
						search_radius=rout,
						model_res=model_res,
						outer_wind_speed_threshold=outer_wind_speed_threshold,
						)

				roci,closedp=compute_roci_RU(latc=nlatc,
					lonc=nlonc,
					lats=lats,
					lons=lons,
					mslp=mslp,
					pmin=npmin,
					dang=d_ang,
					dradius=dr_res,
					search_radius=rout,
					model_res=model_res)
			
			
				
				if np.isnan(roci):roci=0
				if np.isnan(closedp):closedp=0
			
			

				if roci>=critical_outer_radius or closedp>=0:
					center_lats=np.append(center_lats,nlatc)
					center_lons=np.append(center_lons,nlonc)
					outer_r=np.append(outer_r,outer_size)
					mcp=np.append(mcp,npmin)
					mws=np.append(mws,fwind_speed[0])
					outer_p=np.append(outer_p,closedp)
					croci=np.append(croci,roci)


	return center_lats,center_lons,outer_r,mcp,mws,outer_p,croci


def haversine(lon1, lat1, lon2, lat2):
	# convert decimal degrees to radians 
	lon1 = np.deg2rad(lon1)
	lon2 = np.deg2rad(lon2)
	lat1 = np.deg2rad(lat1)
	lat2 = np.deg2rad(lat2)

	# haversine formula 
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
	c = 2 * np.arcsin(np.sqrt(a)) 
	r = 6371
	return c * r



def get_GCD(lon1, lat1, lonc, latc):	

	r = 6371
	radlatc=latc*np.pi/180
	radlonc=lonc*np.pi/180

	radlat1=lat1*np.pi/180
	radlon1=lon1*np.pi/100

	a=np.sin(radlatc)*np.sin(radlat1) + np.cos(radlatc)*np.cos(radlat1)*np.cos(radlonc-radlon1)

	gcd=r*np.arccos(a)

	return gcd


#def get_cyclone_phase(VTU=[None],VTL=[None],B=[None]):
	
	#cyclone_phase=[]
	
	#for i in range(0,len(VTU)-1):
	
		#if VTU[i]<0 and VTL[i]<0 :
			#cyclone_phase=np.append(cyclone_phase,10)
		#elif VTU[i]>0 and VTL[i]>0:
			#cyclone_phase=np.append(cyclone_phase,20)
		#elif VTU[i]<0 and VTL[i]>0:
			#cyclone_phase=np.append(cyclone_phase,30)
		#else:
			#cyclone_phase=np.append(cyclone_phase,-99)
		
	#return cyclone_phase


def get_cyclone_phase(VTU=[None],VTL=[None],B=[None]):
	
	cyclone_phase=[]
	
	for i in range(0,len(VTU)-1):
	
		if VTU[i]>0 and VTL[i]>0  and np.abs(B[i])<=10:
			cyclone_phase=np.append(cyclone_phase,"SDWC")
			
		elif VTU[i]<0 and VTL[i]<0  and np.abs(B[i])<=10:
			cyclone_phase=np.append(cyclone_phase,"SDCC")
		
		elif VTU[i]<0 and VTL[i]>0  and np.abs(B[i])<=10:
			cyclone_phase=np.append(cyclone_phase,"SLWC")
		
		elif VTU[i]>0 and VTL[i]<0  and np.abs(B[i])<=10:
			cyclone_phase=np.append(cyclone_phase,"SLCC")
		
		
		elif VTU[i]>0 and VTL[i]>0  and np.abs(B[i])>10:
			cyclone_phase=np.append(cyclone_phase,"ADWC")
		
		elif VTU[i]<0 and VTL[i]<0  and np.abs(B[i])>10:
			cyclone_phase=np.append(cyclone_phase,"ADCC")
		
		elif VTU[i]<0 and VTL[i]>0  and np.abs(B[i])>10:
			cyclone_phase=np.append(cyclone_phase,"ALWC")
		
		elif VTU[i]>0 and VTL[i]<0  and np.abs(B[i])>10:
			cyclone_phase=np.append(cyclone_phase,"ALCC")
		else:
			cyclone_phase=np.append(cyclone_phase,"UDCC")
		
		
	cyclone_phase=np.append(cyclone_phase,"UDCC")	
	return cyclone_phase



def get_cyclone_type(VTU=None,
					 VTL=None, 
					 B=[None], 
					 cyclone_type="None",
					 core_criteria_length=0,
					 dt_h=6,
					 VTL_threshold=0,
					 VTU_threshold=0,
					 Bhart_threshold=10,):
	
	if cyclone_type.upper()=="TC":
		cont=0

		for i in range(0, len(VTU)-1):
			if  VTU[i]>VTU_threshold and VTL[i]>VTL_threshold and np.abs(B[i])<Bhart_threshold:
				cont=cont+1

				
		if core_criteria_length==-99:
			if cont == (len(VTU)-1):
				return True
			else:
				return False
		elif core_criteria_length>=0:
			if cont >= core_criteria_length:
				return True
			else:
				return False

				
	elif cyclone_type.upper()=="TLC":
		cont=0

		for i in range(0, len(VTU)-1):
			if  VTU[i]>VTU_threshold and VTL[i]>VTL_threshold and np.abs(B[i])<Bhart_threshold:
				cont=cont+1
		
		if core_criteria_length==-99:
			if cont == (len(VTU)-1):
				return True
			else:
				return False
		elif core_criteria_length>=0:
			if cont >= core_criteria_length:
				return True
			else:
				return False
		
		
	
	elif cyclone_type.upper()=="EC":
		cont=0

		for i in range(0, len(VTU)-1):
			if  VTU[i]<VTU_threshold and VTL[i]<VTL_threshold and np.abs(B[i])>Bhart_threshold:
				cont=cont+1
		
		
		if core_criteria_length==-99:
			if cont == (len(VTU)-1):
				return True
			else:
				return False
		elif core_criteria_length>=0:
			if cont >= core_criteria_length:
				return True
			else:
				return False
		
		
	
	elif cyclone_type.upper()=="MC":
		return True
	
	
	elif cyclone_type.upper()=="SC":
		cont=0

		for i in range(0, len(VTU)-1):
			if  VTU[i]<VTU_threshold and VTL[i]>VTL_threshold and np.abs(B[i])<Bhart_threshold:
				cont=cont+1

		if core_criteria_length==-99:
			if cont == (len(VTU)-1):
				return True
			else:
				return False
		elif core_criteria_length>0:
			if cont >= core_criteria_length:
				return True
			else:
				return False
		


def compute_dmslp(latc=None,
				lonc=None,
				lats=np.array(None),
				lons=np.array(None),
				mslp=np.array(None),
				pmin=None,dang=10,
				dradius=50,
				search_radius=750,
				model_res=20,
				great_circle_distance=650,
				filter_center_threshold=350,
				search_limits=[None,None,None,None]
				):


	distance=haversine(lon1=lons, lat1=lats, lon2=lonc, lat2=latc)


	
	#distance=get_GCD(lon1=lons, lat1=lats, lonc=lonc, latc=latc)
	
	
	aux_mslp=np.copy(mslp)
	
	aux_mslp[distance<great_circle_distance]=np.nan
	aux_mslp[distance>search_radius]=np.nan

	dmslp=mslp[np.isfinite(aux_mslp)]


	storm_centre=False
	while storm_centre==False:
		latp,lonp,radius,theta=polar_cords(latc=latc,
					lonc=lonc,
					dth=math.radians(dang),
					dr=dradius,
					search_radius=search_radius)
		pointint=np.empty((len(latp.flatten()),2))
		pointint[:,0]=lonp.flatten()
		pointint[:,1]=latp.flatten()

		pointdata=np.empty((len(lats.flatten()),2))
		pointdata[:,0]=lons.flatten()
		pointdata[:,1]=lats.flatten()

		pminp=griddata(pointdata, mslp.flatten(), pointint, method='nearest')
		pminp=pminp.reshape(int(latp.shape[0]),int(latp.shape[1]))


		nradius,ntheta=np.meshgrid(radius,theta)
		npminp=np.copy(pminp)
		npminp[nradius>filter_center_threshold]=10000

		ii,jj=np.where(npminp==npminp.min())
		ii=int(ii[0])
		jj=int(jj[0])
		npmin=npminp[ii,jj]
		nlatc=latp[ii,jj]
		nlonc=lonp[ii,jj]

		if  npmin<pmin and search_limits[1]<nlatc<search_limits[3] and search_limits[0]<nlonc<search_limits[2]:
			pmin=npmin
			latc=nlatc
			lonc=nlonc
			storm_centre=True
			dmslp,latc, lonc, pmin=compute_dmslp(latc=latc,
												lonc=lonc,
												lats=lats,
												lons=lons,
												mslp=mslp,
												pmin=pmin,
												dang=dang,
												dradius=model_res,
												search_radius=search_radius,
												model_res=model_res,
												great_circle_distance=great_circle_distance,
												filter_center_threshold=filter_center_threshold,
												search_limits=search_limits
							)
		else:

			storm_centre=True
			dmslp=pminp[nradius>great_circle_distance]

	return np.mean(dmslp), latc, lonc, pmin


def relocate_critical_centres(latc=None,
							lonc=None,
							lats=np.array(None),
							lons=np.array(None),
							mslp=np.array(None),
							pmin=None,
							dang=10,
							dradius=50,
							search_radius=750,
							model_res=20,
							great_circle_distance=650,
							filter_center_threshold=350,
							search_limits=[None,None,None,None]):

	
	latp,lonp,radius,theta=polar_cords(latc=latc,
					lonc=lonc,
					dth=math.radians(dang),
					dr=dradius,
					search_radius=search_radius)
	pointint=np.empty((len(latp.flatten()),2))
	pointint[:,0]=lonp.flatten()
	pointint[:,1]=latp.flatten()

	pointdata=np.empty((len(lats.flatten()),2))
	pointdata[:,0]=lons.flatten()
	pointdata[:,1]=lats.flatten()

	pminp=griddata(pointdata, mslp.flatten(), pointint, method='nearest')
	pminp=pminp.reshape(int(latp.shape[0]),int(latp.shape[1]))


	nradius,ntheta=np.meshgrid(radius,theta)
	npminp=np.copy(pminp)
	npminp[nradius>filter_center_threshold]=np.nan
	
	for i in range(0, nradius.shape[0]):
		for j in range(0, nradius.shape[1]):
			if np.isfinite(npminp[i,j]) and  npminp[i,j]<pmin and search_limits[1]<latp[i,j]<search_limits[3] and search_limits[0]<lonp[i,j]<search_limits[2]:
				pmin=npminp[i,j]
				latc=latp[i,j]
				lonc=lonp[i,j]
	

	return pmin, latc, lonc
	


def compute_TC_size(latc=None,
		lonc=None,
		lats=np.array(None),
		lons=np.array(None),
		u=np.array(None),
		v=np.array(None),
		dradius=100,
		dang=5,
		search_radius=2000,
		model_res=20,
		outer_wind_speed_threshold=2.5
		):
	
	latp,lonp,radius,theta=polar_cords(latc=latc,
					lonc=lonc,
					dth=math.radians(dang),
					dr=dradius,
					search_radius=search_radius)
	pointint=np.empty((len(latp.flatten()),2))
	pointint[:,0]=lonp.flatten()
	pointint[:,1]=latp.flatten()

	pointdata=np.empty((len(lats.flatten()),2))
	pointdata[:,0]=lons.flatten()
	pointdata[:,1]=lats.flatten()

	uint=griddata(pointdata, u.flatten(), pointint, method='nearest')
	uint=uint.reshape(int(latp.shape[0]),int(latp.shape[1]))

	vint=griddata(pointdata, v.flatten(), pointint, method='nearest')
	vint=vint.reshape(int(latp.shape[0]),int(latp.shape[1]))

	nradius,ntheta=np.meshgrid(radius,theta)

	ffint=np.sqrt(uint**2+vint**2)


	va=np.empty_like(vint)
	for i in range(0,uint.shape[0]):
		for j in range(0,uint.shape[1]):
			wind_vec_dir=math.atan(uint[i,j]/vint[i,j])
			va[i,j]=ffint[i,j]*math.cos(wind_vec_dir)




	radial_i=[]
	for i in range(0,nradius.shape[0]):
		leg=va[i,:]
		
		
		imax=np.where(leg==leg.max())
		imax=int(imax[0][0])
	
		nav=leg[imax:]
		radii=nradius[i,imax:]
	
		if nav.min()<=outer_wind_speed_threshold<=nav.max():
			finterpolate = interpolate.interp1d(nav, radii)
			radius_out=finterpolate(outer_wind_speed_threshold)
			radial_i=np.append(radial_i,radius_out) 
		elif nav[-1]>outer_wind_speed_threshold:
			radial_i=np.append(radial_i,radii[-1])

		else:
			ii=np.where(nav<=outer_wind_speed_threshold)

			if len(ii[0])>=1:
		
				radial_i=np.append(radial_i,radii[int(ii[0][0])])
		
	outer_size= np.mean(radial_i)
	
	

	if outer_size<=0:
		outer_size=-9999

	return outer_size	



def compute_roci_RU(latc=None,
		lonc=None,
		lats=np.array(None),
		lons=np.array(None),
		mslp=np.array(None),
		pmin=None,dang=10,
		dradius=100,
		search_radius=2000,
		model_res=20):

	
	outerp=900
	outerp_found=False
	while outerp<pmin:
		if search_radius>=1.5*dradius:
			latp,lonp,radius,theta=polar_cords(latc=latc,
						lonc=lonc,
						dth=math.radians(dang),
						dr=dradius,
						search_radius=search_radius)
			pointint=np.empty((len(latp.flatten()),2))
			pointint[:,0]=lonp.flatten()
			pointint[:,1]=latp.flatten()

			pointdata=np.empty((len(lats.flatten()),2))
			pointdata[:,0]=lons.flatten()
			pointdata[:,1]=lats.flatten()

			pminp=griddata(pointdata, mslp.flatten(), pointint, method='nearest')
			pminp=pminp.reshape(int(latp.shape[0]),int(latp.shape[1]))
		
			outsp=[]
			for i in range(0,pminp.shape[0]):
				leg=pminp[i,:]
				if len(leg)>2:
					outsp=np.append(outsp,leg[-1])
					outerp_found=True
			outerp=np.min(outsp)
			search_radius=search_radius-dradius
		else:
			outerp_found=False
			break

	
	if outerp>pmin and outerp_found==True:
		nradius,ntheta=np.meshgrid(radius,theta)
		
			
		critical_rad=[]
		critical_p=[]
		if model_res>50:
			bint=1
		else:
			bint=3

		for i in range(1,nradius.shape[0]):
			leg=pminp[i,:]
			legdr=nradius[i,:]
			check=False
			for j in range(bint,len(leg)-1):
				if  0<=leg[j]-leg[j-1]<=0.00001 and check==False:
					critical_rad=np.append(critical_rad,legdr[j])
					critical_p=np.append(critical_p,leg[j])
					check=True
			if check==False:
				critical_rad=np.append(critical_rad,legdr[-1])
				critical_p=np.append(critical_p,leg[-1])

		closedp=critical_p.min()
		radius_i=[]
		critical_lat=[]
		critical_lon=[]

		for i in range(0,nradius.shape[0]):
			leg=pminp[i,:]
			legdr=nradius[i,:]
			if closedp<leg.min():
				contour_p=leg.min()
			elif closedp>leg.max():
				contour_p=leg.max()
			else:
				contour_p=closedp
		
			if len(leg)<2:
				return np.nan,np.nan
			else:
				finterpolate = interpolate.interp1d(leg, legdr)
				radius_closedp=finterpolate(contour_p)
				radius_i=np.append(radius_i,radius_closedp)
			
			dp=np.abs(leg-closedp)
			ind=np.where(dp==dp.min())
			if len(ind[0])>1:
				ind=int(ind[0][-1])
			else:
				ind=int(ind[0])
			critical_lat=np.append(critical_lat,latp[i,ind])
			critical_lon=np.append(critical_lon,lonp[i,ind])

		SmArea=[]
		for i in range(1,len(critical_rad)):
			SmArea=np.append(SmArea,calc_area(radius_i[i-1],radius_i[i],theta[i]-theta[i-1]))

		sumaArea=np.sum(SmArea)
		roci=np.sqrt(sumaArea/np.pi)
	else:
		roci=0
		closedp=0
	if closedp<pmin:
		roci=0
		closedp=0
	
	return roci, closedp


def polar_cords(latc=None,lonc=None,dth=np.pi/256,dr=0.5,search_radius=2000):
	theta=np.arange(0,2*np.pi+dth,dth)
	radius=np.arange(0,search_radius+dr,dr)
	lat=np.empty((len(theta),len(radius)))
	lon=np.empty((len(theta),len(radius)))
	for k in range(0,len(theta)):

		for j in range(0,len(radius)):
				radii=radius[j]/111
				lon[k,j]=radii*np.cos(theta[k])+lonc
				lat[k,j]=radii*np.sin(theta[k])+latc

	return lat,lon,radius,theta

def calc_area(a=0,b=0,ang=0):
	A=(a*b*np.sin(ang))/2

	return A


def geo_distance(lat1=None, lon1=None, lat2=None, lon2=None):
	
	EARTHRADIUS = 6371 
	lat1 = 90 - lat1
	lat2 = 90 - lat2
	lon1 = 180 + lon1
	lon2 = 180 + lon2

	la1 = (np.pi / 180) * lat1
	la2 = (np.pi / 180) * lat2
	lo1 = (np.pi / 180) * lon1
	lo2 = (np.pi / 180) * lon2
	
	x1 = np.sin(la1) * np.cos(lo1)
	y1 = np.sin(la1) * np.sin(lo1)
	z1= np.cos(la1)

	x2 = np.sin(la2) * np.cos(lo2)
	y2 = np.sin(la2) * np.sin(lo2)
	z2 = np.cos(la2)

	dotp = x1 * x2 + y1 * y2 + z1 * z2
	if dotp>1:
		dotp=1
	angle = math.acos(dotp)
	dist = angle * EARTHRADIUS

	return dist

def tmp_data_loadt(data_file):
	data=np.loadtxt(data_file)
	
	tdata=open(data_file)
	tdata=tdata.readlines()
	
	if len(tdata)==1:
		data=data.reshape(1,len(data))
		
	return data

def get_bearing(lat1, lon1, lat2, lon2):
	dLon = (lon2 - lon1)
	x = math.cos(math.radians(lat2)) * math.sin(math.radians(dLon))
	y = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dLon))
	brng = np.arctan2(x,y)
	brng = np.degrees(brng)
	brng = (brng +360) % 360
	return brng


def compute_VT_series(dates=np.array([None]),
					hours=np.array([None]),
					listlev1=np.array([None]),
					listlev2=np.array([None]),
					liste_lat=np.array([None]),
					liste_lon=np.array([None]),
					max_dist=500,
					lats=np.array([None]),
					lons=np.array([None]),
					levels=np.array([None]),
					Zvar=np.array([None]),
					vtl_vtu_lr=False,
					):
	
	listlev11=[n * 100 for n in listlev1]
	lnP1=np.log(listlev11)
	listlev22=[n * 100 for n in listlev2]
	lnP2=np.log(listlev22)
	
	VTL_seriesb=[]
	VTU_seriesb=[]
	
	VTL_series=[]
	VTU_series=[]

	for ix in range(0,lons.shape[0]):
		for jx in range(0,lons.shape[1]):
			if lons[ix,jx]>=180:
				lons[ix,jx]=lons[ix,jx]-360


	for i in range(0,len(dates)):

		clat = liste_lat[i]
		clon = liste_lon[i]
		
		if clon>=180:
			clon=clon-360
		
		
		
		deltaZ1=[]
		deltaZ2=[]

		distance=haversine(lons, lats, clon, clat)
		for plev in listlev1:

			ilev = list(levels).index(plev)
			
			data_z1 = Zvar[i,ilev,:,:]
			data_z1_m = np.ma.masked_where(distance > max_dist, data_z1)
			zmin=np.min(data_z1_m)
			zmax=np.max(data_z1_m)
			delta_z1=zmax-zmin
			deltaZ1.append(int(delta_z1))

		for plev in listlev2:

			ilev = list(levels).index(plev)
			data_z2 = Zvar[i,ilev,:,:]

			data_z2_m = np.ma.masked_where(distance > max_dist, data_z2)
			zmin=np.min(data_z2_m)
			zmax=np.max(data_z2_m)
			delta_z2=zmax-zmin
			deltaZ2.append(int(delta_z2))


		coef1b=(deltaZ1[-1]-deltaZ1[0])/(lnP1[-1]-lnP1[0])
		VTL_seriesb.append(int(coef1b))

		coef2b=(deltaZ2[-1]-deltaZ2[0])/(lnP2[-1]-lnP2[0])
		VTU_seriesb.append(int(coef2b))
		
		X1 = np.reshape(lnP1, (len(lnP1), 1))
		y1 =  deltaZ1
		model = LinearRegression().fit(X1, y1)
		model.fit(X1, y1)
		yhat1 = model.predict(X1)
		reg = LinearRegression().fit(X1, y1)
		coef1=reg.coef_
		VTL_series.append(int(coef1))
		
		
		X2 = np.reshape(lnP2, (len(lnP2), 1))
		y2 =  deltaZ2
		model = LinearRegression().fit(X2, y2)
		model.fit(X2, y2)
		yhat2 = model.predict(X2)
		reg = LinearRegression().fit(X2, y2)
		coef2=reg.coef_
		VTU_series.append(int(coef2))
		
	if vtl_vtu_lr:
		return VTL_series, VTU_series
	else:
		return VTL_seriesb, VTU_seriesb


def compute_Bparameter_hart(cenlats=np.array([None]),
			cenlons=np.array([None]),
			max_dist=500,
			dates=np.array([None]),
			hours=np.array([None]),
			lats=np.array([None]),
			lons=np.array([None]),
			levels=np.array([None]),
			Zvar=np.array([None])
			):

	for ix in range(0,lons.shape[0]):
		for jx in range(0,lons.shape[1]):
			if lons[ix,jx]>=180:
				lons[ix,jx]=lons[ix,jx]-360


	ilev_top = list(levels).index(600)
	ilev_bot = list(levels).index(900)

	B_series=[]
	for i in range(0,len(cenlats)-1):
		
		for ix in range(0,lons.shape[0]):
			for jx in range(0,lons.shape[1]):
				if lons[ix,jx]>=180:
					lons[ix,jx]=lons[ix,jx]-360
		
		ilev_top = list(levels).index(600)
		ilev_bot = list(levels).index(900)
		
		#thickness=zvar[ilev_top,:]-zvar[ilev_bot,:]
		
		
		
		data=Zvar[i,ilev_top,:]-Zvar[i,ilev_bot,:]
	
		clat = cenlats[i]
		clon = cenlons[i]
		
			
		if clon>=180:
			clon=clon-360
		
		qq_ang_all=np.zeros([lons.shape[0],lons.shape[1]],dtype='f')
		
		for x in range(0,lons.shape[0]):
			for y in range(0,lons.shape[1]):
				qq_ang_all[x,y]=get_bearing(clat, clon, lats[x,y], lons[x,y])
				
		distance=haversine(lons, lats, clon, clat)
		ang=get_bearing(clat, clon, cenlats[i+1], cenlons[i+1])
		Zl=np.zeros([lons.shape[0],lons.shape[1]],dtype='f')
		Zr=np.zeros([lons.shape[0],lons.shape[1]],dtype='f')

		for jlat in range(0,lons.shape[0]-1):
			for jlon in range(0,lons.shape[1]-1):
				if qq_ang_all[jlat,jlon] == ang:
					Zl[jlat,jlon] = 0
					Zr[jlat,jlon] = 0
				elif (ang >= 0 and ang < 180):
					if (qq_ang_all[jlat,jlon] > ang and qq_ang_all[jlat,jlon]  < ang+180):
						Zl[jlat,jlon] = 0
						Zr[jlat,jlon] = data[jlat,jlon]
					else:
						Zr[jlat,jlon] = 0
						Zl[jlat,jlon] = data[jlat,jlon]
				elif (ang >= 180 and ang < 360):
					if (qq_ang_all[jlat,jlon] > ang-180 and qq_ang_all[jlat,jlon] < ang):
						Zr[jlat,jlon] = 0
						Zl[jlat,jlon] = data[jlat,jlon]
					else:
						Zl[jlat,jlon] = 0
						Zr[jlat,jlon] = data[jlat,jlon]


		
		Zr_m = np.ma.masked_where(distance > max_dist, Zr)
		Zl_m = np.ma.masked_where(distance > max_dist, Zl)
		
		Zr_m[Zr_m == 0] = np.nan
		Zl_m[Zl_m == 0] = np.nan
		
		#Zr_m = Zr[distance<max_dist]
		#Zl_m = Zl[distance<max_dist]

		#Zr_m=Zr_m[Zr_m != 0]
		#Zl_m=Zl_m[Zl_m != 0]

		#Zr_mean=np.nanmean(Zr_m)
		#Zl_mean=np.nanmean(Zl_m)
		#B=Zr_mean-Zl_mean


		Zr_m = xr.DataArray(Zr_m)
		Zl_m = xr.DataArray(Zl_m)
		weights = xr.DataArray(np.cos(np.deg2rad(Zr_m.shape[0])))
		Zr_weighted = Zr_m.weighted(weights)
		Zl_weighted = Zl_m.weighted(weights)
		Zr_weighted_mean=Zr_weighted.mean()
		Zl_weighted_mean=Zl_weighted.mean()
		
			
		B=(Zr_weighted_mean-Zl_weighted_mean)
		
		if not np.isnan(B):
			B=int(B)
	
		if clat<0:
			B=-B
		B_series=np.append(B_series,B)

	return B_series

def get_VTU_VTL_B_series(cenlats=np.array([None]),
						cenlons=np.array([None]),
						dates=np.array([None]),
						hours=np.array([None]),
						idir_upper="./",
						source_upperprefix="",
						source="",
						era_date_file_name="",
						search_limits=[None,None,None,None],
						search_region="",
						max_dist=500,
						vtl_vtu_lr=False,
						custom_geopotential_var_name="",
						custom_upper_level_variable_name="",
						varlat="",
						varlon="",
						custom_date_file_name=""
						):
	
	if vtl_vtu_lr:
		listlev1=[900,850,800,750,700,650,600]
		listlev2=[600,550,500,450,400,350,300]
	else:
		listlev1=[900,600]
		listlev2=[600,300]
	
	
	if source.upper()=="ERA5":
		upperfiles=get_era5_files(dates=dates,hours=hours,era_file_prefix=source_upperprefix,era_date_file_name=era_date_file_name)

		sourceZ, sourcelats,sourcelons, source_levels=get_era5_3dvar(idir=idir_upper,
							erafile=upperfiles[0],
							svariable='z',
							varlevel="level",
							search_limits=search_limits,
							search_region=search_region,
							fdate=dates[0]+hours[0],
							full=False
							)


		Zvar=np.empty((len(dates), sourceZ.shape[0], sourceZ.shape[1], sourceZ.shape[2]))
		for idx in range(0, len(dates)):
			fdate=dates[idx]+hours[idx]
			sourceZ, sourcelats,sourcelons, source_levels=get_era5_3dvar(idir=idir_upper,
									erafile=upperfiles[idx],
									svariable='z',
									varlevel="level",
									search_limits=search_limits,
									search_region=search_region,
									fdate=fdate,
									full=False
									 )
		

			Zvar[idx,:,:,:]=sourceZ
			
				
		checklevs=[]
		for plev in listlev1:
			checklevs=np.append(checklevs,plev)
		for plev in listlev2:
			checklevs=np.append(checklevs,plev)

		for ulev in checklevs:
			if not int(ulev) in source_levels:
				if vtl_vtu_lr:
					print_error_message(" Mandatory pressure level = " + str(int(ulev)) +" hPa is missing in input file.\nIf checking_upper_levels_parameters = True and vtl_vtu_lr=True, mandatory pressure levels are from 900 to 300 hPa every 50 hPa\nPlease change or remove your input files directory.\nCyTRACK will download upper levels input files for you, or set checking_upper_levels_parameters = False")
				else:	
					print_error_message(" Mandatory pressure level = " + str(int(ulev)) +" hPa is missing in input file.\nIf checking_upper_levels_parameters = True, mandatory pressure levels are [900,600,300] hPa\nPlease change or remove your input files directory.\nCyTRACK will download upper levels input files for you, or set checking_upper_levels_parameters = False")
							
				
	elif source.upper()=="WRF":
		plevels=np.arange(1000,150,-50)
		upperfiles=get_wrf_files(dates=dates,hours=hours,wrfprefix=source_upperprefix)
		sourcelats,sourcelons,sourceZ, source_levels = get_wrf_uppervar(idir=idir_upper,
																wrffile=upperfiles[0],
																variables=["PH"],
																plevels=plevels)
	
		Zvar=np.empty((len(dates), sourceZ.shape[0], sourceZ.shape[1], sourceZ.shape[2]))
		for idx in range(0, len(dates)):
			sourcelats,sourcelons,sourceZ, source_levels = get_wrf_uppervar(idir=idir_upper,
																wrffile=upperfiles[idx],
																variables=["PH"],
																plevels=plevels)
			Zvar[idx,:,:,:]=sourceZ
					
		
	elif  source.upper()=="CUSTOM":
		
		upperfiles=get_custom_files(dates=dates,hours=hours,custom_file_prefix=source_upperprefix,custom_date_file_name=custom_date_file_name)

		sourceZ, sourcelats,sourcelons, source_levels=get_custom_3dvar(idir=idir_upper,
							customfile=upperfiles[0],
							svariable=custom_geopotential_var_name,
							varlevel=custom_upper_level_variable_name,
							search_limits=search_limits,
							search_region=search_region,
							fdate=dates[0]+hours[0],
							full=False,
							varlat=varlat,
							varlon=varlon,
							)


		Zvar=np.empty((len(dates), sourceZ.shape[0], sourceZ.shape[1], sourceZ.shape[2]))
		for idx in range(0, len(dates)):
			fdate=dates[idx]+hours[idx]
			sourceZ, sourcelats,sourcelons, source_levels=get_custom_3dvar(idir=idir_upper,
									customfile=upperfiles[idx],
									svariable=custom_geopotential_var_name,
									varlevel=custom_upper_level_variable_name,
									search_limits=search_limits,
									search_region=search_region,
									fdate=fdate,
									full=False,
									varlat=varlat,
									varlon=varlon,
									 )
		

			Zvar[idx,:,:,:]=sourceZ
			
				
		checklevs=[]
		for plev in listlev1:
			checklevs=np.append(checklevs,plev)
		for plev in listlev2:
			checklevs=np.append(checklevs,plev)

		for ulev in checklevs:
			if not int(ulev) in source_levels:
				if vtl_vtu_lr:
					print_error_message(" Mandatory pressure level = " + str(int(ulev)) +" hPa is missing in input file.\nIf checking_upper_levels_parameters = True and vtl_vtu_lr=True, mandatory pressure levels are from 900 to 300 hPa every 50 hPa\nPlease change or remove your input files directory.")
				else:	
					print_error_message(" Mandatory pressure level = " + str(int(ulev)) +" hPa is missing in input file.\nIf checking_upper_levels_parameters = True, mandatory pressure levels are [900,600,300] hPa\nPlease change or remove your input files directory.")
	
	
		
		
	
	VTL, VTU = compute_VT_series(dates=dates,
				hours=hours,
				listlev1=listlev1,
				listlev2=listlev2,
				liste_lat=cenlats,
				liste_lon=cenlons,
				max_dist=max_dist,
				lats=sourcelats,
				lons=sourcelons,
				levels=source_levels,
				Zvar=Zvar,
				vtl_vtu_lr=vtl_vtu_lr
				)
	
	Bhart=compute_Bparameter_hart(cenlats=cenlats,
			cenlons=cenlons,
			max_dist=max_dist,
			dates=dates,
			hours=hours,
			lats=sourcelats,
			lons=sourcelons,
			levels=source_levels,
			Zvar=Zvar
			)
	#plot_VTL(dates, hours, VTL[:-1], Bhart,path="./", fname="VTL-B", dpi=600)
	#plot_VTU(dates, hours, VTL[:-1], VTU[:-1], Bhart, path="./", fname="VTU-VTL", dpi=600)
	
	
	return VTU, VTL, Bhart

	
def plot_VTL(dates,hours, VTL_series, B_series,path="./", fname="900-600hPa_Thermal_wind,png", dpi=600):
	fig = plt.figure(figsize=(18., 12.))
	ax = fig.add_subplot(111)
	ax.set_title(dates[0][0:13]+'-'+dates[-1][0:13],loc='right',fontsize=25)
	plt.axhline(y=10, color='grey', linewidth=10,zorder=1)
	plt.axvline(x=0, color='grey', linewidth=10,zorder=1)
	plt.plot(VTL_series, B_series,linewidth=2, color='k')
	plt.scatter(VTL_series,  B_series, c='k',zorder=2.5)
	plt.scatter(VTL_series[0], B_series[0], c='green',zorder=40.5)
	plt.scatter(VTL_series[-1], B_series[-1], c='red',zorder=40.5)
	plt.xlabel("$-V_T^L$ [900-600hPa Thermal wind]", fontsize = 25)
	plt.ylabel("B [900-600hPa thickness symmetry]", fontsize = 25)

	plt.text(-660,-15,'Symmetric', rotation=90., color='red', fontsize = 25)
	plt.text(-660, 60,'Asymmetric', rotation=90., color='red', fontsize = 25)
	plt.text(-450, -30,'Cold core', rotation=0., color='blue', fontsize = 25)
	plt.text(280, -30,'Warm core', rotation=0., color='red', fontsize = 25)
	plt.xlim(-600, 600)
	plt.ylim(-20, 125)
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)

	
	plt.text(VTL_series[0], B_series[0], dates[0] + " "+hours[0], fontsize=20)
	plt.text(VTL_series[-1], B_series[-1], dates[-1] + " "+hours[-1], fontsize=20)
		
	xrange = [(-600, 600)]
	yrange1 = (-20, 30)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='blue', alpha=0.2)
	ax.add_collection(c1)
	xrange = [(-600, 600)]
	yrange1 = (10, 130)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='blue', alpha=0.2)
	ax.add_collection(c1)
	xrange = [(0, 600)]
	yrange1 = (-20, 30)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='red', alpha=0.2)
	ax.add_collection(c1)
	xrange = [(0, 600)]
	yrange1 = (10, 130)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='red', alpha=0.2)
	ax.add_collection(c1)
	plt.tight_layout()
	plt.savefig(path+"/"+fname,dpi=dpi)
	
	
def plot_VTU(dates, hours, VTL_series, VTU_series, B_series, path="./", fname="900-600hPa_Thermal_wind,png", dpi=600):
	fig = plt.figure(figsize=(18., 12.))
	ax = fig.add_subplot(111)
	ax.set_title(dates[0][0:13]+'-'+dates[-1][0:13],loc='right',fontsize=25)
	plt.axhline(y=0, color='grey', linewidth=10,zorder=1)
	plt.axvline(x=0, color='grey', linewidth=10,zorder=1)
	plt.plot(VTL_series, VTU_series, linewidth=2, color='k')
	plt.scatter(VTL_series, VTU_series, c='k',zorder=2.5)
	plt.scatter(VTL_series[0], VTU_series[0], c='green',zorder=40.5)
	plt.scatter(VTL_series[-1], VTU_series[-1], c='red',zorder=40.5)


	plt.xlabel("$-V_T^L$ [900-600hPa Thermal wind]", fontsize = 25)
	plt.ylabel("$-V_T^U$ [600-300hPa Thermal wind]", fontsize = 25)
	plt.text(-350, -580,'Deep cold core', rotation=0., color='blue', fontsize =25, fontweight="bold")
	plt.text(210, 570,'Deep warm core', rotation=0., color='red', fontsize = 25, fontweight="bold")
	plt.text(210, -580,'Shallow warm core', rotation=0., color='orange', fontsize = 25, fontweight="bold")

	plt.xlim(-600, 600)
	plt.ylim(-600, 600)
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)

	plt.text(VTL_series[0], VTU_series[0], dates[0] + " "+hours[0], fontsize=20)
	plt.text(VTL_series[-1], VTU_series[-1], dates[-1] + " "+hours[-1], fontsize=20)
	
	xrange = [(0, 600)]
	yrange1 = (0, 600)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='r', alpha=0.2)
	ax.add_collection(c1)
	xrange = [(-600, 600)]
	yrange1 = (-600, 600)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='b', alpha=0.2)
	ax.add_collection(c1)
	xrange = [(-600, 600)]
	yrange1 = (0, 600)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='orange', alpha=0.2)
	ax.add_collection(c1)
	xrange = [(0, 600)]
	yrange1 = (-600, 600)
	c1 = collections.BrokenBarHCollection(xrange, yrange1, facecolor='orange', alpha=0.2)
	ax.add_collection(c1)

	plt.tight_layout()
	plt.savefig(path+"/"+fname,dpi=dpi)



def paring_centers(cyclone_type="",
		dates=[None], 
		hours=[None],
		dist_threshold=1000,
		tmpdir="."+program_name()+"_tmpdir/", 
		pathoutput="./",
		search_region="NA",
		dt_h=6,
		source="",
		dt_lifetime=24,
		minimum_distance_travelled=1000,
		intensity_threshold=17.5,
		checking_upper_levels_parameters=False,
		idir_upper="./",
		source_upperprefix="",
		era_date_file_name="",
		search_limits=[None,None,None,None],
		max_dist=500,
		vtl_vtu_lr=False,
		core_criteria_length=0,
		VTL_threshold=0,
		VTU_threshold=0,
		Bhart_threshold=10,
		custom_geopotential_var_name="",
		custom_upper_level_variable_name="",
		custom_varlat="",
		custom_varlon="",
		custom_date_file_name=""
		):
	
	
	
	
	search_region=search_region[0:2]
	sys_id=0
	fwrite=open(pathoutput+"/"+program_name()+"_"+search_region+"_"+dates[0]+hours[0]+"-"+dates[-1]+hours[-1]+"_"+source+"_"+cyclone_type+".dat","w")
	print("")
	ProgressBar(0, len(dates), prefix = '	Progress', suffix = '', decimals = 1, length = 80, fill = '+', printEnd = "\r")
	for index in range(0,len(dates)):
		fdate=dates[index]+hours[index]
		if not os.stat(tmpdir+"/critical_centers_"+fdate+".dat").st_size == 0:
			data=tmp_data_loadt(tmpdir+"/critical_centers_"+fdate+".dat")
					
			for i in range(0,data.shape[0]):
				latc=data[i,0]
				lonc=data[i,1]
				
				if lonc>180:
					lonc=lonc-360
				
				if data[i,-1]==1 and data[i,4]!=0:
					nlats=[latc]
					nlons=[lonc]
					nroci=[data[i,5]]
					npmin=[data[i,2]]
					nmws=[data[i,3]]
					nclosedp=[data[i,4]]
					ndate=[dates[index]]
					nhour=[hours[index]]
					nouter_r=[data[i,6]]
					j=index+1

					track_end=False
					check_j=np.copy(j)
					while j < len(dates):
						
						nfdate=dates[j]+hours[j]
						if os.stat(tmpdir+"/critical_centers_"+nfdate+".dat").st_size == 0:
							if (len(nlats)-1)*dt_h>=dt_lifetime and nmws.max()>intensity_threshold:
								
																
								if  checking_upper_levels_parameters:
									VTU, VTL, Bhart=get_VTU_VTL_B_series(cenlats=nlats,
																			cenlons=nlons,
																			dates=ndate,
																			hours=nhour,
																			idir_upper=idir_upper,
																			source_upperprefix=source_upperprefix,
																			source=source,
																			era_date_file_name=era_date_file_name,
																			search_limits=search_limits,
																			search_region=search_region,
																			max_dist=max_dist,
																			vtl_vtu_lr=vtl_vtu_lr,
																			custom_geopotential_var_name=custom_geopotential_var_name,
																			custom_upper_level_variable_name=custom_upper_level_variable_name,
																			varlat=custom_varlat,
																			varlon=custom_varlon,
																			custom_date_file_name=custom_date_file_name
																			)
									Bhart=np.append(Bhart,-99999)
									cyclone_phase=get_cyclone_phase(VTU=VTU,VTL=VTL, B=Bhart)
									
									check_cyclone=get_cyclone_type(VTU=VTU,
																	VTL=VTL, 
																	B=Bhart, 
																	cyclone_type=cyclone_type,
																	core_criteria_length=core_criteria_length,
																	dt_h=dt_h,
																	VTL_threshold=VTL_threshold,
																	VTU_threshold=VTU_threshold,
																	Bhart_threshold=Bhart_threshold,)
								
									
																	
									
									if check_cyclone:
																		
										sys_id=sys_id+1
										write_tracker(cyclone_type=cyclone_type,nlats=nlats,
												nlons=nlons,
												npmin=npmin,
												nmws=nmws,
												nclosedp=nclosedp,
												nroci=nroci,
												nouter_r=nouter_r,
												nctype=cyclone_phase,
												VTU=VTU,
												VTL=VTL,
												Bhart=Bhart,
												search_region=search_region,
												sys_id=sys_id,
												fwrite=fwrite,
												ndates=ndate,
												nhours=nhour,
												)
										nlats=[]
										nlons=[]
										nroci=[]
										npmin=[]
										nmws=[]
										nclosedp=[]
										nouter_r=[]
										nctype=[]
										ndate=[]
										nhour=[]
										track_end=True
																			
									
									
								else:
									sys_id=sys_id+1
									
									nctype=np.empty(len(nlats))
									nctype[:]=-99999
									
									VTL=np.empty_like(nlats)
									VTL[:]=-99999
									
									VTU=np.empty_like(nlats)
									VTU[:]=-99999
									
									Bhart=np.empty_like(nlats)
									Bhart[:]=-99999
									
									
									write_tracker(cyclone_type=cyclone_type,
											nlats=nlats,
											nlons=nlons,
											npmin=npmin,
											nmws=nmws,
											nclosedp=nclosedp,
											nroci=nroci,
											nouter_r=nouter_r,
											nctype=nctype,
											VTU=VTU,
											VTL=VTL,
											Bhart=Bhart,
											search_region=search_region,
											sys_id=sys_id,
											fwrite=fwrite,
											ndates=ndate,
											nhours=nhour,
											)
								nlats=[]
								nlons=[]
								nroci=[]
								npmin=[]
								nmws=[]
								nclosedp=[]
								nouter_r=[]
								nctype=[]
								ndate=[]
								nhour=[]
								track_end=True
						else:
							ndata=tmp_data_loadt(tmpdir+"/critical_centers_"+nfdate+".dat")
							k_paring=None
							back_dist=10000
							for k in range(0,ndata.shape[0]):
								if ndata[k,1]>180:
									ndata[k,1]=ndata[k,1]-360
								if ndata[k,-1]==1 and data[i,4]!=0:
									pdist=geo_distance(lat1=latc, lon1=lonc, lat2=ndata[k,0], lon2=ndata[k,1])
									if pdist<=dist_threshold and pdist<back_dist:
										k_paring=k
										back_dist=np.copy(pdist)
							
							if k_paring!=None:
								nlats=np.append(nlats, ndata[k_paring,0])
								nlons=np.append(nlons, ndata[k_paring,1])
								nroci=np.append(nroci, ndata[k_paring,5])
								npmin=np.append(npmin, ndata[k_paring,2])
								nmws=np.append(nmws, ndata[k_paring,3])
								nclosedp=np.append(nclosedp, ndata[k_paring,4])
								nouter_r=np.append(nouter_r,ndata[k_paring,6])
								ndate=np.append(ndate,dates[j])
								nhour=np.append(nhour,hours[j])
								latc=ndata[k_paring,0]
								lonc=ndata[k_paring,1]
								ndata[k_paring,-1]=0
								check_j=np.copy(j)
								np.savetxt(tmpdir+"/critical_centers_"+nfdate+".dat",ndata)
								track_end=False
							elif j<len(dates)-1 and j-check_j==1:
								nfdate2=dates[j+1]+hours[j+1]
								ndata2=tmp_data_loadt(tmpdir+"/critical_centers_"+nfdate2+".dat")
								k_paring2=None
								back_dist2=10000
								latc=nlats[-1]
								lonc=nlons[-1]
								for k in range(0,ndata2.shape[0]):
									if ndata2[k,1]>180:
										ndata2[k,1]=ndata2[k,1]-360
									if ndata2[k,-1]==1 and data[i,4]!=0:
										pdist2=geo_distance(lat1=latc, lon1=lonc, lat2=ndata2[k,0], lon2=ndata2[k,1])
										if pdist2<=2*dist_threshold and pdist2<back_dist2:
											k_paring2=k
											back_dist2=np.copy(pdist2)

								if k_paring2!=None and len(nlons)>=2:
									m1=(nlons[-2]-nlons[-1])/((nlons[-2]-nlons[-1]))
									m2=(ndata2[k_paring2,1] - nlons[-1])/(ndata2[k_paring2,0] - nlats[1])
									angle=np.arctan(np.abs((m2-m1)/(1+m1*m2)))
								else:
									angle=np.pi/2


								if angle<0.01:
									nextlat=ndata2[k_paring2,0]
									nextlon=ndata2[k_paring2,1]
									nextroci=ndata2[k_paring2,5]
									nextpmin=ndata2[k_paring2,2]
									nextmws=ndata2[k_paring2,3]
									nextclosedp=ndata2[k_paring2,4]
									nextouter_r=ndata2[k_paring2,6]

									nlats=np.append(nlats, (nlats[-1]+nextlat)/2)
									nlons=np.append(nlons, (nlons[-1]+nextlon)/2)
									nroci=np.append(nroci, (nroci[-1]+nextroci)/2)
									npmin=np.append(npmin, (npmin[-1]+nextpmin)/2)
									nmws=np.append(nmws, (nmws[-1]+nextmws)/2)
									nclosedp=np.append(nclosedp, (nclosedp[-1]+nextclosedp)/2)
									nouter_r=np.append(nouter_r,(nouter_r[-1]+nextouter_r)/2)
									ndate=np.append(ndate,dates[j])
									nhour=np.append(nhour,hours[j])
									latc=(nlats[-1]+nextlat)/2
									lonc=(nlons[-1]+nextlon)/2
									track_end=False
							elif (len(nlats)-1)*dt_h>=dt_lifetime and nmws.max()>intensity_threshold:
								check_dist=[]
								for ilat in range(1,len(nlats)):
									aux_dist=geo_distance(lat1=nlats[ilat-1],
											lon1=nlons[ilat-1],
											lat2=nlats[ilat],
											lon2=nlons[ilat])
									check_dist=np.append(check_dist,aux_dist)
								if np.sum(check_dist)>=minimum_distance_travelled:
									
									if  checking_upper_levels_parameters: 
										
										VTU, VTL, Bhart=get_VTU_VTL_B_series(cenlats=nlats,
																			cenlons=nlons,
																			dates=ndate,
																			hours=nhour,
																			idir_upper=idir_upper,
																			source_upperprefix=source_upperprefix,
																			source=source,
																			era_date_file_name=era_date_file_name,
																			search_limits=search_limits,
																			search_region=search_region,
																			max_dist=max_dist,
																			vtl_vtu_lr=vtl_vtu_lr,
																			custom_geopotential_var_name=custom_geopotential_var_name,
																			custom_upper_level_variable_name=custom_upper_level_variable_name,
																			varlat=custom_varlat,
																			varlon=custom_varlon,
																			custom_date_file_name=custom_date_file_name
																			)
										Bhart=np.append(Bhart,-99999)
										cyclone_phase=get_cyclone_phase(VTU=VTU,VTL=VTL,B=Bhart)
										check_cyclone=get_cyclone_type(VTU=VTU,
																		VTL=VTL, 
																		B=Bhart, 
																		cyclone_type=cyclone_type,
																		core_criteria_length=core_criteria_length,
																		dt_h=dt_h,
																		VTL_threshold=VTL_threshold,
																		VTU_threshold=VTU_threshold,
																		Bhart_threshold=Bhart_threshold,)
																		
										if check_cyclone:
																		
											sys_id=sys_id+1
											write_tracker(cyclone_type=cyclone_type,nlats=nlats,
													nlons=nlons,
													npmin=npmin,
													nmws=nmws,
													nclosedp=nclosedp,
													nroci=nroci,
													nouter_r=nouter_r,
													nctype=cyclone_phase,
													VTU=VTU,
													VTL=VTL,
													Bhart=Bhart,
													search_region=search_region,
													sys_id=sys_id,
													fwrite=fwrite,
													ndates=ndate,
													nhours=nhour,
													)
											track_end=True
																		
										
									else:
										sys_id=sys_id+1
										
										nctype=np.empty(len(nlats))
										nctype[:]=-99999
										
										VTL=np.empty_like(nlats)
										VTL[:]=-99999
										
										VTU=np.empty_like(nlats)
										VTU[:]=-99999
										
										Bhart=np.empty_like(nlats)
										Bhart[:]=-99999
										
										
										write_tracker(cyclone_type=cyclone_type,nlats=nlats,
												nlons=nlons,
												npmin=npmin,
												nmws=nmws,
												nclosedp=nclosedp,
												nroci=nroci,
												nouter_r=nouter_r,
												nctype=nctype,
												VTU=VTU,
												VTL=VTL,
												Bhart=Bhart,
												search_region=search_region,
												sys_id=sys_id,
												fwrite=fwrite,
												ndates=ndate,
												nhours=nhour,
												)
									track_end=True
								else:
									track_end=True
							else:
								track_end=True


						if track_end:
							j=len(dates)
							nlats=[]
							nlats=[]
							nroci=[]
							npmin=[]
							nmws=[]
							nouter_r=[]
							nclosedp=[]
							ndate=[]
							nhour=[]
							nctype=[]
						else:
							j=j+1
							
		time.sleep(0.2)
		ProgressBar(index+1, len(dates), prefix = '	Progress', suffix = '', decimals = 1, length = 80, fill = '+', printEnd = "\r")
	print("")
	return sys_id



def write_tracker(cyclone_type,nlats,nlons,npmin,nmws,nclosedp,nroci,nouter_r,nctype,VTU,VTL, Bhart, search_region,sys_id,fwrite,ndates, nhours):
	sp="          "

	nnctype=[]	
	for i in range(0, len(nmws)):
		if nmws[i]>0:
			nmws[i]=nmws[i]*3.6
		if 	nctype[i]==-99999:
			nnctype=np.append(nnctype, "UDCC")
		else:
			nnctype=np.append(nnctype, nctype[i])
			
	nclosedp[nclosedp==0]=-99999
	nroci[nroci==0]=-99999
	Bhart[np.isnan(Bhart)]==-99999
	
	if cyclone_type.upper() in ("EC"):
		line="Cy"+search_region+str(sys_id).zfill(4)+ndates[0][0:4]+", "+str(len(nlats))+","
		fwrite.write(line+"\n")
		outer_r=np.empty_like(nroci)
		nouter_r[:]=-99999
		for k in range(0,len(nlats)):
			
			#line=ndates[k]+", "+nhours[k].zfill(2)+","+sp[:-len(str(nlats[k])[0:5])]+ str(nlats[k])[0:5]+","+sp[:-len(str(nlons[k])[0:6])]+ str(nlons[k])[0:6]+","+sp[:-len(str(npmin[k])[0:6])]+ str(npmin[k])[0:7]+","+sp[:-len(str(nmws[k])[0:5])]+ str(nmws[k])[0:5]+","+sp[:-len(str(nroci[k])[0:7])]+ str(nroci[k])[0:7]+","+sp[:-len(str(nclosedp[k])[0:7])]+ str(nclosedp[k])[0:7]+",   " + nnctype[k] +","+ sp[:-len(str(VTU[k])[0:5])]+ str(VTU[k])[0:5]+","+ sp[:-len(str(VTL[k])[0:5])]+ str(VTL[k])[0:5]+","+ sp[:-len(str(Bhart[k])[0:5])]+ str(Bhart[k])[0:5]+","
			line=ndates[k]+", "+nhours[k].zfill(2)+","+sp[:-len(str(nlats[k])[0:5])]+ str(nlats[k])[0:5]+","+sp[:-len(str(nlons[k])[0:6])]+ str(nlons[k])[0:6]+","+sp[:-len(str(npmin[k])[0:6])]+ str(npmin[k])[0:7]+","+sp[:-len(str(nmws[k])[0:5])]+ str(nmws[k])[0:5]+","+sp[:-len(str(nouter_r[k])[0:7])]+ str(nouter_r[k])[0:7]+","+sp[:-len(str(nclosedp[k])[0:7])]+ str(nclosedp[k])[0:7]+","+sp[:-len(str(nroci[k])[0:7])]+ str(nroci[k])[0:7]+",   " + str(nnctype[k]) +","+ sp[:-len(str(VTU[k])[0:6])]+ str(VTU[k])[0:6]+","+ sp[:-len(str(VTL[k])[0:6])]+ str(VTL[k])[0:6]+","+ sp[:-len(str(Bhart[k])[0:6])]+ str(Bhart[k])[0:6]+","
			fwrite.write(line+"\n")
	elif cyclone_type.upper() in ("TC","TLC","SC","MC"):
		line="Cy"+search_region+str(sys_id).zfill(4)+ndates[0][0:4]+", "+str(len(nlats))+","
		fwrite.write(line+"\n")
		
		for k in range(0,len(nlats)):
			line=ndates[k]+", "+nhours[k].zfill(2)+","+sp[:-len(str(nlats[k])[0:5])]+ str(nlats[k])[0:5]+","+sp[:-len(str(nlons[k])[0:6])]+ str(nlons[k])[0:6]+","+sp[:-len(str(npmin[k])[0:6])]+ str(npmin[k])[0:7]+","+sp[:-len(str(nmws[k])[0:5])]+ str(nmws[k])[0:5]+","+sp[:-len(str(nouter_r[k])[0:7])]+ str(nouter_r[k])[0:7]+","+sp[:-len(str(nclosedp[k])[0:7])]+ str(nclosedp[k])[0:7]+","+sp[:-len(str(nroci[k])[0:7])]+ str(nroci[k])[0:7]+",   " + str(nnctype[k]) +","+ sp[:-len(str(VTU[k])[0:6])]+ str(VTU[k])[0:6]+","+ sp[:-len(str(VTL[k])[0:6])]+ str(VTL[k])[0:6]+","+ sp[:-len(str(Bhart[k])[0:6])]+ str(Bhart[k])[0:6]+","
			fwrite.write(line+"\n")



