import numpy as np
import matplotlib.pylab as plt
import sys
import alarconpy as al
import os
from shapely.geometry import Polygon, MultiPolygon, Point 
import requests



def create_map(search_limits=[None, None, None, None]):
	from cartopy import config
	from cartopy.util import add_cyclic_point
	import cartopy.feature as cfeature
	import cartopy.crs as ccrs
	from cartopy.mpl.geoaxes import GeoAxes
	from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
	import matplotlib.ticker as mticker	
	import math	
	
	min_lon,min_lat=(search_limits[0],search_limits[1])
	max_lon,max_lat=(search_limits[2],search_limits[3])
	paso_h=25
	

	
	crs = ccrs.PlateCarree()
	mapa=plt.subplot(1,1,1,projection=ccrs.PlateCarree(0) )
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
	gl.xlabel_style = {'size': 30, 'color': 'black'}
	gl.ylabel_style = {'size': 30,'color': 'black'}


	return mapa



def get_basin_limits():
	NATL = Polygon(((260, 40), (345, 40), (345, 0), (295, 0), (260, 20)))
	return NATL

def index_row(myList, v):
	j=[]
	for i, x in enumerate(myList):
     
		if v in x:
			j.append(i)
	return j



plt.figure(figsize=(18,12))
mapa=create_map(search_limits=[-100, -5, 0, 60])


basin_limits=get_basin_limits()

cfile=open("CyTRACK_output/CyTRACK_AL_2018090100-2018093018_ERA5_TC.dat")
cfile=cfile.readlines()

lenght=len(cfile)

index=0
cont=0
while index<lenght:
	line_data=cfile[index].split(",")
	
	diff=int(line_data[1].split("\n")[0])
	lats=[]
	lons=[]
	
	for i in range(index+1,index+diff+1):
		line_data=cfile[i]
		n_data=line_data.split(",")
		
		lats=np.append(lats,float(n_data[2]))
		lons=np.append(lons,float(n_data[3]))
		
	
	latg=lats[0]
	long=lons[0]
	if long<0:
		long=long+360

	point=Point(long, latg)
	if point.within(basin_limits):
		mapa.plot(lons,lats, color="r", linewidth=2.5, marker='o')
		cont=cont+1
		
		
	index=index+diff+1


r = requests.get("https://www.nhc.noaa.gov/data/hurdat/hurdat2-atl-02052024.txt", allow_redirects=True)
open('hurdat2-atl-02052024.txt', 'wb').write(r.content)


hdata=open("hurdat2-atl-02052024.txt")
hdata=hdata.readlines()


indexa=index_row(hdata,"AL062018")
indexa=indexa[0]


indexb=index_row(hdata,"AL132018")
indexb=indexb[0]

index=indexa
while index <indexb:
	line_data=hdata[index]
	
	data=line_data.split(",")

	diff=int(data[2])

	lats=[]
	lons=[]
	for j in range(index+1, index+diff+1):
		line_data=hdata[j]
		n_data=line_data.split(",")
		lat=float(n_data[4][0:5])
		lon=float(n_data[5][0:6])
		ch_lon=n_data[5][6:7]

		if ch_lon=="W":
			lon=lon*(-1)

		date=n_data[0]
		if date[4:6]=="09":
			lats=np.append(lats, lat)
			lons=np.append(lons, lon)

	if len(lats)>7:
		mapa.plot(lons,lats, color="k", linewidth=1.5, marker='s')


	index=index+diff+1
plt.savefig("CyTRACK_testing_tracks.png",bbox_inches="tight")
