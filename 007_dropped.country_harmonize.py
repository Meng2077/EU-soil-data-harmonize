import numpy as np
import matplotlib.pyplot as plt
import multiprocess as mp
import glob
import time
from tqdm import tqdm
import os
import sys
import pandas as pd
import dbfread 
import geopandas as gpd
import warnings
import matplotlib.pyplot as plt
import csv
import pyproj

# +belgium
# read in 2 sites
belgium_p = pd.read_csv('/mnt/diskstation/data/soil_points/Belgium/Vlaanderen/Aardewerk-Vlaanderen-2010_Profiel.csv')
belgium_h = pd.read_csv('/mnt/diskstation/data/soil_points/Belgium/Vlaanderen/Aardewerk-Vlaanderen-2010_Horizont.csv',low_memory=False,encoding = "ISO-8859-1")
# merge 2 sites
belgium_p = belgium_p.rename(columns={'ID': 'Profiel_ID'}) 
belgium = belgium_h.merge(belgium_p, on="Profiel_ID", how="inner")
# Define the coordinate systems
lambert72 = pyproj.CRS.from_epsg(31370)  # Lambert72 CRS
wgs84 = pyproj.CRS.from_epsg(4326)  # WGS84 CRS (GPS)
transformer = pyproj.Transformer.from_crs(lambert72, wgs84)
belgium['lat'], belgium['lon'] = transformer.transform(belgium['Coordinaat_Lambert72_X'], belgium['Coordinaat_Lambert72_Y'])
# belgium['Y'], belgium['X'] = transformer.transform(belgium['Coordinaat_Bonne_E'], belgium['Coordinaat_Bonne_N'])

# convert humus to oc
belgium.loc[belgium['Humus_koolstof_nieuwe_formule']==0, 'Humus'] = belgium.loc[belgium['Humus_koolstof_nieuwe_formule']==0, 'Humus']*4/3/1.724
belgium.loc[belgium['Humus_koolstof_nieuwe_formule']==1, 'Humus'] = belgium.loc[belgium['Humus_koolstof_nieuwe_formule']==1, 'Humus']*4/3/2 # new scaler

# extract time info 
belgium['Profilering_Datum'] = belgium['Profilering_Datum'].str.split(' ').str[0]
belgium['Profilering_Datum'] = belgium['Profilering_Datum'].str.split('-').str[-1].astype(float)
belgium.loc[belgium['Profilering_Datum'] >2020, 'Profilering_Datum'] = np.nan

# extract depth info
belgium['hzn_top'] = np.nanmin(belgium[['Diepte_grens_boven1', 'Diepte_grens_boven2']], axis=1)
belgium['hzn_btm'] = np.nanmax(belgium[['Diepte_grens_onder1', 'Diepte_grens_onder2']], axis=1)
belgium.loc[belgium['hzn_top'] > belgium['hzn_btm'],['hzn_top','hzn_btm']] = np.nan

column_names = ['lat','lon','time','hzn_top','hzn_btm','ref']
temp = pd.DataFrame(columns=column_names)
temp['time'] = belgium['Profilering_Datum']
temp['hzn_top'] = belgium['Diepte_grens_boven1']
temp['hzn_btm'] = belgium['Diepte_grens_onder1']
temp.loc[temp['hzn_top'].isna(),'hzn_top'] = belgium.loc[temp['hzn_top'].isna(),'Diepte_grens_boven2']
temp.loc[temp['hzn_btm'].isna(),'hzn_btm'] = belgium.loc[temp['hzn_btm'].isna(),'Diepte_grens_onder2']
temp['hzn_top'] = belgium['hzn_top'] 
temp['hzn_btm'] = belgium['hzn_btm']   
temp['lat'] = belgium['lat']
temp['lon'] = belgium['lon']
temp['oc'] = belgium['Humus']*10
temp['caco3'] = belgium['Calciumcarbonaatgehalte']*10 # %->g/kg
temp['N'] = np.nan
temp['ph_kcl'] = belgium['pH_KCl']
temp['ph_h2o'] = belgium['pH_H2O']
temp['ph_cacl2'] = np.nan
temp['bulk_density'] = np.nan
temp['clay'] = belgium['T0_2']
temp['silt'] = belgium['T2_10']+belgium['T10_20']+belgium['T20_50']
temp['sand'] = belgium['T50_100']+belgium['T100_200']+belgium['T200_500']+belgium['T500_1000']+belgium['T1000_2000']
temp['K'] = np.nan
temp['P'] = np.nan
temp['ref'] = 'vlaanderen.belgium'
temp['nuts0'] = 'BE'

# possible filter
na = temp['time'].isna().sum()
print(f'{na} data with no time info')

na = len(temp[temp['hzn_btm'].isna() | temp['hzn_top'].isna()])
print(f'{na} data with no depth info')

na = len(temp[temp['lat'].isna() | temp['lon'].isna()])
print(f'{na} data with no coordinate info')

print(f'{len(temp)} in total')
# temp.to_csv('/mnt/primus/xuemeng_tmp_harbour/soc_eu/data/belgium_harmonized_v1.csv',index=False)



# Ireland
deep = gpd.read_file('/home/opengeohub/xuemeng/work_xuemeng/spatiotemporal-soc-eu/IE_GSI_Geochemistry_Deeper_Topsoil_S_pH_LOI_IE26_ITM/IE_GSI_Geochemistry_Deeper_Topsoil_S_pH_LOI_Sample_IE26_ITM.shp')
shallow = gpd.read_file('/home/opengeohub/xuemeng/work_xuemeng/spatiotemporal-soc-eu/IE_GSI_GSNI_Geochemistry_Shallow_Topsoil_A_pH_LOI_IE32_ITM/IE_GSI_GSNI_Geochemistry_Shallow_Topsoil_A_pH_LOI_Sample_IE32_ITM.shp')
# deep['hzn_top'] = 35
# deep['hzn_btm'] = 50
# shallow['hzn_top'] = 5
# shallow['hzn_btm'] = 20
# ireland = pd.concat([shallow,deep])   

# gdf_4326 = ireland.to_crs("EPSG:4326")
# temp = pd.DataFrame()
# temp['lat'] = gdf_4326['geometry'].y
# temp['lon'] = gdf_4326['geometry'].x
# temp['hzn_top'] = gdf_4326['hzn_top']
# temp['hzn_btm'] = gdf_4326['hzn_btm']
# temp['nuts0'] = 'IE'
# temp['ref'] = 'ireland.tellus-' + ireland['DATASOURCE']
# temp['ph_cacl2'] = gdf_4326['PH_CACL2']
# temp.to_csv('/home/opengeohub/xuemeng/work_xuemeng/spatiotemporal-soc-eu/ireland_harmonized_v1.csv',index=False)

# scotland
scotland = pd.read_excel('/mnt/diskstation/data/soil_points/Scotland/NSIS_1_10km_grid_gh.xlsx', sheet_name='NSIS1_10km')
osgb36 = pyproj.CRS.from_string("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 +units=m +no_defs")
wgs84 = pyproj.CRS.from_epsg(4326)
transformer = pyproj.Transformer.from_crs(osgb36, wgs84)
scotland['lat'], scotland['lon'] = transformer.transform(scotland['EASTING'], scotland['NORTHING'])

column_names = ['lat','lon','time','hzn_top','hzn_btm','ref']
temp = pd.DataFrame(columns=column_names)
temp['lat'] = scotland['lat']
temp['lon'] = scotland['lon']
temp['nuts0'] = 'UK-scotland'
temp['time'] = scotland['PROFILE_DATE'].astype(str).str[-4:]
temp['hzn_top'] = scotland['HORZ_TOP']
temp['hzn_btm'] = scotland['HORZ_BOTTOM']
temp['ref'] = 'scotland.NSIS1-hutton.ac.uk'

temp['oc'] = scotland['DP1971_ORGANIC_MATTER']*10/1.72
temp['N'] = scotland['DP1971_NITROGEN']*10 # % -> g/kg
temp['caco3'] = np.nan
temp['bulk_density'] = np.nan
temp['ph_kcl'] = np.nan
temp['ph_h2o'] = scotland['DP1971_PH_H2O']
temp['ph_cacl2'] = scotland['DP1971_PH_CACL2']
temp['clay'] = scotland['DP1971_CLAY']
temp['silt'] = scotland['DP1971_UBSILT']
temp['sand'] = scotland['DP1971_UBSAND']
temp['K'] = scotland['NIPAQUA_POTASSIUM'] # ppm = mg/kg
temp['P'] = scotland['NIPAQUA_PHOSPHORUS'] # ppm = mg/kg

# possible filter
na = temp['time'].isna().sum()
print(f'{na} data with no time info')

na = len(temp[temp['hzn_btm'].isna() | temp['hzn_top'].isna()])
print(f'{na} data with no depth info')

na = len(temp[temp['lat'].isna() | temp['lon'].isna()])
print(f'{na} data with no coordinate info')

print(f'{len(temp)} in total')

# temp.to_csv('/mnt/primus/xuemeng_tmp_harbour/soc_eu/data/scotland_harmonized_v1.csv',index=False)

# france
france = pd.read_csv(f'{input_path}/France/RMQS1_analyses_composites_18_11_2021_virgule.csv')
rgf93 = pyproj.CRS.from_epsg(2154)  # RGF 93 coordinate system
wgs84 = pyproj.CRS.from_epsg(4326)  # WGS84 (GPS) coordinate system
transformer = pyproj.Transformer.from_crs(rgf93, wgs84)
france['lat'],france['lon'] = transformer.transform(france['x_theo'], france['y_theo'])

column_names = ['lat','lon','time','hzn_top','hzn_btm','ref']
temp = pd.DataFrame(columns=column_names)
temp['lat'] = france['lat']
temp['lon'] = france['lon']
temp['nuts0'] = 'FR'
temp['time'] = france['date_complete'].str[0:4].astype(float)
france = france.apply(pd.to_numeric, errors='coerce')
temp['hzn_top'] = france['profondeur_hz_sup']
temp['hzn_btm'] = france['profondeur_hz_inf']
temp['ref'] = 'france.RMQS'
temp['oc'] = france['carbone_16_5_1']
temp['N'] = france['n_tot_31_1'] 
temp['caco3'] = france['calc_tot_2_1_2']
temp['bulk_density'] = np.nan
temp['ph_kcl'] = np.nan
temp['ph_h2o'] = france['ph_eau_6_1']
temp['ph_cacl2'] = np.nan
temp['clay'] = france['argile']/10 # g/kg -> %
temp['silt'] = (france['limon_fin']+france['limon_grossier'])/10 # g/kg -> %
temp['sand'] = (france['sable_fin']+france['sable_grossier'])/10 # g/kg -> %
temp['K'] = france['k_tot_hf']*10000 # g/100g -> mg/kg
temp['P'] = france['p_ass_81_1']*1000 # olsen  g/kg -> mg/kg
temp['CEC'] = france['cec_40_1']
# possible filter
na = temp['time'].isna().sum()
print(f'{na} data with no time info')

na = len(temp[temp['hzn_btm'].isna() | temp['hzn_top'].isna()])
print(f'{na} data with no depth info')

na = len(temp[temp['lat'].isna() | temp['lon'].isna()])
print(f'{na} data with no coordinate info')

print(f'{len(temp)} in total')

temp.to_csv(f'{output_path}/france_harmonized_v1.csv',index=False)