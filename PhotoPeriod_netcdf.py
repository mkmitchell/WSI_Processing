import numpy as np
from numpy import *
import netCDF4, os, re, datetime, multiprocessing, time, timeit
from multiprocessing import Lock
from astral import LocationInfo
from astral.sun import sun
from sympy import *
from datetime import date


#Return date object
def getDate(day, year):
        datereturn = datetime.datetime(year, 1, 1) + datetime.timedelta(day)
        datereturn = datereturn.timetuple()
        return datereturn


def doCalc(year):            
        path = '\\\\10.0.0.2\\work\\gis\\projects\\WSI\\data\\'#'D:\\GIS\\projects\\LCC_WSI_Climate\\python code\\wsi\\'
        os.chdir(path)

        #Variables within netCDF that are needed for calculations
        air = 'air.2m.'
        lat = 'lat.nc'
        lon = 'lon.nc'
        #/Variables 

        #Setting up code variables
        goodmonth = (1,2,3,9,10,11,12)
        #Loop that processes every year within yearlist

        yearstr = str(year)

        #Imports the air netcdf file
        f = netCDF4.MFDataset(air + yearstr + '.nc')
        atemp = f.variables['air']
        print('Current year: ' + yearstr)
        ntimes, ny, nx = atemp.shape
        photo_val = np.zeros((ntimes, ny, nx), dtype=float)

        #Lat and lon
        latin = netCDF4.Dataset(lat)
        latintemp = latin.variables['lat']
        lonin = netCDF4.Dataset(lon)
        lonintemp = lonin.variables['lon']
                
        #Looping through time variable
        for i in range(ntimes): #ntimes

                datereturn = getDate(i, year)
                day = datereturn[2]
                month = datereturn[1]

                if month not in goodmonth:
                        continue
                if yearstr == '2019' and month not in (1,2,3,4):
                        continue

                print
                print('Processing '  + str(month) + '/' + str(day) + '/' + yearstr)                            
                #Process only the area of interest
                for b in range(50, 180): #ny 50, 180
                        for c in range(150, 280): #nx 150, 280
                                #need photo period
                                d = date(year, month, day)
                                city = LocationInfo('Atlanta')
                                s = sun(city.observer, date=d)   
                                photo = s['sunset'] - s['sunrise']
                                photosec = photo.total_seconds()
                                photomin = photosec / 60
                                photo_val[i,b,c] = photomin
                                 
                                        
        # create NetCDF file
        print('Creating variables for ', yearstr)
        model_val_year = ('\\\\10.0.0.2\\work\\gis\\projects\\WSI\\data\\Photo_period_' + yearstr + '.nc')
        nco = netCDF4.Dataset(model_val_year,'w',clobber=True)
        nco.createDimension('time', None)
        nco.createDimension('x',nx)
        nco.createDimension('y',ny)

        timeo = nco.createVariable('time', 'f8', ('time'))
        lono = nco.createVariable('lon','f4',('y','x'))
        lato = nco.createVariable('lat','f4',('y','x'))
        xo = nco.createVariable('x','f4',('x'))
        yo = nco.createVariable('y','f4',('y'))
        lco = nco.createVariable('Lambert_Conformal','i4')
        

        photo_val_v = nco.createVariable('photo_period', 'f4',  ('time', 'y', 'x'))
        photo_val_v.units='minutes'
        photo_val_v.long_name='Photo period in minutes'
        photo_val_v.grid_mapping = 'Lambert_Conformal'

        # copy all the variable attributes from original file
        for var in ['time', 'lon','lat','x','y', 'Lambert_Conformal']:
                for att in f.variables[var].ncattrs():
                        setattr(nco.variables[var],att,getattr(f.variables[var],att))

        # copy variable data for lon,lat,x and y
        timeo[:]=f.variables['time'][:]
        lono[:]=f.variables['lon'][:]
        lato[:]=f.variables['lat'][:]
        xo[:]=f.variables['x'][:]
        yo[:]=f.variables['y'][:]

        #  write the model_val data
        photo_val_v[:,:,:] = photo_val


        # copy Global attributes from original file
        for att in f.ncattrs():
                setattr(nco,att,getattr(f,att))

        nco.Conventions='CF-1.6'
        nco.close()


####
#  MAIN
####
if __name__ == '__main__':
        print('Starting script')
        #Start timing entire process
        startfirst = timeit.default_timer()
        
        #Call doCalc
        yearlist = list(range(1979, 2022))
        pool = multiprocessing.Pool(4)
        e = pool.map(doCalc, yearlist)
        pool.close()
        #print("Joining pool")
        pool.join()
        print(e)
        stoplast = timeit.default_timer()
        print(stoplast - startfirst)     
        print("Complete")
