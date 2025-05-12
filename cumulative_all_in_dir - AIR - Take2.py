from numpy import *
import netCDF4
import os
import re
import datetime
import calendar

path = 'D:\\GIS\\projects\\WSI\\data'
os.chdir(path)
yearinput = []
goodmonth = (1, 2, 3, 9, 10, 11, 12)

def getDate(day, year):
	date = datetime.datetime(year, 1, 1) + datetime.timedelta(day)
	date = date.timetuple()
	return date

print('Calculating air accumulation')

for a in os.listdir(path):
        print(datetime.datetime.now())
        f = netCDF4.Dataset(a)
        yearsearch = re.search(r'air.2m.(\d\d\d\d).nc', a)
        atemp = f.variables['air']
        yearstr = yearsearch.group(1)
        print('Current year: ' + yearstr)
        ntimes, ny, nx = atemp.shape
        cold_days = zeros((ntimes, ny, nx), dtype=int)
        
        if yearstr == '2012':
                print('Looping through netCDF')
                for i in range(243,ntimes):
                        print(i)
                        date = getDate(i, int(yearstr))
                        print(date[1], date[2], int(yearstr))
                        for b in range(50, 180): #50, 180
                                for c in range(150,280): #150,280                     
                                        if i == 243:
                                                if atemp[i,b,c] < 273.15:
                                                        cold_days[i,b,c] = 1
                                                else:
                                                        cold_days[i,b,c] = 0
                                        else:
                                                if atemp[i,b,c] < 273.15:
                                                        cold_days[i,b,c] = cold_days[i-1,b,c] + 1
                                                else:
                                                        cold_days[i,b,c] = 0
        else:
                previous = len(yearinput)
                previousyear = int(yearstr) -1
                        
                print('Looping through netCDF')
                for i in range(ntimes): #ntimes
                        print(i)
                        date = getDate(i, int(yearstr))
                        print(date[1], date[2], int(yearstr))
                        month = date[1]
                        day = date[2]
                        if month not in goodmonth:
                                       continue
                        for b in range(50, 180):
                                for c in range(150, 280):                     
                                        if i == 0:
                                                old = netCDF4.Dataset('air.airacc.' + str(previousyear) + '.nc')
                                                oldtimes = old.variables['air_acc']
                                                if atemp[i,b,c] < 273.15:
                                                        if calendar.isleap(previousyear):
                                                                cold_days[i,b,c] = oldtimes[365, b,c] + 1
                                                        else:
                                                                cold_days[i,b,c] = oldtimes[364, b,c] + 1
                                                else:
                                                        cold_days[i,b,c] = 0
                                        else:
                                                if month == 9 and day == 1:
                                                        if atemp[i,b,c] < 273.15:
                                                                cold_days[i,b,c] = 1
                                                        else:
                                                                cold_days[i,b,c] = 0

                                                elif atemp[i,b,c] < 273.15:
                                                        cold_days[i,b,c] = cold_days[i-1,b,c] + 1
                                                else:
                                                        cold_days[i,b,c] = 0

        # create NetCDF file
        cold_days_year = ('air.airacc.' + yearstr + '.nc')
        nco = netCDF4.Dataset(cold_days_year,'w',clobber=True)
        nco.createDimension('time', None)
        nco.createDimension('x',nx)
        nco.createDimension('y',ny)

        timeo = nco.createVariable('time', 'f8', ('time'))
        lono = nco.createVariable('lon','f4',('y','x'))
        lato = nco.createVariable('lat','f4',('y','x'))
        xo = nco.createVariable('x','f4',('x'))
        yo = nco.createVariable('y','f4',('y'))
        lco = nco.createVariable('Lambert_Conformal','i4')
        

        cold_days_v = nco.createVariable('air_acc', 'i4',  ('time', 'y', 'x'))
        cold_days_v.units='days'
        cold_days_v.long_name='total number of days below 0 degC'
        cold_days_v.grid_mapping = 'Lambert_Conformal'

        print('Copying variables')
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

        print('Writing data')
        #  write the cold_days data
        cold_days_v[:,:,:] = cold_days

        # copy Global attributes from original file
        for att in f.ncattrs():
                setattr(nco,att,getattr(f,att))

        nco.Conventions='CF-1.6'
        nco.close()
        yearinput.append(nco)

