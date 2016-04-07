#########################
# WSI calculation following equations in Schummer et al. (2010) and Vanden Elsen et al. (2016) publications
# NARR data variables http://www.emc.ncep.noaa.gov/mmb/rreanl/narr_archive_contents.pdf
#
# Written in Python 2.7
# Mike Mitchell
# mmitchell@ducks.org
#Date 1-15-2016
#########################

import numpy as np
from numpy import *
import netCDF4, os, re, datetime, multiprocessing, time, timeit
from multiprocessing import Lock
from astral import *
from sympy import *
from datetime import date


#bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
#returns rate of change
def getRate(a,b,c,d,e,f,x,y):
        ratechange = np.multiply(np.multiply(x, x), a)
        ratechange = np.add(ratechange, np.multiply(np.multiply(y, y), b))
        ratechange = np.add(ratechange, np.multiply(np.multiply(c, x),y))
        ratechange = np.add(ratechange, np.multiply(x, d))
        ratechange = np.add(ratechange, np.multiply(y, e))
        ratechange = np.add(ratechange, f)
        return ratechange
#discriminant = b**2 - 4*a*c
#return the discriminant
def getDisc(a,b,c,d,e,f,y,ntimes, ny, nx):
        b = np.multiply(np.multiply(y,y), b)
        e = np.multiply(y,e) #y
        db = np.add(d, np.multiply(c,y))
        da = a
        dc = np.zeros((ntimes, ny, nx), dtype=float)
        dc = np.add(np.add(b,e),f)
        disc = np.zeros((ntimes, ny, nx), dtype=float)
        disc = np.add(disc, np.multiply(db,db))
        dc = dc*da*4
        disc = disc - dc
        return disc

#Display values.  Model value is the wsi value
# 1 = Few to no migrants
# 2 = Increasing abundance
# 3 = Decreasing abundance
# 4 = Can not compute
def getDisplay(model, lowerlim,threshold,ntimes,ny,nx):
        display = np.zeros((ntimes, ny, nx), dtype=int)
        disp1 = np.less(model, lowerlim)*1
        disp2 = np.logical_and(model >= lowerlim, model <= threshold)*2
        disp3 = np.greater(model, threshold)*3
        display = np.add(disp1, disp2)
        display = np.add(display, disp3)
        disp4 = np.equal(display, 0)*4
        display = np.add(display, disp4)
        return display
        
##############################################
# Where all the organizing takes place
def doCalc(species, lock=Lock()):
        try:
                tempmean = []
                path = 'R:\\Personal\\Mitchell\\Projects\\LCC_WSI\\wsi\\'#'D:\\GIS\\projects\\LCC_WSI_Climate\\python code\\wsi\\'
                os.chdir(path)

                #Sets yearlist = every year between 1979 and 2013
                yearlist = list(range(1979, 2014))
                ############## TESTING
                #yearlist = list(range(1980, 1981))
                ############## TESTING

                #Variables within netCDF that are needed for calculations
                air = 'air.2m.'
                airaccstr = 'air.airacc.'
                snod = 'snod.'
                snodaccstr = 'snod.acc.'
                photoperiod = 'Photo_period_'
                lat = 'lat.nc'
                lon = 'lon.nc'
                #/Variables 
                for year in yearlist:
                        yearstr = str(year)
                #Loop that process
                #increments counter so check if it's the first year in a process.  I created this so I could run species on different machines
                #and still have the additive features work.  For example if I'm running 1980 and 1979 wasn't run in the same process adds wouldn't work
                #so I have it setup that if you want 1980 you actually start with 1979 and it will run the last few dayd of that year before moving onto 1980

                        #Imports the air netcdf file
                        f = netCDF4.MFDataset(air + yearstr + '.nc')
                        atemp = f.variables['air']
                        with lock:
                                print('Current Species: ' + species, ' year: ', yearstr)
                        ntimes, ny, nx = atemp.shape
                        model_val = np.zeros((ntimes, ny, nx), dtype=float)
                        threshold_val = np.zeros((ntimes, ny, nx), dtype=float)
                        display_val = np.zeros((ntimes, ny, nx), dtype=int)
                        lowerlim_val = np.zeros((ntimes, ny, nx), dtype=float)
                        ratechange_val = np.zeros((ntimes, ny, nx), dtype=float)
                        wsimean = np.zeros((ntimes, ny, nx), dtype=float)
                        sum_val = np.zeros((ntimes, ny, nx), dtype=float)

                        #Air accumulation
                        airaccin = None
                        airaccin = netCDF4.Dataset(airaccstr + yearstr + '.nc')
                        airacctemp = airaccin.variables['air_acc']

                        #Snow depth
                        snodin = None
                        snodin = netCDF4.Dataset(snod + yearstr + '.nc')
                        snodtemp = snodin.variables['snod']

                        #Snow depth accumulation
                        snodaccin = None
                        snodaccin = netCDF4.Dataset(snodaccstr + yearstr + '.nc')
                        snodacctemp = snodaccin.variables['snod_acc']

                        #Photo period
                        photoin = None
                        photoin = netCDF4.Dataset(photoperiod + yearstr + '.nc')
                        phototemp = photoin.variables['photo_period']
                        
                        #Lat and lon
                        latin = netCDF4.Dataset(lat)
                        latintemp = latin.variables['lat']
                        lonin = netCDF4.Dataset(lon)
                        lonintemp = lonin.variables['lon']

                        #Convert air temp(atemp) from K to C and invert
                        atemp = np.subtract(atemp, 273.15)
                        atemp = np.multiply(atemp, -1)
                        
                        #convert snow depth(snod) from m to in
                        snodtemp = np.multiply(snodtemp, 39.3701)
                        if species in ('amwi','gadw', 'gwte', 'nsho'):
                                if year <> 1979:
                                        oldwsi = None
                                        oldwsi = netCDF4.Dataset('G:\\WSI data verification\\dataverification\\WSI_' + species + '.' + str(year-1) + '.nc')
                                        oldwsitemp = oldwsi.variables['wsi']
                                        wsimean[0] = np.mean(oldwsitemp[-7:],axis=0, dtype=float)
                                        for i in range(1,7):
                                                wsimean[i] = np.mean(np.vstack((oldwsitemp[-7+i:],wsimean[:i])), axis=0)
                                                wsimean[i] = np.add(wsimean[i], np.add(airacctemp[i], np.add(snodtemp[i], snodacctemp[i])))
                                        wsimean[6:] = np.reshape([np.mean(atemp[max(i-6, 0):i+1,:], axis=0, dtype=float) for i in range(6,len(atemp))],(ntimes-6, ny, nx))
                                        wsimean[6:] = np.add(wsimean[6:], np.add(airacctemp[6:], np.add(snodtemp[6:], snodacctemp[6:])))
                                else:
                                        wsimean = np.reshape([np.mean(atemp[max(i-6, 0):i+1,:], axis=0, dtype=float) for i in xrange(len(atemp))],(atemp.shape))
                                        wsimean = np.add(wsimean, np.add(airacctemp, np.add(snodtemp, snodacctemp)))

                ##############################################
                        if species == 'abdu':
                                #ABDU
                                #Set X variable
                                xvar = np.zeros((ntimes, ny, nx), dtype=float)
                                xvar = np.multiply(atemp, 0.965975)
                                xvar = np.add(xvar, np.multiply(airacctemp, 0.170796))
                                xvar = np.add(xvar, np.multiply(snodtemp, 0.030813))
                                xvar = np.add(xvar, np.multiply(snodacctemp, 0.196729))
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = -0.0119565155664125 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = -0.00542764724269596 #y**2
                                cvar = 0 #x*y
                                dvar = -0.13440876381107 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = 0.441504992068354 #y
                                fvar = -7.89976167023174                
                ##############################################
                        elif species == 'amwi':
                                #AMWI
                                xvar = wsimean
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = 0.00234616777805522 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = 0.0430177393856143 #y**2
                                cvar = -0.00517629478854324 #x*y
                                dvar = -0.0490211059774208 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = -3.59873018356642 #y
                                fvar = 72.2183219946714                
                #############################################
                        elif species == 'bwte':
                                #BWTE
                                xvar = phototemp
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = -0.00000438137788900899  #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = 0 #y**2
                                cvar = 0.000250098154051365 #x*y
                                dvar = 0.0293178916773958 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = -0.401557844170734 #y
                                fvar = - 10.0987389551483                
                ##############################################
                        elif species == 'gadw':
                                #GADW
                                #********************
                                #I want to make sure I need the invert of temp when calculating mean temp
                                #********************
                                xvar = wsimean
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = 0.00195591282867461 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = 0.0186000710660986 #y**2
                                cvar = -0.010619043348857 #x*y
                                dvar = 0.174605989647859 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = -1.47437807674747 #y
                                fvar = 26.5751850942147                
                ##############################################
                        elif species == 'gwte':
                                #GWTE
                                #********************
                                #I want to make sure I need the invert of temp when calculating mean temp
                                #********************
                                xvar = wsimean
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = -0.00209776377882027 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = 0.00680224416109172 #y**2
                                cvar = -0.0167660524935955 #x*y
                                dvar = 0.395882887147752 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = -0.810374939697731 #y
                                fvar = 19.0499911026847                
                ##############################################
                        elif species == 'mall':
                                #MALL
                                #Set X variable
                                xvar = np.zeros((ntimes, ny, nx), dtype=float)
                                xvar = np.multiply(atemp, 0.932168)
                                xvar = np.add(xvar, np.multiply(airacctemp, 0.235422))
                                xvar = np.add(xvar, np.multiply(snodtemp, 0.051006))
                                xvar = np.add(xvar, np.multiply(snodacctemp, 0.270255))
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = -0.00777604106504123 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = -0.00215591914154136 #y**2
                                cvar = 0.00787506223870345 #x*y
                                dvar = -0.409680237382287 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = 0.192569777713556 #y
                                fvar = -3.63665413633552                
                ##############################################
                        elif species == 'nopi':
                                #NOPI
                                accumwsi = atemp + snodtemp + airacctemp + snodacctemp
                                #Set X variable
                                xvar = accumwsi
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = -0.000516566805091564 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = 0.00614998895508152 #y**2
                                cvar = -0.106466204730593 #x*y
                                dvar = 0.692682347334105 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = -0.618019252861285 #y
                                fvar = 13.3113386287333                
                ##############################################
                        elif species == 'nsho':
                                #NSHO
                                #********************
                                #I want to make sure I need the invert of temp when calculating mean temp
                                #********************
                                xvar = wsimean
                                #WSI = xvar
                                model_val = xvar
                                #Set Y variable
                                yvar = latintemp
                                
                                #Quadratic equation
                                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                                avar = -0.000807498609214904 #x**2
                                bvar = np.zeros((ntimes, ny, nx), dtype=float)
                                bvar = 0.0222284071975206 #y**2
                                cvar = -0.00080088487258268 #x*y
                                dvar = -0.187790381664416 #x
                                evar = np.zeros((ntimes, ny, nx), dtype=float)
                                evar = -1.80451013206125 #y
                                fvar = 34.0866234887647                
                ##############################################
                #Processing for all species 

                        ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                        disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar,ntimes, ny, nx)
                        #calculate out known y to get ax**2 + bx + c
                        # #discriminant = b**2 - 4*a*c
                        disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                        
                        #roots = (-b +/- sq.rt(discriminant))/(2*a)
                        lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                        threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                        
                        model_val = xvar
                        display_val = getDisplay(model_val, lowerlim_val, threshold_val, ntimes, ny, nx)
                        sum_val = np.equal(display_val, 3)
                ##############################################                        
                #Create NetCDF file

                        print('Creating variables for ' + species)
                        model_val_year = ('G:\\WSI data verification\\dataverification\\WSI_' + species + '.' + yearstr + '.nc')
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
                        

                        model_val_v = nco.createVariable('wsi', 'f4',  ('time', 'y', 'x'))
                        model_val_v.units='wsi value'
                        model_val_v.long_name='WSI Value'
                        model_val_v.grid_mapping = 'Lambert_Conformal'

                        model_val_t = nco.createVariable('threshold', 'f4',  ('time', 'y', 'x'))
                        model_val_t.units='threshold value'
                        model_val_t.long_name='Threshold Value'
                        model_val_t.grid_mapping = 'Lambert_Conformal'

                        model_val_d = nco.createVariable('display', 'f4',  ('time', 'y', 'x'))
                        model_val_d.units='display valye'
                        model_val_d.long_name='Display Value'
                        model_val_d.grid_mapping = 'Lambert_Conformal'

                        model_val_ll = nco.createVariable('lower_limit', 'f4',  ('time', 'y', 'x'))
                        model_val_ll.units='lower limit valye'
                        model_val_ll.long_name='Lower Limit Value'
                        model_val_ll.grid_mapping = 'Lambert_Conformal'

                        model_val_cr = nco.createVariable('abundance_change_rate', 'f4',  ('time', 'y', 'x'))
                        model_val_cr.units='LN change rate'
                        model_val_cr.long_name='Rate of change in relative abundance on a log scale'
                        model_val_cr.grid_mapping = 'Lambert_Conformal'

                        model_val_s = nco.createVariable('forcalc', 'f4',  ('time', 'y', 'x'))
                        model_val_s.units='sum valye'
                        model_val_s.long_name='Sum Value'
                        model_val_s.grid_mapping = 'Lambert_Conformal'                

                        with lock:
                                print('Copying variables for ' + species)
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

                        with lock:
                                print('Writing data for ' + species)
                        #  write the model_val data
                        model_val_v[:,:,:] = model_val
                        model_val_t[:,:,:] = threshold_val
                        model_val_d[:,:,:] = display_val
                        model_val_ll[:,:,:] = lowerlim_val
                        model_val_cr[:,:,:] = ratechange_val
                        model_val_s[:,:,:] = sum_val

                        # copy Global attributes from original file
                        for att in f.ncattrs():
                                setattr(nco,att,getattr(f,att))

                        nco.Conventions='CF-1.6'
                        nco.close()
                        with lock:
                                print('Done with ' + species)
        except Exception as e:
                return e
####################################
#  MAIN
####################################
if __name__ == '__main__':
        print('Starting script')
        #Start timing entire process
        startfirst = timeit.default_timer()
        
        #Create pool
        print("Creating pool")
        pool = multiprocessing.Pool(8)
        #Call doCalc and pass the year list
        #specieslist = ('abdu', 'amwi', 'bwte', 'gadw', 'gwte', 'mall', 'nopi', 'nsho')
        specieslist = ('abdu', 'amwi', 'bwte', 'gadw', 'gwte', 'mall', 'nopi', 'nsho')
        e, species, yearstr = pool.map(doCalc, specieslist)
        for result in e:
                if isinstance(result, Exception):
                        print "Error: %s" % result
                        print "Species: %s / year: %s" %(species[count], yearstr[count])
                else:
                        print result
        pool.close()
        pool.join()
        stoplast = timeit.default_timer()
        print(stoplast - startfirst)     
        print("Complete")
