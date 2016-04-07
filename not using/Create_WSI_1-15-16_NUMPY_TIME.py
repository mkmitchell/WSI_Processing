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
        dc = np.zeros(( ntimes, ny, nx), dtype=float)
        dc = np.add(np.add(b,e),f)
        disc = np.zeros(( ntimes, ny, nx), dtype=float)
        disc = np.add(disc, np.multiply(db,db))
        dc = dc*da*4
        disc = disc - dc
        return disc

def getDisplay(model, lowerlim,threshold, ntimes,ny,nx):
        display = np.zeros(( ntimes, ny, nx), dtype=int)
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
        tempmean = []
        path = 'R:\\Personal\\Mitchell\\Projects\\LCC_WSI\\wsi\\'#'D:\\GIS\\projects\\LCC_WSI_Climate\\python code\\wsi\\'
        os.chdir(path)

        #Sets yearlist = every year between 1979 and 2013
        #yearlist = list(range(1979, 2014))
        ############## TESTING
        #yearlist = list(range(1980, 1981))
        ############## TESTING

        #Variables within netCDF that are needed for calculations
        air = 'air.2m.1979_2013_ncrcat.nc'
        airaccstr = 'air.airacc.1979_2013_ncrcat.nc'
        snod = 'snod.1979_2013_ncrcat.nc'
        snodaccstr = 'snod.acc.1979_2013_ncrcat.nc'
        photoperiod = 'Photo_period_1979_2013_ncrcat.nc' #opened when needed then closed
        lat = 'lat.nc'
        lon = 'lon.nc'
        #/Variables 

        #Imports the air netcdf file
        f = netCDF4.MFDataset(air)
        atemp = f.variables['air']
        with lock:
                print('Current Species: ' + species)
        ntimes, ny, nx = atemp.shape
        model_val = np.zeros(( ntimes, ny, nx), dtype=float)
        threshold_val = np.zeros(( ntimes, ny, nx), dtype=float)
        display_val = np.zeros(( ntimes, ny, nx), dtype=int)
        lowerlim_val = np.zeros(( ntimes, ny, nx), dtype=float)
        ratechange_val = np.zeros(( ntimes, ny, nx), dtype=float)

        #Air accumulation
        airaccin = None
        airaccin = netCDF4.Dataset(airaccstr)
        airacctemp = airaccin.variables['air_acc']
        #Snow depth
        snodin = None
        snodin = netCDF4.Dataset(snod)
        snodtemp = snodin.variables['snod']
        #Snow depth accumulation
        snodaccin = None
        snodaccin = netCDF4.Dataset(snodaccstr)
        snodacctemp = snodaccin.variables['snod_acc']
        #Photo period opened when needed then closed

        
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
        wsimean = np.reshape([np.mean(atemp[max(i-6, 0):i+1,:], axis=0) for i in xrange(len(atemp))],(atemp.shape))
        wsimean = np.add(wsimean, airacctemp)
        wsimean = np.add(wsimean, snodtemp)
        wsimean = np.add(wsimean, snodacctemp)
        accumwsi = atemp + snodtemp + airacctemp + snodacctemp
##############################################
        if species == 'abdu':
                #ABDU
                #Set X variable
                xvar = np.zeros(( ntimes, ny, nx), dtype=float)
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
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = -0.00542764724269596 #y**2
                cvar = 0 #x*y
                dvar = -0.13440876381107 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = 0.441504992068354 #y
                fvar = -7.89976167023174                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)
##############################################
        elif species == 'amwi':
                #AMWI
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
                avar = 0.00234616777805522 #x**2
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = 0.0430177393856143 #y**2
                cvar = -0.00517629478854324 #x*y
                dvar = -0.0490211059774208 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = -3.59873018356642 #y
                fvar = 72.2183219946714                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)
#############################################
        elif species == 'bwte':
                photoin = None
                photoin = netCDF4.Dataset(photoperiod)
                phototemp = photoin.variables['photo_period']
                #BWTE
                xvar = phototemp
                #WSI = xvar
                model_val = xvar
                #Set Y variable
                yvar = latintemp
                
                #Quadratic equation
                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                avar = -0.00000438137788900899  #x**2
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = 0 #y**2
                cvar = 0.000250098154051365 #x*y
                dvar = 0.0293178916773958 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = -0.401557844170734 #y
                fvar = - 10.0987389551483                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)
                photoin.close()
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
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = 0.0186000710660986 #y**2
                cvar = -0.010619043348857 #x*y
                dvar = 0.174605989647859 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = -1.47437807674747 #y
                fvar = 26.5751850942147                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)
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
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = 0.00680224416109172 #y**2
                cvar = -0.0167660524935955 #x*y
                dvar = 0.395882887147752 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = -0.810374939697731 #y
                fvar = 19.0499911026847                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)                 
##############################################
        elif species == 'mall':
                #MALL
                #Set X variable
                xvar = np.zeros(( ntimes, ny, nx), dtype=float)
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
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = -0.00215591914154136 #y**2
                cvar = 0.00787506223870345 #x*y
                dvar = -0.409680237382287 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = 0.192569777713556 #y
                fvar = -3.63665413633552                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)

##############################################
        elif species == 'nopi':
                #NOPI
                #Set X variable
                xvar = accumwsi
                #WSI = xvar
                model_val = xvar
                #Set Y variable
                yvar = latintemp
                
                #Quadratic equation
                #bivariate case: f(x,y) = ax**2 + by**2 + cxy + dx + ey + f
                avar = -0.000516566805091564 #x**2
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = 0.00614998895508152 #y**2
                cvar = -0.106466204730593 #x*y
                dvar = 0.692682347334105 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = -0.618019252861285 #y
                fvar = 13.3113386287333                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)

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
                bvar = np.zeros(( ntimes, ny, nx), dtype=float)
                bvar = 0.0222284071975206 #y**2
                cvar = -0.00080088487258268 #x*y
                dvar = -0.187790381664416 #x
                evar = np.zeros(( ntimes, ny, nx), dtype=float)
                evar = -1.80451013206125 #y
                fvar = 34.0866234887647                
                
                ratechange_val = getRate(avar, bvar, cvar, dvar, evar, fvar, xvar, yvar)
                disc = getDisc(avar, bvar, cvar, dvar, evar, fvar, yvar, ntimes, ny, nx)
                #calculate out known y to get ax**2 + bx + c
                # #discriminant = b**2 - 4*a*c
                disc_bnew = np.add(dvar, np.multiply(cvar,yvar))
                
                #roots = (-b +/- sq.rt(discriminant))/(2*a)
                lowerlim_val = (disc_bnew*-1 + np.sqrt(disc))/(2*avar)
                threshold_val = (disc_bnew* -1 - np.sqrt(disc))/(2*avar)                
                
                model_val = xvar
                display_val = getDisplay(model_val, lowerlim_val, threshold_val,  ntimes, ny, nx)  
                        
################ create NetCDF file

                print('Creating variables for ' + species)
                model_val_year = ('R:\\Personal\\Mitchell\\Projects\\LCC_WSI\\wsi\\big_test_' + species + '_wsi.nc')
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

                # copy Global attributes from original file
                for att in f.ncattrs():
                        setattr(nco,att,getattr(f,att))

                nco.Conventions='CF-1.6'
                nco.close()
                with lock:
                        print('Done with ' + species)

####################################
#  MAIN
####################################
if __name__ == '__main__':
        print('Starting script')
        #Start timing entire process
        startfirst = timeit.default_timer()
        
        #Create pool
        print("Creating pool")
        #pool = multiprocessing.Pool(8)
        #Call doCalc and pass the year list
        specieslist = ('abdu','')# 'bwte', 'mall', 'nsho', 'nopi', 'amwi', 'gadw', 'gwte')
        for spec in specieslist:
            doCalc(spec)
            
        # try:
                # print('trying nothing')
                # e = pool.map(doCalc, specieslist)
        # except:
                # print("Error", e)
        # print("Closing pool")
        # pool.close()
        #print("Joining pool")
        #pool.join()
        #print(e)
        stoplast = timeit.default_timer()
        print(stoplast - startfirst)     
        print("Complete")
