import os.path, datetime, multiprocessing
from ocgis import OcgOperations, RequestDataset, RequestDatasetCollection, env
def doCalc(species):
    print 'Working on %s' %(species)
    # Directory holding climate data.
    DATA_DIR = 'G:/WSI data verification/dataverification'
    # Data returns will overwrite in this case. Use with caution!!
    env.OVERWRITE = True
    env.DIR_SHPCABINET = DATA_DIR
    env.DIR_OUTPUT = DATA_DIR
    # Always start with a snippet (if there are no calculations!).
    SNIPPET = False
    #yearstr = str(year)
    # Filename to variable name mapping.
    uri = 'G:/WSI data verification/dataverification/'+ species + '/WSI_OCGIS_'+species+'.1979_2013.nc'
    shp = 'G:/WSI data verification/dataverification/duckzone.shp'
    # RequestDatasetCollection #######################################################
    rdc = RequestDataset(uri, 'forcalc')
    # Return daily sum
    calc = [{'func': 'sum', 'name': 'sum'}]
    ### Write to Shapefile ###########################################################
    prefix = 'WSI_DZ_' + species  
    #print('returning shapefile for ' + species)
    ops = OcgOperations(dataset=rdc, output_format='shp', time_region={'month': [1,2,3,4,9,10,11,12]}, spatial_operation='clip', geom=shp, calc=calc,
                    calc_raw=True, aggregate=True, calc_grouping=['day', 'month', 'year'], prefix=prefix)
    ops.execute()

if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    specieslist = ('abdu','bwte', 'mall', 'nsho', 'nopi', 'gadw', 'gwte', 'amwi')
    e = pool.map(doCalc, specieslist)
    for result in e:
            if isinstance(result, Exception):
                    print "Error: %s" % result
            else:
                    print result
    pool.close()
    pool.join()
    print e
    print("Done")
