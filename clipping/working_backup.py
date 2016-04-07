import os.path, datetime, multiprocessing
from ocgis import OcgOperations, RequestDataset, RequestDatasetCollection, env
yearlist = list(range(1979,1980))
def doCalc(species):

    for year in yearlist:
        yearstr = str(year)
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
        uri = 'G:/WSI data verification/dataverification/WSI_'+species+'.' + yearstr + '.nc'
        shp = 'G:/WSI data verification/dataverification/state.shp'
        # RequestDatasetCollection #######################################################
        rdc = RequestDataset(uri, 'forcalc')
        # Return daily sum
        calc = [{'func': 'sum', 'name': 'sum'}]
        prefix = 'WSI_' + species + '_sum_' + yearstr
        ### Write to Shapefile ###########################################################
        print('returning shapefile...')
        ops = OcgOperations(dataset=rdc, output_format='shp', time_region={'month': [1]}, spatial_operation='clip', geom=shp, calc=calc,
                        calc_raw=True, aggregate=True, calc_grouping='day', prefix='calc')
        ret = ops.execute()

        
if __name__ == '__main__':
    #pool = multiprocessing.Pool(1)
    specieslist = ('abdu','bwte')#, 'mall', 'nsho', 'nopi', 'gadw', 'gwte', 'amwi')
    #e = pool.map(doCalc, specieslist)
    e = doCalc('abdu')
    for result in e:
            if isinstance(result, Exception):
                    print "Error: %s" % result
            else:
                    print result
    pool.close()
    pool.join()
    print e
    print("Done")