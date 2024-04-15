import os.path, datetime, multiprocessing
#from osgeo import gdal, osr
from ocgis import OcgOperations, RequestDataset, env, GeomCabinet


def doCalc(species):
        print('Working on %s' %(species))
        # Directory holding climate data.
        DATA_DIR = '/mnt/d/GIS/projects/WSI/data/'
        # Data returns will overwrite in this case. Use with caution!!
        env.OVERWRITE = True
        env.DIR_SHPCABINET = DATA_DIR
        env.DIR_OUTPUT = DATA_DIR
        # Always start with a snippet (if there are no calculations!).
        SNIPPET = False
        #yearstr = str(year)
        # Filename to variable name mapping.
        aoi = 'flyway'
        uri = '/mnt/d/GIS/projects/WSI/data/WSI_OCGIS_'+ species +'.2012_2021.nc'
        #uri = 'D:/GIS/projects/WSI/data/snod.2013.nc'
        shp = '/mnt/d/GIS/projects/WSI/data/'+ aoi + '.shp'
        env.DIR_GEOMCABINET = DATA_DIR
        sc = GeomCabinet()
        # List the shapefiles available.
        sc.keys()
        # Load geometries from the shapefile.
        geoms = sc.get_geoms('flyway')
        # RequestDatasetCollection #######################################################
        rdc = RequestDataset(uri, 'forcalc')
        rdc.dimension_map.set_variable('x', 'x')
        rdc.dimension_map.set_variable('y', 'y')
        # Return daily sum
        calc = [{'func': 'sum', 'name': 'sum'}]
        ### Write to Shapefile ###########################################################
        prefix = 'WSI_' + aoi + '_' + species  
        #print('returning shapefile for ' + species)
        ops = OcgOperations(dataset=rdc, output_format='shp', time_region={'month': [1,2,3,4,9,10,11,12]}, spatial_operation='clip', geom=geoms, calc=calc,
                        calc_raw=True, aggregate=True, calc_grouping=['day', 'month', 'year'], prefix=prefix)
        ops.execute()

if __name__ == '__main__':
    pool = multiprocessing.Pool()
    specieslist = ('abdu', 'amwi', 'bwte', 'gadw', 'gwte', 'mall', 'nsho', 'nopi')
    #specieslist = ('abdu', 'amwi')
    #doCalc('amwi')
    e = pool.map(doCalc, specieslist)
    for result in e:
            if isinstance(result, Exception):
                    print("Error: %s" % result)
            else:
                    print(result)
    pool.close()
    pool.join()
    print(e)
    print("Done")
