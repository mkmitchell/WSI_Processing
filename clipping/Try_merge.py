import os.path, datetime, multiprocessing
import arcpy

def doCalc(species):
    mergethese = []
    try:
        yearlist = list(range(1979, 2014))
        # Directory holding climate data.
        DATA_DIR = 'G:/WSI data verification/dataverification'
        shapefile = DATA_DIR + '/WSI_' + species + '_big/WSI_' + species + '_big.shp'
        mergethese.append(shapefile)
        arcpy.AddField_management(in_table=shapefile, field_name="SPECIES", field_type="STRING", field_precision="", field_scale="", field_length="",
                                  field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
        arcpy.CalculateField_management(in_table=shapefile, field="SPECIES", expression='"' + species.upper() + '"', expression_type="PYTHON_9.3", code_block="")
                         
    except Exception as e:
        return e



####################################
#  MAIN
####################################
if __name__ == '__main__':
    pool = multiprocessing.Pool(4)
    specieslist = ('gadw', 'gwte', 'amwi')
    e = pool.map(doCalc, specieslist)
    print 'Finished'
    for result in e:
            if isinstance(result, Exception):
                    print "Error: %s" % result
            else:
                    print result
    pool.close()
    pool.join()
    print("Done")
