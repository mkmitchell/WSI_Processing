import os.path, datetime, multiprocessing
import arcpy

def doCalc(species):
    try:
        yearlist = list(range(1979, 2014))
        # Directory holding climate data.
        DATA_DIR = 'G:/WSI data verification/dataverification'
        shapefile = DATA_DIR + '/Merge_' + species + '.shp'
species = 'gwte'
shapefile = 'Merge_' + species.lower()
wsi = 'WSI_' + species.upper()
arcpy.AddField_management(in_table=shapefile, field_name="DT", field_type="DATE", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shapefile, field_name="SQKM", field_type="DOUBLE", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shapefile, field_name="STATE", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shapefile, field_name="CATEGORY", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shapefile, field_name="SUBCAT", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table=shapefile, field_name="SPECIES", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.CalculateField_management(in_table=shapefile, field="DT", expression="Calc( !TIME! )", expression_type="PYTHON_9.3",
                                code_block="from dateutil.parser import parse\nfrom datetime import datetime\ndef Calc(tm):\n  tm = parse(tm)\n  return tm")
arcpy.CalculateField_management(in_table=shapefile, field="SQKM", expression="!SUM! * 1024", expression_type="PYTHON_9.3", code_block="")
arcpy.DeleteField_management(in_table=shapefile, drop_field=["TIME", "LB_TIME", "UB_TIME", "YEAR", "MONTH", "DAY", "SUM", "TID"])



arcpy.CalculateField_management(in_table=shapefile, field='"' + shapefile + '.SUBCAT"', expression="[WSI_abdu_sum_1979_1_ugid.STATE_NAME]", expression_type="VB", code_block="")


arcpy.CalculateField_management(in_table=shapefile, field="SPECIES", expression='"' + species.upper() + '"', expression_type="PYTHON_9.3", code_block="")
arcpy.CalculateField_management(in_table=shapefile, field="CATEGORY", expression='"STATE"', expression_type="PYTHON_9.3", code_block="")
arcpy.AddField_management(in_table=shapefile, field_name="DATE", field_type="DATE", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.CalculateField_management(in_table=shapefile, field="DATE", expression="!DT!", expression_type="PYTHON_9.3", code_block="")
arcpy.DeleteField_management(in_table=shapefile, drop_field=["DT","UGID"])

                                
    except Exception as e:
        return e

####################################
#  MAIN
####################################
if __name__ == '__main__':
    pool = multiprocessing.Pool(4)
    specieslist = ('bwte', 'mall', 'nsho', 'nopi', 'gadw', 'gwte', 'abdu', 'amwi')
    e = doCalc('abdu')
    #e = pool.map(doCalc, specieslist)
    for result in e:
            if isinstance(result, Exception):
                    print "Error: %s" % result
            else:
                    print result
    pool.close()
    pool.join()
    print("Done")
