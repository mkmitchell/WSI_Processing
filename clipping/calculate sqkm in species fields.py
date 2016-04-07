for i in ['ABDU', 'AMWI', 'BWTE', 'GADW', 'GWTE', 'MALL', 'NOPI', 'NSHO']:
        arcpy.AddField_management(in_table="WSI_JV_" + i.lower(), field_name="SPECIES", field_type="STRING", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
        arcpy.CalculateField_management(in_table="WSI_JV_" + i.lower(), field="SPECIES", expression='"' + i + '"', expression_type="PYTHON_9.3", code_block="")

############ RUN MERGE HERE.  Could probably setup a script
############ Rename to Merge_all
for i in ['ABDU', 'AMWI', 'BWTE', 'GADW', 'GWTE', 'MALL', 'NOPI', 'NSHO']:
        arcpy.AddField_management(in_table="Merge_all", field_name=i + "SQKM", field_type="DOUBLE", field_precision="", field_scale="", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
for i in ['ABDU', 'AMWI', 'BWTE', 'GADW', 'GWTE', 'MALL', 'NOPI', 'NSHO']:
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="Merge_all", selection_type="NEW_SELECTION", where_clause='"SPECIES" = \'' + i + '\'')
        arcpy.CalculateField_management(in_table="Merge_all", field=i + "SQKM", expression="Calc( !SPECIES!, !SUM!)", expression_type="PYTHON_9.3", code_block="def Calc(sp, sum):\n  return sum*1024\n\n")

# CLEAR SELECTION

arcpy.SelectLayerByAttribute_management(in_layer_or_view="Merge_all", selection_type="NEW_SELECTION", where_clause='')
arcpy.AddField_management(in_table="Merge_all", field_name="STATE", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table="Merge_all", field_name="CATEGORY", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table="Merge_all", field_name="SUBCAT", field_type="STRING", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.AddField_management(in_table="Merge_all", field_name="DT", field_type="DATE", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
arcpy.CalculateField_management(in_table="Merge_all", field="DT", expression="Calc( !TIME! )", expression_type="PYTHON_9.3",
                                code_block="from dateutil.parser import parse\nfrom datetime import datetime\ndef Calc(tm):\n  tm = parse(tm)\n  return tm")
arcpy.DeleteField_management(in_table="Merge_all", drop_field=["TIME", "LB_TIME", "UB_TIME", "YEAR", "MONTH", "DAY", "SUM", "TID"])
arcpy.AddField_management(in_table="Merge_all", field_name="DATE", field_type="DATE", field_precision="", field_scale="", field_length="",
                          field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")

