import arcpy
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")

arcpy.env.overwriteOutput = True
name = r'Test.gdb'
folder = r'C:\HiWi_Hydro-GIS'
s_init = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L')
cellsize = arcpy.Describe(s_init).children[0].meanCellHeight

#arcpy.CreateFileGDB_management(folder, name)

basin = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\EZG_Schotten2_Vektor'
input_raster = Raster(r'C:\HiWi_Hydro-GIS\PY_SWM_test.gdb\R_20030101')
# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "EZG_Schotten2_Vektor"

raster = input_raster * 10
raster_clip = Raster(arcpy.Clip_management(in_raster=raster,
                                           rectangle="3507687,5 5592837,5 3512037,5 5595787,5",
                                           out_raster="C:/HiWi_Hydro-GIS/Test.gdb/R_clip_20030101.py",
                                           in_template_dataset=basin,
                                           nodata_value="-1,797693e+308",
                                           clipping_geometry="ClippingGeometry",
                                           maintain_clipping_extent="NO_MAINTAIN_EXTENT"))
array = arcpy.RasterToNumPyArray(raster_clip, nodata_to_value=0)


sum_array = array.sum()

print(sum_array)

#Ram angucken

arcpy.gp.RasterCalculator_sa('"AET_20030101" *0 + 1', "c:/Users/Florian/documents/ArcGIS/Default.gdb/aet_20030101_clip_test_3")
