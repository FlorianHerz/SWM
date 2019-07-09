# -*- coding: cp1252 -*-
###################################################################################################
#
# SWM_adapdet by Florian Herz
# Berechnet für einen selbst geählten Zeitraum die PET und speichert die verwendeten Klimadaten
# in einer Tabelle in einer GDB ab.
#
#Erweiterung um die AET (allerdings fehlt die Veränderung des Bodenwasserspeichers des Vortages
#
###################################################################################################


#Import der Systemmodule
import os
import sys
import arcpy
from arcpy.sa import *

#Einstellungen und Erweiterungen laden

arcpy.CheckOutExtension("Spatial")
##arcpy.env.workspace = r'C:\HiWi_Hydro-GIS' #Workingdirectory
##arcpy.env.overwriteOutpt = True

print("Systemmodule geladen.")

#Definition und Erstellen des Dateipfad der Basisdaten und der Ergebnisdatenbank

##name = arcpy.GetParameterAsText(1)
##ordner = arcpy.GetParameterAsText(2)
name = r'SWM_Ergebnisdatenbank_AET_PET_2806-04072004_I.gdb'
ordner = r'C:\HiWi_Hydro-GIS'
ergebnisse = arcpy.CreateFileGDB_management(ordner, name)
arcpy.env.workspace = ordner+r'\{}'.format(name)

#Zeitpunkt angeben

#start = arcpy.GetParameterAsText(3)
#ende = arcpy.GetParameterAsText(4)

start = 20040628
ende  = 20040704

#Abfrage welche Daten gespeichert werden sollen

check_PET = 1
check_AET = 1


#wird nur benötigt um Klimatabelle abzuspeichern
ID = []
temp = []
jahr = []
monat = []
tag = []
relf = []


tftable = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\TempFeuchte' #Grundlegende Klimadaten

#Zuordnung der Haude Rasterdateien zu den Monatszahlen in einem dictionary

print("Lade Rasterdateien des Haudefaktors...")
haude_dic = { 1 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_1'), \
              2 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_2'), \
              3 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_3'), \
              4 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_4'), \
              5 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_5'), \
              6 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_6'), \
              7 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_7'), \
              8 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_8'), \
              9 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_9'), \
              10 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_10'), \
              11 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_11'), \
              12 : Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Haude_12')}
print("...geladen.")

#Variablendefinition
RP_faktor = 0.85 #muss von Benutzer verändert werden können !!!

print("Berechnung der Rasterdatensätze...")
# Laden der Rasterdatensätze
FK = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L')
RP = FK * RP_faktor
WP = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\WP')
gewaesser = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\Gewaesser')
RPdiffWP = RP - WP
Svortag = FK
print("...berechnet.")

#Iteration durch die ausgewählten Klimadaten

with arcpy.da.SearchCursor(tftable, ['TagesID', 'Jahr', 'Monat', 'Tag', 'RelFeu', 'Temp'], \
                           "TagesID >= {0} AND TagesID <= {1}".format(start, ende)) as klimacursor:
    for row in klimacursor:
##        ID.append(row[0])
##        jahr.append(row[1])
##        monat.append(row[2])
##        tag.append(row[3])
##        relf.append(row[4])
##        temp.append(row[5])
        
##        ID = row[0]
##        jahr = row[1]
##        monat = row[2]
##        tag = row[3]
##        relf = row[4]
##        temp = row[5]

        
        #Berechnung der PET
        
        Es = 6.1 * 10**((7.5 * (row[5])/10.0) / ((row[5])/10.0 + 237.2)) #Berechnung des Es
        PET = haude_dic[row[2]] * Es * (1.0 - row[4] / 100.0) #Berechnung der PET
        if check_PET ==1:
            PET.save("PET_{}".format(row[0])) #Speichern der einzelnen PET-Rasterdatensätze
        
        #Berechnung der AET

        AET = Con(gewaesser == 1, PET, Con(Svortag >= RP, PET, Con(RPdiffWP == 0, 0, ((Svortag - WP)/RPdiffWP)*PET)))

        if check_AET == 1:
            AET.save("AET_{}".format(row[0])) # Speichern der einzelnen AET-Rasterdatensätze

        
        print(row[3], row[2], row[1])
        

#Speichern der Klimadaten in einer Tabelle

##arcpy.CreateTable_management(ergebnisse, "Klimatabelle", tftable)
##        
####arcpy.AddField_management("Klimatabelle", "Monat", "DOUBLE")
####arcpy.AddField_management("Klimatabelle", "Monat", "DOUBLE")
####arcpy.AddField_management("Klimatabelle", "rel. Luftfeuchte", "DOUBLE")
####arcpy.AddField_management("Klimatabelle", "Temperatur", "DOUBLE")


##cursor = arcpy.da.InsertCursor("Klimatabelle", "*")
##
##for i in range(len(temp)):
##    cursor.insertRow([i+1, 2634, jahr[i], monat[i], tag[i], relf[i], temp[i]/10, ID[i]])
##
### Delete cursor objects
##del klimacursor
##del cursor
###arcpy.Delete_management("Klimatabelle")






