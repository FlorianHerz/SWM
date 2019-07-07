# -*- coding: cp1252 -*-

""" SWM für ArcGIS """
""" Soil Water Model in Python für ArcGIS für die Lehrveranstaltung "GIS für Hydrologische \
    Fragestellungen" des Fachbereichs 11 Physische Geographie der Goethe Universität Frankfurt """

__author__ = "Florian Herz"
__copyright__ = "Copyright 2019, FH"
__credits__ = ["Florian Herz", "Dr. Hannes Müller Schmied", "Dr. Irene Marzolff"]
__version__ = "1.0"
__maintainer__ = "Florian Herz"
__email__ = "florian13.herz@googlemail.com"
__status__ = "Development"


#Import der Systemmodule
import os
import sys
import arcpy
from arcpy.sa import *

#Einstellungen und Erweiterungen laden

arcpy.CheckOutExtension("Spatial")

arcpy.AddMessage("Systemmodule wurden geladen.")

#Benutzereingaben

basisdaten = arcpy.GetParameterAsText(0)        #C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb
name = arcpy.GetParameterAsText(1)              #SWM_Ergebnisdatenbank_2003-2004_I.gdb
ordner = arcpy.GetParameterAsText(2)            #C:\HiWi_Hydro-GIS
start = int(arcpy.GetParameterAsText(3))        #20030101
ende = int(arcpy.GetParameterAsText(4))         #20041231
RP_faktor = float(arcpy.GetParameterAsText(5))  #0.85
#Abfrage welche Daten gespeichert werden sollen
check_PET = arcpy.GetParameterAsText(6)         #boolean
check_AET = arcpy.GetParameterAsText(7)         #boolean

#Erstellen des Dateipfad der Basisdaten und der Ergebnisdatenbank

ergebnisse = arcpy.CreateFileGDB_management(ordner, name)
arcpy.env.workspace = ordner+r'\{}'.format(name)

arcpy.AddMessage("Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(ordner))

#Zuordnung der Haude Rasterdateien zu den Monatszahlen in einem dictionary

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

# Laden der Rasterdatensätze

tftable = r'{}\TempFeuchte'.format(basisdaten) #Grundlegende Klimadaten
FK = Raster(r'{}\FK_von_L'.format(basisdaten))
RP = FK * RP_faktor
WP = Raster(r'{}\WP'.format(basisdaten))
gewaesser = Raster(r'{}\Gewaesser'.format(basisdaten))
RPdiffWP = RP - WP
Svortag = FK
arcpy.AddMessage("Berechnung und Laden der Rasterdatensätze war erfolgreich.")

#Iteration durch die ausgewählten Klimadaten

with arcpy.da.SearchCursor(tftable, ['TagesID', 'Jahr', 'Monat', 'Tag', 'RelFeu', 'Temp'], \
                           "TagesID >= {0} AND TagesID <= {1}".format(start, ende)) as klimacursor:
    for row in klimacursor:
        
        ID = int(row[0])
        jahr = int(row[1])
        monat = int(row[2])
        tag = int(row[3])
        relf = float(row[4])
        temp = float(row[5])

        
        #Berechnung der PET
        
        Es = 6.1 * 10**((7.5 * (temp)/10.0) / ((temp)/10.0 + 237.2)) #Berechnung des Es
        PET = haude_dic[monat] * Es * (1.0 - relf / 100.0) #Berechnung der PET
        if str(check_PET) == 'true':
            PET.save("PET_{}".format(ID)) #Speichern der einzelnen PET-Rasterdatensätze
        
        #Berechnung der AET

        AET = Con(gewaesser == 1, PET, Con(Svortag >= RP, PET, Con(RPdiffWP == 0, 0, ((Svortag - WP)/RPdiffWP)*PET)))

        if str(check_AET) == 'true':
            AET.save("AET_{}".format(ID)) # Speichern der einzelnen AET-Rasterdatensätze

        
        arcpy.AddMessage("Fertig mit der Berechnung des {0}.{1}.{2}".format(tag, monat, jahr))
        


# Delete cursor objects
del klimacursor






