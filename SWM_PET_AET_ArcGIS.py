# -*- coding: cp1252 -*-

"""Soil Water Model in Python für ArcGIS
Lehrveranstaltung "GIS für hydrologische Fragestellungen" des Fachbereich 11
Institut für Physische Geographie der Goethe Universität Frankfurt am Main.

Dieses einfache Bodenwasser-Modell berechnet für jede Rasterzelle eines
Einzugsgebietes die Boden-Wasser-Bilanz. Ausgabe des Modells ist eine Tabelle
mit täglichem Abflussvolumen in m³ für das Einzugsgebiet, sowie auf Wunsch die
berechneten Rasterdatensätze verschiedener Modellparameter.
Voraussetzung ist das Vorhandensein der folgenden Datensätze in einer
gemeinsamen File Geodatabase (kursiveNamen werden abgefragt, fettemüssen so
vorhanden sein):
Einzugsgebiet(e)- Vektor (Polygon) (Inhalt und Struktur siehe Eingabedatensatz
Einzugsgebiet)
Klimadaten-Tabelle (Inhalt und Struktur siehe Eingabedatensatz Klimadaten-
Tabelle)
FK_von_L- Raster (Feldkapazität in der effektiven Wurzelzone, in mm)
L_in_metern- Raster (effektive Wurzelzone, in m)
WP - Raster (Welkepunkt, in mm)
Gewaesser- Raster (Gewässerflächenmaske mit Gewässer = 1 und nicht-Gewässer 0)
N_Zeitreihen-Tabelle mit Niederschlagsdaten für jeden Tag und jede Messstation.
Die Attributetabelle muss die Felder Stationsnummer, Tagessumme_mm und
TagesID enthalten.
N_Messstationen- Vektor (Punkte); Messstationen des Niederschlags.
Die Attributetabelle muss das Feld Stationsnummer enthalten. 
"""

__author__ = "Florian Herz"
__copyright__ = "Copyright 2019, FH"
__credits__ = ["Florian Herz", "Dr. Hannes Müller Schmied",
               "Dr. Irene Marzolff"]
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

basisdaten = arcpy.GetParameterAsText(0) #Datatype: workspace
#Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb

name = arcpy.GetParameterAsText(1) #Datatype: string
#Default: SWM_Ergebnisdatenbank_2003-2004_I.gdb

ordner = arcpy.GetParameterAsText(2) #Datatype: folder
#Default: C:\HiWi_Hydro-GIS

start = int(arcpy.GetParameterAsText(3)) #Datatype: long
#Default: 20030101

ende = int(arcpy.GetParameterAsText(4)) #Datatype: long
#Default: 20041231

Sinit = arcpy.GetParameterAsText(5) #Datatype: geodataset
#Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L

RP_faktor = float(arcpy.GetParameterAsText(6))#Datatype: double
#Default: 0.85

#Abfrage welche Daten gespeichert werden sollen
check_PET = str(arcpy.GetParameterAsText(7)) #Datatype: boolean
#Default: false

check_AET = str(arcpy.GetParameterAsText(8)) #Datatype: boolean
#Default: false

#Erstellen des Dateipfad der Basisdaten und der Ergebnisdatenbank

ergebnisse = arcpy.CreateFileGDB_management(ordner, name)
arcpy.env.workspace = ordner+r'\{}'.format(name)

arcpy.AddMessage("Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt." \
                 .format(ordner))

#Zuweisung der Rasterdatensätze

haude_dic = { 1 : Raster(r'{}\Haude_1'.format(basisdaten)), \
              2 : Raster(r'{}\Haude_2'.format(basisdaten)), \
              3 : Raster(r'{}\Haude_3'.format(basisdaten)), \
              4 : Raster(r'{}\Haude_4'.format(basisdaten)), \
              5 : Raster(r'{}\Haude_5'.format(basisdaten)), \
              6 : Raster(r'{}\Haude_6'.format(basisdaten)), \
              7 : Raster(r'{}\Haude_7'.format(basisdaten)), \
              8 : Raster(r'{}\Haude_8'.format(basisdaten)), \
              9 : Raster(r'{}\Haude_9'.format(basisdaten)), \
              10 : Raster(r'{}\Haude_10'.format(basisdaten)), \
              11 : Raster(r'{}\Haude_11'.format(basisdaten)), \
              12 : Raster(r'{}\Haude_12'.format(basisdaten))}

tftable = r'{}\TempFeuchte'.format(basisdaten) #Grundlegende Klimadatentabelle
FK = Raster(r'{}\FK_von_L'.format(basisdaten))
RP = FK * RP_faktor
WP = Raster(r'{}\WP'.format(basisdaten))
gewaesser = Raster(r'{}\Gewaesser'.format(basisdaten))
RPdiffWP = RP - WP
Svortag = Sinit
arcpy.AddMessage("Berechnung der Rasterdatensätze war erfolgreich.")

#Iteration durch die ausgewählten Klimadaten

with arcpy.da.SearchCursor(tftable, ['TagesID', 'Jahr', 'Monat', 'Tag',
                                     'RelFeu', 'Temp'], "TagesID >= {0} AND \
                                     TagesID <= {1}".format(start, ende)) \
                                     as klimacursor:
    for row in klimacursor:
        
        ID = int(row[0])
        jahr = int(row[1])
        monat = int(row[2])
        tag = int(row[3])
        relf = float(row[4])
        temp = float(row[5])

        
        #Berechnung der PET
        
        Es = 6.1 * 10**((7.5 * (temp)/10.0) / ((temp)/10.0 + 237.2))
        PET = haude_dic[monat] * Es * (1.0 - relf / 100.0) #Berechnung der PET

        #Speichern der täglichen PET Rasterdateien
        if check_PET == 'true':
            PET.save("PET_{}".format(ID))
        
        #Berechnung der AET

        AET = Con(gewaesser == 1, PET, Con(Svortag >= RP, PET, \
                Con(RPdiffWP == 0, 0, ((Svortag - WP)/RPdiffWP)*PET)))

        #Speichern der täglichen AET Rasterdateien
        if check_AET == 'true':
            AET.save("AET_{}".format(ID))

        
        arcpy.AddMessage("Fertig mit der Berechnung des {0}.{1}.{2}" \
                         .format(tag, monat, jahr))
        
# Delete cursor objects
del klimacursor
