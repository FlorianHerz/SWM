# -*- coding: cp1252 -*-

"""Soil Water Model in Python für ArcGIS
Lehrveranstaltung "GIS für hydrologische Fragestellungen" des Fachbereich 11 Institut für Physische Geographie der
Johann Wolfgang Goethe Universität Frankfurt am Main.

Dieses einfache Bodenwasser-Modell berechnet für jede Rasterzelle eines Einzugsgebietes die Boden-Wasser-Bilanz.
Ausgabe des Modells ist eine Tabelle mit täglichem Abflussvolumen in m³ für das Einzugsgebiet, sowie auf Wunsch die
berechneten Rasterdatensätze verschiedener Modellparameter. Voraussetzung ist das Vorhandensein der folgenden Datensätze
in einer gemeinsamen File Geodatabase (kursiveNamen werden abgefragt, fettemüssen so vorhanden sein):
Einzugsgebiet(e)- Vektor (Polygon) (Inhalt und Struktur siehe Eingabedatensatz Einzugsgebiet)
Klimadaten-Tabelle (Inhalt und Struktur siehe Eingabedatensatz Klimadaten-Tabelle)
FK_von_L- Raster (Feldkapazität in der effektiven Wurzelzone, in mm)
L_in_metern- Raster (effektive Wurzelzone, in m)
WP - Raster (Welkepunkt, in mm)
Gewaesser- Raster (Gewässerflächenmaske mit Gewässer = 1 und nicht-Gewässer 0)
N_Zeitreihen-Tabelle mit Niederschlagsdaten für jeden Tag und jede Messstation. Die Attributetabelle muss die Felder
Stationsnummer, Tagessumme_mm und TagesID enthalten.
N_Messstationen- Vektor (Punkte); Messstationen des Niederschlags. Die Attributetabelle muss das Feld Stationsnummer
enthalten.
"""

__author__ = "Florian Herz"
__copyright__ = "Copyright 2019, FH"
__credits__ = ["Florian Herz", "Dr. Hannes Müller Schmied",
               "Dr. Irene Marzolff"]
__version__ = "1.0"
__maintainer__ = "Florian Herz"
__email__ = "florian13.herz@googlemail.com"
__status__ = "Development"

#  Import der Systemmodule
import arcpy
from arcpy.sa import *

#  Einstellungen und Erweiterungen laden

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

data = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb'
name = r'SWM_Ergebnisdatenbank_p.gdb'
folder = r'C:\HiWi_Hydro-GIS'
output = arcpy.CreateFileGDB_management(folder, name)
ordner = arcpy.env.workspace = folder+r'\{}'.format(name)

s_init = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L')

date = 20040714
idw_pow = "1"

p_data = r'{}\N_Zeitreihen'.format(data)
p_stations = r'{}\N_Messstationen'.format(data)
p_temp = arcpy.Copy_management(p_stations, "p_temp")

cellsize = arcpy.GetRasterProperties_management(s_init, "CELLSIZEX")
station = []
p_sum = []

arcpy.AddField_management(p_temp, "Tagessumme_mm", "DOUBLE")

with arcpy.da.SearchCursor(p_data, ["Stationsnummer", "Tagessumme_mm", "TagesID"], "TagesID = {}".format(date)) as cursor:
    for row in cursor:
        station.append(row[0])
        p_sum.append(row[1])

for i in range(len(station)):
    with arcpy.da.UpdateCursor(p_temp, ["Stationsnummer", "Tagessumme_mm"]) as p_cursor:
        for row in p_cursor:
            if row[0] == station[i]:
                row[1] = p_sum[i]

            p_cursor.updateRow(row)

idw = Idw(p_temp, "Tagessumme_mm", cellsize, idw_pow, "FIXED 20000 5", "")

idw.save("IDW_{}".format(date))

# Delete cursor objects
del cursor
del p_cursor
arcpy.Delete_management("p_temp")
