# -*- coding: cp1252 -*-

"""Soil Water Model in Python f�r ArcGIS
Lehrveranstaltung "GIS f�r hydrologische Fragestellungen" des Fachbereich 11 Institut f�r Physische Geographie der
Johann Wolfgang Goethe Universit�t Frankfurt am Main.

Dieses einfache Bodenwasser-Modell berechnet f�r jede Rasterzelle eines Einzugsgebietes die Boden-Wasser-Bilanz.
Ausgabe des Modells ist eine Tabelle mit t�glichem Abflussvolumen in m� f�r das Einzugsgebiet, sowie auf Wunsch die
berechneten Rasterdatens�tze verschiedener Modellparameter. Voraussetzung ist das Vorhandensein der folgenden Datens�tze
in einer gemeinsamen File Geodatabase (kursiveNamen werden abgefragt, fettem�ssen so vorhanden sein):
Einzugsgebiet(e)- Vektor (Polygon) (Inhalt und Struktur siehe Eingabedatensatz Einzugsgebiet)
Klimadaten-Tabelle (Inhalt und Struktur siehe Eingabedatensatz Klimadaten-Tabelle)
FK_von_L- Raster (Feldkapazit�t in der effektiven Wurzelzone, in mm)
L_in_metern- Raster (effektive Wurzelzone, in m)
WP - Raster (Welkepunkt, in mm)
Gewaesser- Raster (Gew�sserfl�chenmaske mit Gew�sser = 1 und nicht-Gew�sser 0)
N_Zeitreihen-Tabelle mit Niederschlagsdaten f�r jeden Tag und jede Messstation. Die Attributetabelle muss die Felder
Stationsnummer, Tagessumme_mm und TagesID enthalten.
N_Messstationen- Vektor (Punkte); Messstationen des Niederschlags. Die Attributetabelle muss das Feld Stationsnummer
enthalten.
"""

__author__ = "Florian Herz"
__copyright__ = "Copyright 2019, FH"
__credits__ = ["Florian Herz", "Dr. Hannes M�ller Schmied",
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

data = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb'
name = r'SWM_Ergebnisdatenbank_XVI.gdb'
folder = r'C:\HiWi_Hydro-GIS'
output = arcpy.CreateFileGDB_management(folder, name)
ordner = arcpy.env.workspace = folder+r'\{}'.format(name)

s_init = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L')

date = 20040609
idw_pow = "1"

p_data = r'{}\N_Zeitreihen'.format(data)
p_stations = r'{}\N_Messstationen'.format(data)
cellsize = arcpy.GetRasterProperties_management(s_init, "CELLSIZEX")

arcpy.MakeQueryTable_management(in_table=[p_data, p_stations], out_table="p_temp", in_key_field_option="USE_KEY_FIELDS",
                                in_key_field="N_Messstationen.Stationsnummer;N_Zeitreihen.Stationsnummer", in_field="",
                                where_clause="N_Zeitreihen.Stationsnummer = N_Messstationen.Stationsnummer \
                                AND N_Zeitreihen.TagesID = {}".format(date))

arcpy.gp.Idw_sa("p_temp", "N_Zeitreihen.Tagessumme_mm", "{0}/Idw_p_{1}".format(arcpy.env.workspace, date),
                cellsize, idw_pow, "FIXED 20000 5", "")

arcpy.Delete_management("p_temp")