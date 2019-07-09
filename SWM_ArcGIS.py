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

arcpy.AddMessage("Systemmodule wurden geladen.")

#  Benutzereingaben

data = arcpy.GetParameterAsText(0)  # Datatype: workspace #  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb

name = arcpy.GetParameterAsText(1)  # Datatype: string #  Default: SWM_Ergebnisdatenbank_2003-2004_I.gdb

folder = arcpy.GetParameterAsText(2)  # Datatype: folder #  Default: C:\HiWi_Hydro-GIS

start = int(arcpy.GetParameterAsText(3))  # Datatype: long #  Default: 20030101

end = int(arcpy.GetParameterAsText(4))  # Datatype: long #  Default: 20041231

s_init = Raster(arcpy.GetParameterAsText(5))  # Datatype: geodataset
#  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L

rp_factor = float(arcpy.GetParameterAsText(6))  # Datatype: double #  Default: 0.85

#  Abfrage welche Daten gespeichert werden sollen
check_pet = str(arcpy.GetParameterAsText(7))  # Datatype: boolean #  Default: false

check_aet = str(arcpy.GetParameterAsText(8))  # Datatype: boolean #  Default: false

#  Erstellen der Ergebnisdatenbank und Wahl des Arbeitsverzeichnisses

arcpy.CreateFileGDB_management(folder, name)
arcpy.env.workspace = r'{0}\{1}'.format(folder, name)
arcpy.env.scratchWorkspace = arcpy.env.workspace

arcpy.AddMessage("Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(folder))

#  Zuweisung der Rasterdatensätze

haude_dic = {1: Raster(r'{}\Haude_1'.format(data)),
             2: Raster(r'{}\Haude_2'.format(data)),
             3: Raster(r'{}\Haude_3'.format(data)),
             4: Raster(r'{}\Haude_4'.format(data)),
             5: Raster(r'{}\Haude_5'.format(data)),
             6: Raster(r'{}\Haude_6'.format(data)),
             7: Raster(r'{}\Haude_7'.format(data)),
             8: Raster(r'{}\Haude_8'.format(data)),
             9: Raster(r'{}\Haude_9'.format(data)),
             10: Raster(r'{}\Haude_10'.format(data)),
             11: Raster(r'{}\Haude_11'.format(data)),
             12: Raster(r'{}\Haude_12'.format(data))}

climatedata = r'{}\TempFeuchte'.format(data)  # Grundlegende Klimadatentabelle
fk = Raster(r'{}\fk_von_L'.format(data))
rp = fk * rp_factor
wp = Raster(r'{}\wp'.format(data))
gewaesser = Raster(r'{}\Gewaesser'.format(data))
rpwp_dif = rp - wp
s_pre = s_init

arcpy.AddMessage("Berechnung der Rasterdatensätze war erfolgreich.")

#  Iteration durch die ausgewählten Klimadaten

with arcpy.da.SearchCursor(climatedata, ['Tagesid', 'Jahr', 'Monat', 'Tag', 'RelFeu', 'Temp'], "Tagesid >= {0} AND\
                                     Tagesid <= {1}".format(start, end)) as cursor:
    for row in cursor:
        
        id_day = int(row[0])
        year = int(row[1])
        month = int(row[2])
        day = int(row[3])
        humid = float(row[4])
        temp = float(row[5])/10.0

        # Berechnung der PET

        es = 6.1 * 10**((7.5 * temp) / (temp + 237.2))
        pet = haude_dic[month] * es * (1.0 - humid / 100.0)  # Berechnung der PET

        #  Speichern der täglichen PET Rasterdateien
        if check_pet == 'true':
            pet.save("PET_{}".format(id_day))
        
        # Berechnung der AET

        aet = Con(gewaesser == 1, pet, Con(s_pre >= rp, pet, Con(rpwp_dif == 0, 0, ((s_pre - wp)/rpwp_dif)*pet)))

        # Speichern der täglichen AET Rasterdateien
        if check_aet == 'true':
            aet.save("AET_{}".format(id_day))

        arcpy.AddMessage("Fertig mit der Berechnung des {0}.{1}.{2}".format(day, month, year))
        
# Delete cursor objects
del cursor
