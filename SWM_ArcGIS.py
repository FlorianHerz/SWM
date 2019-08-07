# -*- coding: cp1252 -*-

"""Soil Water Model in Python für ArcGIS
Lehrveranstaltung "GIS für hydrologische Fragestellungen" des Fachbereich 11 Institut für Physische Geographie der
Johann Wolfgang Goethe Universität Frankfurt am Main.

Dieses einfache Bodenwasser-Modell berechnet für jede Rasterzelle eines Einzugsgebietes ein Boden-Wasser-Bilanz. Ausgabe
des Modells ist eine Tabelle mit täglichem Abflussvolumen in m³ für das Einzugsgebietsowie aufWunsch die berechneten
Rasterdatensätze verschiedener Modellparameter.
Voraussetzung ist das Vorhandensein der folgenden Datensätze in einer gemeinsamen File Geodatabase:
Einzugsgebiet(e)- Vektor (Polygon) (Inhalt und Struktur siehe Eingabedatensatz Einzugsgebiet)
TempFeuchte - Tabelle mit den Klimadaten. Die Attributetabelle muss die Felder Tagesid, Jahr, Monat, Tag, RelFeu, und
Temp enthalten
fk_von_l - Raster (Feldkapazität in der effektiven Wurzelzone, in mm)
L_in_metern- Raster (effektive Wurzelzone, in m)
wp - Raster (Welkepunkt, in mm)
Gewaesser- Raster (Gewässerflächenmaske mit Gewässer = 1 und nicht-Gewässer 0 )
N_Zeitreihen- Tabelle mit Niederschlagsdaten für jeden Tag und jede Messstation. Die Attributetabelle muss die Felder
Stationsnummer, Tagessumme_mm und TagesID enthalten
N_Messstationen- Vektor (Punkte); Messstationen des Niederschlags. Die Attributetabelle muss das Feld Stationsnummer
enthalten
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
import time

arcpy.CheckOutExtension("Spatial")
arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Systemmodule geladen.")

#  Funktionen


def get_pet(haude_factor, temperature, humidity):
    """Funktion zur Berechnung der PET nach Haude. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen sind die
    Rasterdatensätze zum Haudefaktor, zur Temperatur und zur relativen Luftfeuchtigkeit. Die PET wird anhand einer
    Formel berechnet.
    """
    pet_raster = haude_factor * (6.1 * 10 ** ((7.5 * temperature) / (temperature + 237.2))) * (1.0 - humidity / 100.0)
    return pet_raster


def get_aet(pet_raster, water_raster, s_pre_raster, rp_raster, rpwp_dif_raster, wp_raster):
    """Funktion zur Berechnung der AET. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen sind die Rasterdatensätze
    zur PET, der Gewässerpixel, zum Bodenwasserspeicher des Vortages, zum Reduktionspunkt, zur Reduktionspunkt-Welke-
    punkt-Differenz und zum Welkepunkt. Die AET wird entweder über Abfragen bestimmt, oder über eine Formel berechnet.
    Der Wert der AET entspricht bei Gewässerpixeln, sowie Zellen in denen der Bodenwassersgehalt oberhalb des
    ReduktionspunktesP liegt der PET. Die AET ist 0, wenn der Reduktionspunkt gleich dem Welkepunkt ist. Bei den anderen
    Pixeln wird die AET nach einer Formel kalkuliert.
    """
    aet_raster = Con(water_raster == 1, pet_raster,
                     Con(s_pre_raster >= rp_raster, pet_raster,
                         Con(rpwp_dif_raster == 0, 0, ((s_pre_raster - wp_raster) / rpwp_dif_raster) * pet_raster)))
    return aet_raster


def get_precipitation(timeseries, p_stations_temp, date, idw_pow, rastercellsize):
    """Funktion zur Interpolation des Niederschlags anhand der IDW-Methode. Die Ausgabe ist ein Rasterdatensatz.
    Eingabevariablen sind die Tabellen zu den Niederschlagsmesswerten und Niederschlagsstationen, die TagesID, der
    Exponent der IDW-Methode und die Größe der Rasterzellen. Die Berechnung erfolgt über das IDW-Tool der Erweiterung
    "Spatial Analyst". Zunächst werden die Niederschlagswerte eines Tages je Station in eine Liste geschrieben. Der
    Index der Niederschlagswerte je Station entspricht dem Index des dazugehörigen Stationsnamens in einer zweiten
    Liste. Diese Listen werden mittels eines Searchcursors erstellt. Über einen Updatecursor werden die Niederschlags-
    messwerte in eine Tabelle geschrieben, die Informationen zum Stationsnamen und deren Koordinaten enthalten. Der
    Interpolationsradius ist auf das EZG Schotten 2 angepasst.
    """
    station = []  # Liste der Stationsnamen, um den Index der N-Messerwerte der zweiten Liste den Stationen zuzuordnen
    p_sum = []  # Liste der N-Messwerte eines Tages

    # Speichert alle N-Messwerte mit den dazugehörigen Stationsnamen eines Tages in die Listen
    with arcpy.da.SearchCursor(timeseries, ["Stationsnummer", "Tagessumme_mm", "TagesID"],
                               "TagesID = {}".format(date)) as p_cursor:
        for information in p_cursor:
            station.append(information[0])
            p_sum.append(information[1])

    # Aktualisiert den N-Messwert je Tag in der temporären Niederschlagsdatei
    for i in range(len(station)):
        with arcpy.da.UpdateCursor(p_stations_temp, ["Stationsnummer", "Tagessumme_mm"]) as p_incursor:
            for p_row in p_incursor:
                if p_row[0] == station[i]:
                    p_row[1] = p_sum[i]

                p_incursor.updateRow(p_row)

    # Niederschlagsinterpolation mittels Inverse Distance Weighting (IDW)
    idw = Idw(p_stations_temp, "Tagessumme_mm", rastercellsize, idw_pow, RadiusFixed(20000.00000, 5), "")

    # Löschen der temporären Objekte innerhalb der Funktion
    del p_cursor
    del p_incursor

    return idw


def overflow(p_raster, s_pre_raster, fk_raster):
    """Funktion zur Berechnung des direkten Abflusses des Niederschlags bei gesättigtem Bodenwasserspeicher. Die Ausgabe
    ist ein Rasterdatensatz. Eingabevariablen sind die Rasterdateien zum Niederschlag, zum Bodenwasserspeicher des Vor-
    tages und zur Feldkapazität.
    """
    r_overflow = Con(p_raster + s_pre_raster > fk_raster,
                     p_raster + s_pre_raster - fk_raster, 0)
    return r_overflow


def runoff_land(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster):
    """Funktion zur Berechnung des Abflusses von Landpixeln. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen sind
    die für die Overflow-Funktion benötigten Variablen, sowie die Rasterdatensätze zu den Gewässerpixeln, zum Lamda-Wert
    und zum Welkepunkt.
    """
    r_land = Con(water_raster == 0,
                 lambda_wert * ((s_pre_raster - wp_raster) ** 2) + overflow(p_raster, s_pre_raster, fk_raster), 0)
    return r_land


def get_runoff(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster, pet_raster):
    """Funktion zur Berechnung des Gesamtabflusses. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen sind die für
    die Overflow-Funktion und für die runoff-land-Funktion benötigten Variablen, sowie der Rasterdatensatz der PET. Der
    Gesamtabfluss ist die Summe des Abflusses von Landpixeln und von Gewässerpixeln. Ersterer wird über die runoff-land-
    Funktion berechnet, während letzterer der Differenz aus dem Niederschlag und der PET entspricht.
    """
    r = runoff_land(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster) +\
        water_raster * (p_raster - pet_raster)
    return r


def get_soilwater(s_pre_raster, p_raster, aet_raster, runoff_raster):
    """Funktion zur Berechnung des aktuellen Bodenwasserspeichers. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen
    sind die Rasterdatensätze des Bodenwasserspeichers des Vortages, des Niederschlags, der AET und des Gesamtabflusses.
    Die Berechnung erfolgt anhand der allgemeinen Wasserhaushaltsgleichung.
    """
    soilwater = s_pre_raster + p_raster - aet_raster - runoff_raster
    return soilwater


#  Benutzereingaben

data = arcpy.GetParameterAsText(0)  # Datatype: workspace #  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb
basin = arcpy.GetParameterAsText(1)  # Datatype: Featurelayer
#  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\EZG_Schotten2_Vektor
basin_id = arcpy.GetParameterAsText(2)  # String #  Default: Id
name = arcpy.GetParameterAsText(3)  # Datatype: string #  Default: SWM_Ergebnisdatenbank_2003-2004_I.gdb
folder = arcpy.GetParameterAsText(4)  # Datatype: folder #  Default: C:\HiWi_Hydro-GIS
start = int(arcpy.GetParameterAsText(5))  # Datatype: long #  Default: 20030101
end = int(arcpy.GetParameterAsText(6))  # Datatype: long #  Default: 20041231
s_init = Raster(arcpy.GetParameterAsText(7))  # Datatype: geodataset
#  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L
rp_factor = float(arcpy.GetParameterAsText(8))  # Datatype: double #  Default: 0.85
c = int(arcpy.GetParameterAsText(9))  # Datatype: short #  Default: 150
idw_exponent = arcpy.GetParameterAsText(10)  # Datatype: string #  Default: 1
outname = arcpy.GetParameterAsText(16)  # Datatype: string # Default: Durchfluss

#  Abfrage welche Daten gespeichert werden sollen
check_pet = str(arcpy.GetParameterAsText(11))  # Datatype: boolean #  Default: false
check_aet = str(arcpy.GetParameterAsText(12))  # Datatype: boolean #  Default: false
check_p = str(arcpy.GetParameterAsText(13))  # Datatype: boolean #  Default: false
check_r = str(arcpy.GetParameterAsText(14))  # Datatype: boolean #  Default: false
check_s = str(arcpy.GetParameterAsText(15))  # Datatype: boolean #  Default: false

#  Erstellen der Ergebnisdatenbank und Wahl des Arbeitsverzeichnisses

arcpy.env.overwriteOutput = True
arcpy.env.extent = s_init
arcpy.CreateFileGDB_management(folder, name)
arcpy.env.workspace = r'{0}\{1}'.format(folder, name)
arcpy.env.scratchWorkspace = arcpy.env.workspace

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(folder))

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

climatedata = r'{}\TempFeuchte'.format(data)
fk = Raster(r'{}\fk_von_L'.format(data))
rp = fk * rp_factor
wp = Raster(r'{}\wp'.format(data))
water = Raster(r'{}\Gewaesser'.format(data))
rpwp_dif = rp - wp
s_pre = s_init
p_data = r'{}\N_Zeitreihen'.format(data)
# Kopie der Tabelle der Niederschlagsstationen, wird am Ende des Skripts gelöscht
p_temp = arcpy.Copy_management(r'{}\N_Messstationen'.format(data), "p_temp")
arcpy.AddField_management(p_temp, "Tagessumme_mm", "DOUBLE")  # Hinzufügen eines neuen Attributefelds für die Messwerte
cellsize = arcpy.Describe(s_init).children[0].meanCellHeight
lambda_parameter = (c / (Raster(r'{}\L_in_metern'.format(data)) * 1000)**2)
runoff_daily = []
id_list = []
date_list = []

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Berechnung der Rasterdatensaetze war erfolgreich.")

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

        # Berechnungen je Tag
        pet = get_pet(haude_dic[month], temp, humid)
        aet = get_aet(pet, water, s_pre, rp, rpwp_dif, wp)
        precipitation = get_precipitation(p_data, p_temp, id_day, idw_exponent, cellsize)
        runoff = get_runoff(water, lambda_parameter, wp, precipitation, s_pre, fk, pet)
        s = get_soilwater(s_pre, precipitation, aet, runoff)

        #  Speichern der täglichen Rasterdateien wenn ausgewählt
        if check_pet == 'true':
            pet.save("PET_{}".format(id_day))

        if check_aet == 'true':
            aet.save("AET_{}".format(id_day))

        if check_p == 'true':
            precipitation.save("Niederschlag_{}".format(id_day))

        if check_r == 'true':
            runoff.save("R_{}".format(id_day))

        if check_s == 'true':
            s.save("S_{}".format(id_day))

        s_pre = s  # Überschreiben des Bodenwasserspeichers des Vortages

        # Summierung des Gesamtabflusses des Einzugsgebietes in der Tabelle "Q_day"
        ZonalStatisticsAsTable(basin, basin_id, runoff, "Q_day", "true", "SUM")

        # Berechnung des Durchflusses
        with arcpy.da.SearchCursor("Q_day", "Sum") as r_cursor:
            for r_sum in r_cursor:
                runoff_daily.append(r_sum[0] * 0.001 * cellsize ** 2)
        id_list.append(id_day)
        date_list.append("{0}.{1}.{2}".format(day, month, year))

        arcpy.AddMessage(time.strftime("%H:%M:%S: ") +
                         "Fertig mit der Berechnung des {0}.{1}.{2}".format(day, month, year))

# Erstellen der Tabelle mit den Durchflusswerten des Untersuchungszeitraumes
arcpy.CreateTable_management(arcpy.env.workspace, outname)
arcpy.AddField_management("{0}\{1}".format(arcpy.env.workspace, outname), "Datum", "TEXT")
arcpy.AddField_management("{0}\{1}".format(arcpy.env.workspace, outname), "Q", "DOUBLE")
arcpy.AddField_management("{0}\{1}".format(arcpy.env.workspace, outname), "Tages_ID", "LONG")

q_cursor = arcpy.da.InsertCursor("{0}\{1}".format(arcpy.env.workspace, outname), ["Datum", "Q", "Tages_ID"])

for i in range(len(id_list)):
    row = [date_list[i], runoff_daily[i], id_list[i]]
    q_cursor.insertRow(row)

# Löschen der temporären Objekte
del cursor
del q_cursor
del row
arcpy.Delete_management("Q_day")
arcpy.Delete_management("p_temp")

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Modellierung Abgeschlossen.")
