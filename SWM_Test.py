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
import os
import time

arcpy.CheckOutExtension("Spatial")
arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Systemmodule geladen.")


#  Funktionen


def get_pet(haude_factor, temperature, humidity, bool_pet, date):
    """
    Funktion zur Berechnung der PET nach Haude. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen sind die
    Rasterdatensätze zum Haudefaktor, zur Temperatur und zur relativen Luftfeuchtigkeit. Die PET wird anhand einer
    Formel berechnet.
    :param haude_factor:
    :param temperature:
    :param humidity:
    :param bool_pet:
    :param date:
    :return:
    """
    pet_raster = haude_factor * (6.1 * 10 ** ((7.5 * temperature) / (temperature + 237.2))) * (1.0 - humidity / 100.0)

    if bool_pet:
        pet_raster.save("PET_{}".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "PET ausgeführt.")
    return pet_raster


def get_aet(pet_raster, water_raster, s_pre_raster, rp_raster, rpwp_dif_raster, wp_raster, bool_aet, date):
    """
    Funktion zur Berechnung der AET. Die Ausgabe ist ein Rasterdatensatz. Eingabevariablen sind die Rasterdatensätze
    zur PET, der Gewässerpixel, zum Bodenwasserspeicher des Vortages, zum Reduktionspunkt, zur Reduktionspunkt-Welke-
    punkt-Differenz und zum Welkepunkt. Die AET wird entweder über Abfragen bestimmt, oder über eine Formel berechnet.
    Der Wert der AET entspricht bei Gewässerpixeln, sowie Zellen in denen der Bodenwassersgehalt oberhalb des
    ReduktionspunktesP liegt der PET. Die AET ist 0, wenn der Reduktionspunkt gleich dem Welkepunkt ist. Bei den anderen
    Pixeln wird die AET nach einer Formel kalkuliert.
    :param pet_raster:
    :param water_raster:
    :param s_pre_raster:
    :param rp_raster:
    :param rpwp_dif_raster:
    :param wp_raster:
    :param bool_aet: check_aet
    :param date: TagesID
    :return:
    """
    aet_raster = Con(water_raster == 1, pet_raster,
                     Con(s_pre_raster >= rp_raster, pet_raster,
                         Con(rpwp_dif_raster == 0, 0, ((s_pre_raster - wp_raster) / rpwp_dif_raster) * pet_raster)))

    if bool_aet:
        aet_raster.save("AET_{}".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "AET ausgeführt.")
    return aet_raster


def get_precipitation(timeseries, p_stations_temp, date, idw_pow, rastercellsize, bool_p):
    """
    Funktion zur Interpolation des Niederschlags anhand der IDW-Methode. Die Ausgabe ist ein Rasterdatensatz.
    Eingabevariablen sind die Tabellen zu den Niederschlagsmesswerten und Niederschlagsstationen, die TagesID, der
    Exponent der IDW-Methode und die Größe der Rasterzellen. Die Berechnung erfolgt über das IDW-Tool der Erweiterung
    "Spatial Analyst". Zunächst werden die Niederschlagswerte eines Tages je Station in eine Liste geschrieben. Der
    Index der Niederschlagswerte je Station entspricht dem Index des dazugehörigen Stationsnamens in einer zweiten
    Liste. Diese Listen werden mittels eines Searchcursors erstellt. Über einen Updatecursor werden die Niederschlags-
    messwerte in eine Tabelle geschrieben, die Informationen zum Stationsnamen und deren Koordinaten enthalten. Der
    Interpolationsradius ist auf das EZG Schotten 2 angepasst.
    :param timeseries:
    :param p_stations_temp:
    :param date:
    :param idw_pow:
    :param rastercellsize:
    :param bool_p:
    :return:
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
    del information
    del p_incursor
    del p_row

    if bool_p:
        idw.save("Niederschlag_{}".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "P ausgeführt.")

    return idw


def get_runoff(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster, pet_raster, bool_r, date):
    """
    Funktion zur Berechnung des Gesamtabflusses. Der Gesamtabfluss ist die Summe des Abflusses von Landpixeln und von
    Gewässerpixeln. Ersterer wird über die runoff_land-Funktion berechnet, während letzterer der Differenz aus dem
    Niederschlag und der PET entspricht.
    :param water_raster:
    :param lambda_wert:
    :param wp_raster:
    :param p_raster:
    :param s_pre_raster:
    :param fk_raster:
    :param pet_raster:
    :param bool_r:
    :param date:
    :return:
    """
    r_overflow = Con(p_raster + s_pre_raster > fk_raster, p_raster + s_pre_raster - fk_raster, 0)

    r_land = Con(water_raster == 0, lambda_wert * ((s_pre_raster - wp_raster) ** 2) + r_overflow, 0)

    r = r_land + water_raster * Con(p_raster > pet_raster, p_raster - pet_raster, 0)

    if bool_r:
        r.save("R_{}".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "gesamtabfluss ausgeführt.")
    return r


def get_soilwater(s_pre_raster, p_raster, aet_raster, runoff_raster, bool_s, date):
    """
    Funktion zur Berechnung des aktuellen Bodenwasserspeichers. Die Berechnung erfolgt anhand der allgemeinen Wasser-
    haushaltsgleichung.
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Raster)
    :param p_raster: Ausgaberaster der Niederschlagsinterpolation
    :param aet_raster: Ausgaberaster der AET
    :param runoff_raster: Ausgaberaster des Gesamtabflusses
    :param bool_s: check_s
    :param date: TagesID
    :return: Rasterdatensatz des aktuellen Bodenwasserspeichers
    """
    soilwater = s_pre_raster + p_raster - aet_raster - runoff_raster

    if bool_s:
        soilwater.save("S_{}".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Bodenwasserspeicher ausgeführt.")
    return soilwater


#def get_q_m3(runoff_raster, rastercellsize):
#    """
#    Berechnet den Gesamtabfluss des Einzugsgebietes in m^3.
#    :param runoff_raster:
#    :param rastercellsize: Größe der Rasterzellen
#    :return: täglicher Gesamtabfluss in m^3
#    """
#    import numpy
#    numpy.array = arcpy.RasterToNumPyArray(runoff_raster, nodata_to_value=0)
#    r_sum = numpy.array.sum()
#    r_m3 = (r_sum * 0.001 * rastercellsize ** 2)
#
#    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Q in m3 berechnet.")
#
#    return r_m3


def set_resulttable(workspace, tablename, list_date, list_runoff, list_id):
    """
    Erstellt die Ergebnistabelle in der Ergebnisdatenbank.
    :param workspace: Dateipfad der Ergebnisdatenbank.
    :param tablename: Name der Ergebnistabelle
    :param list_date: Liste der Daten des Untersuchungszeitraums
    :param list_runoff: Liste der täglichen Gesamtabflüsse
    :param list_id: Liste der Tages_IDs des Untersuchungszeitraums
    :return: Ergebnistabelle
    """
    if not(os.path.isdir(workspace + "\\" + tablename)):
        arcpy.CreateTable_management(workspace, tablename)
        arcpy.AddField_management(workspace + "\\" + tablename, "Datum", "TEXT")
        arcpy.AddField_management(workspace + "\\" + tablename, "Q", "DOUBLE")
        arcpy.AddField_management(workspace + "\\" + tablename, "Tages_ID", "LONG")

    q_cursor = arcpy.da.InsertCursor(workspace + "\\" + tablename, ["Datum", "Q", "Tages_ID"])

    for i in range(len(list_id)):
        result = [list_date[i], list_runoff[i], list_id[i]]
        q_cursor.insertRow(result)

    del q_cursor
    del result

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Ergebnistabelle erstellt.")


#  Benutzervorgaben

data = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb'
basin = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\EZG_Schotten2_Vektor'
basin_id = "Id"
name = r'PY_SWM_2004_test.gdb'
folder = r'C:\HiWi_Hydro-GIS'
start = 20040101
end = 20040110
s_init = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L')
rp_factor = 0.85
c = 150
idw_exponent = "1"
check_pet = 0
check_aet = 0
check_p = 0
check_r = 0
check_s = 0
outname = "Ergebnistabelle_{}_{}".format(start, end)
#  Erstellen der Ergebnisdatenbank und Wahl des Arbeitsverzeichnisses

arcpy.env.overwriteOutput = True
arcpy.env.extent = basin
if not(os.path.isdir(folder + "\\" + name)):
    arcpy.CreateFileGDB_management(folder, name)
    arcpy.AddMessage(
        time.strftime("%H:%M:%S: ") + "Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(folder))
arcpy.env.workspace = folder + "\\" + name
arcpy.env.scratchWorkspace = r'C:\HiWi_Hydro-GIS\Scratch'
#  arcpy.env.mask = basin

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Vorhandene Ergebnisdatenbank wird überschrieben.")

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
cellsize = s_init.meanCellHeight
lambda_parameter = (c / (Raster(r'{}\L_in_metern'.format(data)) * 1000) ** 2)
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
        temp = float(row[5]) / 10.0

        # Berechnungen je Tag
        pet = get_pet(haude_dic[month], temp, humid, check_pet, id_day)
        aet = get_aet(pet, water, s_pre, rp, rpwp_dif, wp, check_aet, id_day)
        precipitation = get_precipitation(p_data, p_temp, id_day, idw_exponent, cellsize, check_p)
        #  runoff = get_runoff(water, lambda_parameter, wp, precipitation, s_pre, fk, pet, check_r, id_day)
        r_overflow = Con(precipitation + s_pre > fk, precipitation + s_pre - fk, 0)

        r_land = Con(water == 0, lambda_parameter * ((s_pre - wp) ** 2) + r_overflow, 0)

        r = r_land + water * Con(precipitation > pet, precipitation - pet, 0)
        arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Abfluss berechnet.")
        ZonalStatisticsAsTable(basin, basin_id, r, "Q_day", "true", "SUM")
        #  runoff_m3 = aet.mean * 9309  # * 0.001 * cellsize * cellsize
        arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Abfluss summiert in m3.")
        s = get_soilwater(s_pre, precipitation, aet, r, check_s, id_day)

        s_pre = s  # Überschreiben des Bodenwasserspeichers des Vortages

        #  runoff_daily.append(round(runoff_m3, 4))  # Berechnung des Durchflusses
        with arcpy.da.SearchCursor("Q_day", "Sum") as r_cursor:
            for r_sum in r_cursor:
                runoff_daily.append(r_sum[0] * 0.001 * cellsize ** 2)
        print(runoff_daily)
        id_list.append(id_day)
        date_list.append("{0}.{1}.{2}".format(day, month, year))

        arcpy.AddMessage(time.strftime("%H:%M:%S: ") +
                         "Fertig mit der Berechnung des {0}.{1}.{2}".format(day, month, year))

# Erstellen der Tabelle mit den Durchflusswerten des Untersuchungszeitraumes
set_resulttable(arcpy.env.workspace, outname, date_list, runoff_daily, id_list)

# Löschen der temporären Objekte
del cursor
del row
del r_cursor
arcpy.Delete_management("p_temp")

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Modellierung Abgeschlossen.")
