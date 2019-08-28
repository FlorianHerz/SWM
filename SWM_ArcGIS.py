# -*- coding: cp1252 -*-

"""Soil Water Model in Python f�r ArcGIS
Lehrveranstaltung "GIS f�r hydrologische Fragestellungen" des Fachbereich 11 Institut f�r Physische Geographie der
Johann Wolfgang Goethe Universit�t Frankfurt am Main.

Dieses einfache Bodenwasser-Modell berechnet f�r jede Rasterzelle eines Einzugsgebietes ein Boden-Wasser-Bilanz. Ausgabe
des Modells ist eine Tabelle mit t�glichem Abflussvolumen in m� f�r das Einzugsgebietsowie aufWunsch die berechneten
Rasterdatens�tze verschiedener Modellparameter.
Voraussetzung ist das Vorhandensein der folgenden Datens�tze in einer gemeinsamen File Geodatabase:
Einzugsgebiet(e)- Vektor (Polygon) (Inhalt und Struktur siehe Eingabedatensatz Einzugsgebiet)
TempFeuchte - Tabelle mit den Klimadaten. Die Attributetabelle muss die Felder Tagesid, Jahr, Monat, Tag, RelFeu, und
Temp enthalten
fk_von_l - Raster (Feldkapazit�t in der effektiven Wurzelzone, in mm)
L_in_metern- Raster (effektive Wurzelzone, in m)
wp - Raster (Welkepunkt, in mm)
Gewaesser- Raster (Gew�sserfl�chenmaske mit Gew�sser = 1 und nicht-Gew�sser 0 )
N_Zeitreihen- Tabelle mit Niederschlagsdaten f�r jeden Tag und jede Messstation. Die Attributetabelle muss die Felder
Stationsnummer, Tagessumme_mm und TagesID enthalten
N_Messstationen- Vektor (Punkte); Messstationen des Niederschlags. Die Attributetabelle muss das Feld Stationsnummer
enthalten
"""

__author__ = "Florian Herz"
__copyright__ = "Copyright 2019, FH"
__credits__ = ["Florian Herz", "Dr. Hannes M�ller Schmied",
               "Dr. Irene Marzolff"]
__version__ = "1.0"
__maintainer__ = "Florian Herz"
__email__ = "florian13.herz@googlemail.com"
__status__ = "Development"

########################################################################################################################
#  Import der Systemmodule und Aktivierung der ArcGIS Erweiterungen
########################################################################################################################

import arcpy
from arcpy.sa import *
import time

arcpy.CheckOutExtension("Spatial")
arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Systemmodule geladen.")

########################################################################################################################
#  Definition der Funktionen
########################################################################################################################


def get_pet(haude_factor, temperature, humidity, bool_pet, date):
    """
    Funktion zur Berechnung der potentiellen Evapotranspiration (PET) nach Haude.
    :param haude_factor: Haudefaktor des entsprechenden Monats (Typ: raster)
    :param temperature: Tagestemperatur in Grad Celsius (Typ: float)
    :param humidity: relative Luftfeuchte (Typ: float)
    :param bool_pet: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param date: Tages ID (Typ: integer)
    :return: PET des Tages (Typ: raster)
    """
    pet_raster = haude_factor * (6.1 * 10 ** ((7.5 * temperature) / (temperature + 237.2))) * (1.0 - humidity / 100.0)

    if bool_pet:
        pet_raster.save("PET_{}".format(date))

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "PET ausgef�hrt.")
    return pet_raster


def get_aet(pet_raster, water_raster, s_pre_raster, rp_raster, rpwp_dif_raster, wp_raster, bool_aet, date):
    """
    Funktion zur Berechnung der AET. Die AET wird entweder �ber Abfragen bestimmt, oder �ber eine Formel berechnet.
    Der Wert der AET entspricht bei Gew�sserpixeln, sowie Zellen in denen der Bodenwassersgehalt oberhalb des
    Reduktionspunktes liegt, der PET. Die AET ist 0, wenn der Reduktionspunkt gleich dem Welkepunkt ist. Bei den anderen
    Pixeln wird die AET ueber eine Formel kalkuliert.
    :param pet_raster: Ausgabe der Funktion "get_pet" (Typ: raster)
    :param water_raster: Maske der Gewaesserpixel (Typ: rasater)
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param rp_raster: Reduktionspunkt (rp) (Typ: raster)
    :param rpwp_dif_raster: Differenz zwischen dem Reduktionspunkt (rp) und dem Welkepunkt (wp) (Typ: raster)
    :param wp_raster: Welkepunkt (wp) (Typ: raster)
    :param bool_aet: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param date: Tages ID (Typ: integer)
    :return: AET des Tages (Typ: raster)
    """
    aet_raster = Con(water_raster == 1, pet_raster,
                     Con(s_pre_raster >= rp_raster, pet_raster,
                         Con(rpwp_dif_raster == 0, 0, ((s_pre_raster - wp_raster) / rpwp_dif_raster) * pet_raster)))

    if bool_aet:
        aet_raster.save("AET_{}".format(date))

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "AET ausgef�hrt.")
    return aet_raster


def get_precipitation(timeseries, p_stations_temp, date, idw_pow, rastercellsize, bool_p):
    """
    Funktion zur Interpolation des Niederschlags anhand des Inverse Distance Weighting (IDW). Die Berechnung erfolgt
    �ber das IDW-Tool der Erweiterung "Spatial Analyst". Zun�chst werden die Niederschlagswerte eines Tages je Station
    in eine Liste geschrieben. Der Index der Niederschlagswerte je Station entspricht dem Index des dazugeh�rigen
    Stationsnamens in einer zweiten Liste. Diese Listen werden mittels eines Searchcursors erstellt. �ber einen
    Updatecursor werden die Niederschlagsmesswerte in eine Tabelle geschrieben, die Informationen zum Stationsnamen und
    deren Koordinaten enthalten. Der Interpolationsradius ist auf das EZG Schotten 2 angepasst.
    :param timeseries: Niederschlagswerte je Station in mm (Typ: table)
    :param p_stations_temp: Kopie der Tabelle der Stationsnamen und -koordinaten in der Ergebnisdatenbank (Typ: table)
    :param date: Tages ID (Typ: integer)
    :param idw_pow: Exponent des IDW (Typ: integer)
    :param rastercellsize: Groesse der Rasterzellen (Typ: foat)
    :param bool_p: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :return: Interpolation des Niederschlags (Typ: raster)
    """
    station = []  # Liste der Stationsnamen, um den Index der N-Messerwerte der zweiten Liste den Stationen zuzuordnen
    p_sum = []  # Liste der N-Messwerte eines Tages

    # Speichert alle N-Messwerte mit den dazugeh�rigen Stationsnamen eines Tages in die Listen
    with arcpy.da.SearchCursor(timeseries, ["Stationsnummer", "Tagessumme_mm", "TagesID"],
                               "TagesID = {}".format(date)) as p_cursor:
        for information in p_cursor:
            station.append(information[0])
            p_sum.append(information[1])

    del p_cursor
    del information

    # Aktualisiert den N-Messwert je Tag in der tempor�ren Niederschlagsdatei
    for i in range(len(station)):
        p_incursor = arcpy.da.InsertCursor(p_stations_temp, ["Tagessumme_mm"])
        p_row = p_sum[i]
        p_incursor.updateRow(p_row)

    del p_incursor
    del p_row

    # Niederschlagsinterpolation mittels Inverse Distance Weighting (IDW)
    idw = Idw(p_stations_temp, "Tagessumme_mm", rastercellsize, idw_pow, RadiusFixed(20000.00000, 5), "")

    if bool_p:
        idw.save("Niederschlag_{}".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "P ausgef�hrt.")

    return idw


def get_runoff(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster, pet_raster, bool_r, date):
    """
    Funktion zur Berechnung des Gesamtabflusses. Der Gesamtabfluss ist die Summe des Abflusses von Landpixeln und von
    Gew�sserpixeln. Ersterer wird �ber die runoff_land-Funktion berechnet, w�hrend letzterer der Differenz aus dem
    Niederschlag und der PET entspricht.
    :param water_raster: Maske der Gewaesserpixel (Typ: rasater)
    :param lambda_wert: Lambda Wert
    :param wp_raster: Welkepunkt (wp) (Typ: raster)
    :param p_raster: Ausgabe der Funktion "get_precipitation" (Typ: raster)
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param fk_raster: Feldkapazitaet (fk) (Typ: raster)
    :param pet_raster: Ausgabe der Funktion "get_pet" (Typ: raster)
    :param bool_r: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param date: Tages ID (Typ: integer)
    :return: Gesamtabfluss je Rasterzelle (Typ: raster)
    """
    r_overflow = Con(p_raster + s_pre_raster > fk_raster, p_raster + s_pre_raster - fk_raster, 0)

    r_land = Con(water_raster == 0, lambda_wert * ((s_pre_raster - wp_raster) ** 2) + r_overflow, 0)

    r = r_land + water_raster * Con(p_raster > pet_raster, p_raster - pet_raster, 0)

    if bool_r:
        r.save("R_{}".format(date))

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "gesamtabfluss ausgef�hrt.")

    return r


def get_soilwater(s_pre_raster, p_raster, aet_raster, runoff_raster, bool_s, date):
    """
    Funktion zur Berechnung des aktuellen Bodenwasserspeichers. Die Berechnung erfolgt anhand der allgemeinen Wasser-
    haushaltsgleichung.
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param p_raster: Ausgabe der Funktion "get_precipitation" (Typ: raster)
    :param aet_raster: Ausgabe der Funktion "get_aet" (Typ: raster)
    :param runoff_raster: Ausgabe der Funktion "get_runoff" (Typ: raster)
    :param bool_s: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param date: Tages ID (Typ: integer)
    :return: Aktueller Bodenwasserspeicher (Typ: raster)
    """
    soilwater = s_pre_raster + p_raster - aet_raster - runoff_raster

    if bool_s:
        soilwater.save("S_{}".format(date))

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Bodenwasserspeicher ausgef�hrt.")

    return soilwater


def get_q_m3(runoff_raster, rastercellsize):
    """
    Berechnet den Gesamtabfluss des Einzugsgebietes in m^3.
    :param runoff_raster:  Ausgabe der Funktion "get_runoff" (Typ: raster)
    :param rastercellsize: Groesse der Rasterzellen (Typ: foat)
    :return: Gesamtabfluss des EZGs in Kubikmeter (Typ: float)
    """
    array = arcpy.RasterToNumPyArray(runoff_raster, nodata_to_value=0)
    r_sum = array.sum()
    r_m3 = (r_sum * 0.001 * rastercellsize ** 2)

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Q in m3 berechnet.")

    return r_m3


########################################################################################################################
#  Benutzereingaben in ArcMap
########################################################################################################################

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
outname = arcpy.GetParameterAsText(16)  # Datatype: string # Default: Ergebnistabelle_Gesamtabfluss_m3

#  Abfrage welche Daten gespeichert werden sollen
check_pet = bool(arcpy.GetParameterAsText(11))  # Datatype: boolean #  Default: false
check_aet = bool(arcpy.GetParameterAsText(12))  # Datatype: boolean #  Default: false
check_p = bool(arcpy.GetParameterAsText(13))  # Datatype: boolean #  Default: false
check_r = bool(arcpy.GetParameterAsText(14))  # Datatype: boolean #  Default: false
check_s = bool(arcpy.GetParameterAsText(15))  # Datatype: boolean #  Default: false

########################################################################################################################
#  Erstellen der Ergebnisdatenbank und Bestimmen des Arbeitsvereichnisses
########################################################################################################################

arcpy.env.overwriteOutput = True  # Vorhandene Ergebnisse k�nnen �berschrieben werden
arcpy.env.extent = basin
arcpy.env.mask = basin  # Erstellte Datens�tze werden nur f�r das EZG berechnet

if arcpy.Exists(folder + "\\" + name + ".gdb"):
    arcpy.Delete_management(folder + "\\" + name + ".gdb")
arcpy.CreateFileGDB_management(folder, name)

arcpy.env.workspace = folder + "\\" + name
# arcpy.env.scratchWorkspace = arcpy.env.workspace

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(folder))

########################################################################################################################
#  Zuweisung der Rasterdatens�tze
#  (Die Dateien m�ssen in der Basisdatenbank unter den vorgegebenen Namen gespeichert sein.
########################################################################################################################

# Dictionary zum Haudefaktor
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
cellsize = arcpy.Describe(s_init).children[0].meanCellHeight
lambda_parameter = (c / (Raster(r'{}\L_in_metern'.format(data)) * 1000) ** 2)

# Eine Kopie der Tabelle der Niederschlagsstationen, wird in der Ergebnisdatenbank erstellt
# Der Kopie wird ein neues Attributefeld f�r die t�glichen Niederschagsmesswerte je Station hinzugef�gt
p_temp = arcpy.Copy_management(r'{}\N_Messstationen'.format(data), "p_temp")
arcpy.AddField_management(p_temp, "Tagessumme_mm", "DOUBLE")

# Erstellen der Ergebnistabelle
arcpy.CreateTable_management(arcpy.env.workspace, outname)
arcpy.AddField_management(arcpy.env.workspace + "\\" + outname, "Datum", "TEXT")
arcpy.AddField_management(arcpy.env.workspace + "\\" + outname, "Q", "DOUBLE")
arcpy.AddField_management(arcpy.env.workspace + "\\" + outname, "Tages_ID", "LONG")

# Listen der Ergebnisse f�r die Ergebnistabelle
runoff_daily = []  # t�gliche Gesamtabfluesse des EZG in Kubikmeter
id_list = []  # Tages IDs
date_list = []  # Dati

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Berechnung der Rasterdatensaetze war erfolgreich.")

########################################################################################################################
#  Beginn des Hauptprogramms
#  Start der Modellierung
#  Iteration durch die Klimadaten des Untersuchungszeitraums
########################################################################################################################


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
        runoff = get_runoff(water, lambda_parameter, wp, precipitation, s_pre, fk, pet, check_r, id_day)
        runoff_m3 = get_q_m3(runoff, cellsize)
        s = get_soilwater(s_pre, precipitation, aet, runoff, check_s, id_day)

        s_pre = s  # �berschreiben des Bodenwasserspeichers des Vortages

        # Einf�gen des Ablusses in die Ergebnistabelle
        q_cursor = arcpy.da.InsertCursor(arcpy.env.workspace + "\\" + outname, ["Datum", "Q", "Tages_ID"])
        output_row = ["{0}.{1}.{2}".format(day, month, year), runoff_m3, id_day]
        q_cursor.insertRow(output_row)

        arcpy.AddMessage(time.strftime("%H:%M:%S: ") +
                         "Fertig mit der Berechnung des {0}.{1}.{2}".format(day, month, year))

# L�schen der Cursor und tempor�ren Objekte
del cursor
del q_cursor
arcpy.Delete_management("p_temp")

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Modellierung Abgeschlossen.")
