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


def get_pet(haude_factor, temperature, humidity, date):
    """
    Funktion zur Berechnung der potentiellen Evapotranspiration (PET) nach Haude.
    :param haude_factor: Haudefaktor des entsprechenden Monats (Typ: raster)
    :param temperature: Tagestemperatur in Grad Celsius (Typ: float)
    :param humidity: relative Luftfeuchte (Typ: float)
    :param date: Tages ID (Typ: integer)
    :return: PET des Tages (Typ: raster)
    """
    pet_raster = haude_factor * (6.1 * 10 ** ((7.5 * temperature) / (temperature + 237.2))) * (1.0 - humidity / 100.0)
    pet_raster.save("PET_{}.tif".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "PET ausgeführt.")

    return pet_raster


def get_aet(pet_raster, water_raster, s_pre_raster, rp_raster, rpwp_dif_raster, wp_raster, date):
    """
    Funktion zur Berechnung der AET. Die AET wird entweder über Abfragen bestimmt, oder über eine Formel berechnet.
    Der Wert der AET entspricht bei Gewässerpixeln, sowie Zellen in denen der Bodenwassersgehalt oberhalb des
    Reduktionspunktes liegt, der PET. Die AET ist 0, wenn der Reduktionspunkt gleich dem Welkepunkt ist. Bei den anderen
    Pixeln wird die AET ueber eine Formel kalkuliert.
    :param pet_raster: Ausgabe der Funktion "get_pet" (Typ: raster)
    :param water_raster: Maske der Gewaesserpixel (Typ: rasater)
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param rp_raster: Reduktionspunkt (rp) (Typ: raster)
    :param rpwp_dif_raster: Differenz zwischen dem Reduktionspunkt (rp) und dem Welkepunkt (wp) (Typ: raster)
    :param wp_raster: Welkepunkt (wp) (Typ: raster)
    :param date: Tages ID (Typ: integer)
    :return: AET des Tages (Typ: raster)
    """
    aet_raster = Con(water_raster == 1, pet_raster,
                     Con(s_pre_raster >= rp_raster, pet_raster,
                         Con(rpwp_dif_raster == 0, 0, ((s_pre_raster - wp_raster) / rpwp_dif_raster) * pet_raster)))
    aet_raster.save("AET_{}.tif".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "AET ausgeführt.")

    return aet_raster


def get_precipitation(dataspace, date, idw_pow, rastercellsize):
    """
    Funktion zur Interpolation des Niederschlags anhand des Inverse Distance Weighting (IDW). Die Berechnung erfolgt
    über das IDW-Tool der Erweiterung "Spatial Analyst". Zunächst werden die Niederschlagswerte eines Tages je Station
    in eine Liste geschrieben. Der Index der Niederschlagswerte je Station entspricht dem Index des dazugehörigen
    Stationsnamens in einer zweiten Liste. Diese Listen werden mittels eines Searchcursors erstellt. Über einen
    Updatecursor werden die Niederschlagsmesswerte in eine Tabelle geschrieben, die Informationen zum Stationsnamen und
    deren Koordinaten enthalten. Der Interpolationsradius ist auf das EZG Schotten 2 angepasst.
    :param dataspace: Dateipfad der Basisdaten (Typ: string)
    :param date: Tages ID (Typ: integer)
    :param idw_pow: Exponent des IDW (Typ: integer)
    :param rastercellsize: Groesse der Rasterzellen (Typ: foat)
    :return: Interpolation des Niederschlags (Typ: raster)
    """
    arcpy.MakeQueryTable_management(r'{}\N_Messstationen;'.format(dataspace) +
                                    r'{}\N_Zeitreihen'.format(dataspace), "p_temp",
                                    "USE_KEY_FIELDS", "N_Messstationen.Stationsnummer;N_Zeitreihen.TagesID",
                                    "N_Messstationen.Stationsnummer;N_Messstationen.Stationsname; N_Messstationen.Shape\
                                    ;N_Zeitreihen.Tagessumme_mm;N_Zeitreihen.TagesID", "N_Zeitreihen.Stationsnummer =\
                                    N_Messstationen.Stationsnummer AND N_Zeitreihen.TagesID = {}".format(date))

    # Niederschlagsinterpolation mittels Inverse Distance Weighting (IDW)
    idw = Idw("p_temp", "N_Zeitreihen.Tagessumme_mm", rastercellsize, idw_pow, RadiusFixed(20000.00000, 5), "")
    idw.save("IDW_{}.tif".format(date))
    arcpy.Delete_management("p_temp")
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "P ausgeführt.")

    return idw


def overflow(p_raster, s_pre_raster, fk_raster):
    """
    Berechnet den direkten Abfluss des Niederschlags bei gesättigtem Bodenwasser.
    :param p_raster: Ausgabe der Funktion "get_precipitation" (Typ: raster)
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param fk_raster: Feldkapazitaet (fk) (Typ: raster)
    """
    r_overflow = Con(p_raster + s_pre_raster > fk_raster,
                     p_raster + s_pre - fk_raster, 0)

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Ueberlauf ausgeführt.")
    return r_overflow


def runoff_land(water_raster, lambda_param, wp_raster, p_raster, s_pre_raster, fk_raster):
    """
    Berechnet den Abfluss von Landpixeln.
    :param water_raster: Maske der Gewaesserpixel (Typ: rasater)
    :param lambda_param: Lambda Wert
    :param wp_raster: Welkepunkt (wp) (Typ: raster)
    :param p_raster: Ausgabe der Funktion "get_precipitation" (Typ: raster)
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param fk_raster: Feldkapazitaet (fk) (Typ: raster)
    """
    r_land = Con(water_raster == 0,
                 lambda_param * ((s_pre - wp_raster) ** 2) + overflow(p_raster, s_pre_raster, fk_raster), 0)

    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Landabfluss ausgeführt.")
    return r_land


def get_runoff(water_raster, lambda_param, wp_raster, p_raster, s_pre_raster, fk_raster, pet_raster, date):
    """
    Funktion zur Berechnung des Gesamtabflusses. Der Gesamtabfluss ist die Summe des Abflusses von Landpixeln und von
    Gewässerpixeln. Ersterer wird über die runoff_land-Funktion berechnet, während letzterer der Differenz aus dem
    Niederschlag und der PET entspricht.
    :param water_raster: Maske der Gewaesserpixel (Typ: rasater)
    :param lambda_param: Lambda Wert
    :param wp_raster: Welkepunkt (wp) (Typ: raster)
    :param p_raster: Ausgabe der Funktion "get_precipitation" (Typ: raster)
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param fk_raster: Feldkapazitaet (fk) (Typ: raster)
    :param pet_raster: Ausgabe der Funktion "get_pet" (Typ: raster)
    :param date: Tages ID (Typ: integer)
    :return: Gesamtabfluss je Rasterzelle (Typ: raster)
    """
    r = runoff_land(water_raster, lambda_param, wp_raster, p_raster, s_pre_raster, fk_raster) +\
        water_raster * Con(p_raster > pet_raster, p_raster - pet_raster, 0)
    r.save("R_{}.tif".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Gesamtabfluss ausgeführt.")

    return r


def get_soilwater(s_pre_raster, p_raster, aet_raster, runoff_raster, date):
    """
    Funktion zur Berechnung des aktuellen Bodenwasserspeichers. Die Berechnung erfolgt anhand der allgemeinen Wasser-
    haushaltsgleichung.
    :param s_pre_raster: Bodenwasserspeicher des Vortages (Typ: raster)
    :param p_raster: Ausgabe der Funktion "get_precipitation" (Typ: raster)
    :param aet_raster: Ausgabe der Funktion "get_aet" (Typ: raster)
    :param runoff_raster: Ausgabe der Funktion "get_runoff" (Typ: raster)
    :param date: Tages ID (Typ: integer)
    :return: Aktueller Bodenwasserspeicher (Typ: raster)
    """
    soilwater = s_pre_raster + p_raster - aet_raster - runoff_raster
    soilwater.save("S_{}.tif".format(date))
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Bodenwasserspeicher ausgeführt.")

    return soilwater


def get_q_m3(runoff_raster, rastercellsize):
    """
    Berechnet den Gesamtabfluss des Einzugsgebietes in m^3.
    :param runoff_raster: Ausgabe der Funktion "get_runoff" (Typ: raster)
    :param rastercellsize: Größe der Rasterzellen
    :return: täglicher Gesamtabfluss in m^3
    """
    array = arcpy.RasterToNumPyArray(runoff_raster, nodata_to_value=0)
    r_sum = array.sum()
    r_m3 = (r_sum * 0.001 * rastercellsize ** 2)
    arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Q in m3 berechnet.")

    return r_m3


def delete_raster(yesterday, bool_pet, bool_aet, bool_p, bool_r, bool_s):
    """
    Löscht wenn gewünscht die erstellten Rasterdateien. Es können nur die Dateien des Vortages gelöscht werden. Die des
    ersten und letzten Tages müssen nach dem Cursor des Hauptprogramms separat gelöscht werden.
    :param yesterday: Tages ID des Vortages (Typ: integer)
    :param bool_pet: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param bool_aet: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param bool_p: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param bool_r: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :param bool_s: Angabe ob der Ausgabedatensatz gespeichert werden soll (Typ: boolean)
    :return:
    """
    if not bool_pet:
        arcpy.Delete_management("PET_{}.tif".format(yesterday))
    if not bool_aet:
        arcpy.Delete_management("AET_{}.tif".format(yesterday))
    if not bool_p:
        arcpy.Delete_management("IDW_{}.tif".format(yesterday))
    if not bool_r:
        arcpy.Delete_management("R_{}.tif".format(yesterday))
    if not bool_s:
        arcpy.Delete_management("S_{}.tif".format(yesterday))


def write_to_table(resultspace, tablename, result, date):
    """
    Schreibt einen Ergebniswert und das zugehörige Datum in die Ergebnistabelle.
    :param resultspace: Dateipfad des Workspace (Typ: sting)
    :param tablename: Name der Ergebnistabelle (Typ: string)
    :param result: Gesamtabfluss des EZGs in m^3 (Typ: float)
    :param date: Tages ID (Typ: string)
    :return:
    """
    q_cursor = arcpy.da.InsertCursor(r'{0}\Ergebnistabellen.gdb\{1}'.format(resultspace, tablename),
                                             ["Datum", "Q"])
    output_row = ["{0}.{1}.{2}".format(date[-2:], date[-4:-2], date[:3]), result]
    q_cursor.insertRow(output_row)
    del q_cursor


########################################################################################################################
#  Benutzereingaben in ArcMap
########################################################################################################################

data = arcpy.GetParameterAsText(0)  # Datatype: workspace #  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb
basin = arcpy.GetParameterAsText(1)  # Datatype: Featurelayer
#  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\EZG_Schotten2_Vektor
basin_id = arcpy.GetParameterAsText(2)  # Field #  Default: Id
name = arcpy.GetParameterAsText(3)  # Datatype: string #  Default: SWM_Schotten2_2003-2004
folder = arcpy.GetParameterAsText(4)  # Datatype: folder #  Default: C:\HiWi_Hydro-GIS
id_yesterday = start = int(arcpy.GetParameterAsText(5))  # Datatype: long #  Default: 20030101
end = int(arcpy.GetParameterAsText(6))  # Datatype: long #  Default: 20041231
s_init = ExtractByMask(Raster(arcpy.GetParameterAsText(7)), basin)  # Datatype: geodataset
#  Default: C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L
rp_factor = float(arcpy.GetParameterAsText(8))  # Datatype: string #  Default: 0.85
c = int(arcpy.GetParameterAsText(9))  # Datatype: string #  Default: 150
idw_exponent = arcpy.GetParameterAsText(10)  # Datatype: string #  Default: 1
#  Abfrage welche Daten gespeichert werden sollen
check_pet = arcpy.GetParameterAsText(11)  # Datatype: boolean #  Default: true
if check_pet == 'false':
    check_pet = False
else:
    check_pet = True
check_aet = arcpy.GetParameterAsText(12)  # Datatype: boolean #  Default: true
if check_aet == 'false':
    check_aet = False
else:
    check_aet = True
check_p = arcpy.GetParameterAsText(13)  # Datatype: boolean #  Default: true
if check_p == 'false':
    check_p = False
else:
    check_p = True
check_r = arcpy.GetParameterAsText(14)  # Datatype: boolean #  Default: true
if check_r == 'false':
    check_r = False
else:
    check_r = True
check_s = arcpy.GetParameterAsText(15)  # Datatype: boolean #  Default: true
if check_s == 'false':
    check_s = False
else:
    check_s = True
outname = arcpy.GetParameterAsText(16) + "_{}bis{}".format(start, end)  # Datatype: string # Default: Ergebnistabelle
########################################################################################################################
#  Erstellen der Ergebnisdatenbank und Bestimmen des Arbeitsvereichnisses
########################################################################################################################

arcpy.env.overwriteOutput = True  # Vorhandene Ergebnisse können überschrieben werden
arcpy.env.extent = s_init
workspace = arcpy.env.workspace = r'{}\{}'.format(folder, name)

if arcpy.Exists(workspace):
    arcpy.Delete_management(workspace)
arcpy.CreateFolder_management(folder, name)

arcpy.CreateFileGDB_management(workspace, "Ergebnistabellen.gdb")

arcpy.CreateFolder_management(folder, "Scratch")
arcpy.env.scratchWorkspace = r'{}\Scratch'.format(folder)

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(folder))

########################################################################################################################
#  Zuweisung der Rasterdatensätze
#  (Die Dateien müssen in der Basisdatenbank unter den vorgegebenen Namen gespeichert sein.
########################################################################################################################

# Dictionary zum Haudefaktor
haude_dic = {1: ExtractByMask(Raster(r'{}\Haude_1'.format(data)), basin),
             2: ExtractByMask(Raster(r'{}\Haude_2'.format(data)), basin),
             3: ExtractByMask(Raster(r'{}\Haude_3'.format(data)), basin),
             4: ExtractByMask(Raster(r'{}\Haude_4'.format(data)), basin),
             5: ExtractByMask(Raster(r'{}\Haude_5'.format(data)), basin),
             6: ExtractByMask(Raster(r'{}\Haude_6'.format(data)), basin),
             7: ExtractByMask(Raster(r'{}\Haude_7'.format(data)), basin),
             8: ExtractByMask(Raster(r'{}\Haude_8'.format(data)), basin),
             9: ExtractByMask(Raster(r'{}\Haude_9'.format(data)), basin),
             10: ExtractByMask(Raster(r'{}\Haude_10'.format(data)), basin),
             11: ExtractByMask(Raster(r'{}\Haude_11'.format(data)), basin),
             12: ExtractByMask(Raster(r'{}\Haude_12'.format(data)), basin)}

climatedata = r'{}\TempFeuchte'.format(data)
fk = ExtractByMask(Raster(r'{}\fk_von_L'.format(data)), basin)
rp = fk * rp_factor
wp = ExtractByMask(Raster(r'{}\wp'.format(data)), basin)
water = ExtractByMask(Raster(r'{}\Gewaesser'.format(data)), basin)
rpwp_dif = rp - wp
s_pre = s_init
p_data = r'{}\N_Zeitreihen'.format(data)
cellsize = s_init.meanCellHeight
lambda_parameter = (c / (ExtractByMask(Raster(r'{}\L_in_metern'.format(data)), basin) * 1000) ** 2)

# Erstellen der Ergebnistabelle
result_path = arcpy.CreateTable_management(r'{}\Ergebnistabellen.gdb'.format(workspace), outname)
arcpy.AddField_management(result_path, "Datum", "TEXT")
arcpy.AddField_management(result_path, "Q", "DOUBLE")

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
        pet = get_pet(haude_dic[month], temp, humid, id_day)
        aet = get_aet(pet, water, s_pre, rp, rpwp_dif, wp, id_day)
        precipitation = get_precipitation(data, id_day, idw_exponent, cellsize)
        runoff = get_runoff(water, lambda_parameter, wp, precipitation, s_pre, fk, pet, id_day)
        runoff_m3 = get_q_m3(runoff, cellsize)
        s = get_soilwater(s_pre, precipitation, aet, runoff, id_day)

        s_pre = s  # Überschreiben des Bodenwasserspeichers des Vortages
        if not id_yesterday == start:
            delete_raster(id_yesterday, check_pet, check_aet, check_p, check_r, check_s)
        id_yesterday = id_day
        write_to_table(workspace, outname, runoff_m3, str(id_day))

        arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Fertig mit der Berechnung des {0}.{1}.{2}".format(day, month, year))

del cursor
for i in [start, end]:
    delete_raster(i, check_pet, check_aet, check_p, check_r, check_s)
arcpy.Delete_management(r'{}\Scratch'.format(folder))

arcpy.AddMessage(time.strftime("%H:%M:%S: ") + "Modellierung Abgeschlossen.")
