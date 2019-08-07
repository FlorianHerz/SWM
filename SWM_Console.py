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
import time


print(time.strftime("%H:%M:%S: ") + "Systemmodule geladen.")

#  Funktionen


def get_pet(haude_factor, temperature, humidity):
    """Berechnung eines Rasterdatensatzes der PET nach Haude.
    Es wird der Rasterdatensatz zum  Haude-Faktor mit der Temperatur und der relativen Luftfeuchtigkeit benötigt.
    """
    pet_raster = haude_factor * (6.1 * 10 ** ((7.5 * temperature) / (temperature + 237.2))) * (1.0 - humidity / 100.0)
    return pet_raster


def get_aet(pet_raster, water_raster, s_pre_raster, rp_raster, rpwp_dif_raster, wp_raster):
    """Berechnung der AET je Rasterzelle.
    Der Wert der AET entspricht bei Gewässerpixeln, sowie Zellen in denen der Bodenwassersgehalt oberhalb des RP liegt
    der PET. Die AET ist 0 wenn der RP gleich dem WP ist. Bei den anderen Pixeln wird die AET gemäß des SWMs berechnet.
    """
    aet_raster = Con(water_raster == 1, pet_raster,
                     Con(s_pre_raster >= rp_raster, pet_raster,
                         Con(rpwp_dif_raster == 0, 0, ((s_pre_raster - wp_raster) / rpwp_dif_raster) * pet_raster)))
    return aet_raster


def get_precipitation(timeseries, p_stations_temp, date, idw_pow, rastercellsize):
    """Niederschlagsinterpolation
    """
    station = []
    p_sum = []
    with arcpy.da.SearchCursor(timeseries, ["Stationsnummer", "Tagessumme_mm", "TagesID"],
                               "TagesID = {}".format(date)) as p_cursor:
        for information in p_cursor:
            station.append(information[0])
            p_sum.append(information[1])

    for i in range(len(station)):
        with arcpy.da.UpdateCursor(p_stations_temp, ["Stationsnummer", "Tagessumme_mm"]) as p_incursor:
            for p_row in p_incursor:
                if p_row[0] == station[i]:
                    p_row[1] = p_sum[i]

                p_incursor.updateRow(p_row)

    idw = Idw(p_stations_temp, "Tagessumme_mm", rastercellsize, idw_pow, RadiusFixed(20000.00000, 5), "")

    # Delete temporary objects
    del p_cursor
    del p_incursor

    return idw


def overflow(p_raster, s_pre_raster, fk_raster):
    """Berechnet den direkten Abfluss des Niederschlags bei gesättigtem Bodenwasser
    """
    r_overflow = Con(p_raster + s_pre_raster > fk_raster,
                     p_raster + s_pre - fk_raster, 0)
    return r_overflow


def runoff_land(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster):
    """Berechnet den Abfluss von Landpixeln
    """
    r_land = Con(water_raster == 0,
                 lambda_wert * ((s_pre - wp_raster) ** 2) + overflow(p_raster, s_pre_raster, fk_raster), 0)
    return r_land


def get_runoff(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster, pet_raster):
    """Berechnung des Oberflächenabflusses
    Niederschlag_id_day, s_pre, fk, water, lambda, wp, pet
    """
    r = runoff_land(water_raster, lambda_wert, wp_raster, p_raster, s_pre_raster, fk_raster) +\
        water_raster * (p_raster - pet_raster)
    return r


def get_soilwater(s_pre_raster, p_raster, aet_raster, runoff_raster):
    """Berechnung des aktuellen Bodenwasserspeichers
    """
    soilwater = s_pre_raster + p_raster - aet_raster - runoff_raster
    return soilwater


arcpy.CheckOutExtension("Spatial")

# Vorgabe der Benutzereingaben

data = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb'
basin = r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\EZG_Schotten2_Vektor'
basin_id = "Id"
name = r'SWM_Ergebnisdatenbank.gdb'
folder = r'C:\HiWi_Hydro-GIS'
start = 20040628
end = 20040704
s_init = Raster(r'C:\HiWi_Hydro-GIS\MTP_HydroGIS_Basisdaten.gdb\FK_von_L')
rp_factor = 0.85
c = 150
idw_exponent = "1"
check_pet = "true"
check_aet = "true"
check_p = "true"
check_r = "true"
check_s = "true"

#  Erstellen der Ergebnisdatenbank und Wahl des Arbeitsverzeichnisses

arcpy.env.overwriteOutput = True
arcpy.env.extent = s_init
arcpy.CreateFileGDB_management(folder, name)
arcpy.env.workspace = folder+r'\{}'.format(name)

print(time.strftime("%H:%M:%S: ") + "Die Ergebnisdatenbank wurde im Verzeichnis {} erstellt.".format(folder))

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
water = Raster(r'{}\Gewaesser'.format(data))
rpwp_dif = rp - wp
s_pre = s_init
p_data = r'{}\N_Zeitreihen'.format(data)
p_temp = arcpy.Copy_management(r'{}\N_Messstationen'.format(data), "p_temp")
arcpy.AddField_management(p_temp, "Tagessumme_mm", "DOUBLE")
cellsize = arcpy.Describe(s_init).children[0].meanCellHeight
lambda_parameter = (c / (Raster(r'{}\L_in_metern'.format(data)) * 1000)**2)
runoff_daily = []
id_list = []
date_list = []

print(time.strftime("%H:%M:%S: ") + "Berechnung der Rasterdatensaetze war erfolgreich.")

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

        # Berechnung der PET

        pet = get_pet(haude_dic[month], temp, humid)

        # Berechnung der AET

        aet = get_aet(pet, water, s_pre, rp, rpwp_dif, wp)

        #  Niederschlagsinterpolation

        precipitation = get_precipitation(p_data, p_temp, id_day, idw_exponent, cellsize)

        #  Ablfussberechnung

        runoff = get_runoff(water, lambda_parameter, wp, precipitation, s_pre, fk, pet)

        #  Berechnung des Bodenwasserspeichers

        s = get_soilwater(s_pre, precipitation, aet, runoff)

        #  Speichern der täglichen Rasterdateien
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

        s_pre = s
        ZonalStatisticsAsTable(basin, basin_id, runoff, "Q_day", "true", "SUM")
        with arcpy.da.SearchCursor("Q_day", "Sum") as r_cursor:
            for r_sum in r_cursor:
                runoff_daily.append(r_sum[0] * 0.001 * cellsize ** 2)
        id_list.append(id_day)
        date_list.append("{0}.{1}.{2}".format(day, month, year))

        print(time.strftime("%H:%M:%S: ") + "Fertig mit der Berechnung des {0}.{1}.{2}".format(day, month, year))

arcpy.CreateTable_management(arcpy.env.workspace, "Durchfluss")
arcpy.AddField_management("Durchfluss", "Datum", "STRING")
arcpy.AddField_management("Durchfluss", "Q", "DOUBLE")
arcpy.AddField_management("Durchfluss", "Tages_ID", "LONG")


q_cursor = arcpy.da.InsertCursor("Durchfluss", ["Datum", "Q", "Tages_ID"])

for i in range(len(id_list)):
    row = [date_list[i], runoff_daily[i], id_list[i]]
    q_cursor.insertRow(row)


# Delete cursor objects
del cursor
del q_cursor
del row
arcpy.Delete_management("Q_day")
arcpy.Delete_management("p_temp")

print(time.strftime("%H:%M:%S: ") + "Modellierung Abgeschlossen.")
