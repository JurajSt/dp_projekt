import math
from input_data import *
from azi_ele_coor import *
import ogr
import sys

try:
    import psycopg2
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"
try:
    import xlsxwriter
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    import matplotlib.pyplot as plt
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"


connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

def fFresnelZ(azi, ele):

    cursor.execute("DROP TABLE IF EXISTS FresnelZone")
    cursor.execute("CREATE TABLE FresnelZone (stname CHAR (4), azimuth REAL, elevangle REAL, geom GEOMETRY)")

    c = 0
    while c < len(stanice):
        pi = math.pi
        h = 1.6
        if len(azi) == 0:
            azi = stanice[c][2]
        if len(ele) == 0:
            ele = np.linspace(stanice[c][1][0], stanice[c][1][1], 3)
        n=1
        l = 0.24
        d = n*l/2
        elipsa_list = []


        row = 1
        col = 0
        workbook = xlsxwriter.Workbook(cesta_vystup+'FresnelZone'+ stanice[c][4]+'.'+ulozenie)
        worksheet = workbook.add_worksheet()
        worksheet.write(0, 0, 'C')
        worksheet.write(0, 1, 'x')
        worksheet.write(0, 2, 'y')

        for j in range(len(azi)):

            azi_rad = math.radians(azi[j])

            for i in range(len(ele)):

                ele_rad = math.radians(int(ele[i]))
                line = ogr.Geometry(ogr.wkbLineString)

                R = h/math.tan(ele_rad) + (d/math.sin(ele_rad))/math.tan(ele_rad)
                b = math.sqrt(((2*d*h)/math.sin(ele_rad)) + (math.pow(d/math.sin(ele_rad),2)))
                a = b/math.sin(ele_rad)
                k = 0
                st = math.radians(1)      # stupen na radiany

                x_elipsa = []
                y_elipsa = []
                while k <= 2*pi:
                    xi = a * math.cos(k) + R
                    yi = b * math.sin(k)
                    xe = math.sin(azi_rad)*xi - math.cos(azi_rad)*yi
                    ye = math.sin(azi_rad)*yi + math.cos(azi_rad)*xi
                    k = k+st
                    Lat, Lon, hel = fXYZ_to_LatLonH(telg[0][0], telg[0][1], telg[0][2])
                    WGSmercatorX, WGSmercatorY = fwgs4326towgs3857(Lat, Lon)


                    worksheet.write(row, col, ele[i])
                    worksheet.write(row, col, ele[i])
                    worksheet.write(row, col + 1, WGSmercatorX + xe)
                    worksheet.write(row, col + 2, WGSmercatorY + ye)
                    row = row + 1

                    line.AddPoint(round(WGSmercatorX + xe,3), round(WGSmercatorY + ye,3))

                    x_elipsa.append(stanice[c][0][0] + xe)
                    y_elipsa.append(stanice[c][0][1] + ye)
                    list = [x_elipsa,y_elipsa]
                elipsa_list.append(list)
                wkt_line = line.ExportToWkt()
                #print wkt_line

                cursor.execute("INSERT INTO FresnelZone (stname, azimuth, elevangle, geom)"
                               " VALUES ('%s',  %s, %s, ST_GeometryFromText('%s'))" %(stanice[c][4], azi[j], ele[i], wkt_line))

                #print wkt_line

        plt.plot(elipsa_list[0][0], elipsa_list[0][1], 'r', label='ele uhol'+str(ele[0]))
        plt.plot(elipsa_list[1][0], elipsa_list[1][1], 'g', label='ele uhol'+str(ele[0]))
        plt.plot(elipsa_list[2][0], elipsa_list[2][1], 'b', label='ele uhol'+str(ele[0]))
        plt.show()

        workbook.close()
        c = c+1
