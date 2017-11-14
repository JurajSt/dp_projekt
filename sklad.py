import psycopg2
from geographiclib import geodesic
from input_data import *
import ogr
import azi_ele_coor as aec

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()
# vytvorenie tabulky pre body azimutu v epsg:3857 wgs84 web mercator
cursor.execute('DROP TABLE IF EXISTS kruh')
cursor.execute('CREATE TABLE  kruh (id INT NOT NULL, azimuth REAL, geom GEOMETRY, PRIMARY KEY (id))')
number_of_vertices = 36

radius = 100
vertices = []
latitude, longitude, h  = aec.fXYZ_to_LatLonH(3929181.851, 1455236.510, 4793653.699)

for i in range(number_of_vertices):
    degree = 360.0/number_of_vertices*i
    vertex = geodesic.Geodesic.WGS84.Direct(latitude, longitude, degree, radius)
    vertices.append((vertex['lat2'], vertex['lon2']))
    coor = aec.fLatLonH_to_XYZ(vertex['lat2'], vertex['lon2'], h)
    point = ogr.Geometry(ogr.wkbPoint)
    angle = aec.fComputeAziEle([3929181.851, 1455236.510, 4793653.699], [coor[0], coor[1], coor[2]])
    point.AddPoint(vertex['lon2'], vertex['lat2'])
    wkt = point.ExportToWkt()
    cursor.execute('INSERT INTO kruh (id, azimuth, geom) VALUES (%s, %s, ST_GeometryFromText(%s))',
        (i,angle[0], wkt))
    vertices.append(vertices[0])
#print vertices