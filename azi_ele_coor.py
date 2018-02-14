#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      kac072
#
# Created:     12.07.2017
# Copyright:   (c) kac072 2017
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import math
import numpy

def fwgs4326towgs3857(lat, lon):
    x = lon * 20037508.34 / 180;
    y = math.log(math.tan((90 + lat) * math.pi / 360)) / (math.pi / 180);
    y = y * 20037508.34 / 180;

    return x, y

def fLatLonH_to_XYZ(B,L,H):

    a=6378137.0
    b=6356752.3142
    e_2=((a*a)-(b*b))/(a*a)
    B = math.radians(float(B))
    L = math.radians(float(L))

    N = a/math.sqrt(1-e_2*math.pow(math.sin(B),2))

    X = (N+H)*math.cos(B)*math.cos(L)
    Y = (N+H)*math.cos(B)*math.sin(L)
    Z = (N*(1-e_2)+H)*math.sin(B)

    return X, Y, Z

def fXYZ_to_LatLonH(X,Y,Z):
    X=float(X)
    Y=float(Y)
    Z=float(Z)

    a=6378137.0
    b=6356752.3142
    e_2=((a*a)-(b*b))/(a*a)
    e=math.sqrt(e_2)

    lon = math.atan(Y/X)
    lon_degrees = math.degrees(lon)
    p = math.sqrt(X*X+Y*Y)

    lat0=math.atan((Z/p)*(1/(1-e_2)))
    cos2lat0 = (1+math.cos(2*lat0))/2
    sin2lat0 = (1-math.cos(2*lat0))/2
    N0 = (a*a)/math.sqrt(a*a*cos2lat0 + b*b*sin2lat0)
    h0=(p/math.cos(lat0))-N0
    lat1=math.atan((Z/p)*(1/(1-(e_2)*(N0/(N0+h0)))))

    cos2lat1 = (1+math.cos(2*lat1))/2
    sin2lat1 = (1-math.cos(2*lat1))/2
    N1 = (a*a)/math.sqrt(a*a*cos2lat1 + b*b*sin2lat1)
    h1=(p/math.cos(lat1))-N1
    lat2=math.atan((Z/p)*(1/(1-(e_2)*(N1/(N1+h1)))))
    lat2_degrees = math.degrees(lat2)

    cos2lat2 = (1+math.cos(2*lat2))/2
    sin2lat2 = (1-math.cos(2*lat2))/2
    N2 = (a*a)/math.sqrt(a*a*cos2lat2 + b*b*sin2lat2)
    h2=(p/math.cos(lat2))-N2
    lat3=math.atan((Z/p)*(1/(1-(e_2)*(N2/(N2+h2)))))
    lat3_degrees = math.degrees(lat3)

    return lat3_degrees, lon_degrees, h2

def fComputeAziEle(receiver, satellite): #compute azimuth and elevation angles for given coordinats of receiver and satellite,
    #coordinates are given in XYZ format as a Python list = [X,Y,Z]
    #https://www.mathworks.com/matlabcentral/fileexchange/41364-gps-navigation-toolbox/content/GNT08.1.2/final-GPS/Calc_Azimuth_Elevation.m
    #http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
    vector =[satellite[0]-receiver[0],satellite[1]-receiver[1],satellite[2]-receiver[2]]
    vector_size = math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])
    unit_vector = numpy.matrix([vector[0]/vector_size, vector[1]/vector_size, vector[2]/vector_size])
    unit_vector_T = unit_vector.T

    rec_latlonH = fXYZ_to_LatLonH(receiver[0],receiver[1],receiver[2])
    rec_lat = math.radians(rec_latlonH[0])
    rec_lon = math.radians(rec_latlonH[1])
    ENU = numpy.matrix([[-1*math.sin(rec_lon), math.cos(rec_lon), 0], [-1*math.sin(rec_lat)*math.cos(rec_lon), -1*math.sin(rec_lat)*math.sin(rec_lon), math.cos(rec_lat)],
    [math.cos(rec_lat)*math.cos(rec_lon), math.cos(rec_lat)*math.sin(rec_lon), math.sin(rec_lat)]])

    elevation = math.degrees(math.asin(float(ENU[2]*unit_vector_T)))
    azimuth_first = math.degrees(math.atan(float(ENU[0]*unit_vector_T)/float(ENU[1]*unit_vector_T)))
    #azimuth correction to place it in a correct quadrant
    if (float(ENU[0]*unit_vector_T) > 0 and float(ENU[1]*unit_vector_T) < 0) or (float(ENU[0]*unit_vector_T) < 0 and float(ENU[1]*unit_vector_T) < 0):
        azi_correction = 180
    elif float(ENU[0]*unit_vector_T) < 0 and float(ENU[1]*unit_vector_T) > 0:
        azi_correction = 360
    else:
        azi_correction = 0

    azimuth = azimuth_first + azi_correction
    return [azimuth, elevation]


def fCalculateAzimuth(receiver, satellite):

    #rec_latlonH = fXYZ_to_LatLonH(receiver[0], receiver[1], receiver[2])
    #sat_latlonH = fXYZ_to_LatLonH(satellite[0], satellite[1], satellite[2])
    dX = satellite[0] - receiver[0]
    dY = satellite[1] - receiver[1]
    PI = math.pi
    Azimuth = 0  # Default case, dX = 0 and dY >= 0
    if dX > 0:
        Azimuth = 90 - math.atan(dY / dX) * 180 / PI
    elif dX < 0:
        Azimuth = 270 - math.atan(dY / dX) * 180 / PI
    elif dY < 0:
        Azimuth = 180

    return Azimuth

