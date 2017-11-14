import math

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