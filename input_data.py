
# pripojenie k db
hostname = 'localhost'
username = 'postgres'
password = 'sql'
database = 'dp'
DB = "host=" + hostname + " user=" + username + " password=" + password + " dbname=" + database

# cesty k suborom
nav_path = 'D:/diplomka/dp_projekt/data/test/navig/'#'D:/diplomka/dp_projekt/data/bordel/test/navig/'
obs_path = 'D:/diplomka/dp_projekt/data/test/obs/'  # 'D:/diplomka/telg_2/'
eph_path = 'D:/diplomka/dp_projekt/data/test/eph/'  # igr18644.sp3'
csvpath = 'data/SHMU/HOEFN_15_16.csv'


# datum zaciatku a konca vypoctu
rok_start, mes_start, den_start = 2015, 9, 1
rok_konec, mes_konec, den_konec = 2017, 3, 10

# udaje o referencnej stanici
#nazovStanice = [[x,y,z], [elev. uhol od-do], [azimut od-do], referencna_vyska_stanice, nazov_stanice_rinex, nazov_meteostanice]
# pozn. suradnice vo WGS84, azimuta elev. uhol v stupnoch. Ak neznami azimut tak prazdny seznam,
# azimut zadavany podla wgs84 web mercator!!!
kame = [[3892532.358, 1572220.333, 4785952.565], [5,25], [50, 80], 1.6, 'kame', 'Kamenica nad Cirochou']
telg = [[3947396.223, 1451396.020, 4780197.834], [5,30], [190, 220], 1.6, 'telg', 'Telgart']
hofn = [[2679690.298, -727951.336, 5722789.244], [5,25], [110, 150], 4.45, 'hofn', 'Hoefn']
ganp = [[3929181.9, 1455236.5, 4793653.8]]

stanice = [hofn]
# odhadovana vyska anteny xx, yy, zz>> min a max vyska a krok v metroch
vyska_min = 1
vyska_max = 2
krok = 0.01   # krok po 0.01  tj  1 cm

#['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G10' 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18'
# , 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30', 'G31', 'G32']
satellite = ['G14'] # zvol cislo druzice

# cisla druzic pre ktore bol v testovacich observacnych datach prijimany signal na L5, vychadza to z http://gpsworld.com/the-almanac/druzice block IIF
satellites_S5 = ['G01', 'G03','G06', 'G08', 'G09', 'G10', 'G24','G25', 'G26','G27', 'G30', 'G32']
# signaly do vypoctu brane automaticky SNR1 a SNR2, SNR 5 ak je k dispozicii

pocet_observacii = 70  # minimalny pocet observacii pre analyzu

platnos_spravy = 7200           #7200 sekund = 2 hodiny

# data potrebne k ulozeniu a analyze
analyza1 = 'a'   # a = urobit analyzu pre vyber druzice, n = bez analyzy. Pri jednosekundovych datach bude cas vypoctu pre analyzu dlhsi.
ulozenie = 'csv'   # ulozenie db = do databaze, csv (xlsx) = do csv, xlsx musi byt zadana aj cesta!!!
cesta_vystup = 'D:/diplomka/vystupcsv/'
sursystem = 'xy' # xyz = pravouhle suradnice WGS84 , llh = zemepisne suradnice WGS84(epsg: 4326), xy = WGS84 web mercator (epsg: 3857)
graf = 'a'   # a = vysledne grafi sa ulozia , n = neulozia

# parameter pre fresnelove zony , bude pouzity rovnaky pre vsetky stanice, ak nezadane nic parametre sa preberu  z dat o stanici vyssie
fresnel_azi = []    #[100, 150]
fresnel_ele = []    #[7,15,30]