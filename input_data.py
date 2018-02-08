# pripojenie k db
hostname = 'localhost'
username = 'postgres'
password = 'sql'
database = 'dp'
DB = "host=" + hostname + " user=" + username + " password=" + password + " dbname=" + database

rinexObspath = 'data/bordel/test/obs/'                  #data/ganp3060.16o'              # demo.10o'
rinexNavpath = 'data/bordel/test/navig/'# demo.10n
csvpath = 'data/SHMU/LIE1_TELG_KAME.csv'

# station
ganp = [3929181.851, 1455236.510, 4793653.699]
kame = [3892532.358, 1572220.333, 4785952.565]
trf2 = [4119400.427, 1170248.492, 4712323.807]
sgb2 = [4180931.236, 973735.204, 4703203.291]
vae6 = [3249402.748, 692761.948, 5426399.896]
suld = [3446394.505, 591712.939, 5316383.267]
hofn = [2679690.298, -727951.336, 5722789.244]
nya1 = [1202434.130, 252632.221, 6237772.435]

Lies = []
Telg = []
Kame = []

platnos_spravy = 7200           #7200 sekund = 2 hodiny
interval = 30                   # interval vypoctu polohy druzice v sekundach