.\main.exe SingleGold FourGold Scan -1 2 0 0 -1 -1 | py .\plot2d.py SingleGold FourGold Scan -1 2 0 0 -1 -1 -1
.\main.exe SingleGold FourGold Scan -1 3.5 0 0 -1 -1 | py .\plot2d.py SingleGold FourGold Scan -1 3.5 0 0 -1 -1 -1
.\main.exe SingleGold FourGold Scan -1 4 0 0 -1 -1 | py .\plot2d.py SingleGold FourGold Scan -1 4 0 0 -1 -1 -1
.\main.exe SingleGold FourGold Scan -1 5 0 0 -1 -1 | py .\plot2d.py SingleGold FourGold Scan -1 5 0 0 -1 -1 -1

start run0.bat 3 2
start run1.bat 3 2
start run2.bat 3 2
start run3.bat 3 2
timeout 18000 /nobreak
start run4.bat 3 2
start run5.bat 3 2
start run6.bat 3 2