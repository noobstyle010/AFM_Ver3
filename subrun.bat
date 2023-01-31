for /l %%j in (1,1,25) do (
    for /l %%i in (1,1,100) do (
        rem .\main.exe probe surface scan R H x y l a
        .\main.exe FunctionalizedTip FourGold Scan %1 %2 0 0 %3 %%j | py .\plot2d.py FunctionalizedTip FourGold Scan %1 %2 0 0 %3 %%j %%i
    )
)