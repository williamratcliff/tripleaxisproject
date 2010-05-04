del *.o
del *.mod
g95 -c spg_info.f90  -O -I..\LibG95
g95 *.o -o spg_info  -O -L..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
