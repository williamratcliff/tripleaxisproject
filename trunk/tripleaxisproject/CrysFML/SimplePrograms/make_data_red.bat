del *.o
del *.mod
g95 -c data_red.f95  -O -I..\LibG95
g95 *.o -o data_red  -O -L..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
