del *.o
del *.mod
g95 -c subgroups.f90  -O -I..\LibG95
g95 *.o -o subgroups  -O -L..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
