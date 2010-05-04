del *.o
del *.mod
g95 -c check_group.f90  -O -I..\LibG95
g95 *.o -o check_group  -O -L..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
