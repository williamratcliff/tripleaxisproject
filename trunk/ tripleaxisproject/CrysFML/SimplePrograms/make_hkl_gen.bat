del *.o
del *.mod
g95 -c hkl_gen.f90  -O -I..\LibG95
g95 *.o -o hkl_gen  -O -L..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
