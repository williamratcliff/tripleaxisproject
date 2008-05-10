del *.o
g95 -c sup_exc.f90 -O -I..\LibG95
g95 -c simbo.f90   -O -I..\LibG95
g95 *.o -o simbo_g95 -O -LC:\CrysFML\LibG95 -lgcrysfml -Wl,--heap=0x01000000
