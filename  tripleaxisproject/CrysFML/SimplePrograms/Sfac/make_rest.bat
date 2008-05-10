del *.o
g95 -c observ.f90          -O -I..\..\LibG95
g95 -c cost_functions.f90  -O -I..\..\LibG95
g95 -c Opt_restraints.f90  -O -I..\..\LibG95
g95  *.o -o Opt_restraints -L..\..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
