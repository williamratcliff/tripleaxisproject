del *.o
g95 -c observ.f90          -O -I..\..\LibG95
g95 -c cost_functions.f90  -O -I..\..\LibG95
g95 -c Optim_Sfac.f90      -O -I..\..\LibG95
g95  *.o -o Optim_Sfac     -L..\..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
