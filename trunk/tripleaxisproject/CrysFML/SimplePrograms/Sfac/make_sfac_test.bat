del *.o
g95 -c sfac_test.f90  -O -I..\..\LibG95
g95 *.o -o sfac_test  -O -L..\..\LibG95 -lgcrysfml  -Wl,--heap=0x01000000
