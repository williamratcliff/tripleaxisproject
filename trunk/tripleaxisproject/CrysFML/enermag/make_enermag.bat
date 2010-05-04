del *.obj
lf95 -c sup_exc.f90  -nstchk -nchk -o1 -mod ".;C:\crysFML\LibC"
lf95 -c enermag.f90  -nstchk -nchk -o1 -mod ".;C:\crysFML\LibC"
lf95 *.obj -out enermag -nomap -nstchk -nchk -o1 -mod ".;C:\crysFML\LibC" -lib C:\crysFML\LibC\crysFML
    upx enermag.exe
    copy enermag.exe d:\progs\.
