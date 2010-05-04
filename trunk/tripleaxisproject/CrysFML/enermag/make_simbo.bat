del *.obj
lf95 -c sup_exc.f90 -stchk -chk -o0 -mod ".;c:\crysFML\LibC"
lf95 -c simbo.f90   -stchk -chk -o0 -mod ".;c:\crysFML\LibC"
lf95 *.obj -out simbo -nomap -stchk -chk -o0 -mod ".;c:\crysFML\LibC" -lib c:\crysFML\LibC\crysFML
    upx simbo.exe
    copy simbo.exe d:\progs\.
