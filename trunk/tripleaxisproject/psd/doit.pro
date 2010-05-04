pro doit

;triangulate,xran,yran,triangles
;zgrid=trigrid(xran,yran,zran,triangles,$
;xgrid=xgrid, ygrid=ygrid)
;surface, zgrid, xgrid, ygrid

;file=dialog_pickfile(filter='*.dat',path='c:\sqltest',file='temp.dat')
;FMT = 'F,F,F,F'
;READCOL,'c:\sqltest\temp.dat',F=FMT,h,k,dE,intensity
RDFLOAT,'c:\sqltest\temp.dat',h,k,dE,intensity,skipline=1

;hsize=size(h,/n_elements)
;Earr=floatarr(hsize,hsize)
;Iarr=floatarr(hsize,hsize)
Earr=diag_matrix(dE)
Iarr=diag_matrix(intensity)
s = Size(Iarr, /Dimensions)
Vert_Colors=Reform(BytScl(Iarr), s[0]*s[1])
fsc_surface,Earr,h,k,colortable=13,/elevation_shading,data2=Iarr ; works ok
thispalette = obj_new('idlgrpalette')
thispalette->loadct, 13
;isurface, Earr,h,k,vert_colors=vert_colors,style=0,palette=thispalette,background_color=[0,0,0]
;loadcolors
;bottom=16
;loadct, 13,bottom=bottom
;;loadct, 13
;triangulate, h,k,triangles
;N_Colors=10
;dEgrid=trigrid(h,k,DE,triangles, xgrid=hgrid, ygrid=kgrid)
;Igrid=trigrid(h,k,intensity,triangles, xgrid=xgrid, ygrid=ygrid)
;shades=bytscl(Igrid,top=!d.table_size-1-bottom)+byte(bottom)
;shade_surf, dEgrid,hgrid,kgrid, shades=shades;shades=BYTSCL(Igrid,top=!d.table_size)
end