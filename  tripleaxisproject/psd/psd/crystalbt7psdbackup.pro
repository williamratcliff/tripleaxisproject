;Realse Note
;Author: Jiying Li, version 1.0
;Function:  For the BT7 a-axis powder mode data reduction. It can also be used to view the raw data
;collected under any modes using the psd
;

;***************************************************************
;The following function read the data for every psd channel
;Record the counting unit (time) or monitor
;It reads the independent variable and save it together with psd data
pro load_psd_data, event
	widget_control,event.top,get_uvalue = pstate
	header = ['QX', 'QY', 'QZ', 'E', 'A4', 'Time', 'Temp', 'Monitor', $
			'Detector', ',DFMRot','A2', 'A3', 'A5', 'A6', 'DFM', $
			'A1', 'AnalyzerBlade1', 'AnalyzerBlade2', 'AnalyzerBlade3', 'AnalyzerBlade4', $
			 'AnalyzerBlade5', 'AnalyzerBlade6', 'AnalyzerBlade7', 'AnalyzerBlade8', 'AnalyzerBlade9', $
			 'AnalyzerBlade10', 'AnalyzerBlade11','AnalyzerBlade12', 'AnalyzerBlade13', 'AnalyzerRotation', $
			 'ApertHori', 'ApertVert', 'BkSltHght', 'BkSltWdth', 'DDC0', 'DDC1', 'DDC2', 'DiffDet', 'FLIP', 'FilRot', $
			'FilTilt', 'FilTran', 'FocusCu', 'FocusPG', 'H', 'HKL', 'K', 'L','Monitor2',$
			'MonoBlade01', 'MonoBlade02', 'MonoBlade03', 'MonoBlade04', 'MonoBlade05', $
			'MonoBlade06','MonoBlade07','MonoBlade08', 'MonoBlade09', 'MonoBlade10','MonoElev','MonoTrans',$
			'PSDC00','PSDC01', 'PSDC02', 'PSDC03', 'PSDC04', 'PSDC05', 'PSDC06', 'PSDC07', $
			'PSDC08', 'PSDC09', 'PSDC10', 'PSDC11', 'PSDC12', 'PSDC13', 'PSDC14','PSDC15', $
			'PSDC16', 'PSDC17', 'PSDC18', 'PSDC19', 'PSDC20', 'PSDC21', 'PSDC22','PSDC23', $
			'PSDC24', 'PSDC25', 'PSDC26', 'PSDC27', 'PSDC28', 'PSDC29', 'PSDC30', 'PSDC31',  $
			'PSDC32', 'PSDC33', 'PSDC34', 'PSDC35', 'PSDC36', 'PSDC37', 'PSDC38', 'PSDC39',  $
			'PSDC40', 'PSDC41', 'PSDC42', 'PSDC43', 'PSDC44', 'PSDC45', 'PSDC46', 'PSDC47',  $
			'PSDet', 'PostAnaColl', 'PostMonoColl', 'PreAnaColl', 'PreMonoColl', 'RC', 'SC']
	varname = ['QX', 'QY', 'QZ', 'E', 'A3', 'A4','H','K','L']
	fn = dialog_pickfile(filter='*.bt7',title='Select your data file:')
    if fn ne '' then begin
		nlines = file_lines(fn)
		; creat all the tage names
		name=sindgen(48)
		for n = 0, 47 do begin
			if n lt 10 then name[n]='PSDC0'+string(strtrim(n,2))
		    if n ge 10 then	name[n]='PSDC'+string(strtrim(n,2))
		endfor
		if fn ne '' then begin
    		;CHECK HOW MANY LINES ARE IN THE FILE AND READ THE FILE.
	        s = strarr(nlines)
	        openr,lun,fn,/get_lun
	        readf,lun,s
	        free_lun,lun
	        tagline = 0
	        indepvarline=0
	        referline = 0
	        for i=0,nlines-1 do begin
	        	if stregex(s[i],'#Npoints',/boolean) gt 0 then npoint =i
	        	if stregex(s[i],'#Reference',/boolean) gt 0 then referline =i
	        	if stregex(s[i],'#Scan  ',/boolean) gt 0 then indepvarline = i
	        	if stregex(s[i],'#FixedE  ',/boolean) gt 0 then fixedeline = i
	        	if stregex(s[i],'#Lattice ',/boolean) gt 0 then latticeline = i
	        	if stregex(s[i],'#Orient  ',/boolean) gt 0 then orientline = i
            if stregex(s[i],'#MonoSpacing ',/boolean) gt 0 then monspaceline =i
            if stregex(s[i],'#Columns',/boolean) gt 0 then tagline =i
	        endfor
          npot = strsplit(s[npoint],' ',/extract)
          (*pstate).npoints = fix(npot(1))
          refer = strsplit(s[referline],/extract)
	       	indepvar = strsplit(s[indepvarline],/extract)
	       	segs = strsplit(s[tagline],/extract)
	       	if (*pstate).npoints ge (nlines-tagline) then (*pstate).npoints=(nlines-tagline)-1
	       	if  (*pstate).npoints lt 1 then begin
              void = dialog_message('The Data File is Empty!')
              return
          endif
	       	headnumb = indgen(2)
	       	headnumb(0) = indepvar(1)
          headnumb(1) = refer(1)
          colname = header(headnumb(0)-1)
          varname = [colname,varname]
          refername = header(headnumb(1)-1)
          (*pstate).unitname = refername
          fixe=strsplit(s[fixedeline],' ',/extract)
          sze=size(fixe,/n_elements)
          if sze lt 2 then begin
            (*pstate).ef = 0
            endif else begin
              (*pstate).ef = float(fixe(1))
          endelse
          latt=strsplit(s[latticeline],/extract)
          (*pstate).lattice = latt(1:6)
          ori=strsplit(s[orientline],/extract)
          (*pstate).orient = ori(1:6)
          mond=strsplit(s[monspaceline],/extract)
          (*pstate).monspacing = double(mond[1])
          ntags = n_elements(segs)-1
	        tags = segs[1:ntags]
	        nvar=size(varname)
	        vardata=dblarr(nvar(1),(*pstate).npoints)
	        raw_data= dblarr(48,(*pstate).npoints)
			;look for the reference unit
			wh = where(strlowcase(tags) eq strlowcase(refername),count)
			if count ne 0 then begin
               segs = strsplit(s[tagline+1],/extract)
               (*pstate).unit = double(segs[wh])
            endif else begin
               void = dialog_message('That reference was not found in the file!')
               print,  reference
            endelse

    		;FIND THE POSITION OF THE INPUT COLUMN NAME

			for j = 0, 47 do begin
               	wh = where(strlowcase(tags) eq strlowcase(name(j)),count)
               	if count ne 0 then begin
   					;EXTRACT THE DATA IF THE COLUMN IS AVAILABLE.
                    for i=tagline+1,(tagline+(*pstate).npoints) do begin
                        segs = strsplit(s[i],/extract)
                        raw_data(j,i-tagline-1) = double(segs[wh])
                    endfor;i
               	endif else begin
                   	void = dialog_message('That column name was not found in the file!')
                   	print,  name(j)
           		endelse

 			endfor;j

			;Extarct all the variable for the second data file
			for j = 0, nvar(1)-1 do begin
               	wh = where(strlowcase(tags) eq strlowcase(varname(j)),count)
               	if count ne 0 then begin
   			        for i=tagline+1,(tagline+(*pstate).npoints) do begin
                        segs = strsplit(s[i],/extract)
                        vardata(j,i-tagline-1) = double(segs[wh])
                    endfor;i
               	endif else begin
                   	void = dialog_message('That column name was not found in the file!')
                   	print,  varname(j)
           		endelse

 			endfor;j

		endif
	if colname eq 'A4' || colname eq 'Temp' then daxis_Ch_Eff, event
  if colname eq 'E'  then triaxis_Ch_Eff, event
  if colname eq 'PSDet' then daxis_Ch_Eff, event
 	*(*pstate).data = raw_data
	*(*pstate).vardata = vardata
	(*pstate).colname = colname
	sz = size(raw_data)
	nvars = sz[2]
	widget_control,(*pstate).var_slider,set_slider_max=nvars-1

	endif;fn ne ''

end ;Load_psd_data
;********************************************************
pro Exclude_PSD_Ch, event
	widget_control,event.top,get_uvalue = pstate
	widget_control,(*pstate).BadCh, get_value=BadCh
	if (*pstate).colname eq 'A4' || (*pstate).colname eq 'Temp' then daxis_Ch_Eff, event
  if (*pstate).colname eq 'E'  then triaxis_Ch_Eff, event
  if (*pstate).colname eq 'PSDet' then daxis_Ch_Eff, event
	if strlen(BadCh) gt 0 then begin
	  BadCh = strsplit(BadCh,',',/extract)
    for i=0,size(BadCh,/n_elements)-1 do	(*pstate).Ch_axis_Eff(fix(BadCh(i)))=0
  endif
end;Exclude_PSD_CH
;**********************************************
;this function does the channel efficience correction
pro daxis_Ch_Eff,event
	widget_control,event.top,get_uvalue = pstate
	Ch_Eff = findgen(48)
	mydirectory='c:\psd\'
	fn = strjoin([mydirectory,'PSD_Channeal_Eff.dat'])
	print, 'daxis ',fn
	if fn ne '' then begin
		openr, fun, fn, /get_lun
		readf, fun, Ch_Eff
		free_lun, fun

	endif
   (*pstate).Ch_axis_Eff=Ch_Eff

end; Correct_data

;**********************************************
;this function does the channel efficience correction while using PSD
pro triaxis_Ch_Eff, event
	widget_control,event.top,get_uvalue = pstate
	Ch_dE = findgen(48)
	Ch_Eff = findgen(48)
	Ch_Eff = findgen(48)
	mydirectory='c:\psd\'
	fn = strjoin([mydirectory,'psd_Ef_spacing.dat'])
	if fn ne '' then begin
		openr, fun, fn, /get_lun
		readf, fun, Ch_dE,Ch_eff
		free_lun, fun
	endif
	(*pstate).Ch_axis_Eff=Ch_Eff
	(*pstate).Ch_Eng_Eff=Ch_dE
end; Correct_data


;***************************************************
;this function calculate the hkl for every psd channel in every step
;using the vardata file and the 2theta, E calibration data files
;the new hkl set for every channel is calculated using the algro from Igor
function scalarm, x1, y1, z1, x2, y2, z2,lattemp
	s=x1* x2*lattemp(0)^2+y1*y2*lattemp(1)^2+z1*z2*lattemp(2)^2+$
   		(x1*y2+x2*y1)*lattemp(0)*lattemp(1)*cos(lattemp(5))+$
   		(x1*z2+x2*z1)*lattemp(0)*lattemp(2)*cos(lattemp(4))+$
   		(z1*y2+z2*y1)*lattemp(2)*lattemp(1)*cos(lattemp(3))
	return, s
end;scalar

;calculate the reciprocal lattice parameter
function  star,lat
	V=2*lat(0)*lat(1)*lat(2)*sqrt(sin((lat(3)+lat(4)+lat(5))/2) * $
		sin((-lat(3)+lat(4)+lat(5))/2) * sin((lat(3)-lat(4)+lat(5))/2)*$
		sin((lat(3)+lat(4)-lat(5))/2))
	Vstar = (2*!PI)^3/V
	starlat = findgen(6)
	starlat[0] = 2*!pi*lat(1)*lat(2)*sin(lat(3))/V
	starlat[1] = 2*!pi*lat(0)*lat(2)*sin(lat(4))/V
	starlat[2] = 2*!pi*lat(1)*lat(0)*sin(lat(5))/V
	starlat[3] = acos((cos(lat(4))*cos(lat(5))-cos(lat(3)))/(sin(lat(4))*sin(lat(5))) )
	starlat[4] = acos((cos(lat(3))*cos(lat(5))-cos(lat(4)))/(sin(lat(3))*sin(lat(5))) )
	starlat[5] = acos((cos(lat(3))*cos(lat(4))-cos(lat(5)))/(sin(lat(3))*sin(lat(4))) )

	return, starlat
end; star function

function StandardSystem, lat,orient
	starlat = findgen(6)
	starlat = star(lat)
	orient1=findgen(3)
	orient2=findgen(3)
	orient1=orient(0:2)
	orient2=orient(3:5)
	modx = sqrt(scalarm(orient1(0),orient1(1),orient(2),orient1(0),orient1(1),orient(2),starlat))
	x=orient1/modx

	proj=scalarm(orient2(0),orient2(1), orient(2),x(0),x(1),x(2),starlat)
	y=orient2
	y=y-x*proj
	mody=sqrt(scalarm(y(0),y(1),y(2),y(0),y(1),y(2),starlat))
	if mody le 0 then begin
		void = dialog_message('Fatal error, orienting vectors are colinear!')
		return, 1
	endif
	y=y/mody
	z=y
	z(0)=x(1)*y(2)-y(1)*x(2)
	z(1)=x(2)*y(0)-y(2)*x(0)
	z(2)=x(0)*y(1)-y(0)*x(1)
	proj=scalarm(z(0),z(1),z(2),x(0),x(1),x(2),starlat)
	z(0)=z(0)-x(0)*proj
	z(1)=z(1)-x(1)*proj
	z(2)=z(2)-x(2)*proj
	proj=scalarm(z(0),z(1),z(2),y(0),y(1),y(2),starlat)
	z(0)=z(0)-y(0)*proj
	z(1)=z(1)-y(1)*proj
	z(2)=z(2)-y(2)*proj
	modz=sqrt(scalarm(z(0),z(1),z(2),z(0),z(1),z(2),starlat))
	z=z/modz
	q=[x,y,z]
	return,q
end;standardsystem

;****************************************************************
pro wherespectrometer, event,H=H,K=K,L=L,dE=dE
	widget_control,event.top,get_uvalue = pstate
	widget_control,(*pstate).PsdCenter, get_value=centerch
	centerch=fix(centerch)

	lat=findgen(6)
	orient=findgen(6)
	lat= (*pstate).lattice
	orient= (*pstate).orient
	lat(3:5)=!dtor*lat(3:5)
	Npn=(*pstate).npoints
	anatau=2*!PI/(*pstate).monspacing
	anatau=1.87325
	Ef=(*pstate).Ef
	Ef=14.7
	dE=(*(*pstate).vardata)(4,*)
	A3=!dtor*(*(*pstate).vardata)(5,*)
	A4=!dtor*(*(*pstate).vardata)(6,*)
	Ei=dE+Ef
	Efpsd=findgen(48)
	Ch_eff=findgen(48)
	dA4=findgen(48)
	mydirectory='c:\psd\'
	fn = strjoin([mydirectory,'psd_Ef_spacing.dat'])
	if fn ne '' then begin
		openr, fun, fn, /get_lun
		readf, fun, Efpsd,Ch_eff
		free_lun, fun
	endif
	Ki=sqrt(Ei/2.072142)
	Kfpsd=sqrt(Efpsd/2.072142)

	;calculate the difference between every channel to the center channel in A4
	for i=0,47 do dA4(i)=asin(anatau/2/Kfpsd(i))-asin(anatau/2/Kfpsd(centerch))
	A4=transpose(A4)#(fltarr(48)+1.)+(fltarr(Npn)+1.)#dA4
	A3=transpose(A3)#(fltarr(48)+1.)
	Ki=transpose(Ki)#(fltarr(48)+1)

	Kfpsd=(fltarr(Npn)+1)#Kfpsd
    M2 = 2*asin(anatau/(2*Ki))
   	A2 = 2*asin(anatau/(2*Kfpsd))
   	Q=sqrt(Ki^2+Kfpsd^2-2*Ki*Kfpsd*cos(A4))
	delta=abs(acos((Q^2+Ki^2-Kfpsd^2)/(2*Ki*Q)))
	psi=A3+delta-!pi/2
	qx=Q*cos(psi)
	qy=Q*sin(psi)
	unithkl=StandardSystem(lat,orient)

	unitH=unithkl(0:2)
	unitK=unithkl(3:5)
	unitL=unithkl(6:8)
	H=qx*unitH(0)+qy*unitK(0)
	K=qx*unitH(1)+qy*unitK(1)
	L=qx*unitH(2)+qy*unitK(2)
	Ei=transpose(Ei)#(fltarr(48)+1)
	Ef=(fltarr(Npn)+1)#Efpsd
  	dE=Ei-Ef

end;wherespectro
;********************************************************************
pro xstal3axis, event
	widget_control,event.top,get_uvalue = pstate
	widget_control,(*pstate).PsdCenter, get_value=centerch
	widget_control,(*pstate).PsdLeft, get_value=leftch
	widget_control,(*pstate).PsdRight, get_value=rightch
	wset,(*pstate).winvis
	chrange=indgen(3)
	chrange(0)=leftch
	chrange(1)=centerch
	chrange(2)=rightch
	Npn=(*pstate).npoints
   	data=(*(*pstate).data)
   	data=transpose(data)
   	wherespectrometer, event, H=H, K=K, L=L, dE=dE
   	;Now construct the data to one-dimensional array for using xplot3D
	usedch=abs(chrange(2)-chrange(0)+1)
	xdat=H(*,chrange(0):chrange(2))
	ydat=K(*,chrange(0):chrange(2))
	zdat=dE(*,chrange(0):chrange(2))
	idat=data(*,chrange(0):chrange(2))
	xdat=reform(xdat,usedch*Npn)
	ydat=reform(ydat,usedch*Npn)
	zdat=reform(zdat,usedch*Npn)
    idat=reform(idat,usedch*Npn)
   	device, decomposed=0
	loadct, 13
	zcolors=bytscl(idat,top=!D.N_colors-1)
	surface,dist(5),/Nodat,/save,xrange=[min(H),max(H)],yrange=[min(K),max(K)],zrange=[min(zdat),max(zdat)], $
		xstyle=1, ystyle=1, zstyle=1, charsize=2, color=!D.N_colors-1,xtitle='H',ytitle='K',ztitle='dE'
	plots,xdat,ydat,zdat,psym=4,color=zcolors,symsize=1,/T3D

	(*pstate).xdata = ptr_new(xdat)
	(*pstate).ydata = ptr_new(ydat)
	*(*pstate).stitchdata = zdat
	(*pstate).intdata = ptr_new(idat)


	;Use the code from william to do a object window

;	Earr=diag_matrix(zdat)
;	Iarr=diag_matrix(idat)
;	s = Size(Iarr, /Dimensions)
;	Vert_Colors=Reform(BytScl(Iarr), s[0]*s[1])
;	fsc_surface,Earr,xdat,ydat,colortable=13,/elevation_shading,data2=Iarr ; works ok
;	thispalette = obj_new('idlgrpalette')
;	thispalette->loadct, 13

	;Now construct the data file for the use of surface plot
	;The plot was done using Yiming Qiu's object
	;usedch=abs(chrange(2)-chrange(0)+1)
	;xdat=H(*,chrange(0):chrange(2))
	;ydat=K(*,chrange(0):chrange(2))
	;tempzdat=dE(*,chrange(0):chrange(2))
	;tempidat=data(*,chrange(0):chrange(2))
	;xdat=reform(xdat,usedch*Npn)
	;ydat=reform(ydat,usedch*Npn)
	;zdat=fltarr(usedch*Npn,usedch*Npn)
    ;idat=fltarr(usedch*Npn,usedch*Npn)
    ;for i=0,usedch*Npn-1 do begin
   	;	for j=0,usedch*Npn-1 do begin
   	;		if i eq j then begin
   	;			zdat(i,j)=tempzdat(i)
   	;			idat(i,j)=tempidat(i)
   	;		endif
   	;	endfor
   	;endfor
    ;szidat=size(idat)
   	;szxdat=size(xdat)
   	;tmp = obj_new('dm_plot',xdat,ydat,zdat,idat,/surfplot)
    ;tmp->draw

end; xstal3axis
;*****************************************************************************
function bw2rgb2,grey,noscale=noscale
; grey is an array of gray scale colors
; returns a floating point array, 3 x N
; /noscale indicates grey values are between 0 and 1, and don't need
; to be scaled (default rescales grey array to floats between 0 and 1)

g=bytarr(256)
g(97:162)=255b
g(97-64:96)=indgen(64)*4
g(163:163+63)=252b-indgen(64)*4
r=shift(g,66)
b=shift(g,-66)

if keyword_set(noscale) then begin
	red=interpolate(r,255.99*grey)
	green=interpolate(g,255.99*grey)
	blue=interpolate(b,255.99*grey)
endif else begin
	maxg=max(grey)
	ming=min(grey)
	if maxg ne ming then begin
		scaleg=1.0/(maxg-ming)
	endif else begin
		message,'warning: all the same color, unpredictable results.',/inf
		scaleg=1.0
	endelse
	i=(255.99*(grey-ming)*scaleg)
	red=r(i) & green=g(i) & blue=b(i)
	red=interpolate(r,i)
	green=interpolate(g,i)
	blue=interpolate(b,i)
endelse

; this helps solve problems if 'grey' is higher dimensional
; (although it doesn't actually do what we'd like, which is to
; generate a 3 x whatever x whatever... array.)
result=fltarr(3,n_elements(grey))
result(0,*)=red/255.0
result(1,*)=green/255.0
result(2,*)=blue/255.0

return,result
end
;****************************************************
pro PSD_cleanup,tlb
	widget_control,tlb,get_uvalue = pstate
	wdelete,(*pstate).winpix
	ptr_free, (*pstate).data, (*pstate).stitchdata, (*pstate).vardata
	ptr_free, (*pstate).xdata, (*pstate).ydata, (*pstate).errdata
	ptr_free, (*pstate).image, (*pstate).xptr, (*pstate).yptr
	ptr_free, (*pstate).zoomImage
	heap_free,pstate
end
;**************************************
pro PSD_event,event
	uname = widget_info(event.id,/uname)
	case uname of
		'LoadDat': Load_psd_data, event
		'QUIT':  widget_control,event.top,/destroy
		'PSD_CH': 	PSD_Ch_Plot,event
		'VAR':		PSD_Var_Plot, event
		'ExcludeCh':Exclude_PSD_Ch,event
		'SUM': 	PSD_Sum_Plot, event
		'STICH': PSD_Stitch_Plot, event
		'SaveStich': Save_PSD_Data, event
		'WIN': Psd_plot_refresh, event
	else:
	endcase
end
;*****************************************
;Plot the each channel vs the independent variable
pro PSD_Ch_Plot,event
	widget_control,event.top,get_uvalue = pstate
	Ch_Eff=findgen(48)
	Ch_Eff=(*pstate).Ch_axis_Eff
	datasize = size(*((*pstate).data))
	if datasize(1) le 0 then begin
		void = dialog_message('The Data File is Empty!')
		return
	endif
	widget_control,(*pstate).ch_slider, get_value=slval
	vardata = *((*pstate).vardata)
	final_data = *((*pstate).data)
	x = vardata(0,*)
	y = final_data(slval,*)*Ch_Eff(slval)
	(*pstate).plot_title = 'PSD channel Plot, Channel = '   $
	 		+strtrim(string(slval,format = '(I3)'),2)
	(*pstate).xtitle = 'Independent Varialbe  '+ (*pstate).colname
	(*pstate).ytitle = 'Intensity ( '+ string(strtrim((*pstate).unit, 2))+'  / '+(*pstate).unitname +' )'
	wset,(*pstate).winvis
	(*pstate).xdata = ptr_new(x)
	(*pstate).ydata = ptr_new(y)
	(*pstate).errdata = ptr_new(sqrt(y))
	minx = min(x, max=maxx)
	miny = min(y, max=maxy)

	plot,x,y,psym = 4,thick = 4.0, title = (*pstate).plot_title, xtitle=(*pstate).xtitle,$
		xrange = [minx, maxx], yrange = [miny, maxy], ytitle = (*pstate).ytitle

	errplot, x, (y-sqrt(y)),(y+sqrt(y))
	p = findgen(4)
	sz=size(y)
	Gauss_Init_Guess, x, y, (sz(2)-1), p = p
	PSD_Gauss_Fit, event,x, y, p = p
	image = tvrd(true=1)
	(*pstate).image = ptr_new(image)
	(*pstate).zoomImage=ptr_new(image)
	(*pstate).imageflage = 1
	xrdata=findgen(3,sz(2))
  for k = 0, sz(2)-1 do begin
	  xrdata(0, k) = x(k)
		xrdata(1, k) = y(k)
		xrdata(2, k) = sqrt(y(k))
	endfor
	*(*pstate).stitchdata = xrdata

end; PSD_Ch_Plot
;*******************************************************
;Plot all the psd channel at each indeendent variable step
pro PSD_Var_Plot,event
	widget_control,event.top,get_uvalue = pstate
	Ch_Eff = findgen(48)
	Ch_Eff=(*pstate).Ch_axis_Eff
	datasize = size(*((*pstate).data))
	if datasize(1) le 0 then begin
		void = dialog_message('The Data File is Empty!')
		return
	endif
	widget_control,(*pstate).var_slider, get_value=slval
	widget_control,(*pstate).PsdCenter, get_value=centerch
	if centerch lt 0 || centerch gt 47 then centerch =23
	x = findgen(48)
	final_data = *((*pstate).data)
	vardata = *((*pstate).vardata)
	A4=(*(*pstate).vardata)(6,slval)
	ch_space=findgen(48)
	ch_position=findgen(48)
	mydirectory='c:\psd\'
	fn = strjoin([mydirectory,'PSD_A4_Spacing.dat'])
	if fn ne '' then begin
		openr, fun, fn, /get_lun
		readf, fun, ch_position
		free_lun, fun
	endif
	;Now calculate the spacing between every channel and the specified center channel


	for k = 0, 47 do begin
		ch_space(k) = ch_position(fix(centerch))-ch_position(k)+A4
	endfor

	sz=size(final_data)
	if slval ge sz(2) then slval = sz(2)-1
	y = final_data(0:47,slval)
	wset,(*pstate).winvis
	(*pstate).xdata = ptr_new(ch_space)
	(*pstate).ydata = ptr_new(y)
	(*pstate).errdata = ptr_new(sqrt(y))
	xrdata=findgen(3,48)
    for k = 0, 47 do begin
		xrdata(0, k) = ch_space(k)
		xrdata(1, k) = y(k)
		xrdata(2, k) = sqrt(y(k))
	endfor
	*(*pstate).stitchdata = xrdata
	minx = min(ch_space, max=maxx)
	miny = min(y, max=maxy)
	(*pstate).plot_title = 'PSD Variable Plot,'$
		 +strtrim(string((*pstate).colname), 2) + ' = ' $
	     +strtrim(string(vardata(0,slval),format = '(f10.3)'),2)
	(*pstate).xtitle = 'A4 '
	(*pstate).ytitle ='Neutron ( '+ string(strtrim((*pstate).unit, 2))+'/'+(*pstate).unitname +' )'
	plot,ch_space,y,psym = 4,thick = 4.0, ymargin=[5, 3],subtitle = (*pstate).plot_title, xtitle=(*pstate).xtitle, $
		xrange = [minx, maxx], yrange = [miny, maxy], ytitle=(*pstate).ytitle,xstyle=1
	errplot, ch_space, (y-sqrt(y)),(y+sqrt(y))
	axis, xaxis=1, xrange=[47,0], xtitl='PSD Channel Number',xstyle=1
	p = findgen(4)
	sz=size(y)
	Gauss_Init_Guess, ch_space, y, (sz(1)-1), p = p
	PSD_Gauss_Fit, event,ch_space, y, p=p
	image = tvrd(true=1)
	(*pstate).image = ptr_new(image)
	(*pstate).zoomImage=ptr_new(image)
	(*pstate).imageflage = 1

end;PSD_Var_Plot

;**********************************************************
;Plot the sum of all the psd channel at each independent variable step
pro PSD_Sum_Plot, event
	widget_control,event.top,get_uvalue = pstate
	datasize = size(*((*pstate).data))
	if datasize(1) le 0 then begin
		void = dialog_message('The Data File is Empty!')
		return
	endif
	widget_control,(*pstate).PsdLeft, get_value=leftch
	widget_control,(*pstate).PsdRight, get_value=rightch

	Ch_eff=findgen(48)
	Ch_eff=(*pstate).Ch_2axis_Eff

	minch=leftch
	maxch=rightch
	range=indgen(2)
	range(0) = minch
	range(1) = maxch
	if range(0) lt 0 || range(0) gt 47 then begin
		void = dialog_message('The left of the PSD is out of [0 47], set to default value [0 47]')
		range(0) = 0
		range(1) = 47
		widget_control,(*pstate).PsdLeft,set_value = '0'
		widget_control,(*pstate).PsdRight,set_value = '47'
	endif
	if range(1) lt range(0) || range(1) gt 47 then begin
		void = dialog_message('The left of the PSD is greater right, set default value to [0 47]')
		widget_control,(*pstate).PsdLeft,set_value = '0'
		widget_control,(*pstate).PsdRight,set_value = '47'
		range(0) = 0
		range(1) = 47
	endif
	data=(*(*pstate).data)
	vardata=(*(*pstate).vardata)
	sz=size(data)
	y=findgen(sz(2))
	x = vardata(0,*)
	for i=0, sz(2)-1 do begin
		sum=0
		for j = range(0), range(1) do begin
			sum=sum+data(j,i)
		endfor
		y(i)=sum
	endfor
	wset,(*pstate).winvis
	(*pstate).plot_title = string('Sum of the PSD Channel')
	(*pstate).xtitle = (*pstate).colname
	(*pstate).ytitle = 'Intensity ( '+ string(strtrim((*pstate).unit, 2))+'  / '+(*pstate).unitname +' )'
	minx = min(x, max=maxx)
	miny = min(y, max=maxy)
	(*pstate).xdata = ptr_new(x)
	(*pstate).ydata = ptr_new(y)
	(*pstate).errdata = ptr_new(sqrt(y))
	xrdata=findgen(3,fix(sz(2)))
    for k = 0, sz(2)-1 do begin
		xrdata(0, k) = x(k)
		xrdata(1, k) = y(k)
		xrdata(2, k) = sqrt(y(k))
	endfor
	*(*pstate).stitchdata = xrdata
	plot,x,y,psym = 4,thick = 4.0, title = (*pstate).plot_title, xtitle = (*pstate).xtitle, $
		xrange = [minx, maxx], yrange = [miny, maxy],  $
		ytitle= (*pstate).ytitle
	errplot, x, (y-sqrt(y)),(y+sqrt(y))
	;image = tvrd(true=1)
	(*pstate).image = ptr_new(image)
	(*pstate).zoomImage=ptr_new(image)
	(*pstate).imageflage = 1
end;PSD_Sum_Plot
;******************************************************
; This code will rebin the data regarding to A4. the channel spacing is read out from
; and data file. The output is a equal step 0.1 data file with A4 between the specified
; A4 minimium and A4 maximum. (The data points from the first front half of the psd and
; the last step of rear psd is removed.


pro PSD_Stitch_Plot, event
	widget_control,event.top,get_uvalue = pstate
	datasize = size(*((*pstate).data))
	if datasize(1) le 0 then begin
		void = dialog_message('The Data File is Empty!')
		return
	endif
	widget_control,(*pstate).PsdCenter, get_value=centerch
	widget_control,(*pstate).PsdLeft, get_value=leftch
	widget_control,(*pstate).PsdRight, get_value=rightch

	range=indgen(3)
	range(0) = leftch
	range(1) = centerch
	range(2) = rightch
	if range(1) lt 0 || range(1) gt 47 then begin
		void = dialog_message('The center of the PSD is out of [0 47], set to default value [0 23 47]')
		widget_control,(*pstate).PsdCenter,set_value = '23'
		widget_control,(*pstate).PsdLeft,set_value = '0'
		widget_control,(*pstate).PsdRight,set_value = '47'
		range(0) = 0
		range(1) = 23
		range(2) = 47
	endif
	if range(0) ge range(1) || range(0) ge range(2) || range(0) lt 0 || range(0) gt 47 then begin
		void = dialog_message('The left of the PSD is greater center, set default value to [0 23 47]')
		widget_control,(*pstate).PsdCenter,set_value = '23'
		widget_control,(*pstate).PsdLeft,set_value = '0'
		widget_control,(*pstate).PsdRight,set_value = '47'
		range(0) = 0
		range(1) = 23
		range(2) = 47
	endif
	if range(2) lt range(0) || range(2) lt range(1) || range(2) lt 0 || range(2) gt 47 then begin
		void = dialog_message('The left of the PSD is greater than the right, set default value to 0, 47')
		widget_control,(*pstate).PsdCenter,set_value = '23'
		widget_control,(*pstate).PsdLeft,set_value = '0'
		widget_control,(*pstate).PsdRight,set_value = '47'
		range(0) = 0
		range(1) = 23
		range(2) = 47
	endif

	Ch_eff=findgen(48)
	mydirectory='c:\psd\'
	;fn = strjoin([mydirectory,'PSD_Channeal_Eff.dat'])
	;if fn ne '' then begin
;		openr, fun, fn, /get_lun
;		readf, fun, Ch_eff
;		free_lun, fun
;	endif
	Ch_eff=(*pstate).Ch_axis_Eff
	if (*pstate).colname eq 'A4' then begin
		data = (*(*pstate).data)
		var =(*(*pstate).vardata)
		sz = size(data)
		data_err = dindgen(sz(1), sz(2))
		data_err = sqrt(data)

		A4_begin = var(0, 0)
		A4_end = var(0, sz(2)-1)

		if A4_begin eq A4_end then begin
			void = dialog_message('The step size of A4 is 0, could not stitch !')
			return
		endif
		output_width = 0.1
		output_npt = round(abs(A4_end-A4_begin)/output_width)
		output_data = dblarr(output_npt)
		data_norm=dblarr(output_npt)+1
		output_data_err = dblarr(output_npt)
		output_data_left = findgen(output_npt+1)
		output_data_right = findgen(output_npt)
		;data_plus_count = findgen(output_npt)
		dis = dblarr(48, output_npt)
		frac = dblarr(48, output_npt)
		z_in = dblarr(48)
		dz_in = dblarr(48)
		output_data_left=min([A4_begin,A4_end])+output_data_left*output_width
		;for i = 0, (output_npt-1) do begin
		;	if A4_begin lt A4_end then begin
		;		output_data_left(i) = A4_begin + i * output_width
		;	endif
		;	if A4_begin gt A4_end then begin
		;		output_data_left(i) = A4_begin - i * output_width
		;	endif
		;endfor
		output_data_left(output_npt) = max([A4_begin,A4_end])
		output_data_right = output_data_left(1:output_npt)
		;output_data_left=reverse(output_data_left)
		;Now open the A4 spacing calibration file and readout
		ch_position = dindgen(48)
		ch_space = dindgen(48)
		mydirectory='c:\psd\'
		fn = strjoin([mydirectory,'PSD_A4_Spacing.dat'])
		if fn ne '' then begin
			openr, fun, fn, /get_lun
			readf, fun, ch_position
			free_lun, fun
		endif
		;Now calculate the spacing between every channel and the specified center channel
		;for k = 0, 47 do begin
		;	ch_space(k) = ch_position(range(1))-ch_position(k)
		;endfor
		ch_space=ch_position(range(1))-ch_position
		ch_space_extend=fltarr(49)
		ch_space_extend(0:47)=ch_space
		ch_space_extend(48)=ch_space_extend(47)
		ch_space=ch_space_extend
		; After align the middle of PSD with the every A4 position and then build
		; the left and right coordinate for every psd position
		ch_left = dblarr(49)
		ch_right =dblarr(48)
		;for i = 0, output_npt-1 do begin
		;	data_plus_count(i) = 0.0
		;endfor
		data_plus_count=dblarr(output_npt)
		mon_in=ch_left(0:47)+1.0
		for l=0,47 do begin
			if Ch_eff(l) gt 1E-2 then begin
			mon_in(l)=0
			endif
			endfor ; for
		output_tmp=data_plus_count
		dz_output_tmp=data_plus_count
		output_mon=data_plus_count
		dz_output_mon=data_plus_count
		print,'sz=',sz(2)
		for i = 0, sz(2)-1 do begin
			;for j = 1, 47 do begin
			;	ch_left(j) = ch_space(j) + 0.5*(ch_space(j-1) -ch_space(j))+var(0,i)
			;endfor
			ch_left=ch_space+0.5*(-ch_space+shift(ch_space,1))+var(0,i)
			ch_left(0) = ch_space(0) + 0.5*(ch_space(0) - ch_space(1)) +var(0, i)
			ch_left(48) = ch_space(47) - 1.0*(ch_space(46)-ch_space(47)) + var(0, i)
			ch_right = ch_left(1:48)
			z_in = data(*, i)*Ch_eff
			dz_in = data_err(*,i)*Ch_eff
			;help,ch_left, z_in,dz_in,output_data_left
			ch_left=reverse(ch_left)
			z_in=reverse(z_in)
			dz_in=reverse(dz_in)
			;output_data_left=reverse(output_data_left)
			drebin,ch_left,z_in,dz_in,output_data_left,output_tmp,dz_output_tmp,/histogram,/to_histogram,err=err,emsg=emsg
			print, emsg
			drebin,ch_left,mon_in,mon_in,output_data_left,output_mon,dz_output_mon,/histogram,/to_histogram,err=err,emsg=emsg
			output_data=output_data+output_tmp*output_mon
			output_data_err=output_data_err+dz_output_tmp*output_mon
			data_norm=data_norm+output_mon
			;print, 'emsg2 ',emsg
			;print,'output_mon ',n_elements(output_mon)
			;print,'output_tmp ',n_elements(output_tmp)
			;print,'ch_left',(ch_left)
			;print,'data',n_elements(z_in)
			;drebin_histo,x_in,z_in,dz_in,x_out,z_out,dz_out
			;for j = range(0), range(2) do begin
			;	for k = 0, output_npt-1 do begin
			;		min_dis_righ = min([ch_left(j), output_data_right(k)])
			;		max_dis_left = max([ch_right(j), output_data_left(k)])
			;		dis(j, k) = max([0, (min_dis_righ-max_dis_left)])
			;		frac(j, k) = dis(j, k)/abs(ch_left(j+1)-ch_left(j))
;
;					output_data(k) = output_data(k)+z_in(j)*frac(j, k)
;					output_data_err(k)= output_data_err(k) + (dz_in(j))^2*frac(j,k)
;					if Ch_eff(j) gt 1E-2 then begin
;						;print,j,Ch_eff(j)
;						 data_plus_count(k) = data_plus_count(k) +frac(j, k)
;					endif
;				endfor
;			endfor
		endfor; i
		output_data=output_data/data_norm
		output_data_err=output_data_err/data_norm
		;for k = 0, output_npt-1 do begin
		;	if data_plus_count(k) eq 0 then begin
		;		print, 'Frac(', k,')=', 0
		;		continue
		;	endif
		;	output_data(k) = output_data(k)/data_plus_count(k)
		;	output_data_err(k)=sqrt(output_data_err(k))/data_plus_count(k)
		;endfor
		xrdata = dindgen(3, output_npt)
		for k = 0, output_npt-1 do begin
			xrdata(0, k) = output_data_left(k)
			xrdata(1, k) = output_data(k)
			xrdata(2, k) = output_data_err(k)
		endfor
		*(*pstate).stitchdata = xrdata
		wset,(*pstate).winvis
		minx = min(xrdata(0,*), max=maxx)
		miny = min(xrdata(1,*), max=maxy)
		(*pstate).xdata = ptr_new(xrdata(0,*))
		(*pstate).ydata = ptr_new(xrdata(1,*))
		(*pstate).errdata = ptr_new(xrdata(2,*))
		(*pstate).plot_title = (*pstate).colname +' Stitch Plot'
		(*pstate).xtitle = 'Independent Variable' +	(*pstate).colname
		(*pstate).ytitle = 'Intensity ( '+ string(strtrim((*pstate).unit, 2))+' / '+(*pstate).unitname +' )'
		plot, xrdata(0, *), xrdata(1,*), psym = 4, thick = 4.0, $
			title = (*pstate).plot_title, xtitle = (*pstate).xtitle, $
			xrange = [minx, maxx], yrange = [miny, maxy], $
			ytitle= (*pstate).ytitle
		errplot, xrdata(0, *), (xrdata(1,*)-xrdata(2, *)),(xrdata(1,*)+xrdata(2, *))
		image = tvrd(true=1)
		(*pstate).image = ptr_new(image)
		(*pstate).zoomImage=ptr_new(image)
		(*pstate).imageflage = 1
	endif; for A4

	;now deal with the situation for single crystal Ei scan
	if (*pstate).colname eq 'E' then begin
		xstal3axis, event
		image = tvrd(true=1)
		(*pstate).image = ptr_new(image)
		(*pstate).zoomImage=ptr_new(image)
		(*pstate).imageflage = 1
	endif; for independent variable is E
	;Now deal with the independent is temperature situation
	if (*pstate).colname eq 'Temp' then begin



	endif;while independent variable is temp


end; Rebin_PSD_Data
;*************************************************
; used to save the XRD pattern
pro Save_PSD_Data, event
	widget_control,event.top,get_uvalue = pstate
	datasize = size(*((*pstate).stitchdata))
	if datasize(1) le 0 then begin
		void = dialog_message('The Stitched Data File is Empty!')
		return
	endif
	xrdata = (*(*pstate).stitchdata)
	fn = dialog_pickfile(default_extension='dat',/write, title='Select the data file to save XRD pattern:')
   	if (*pstate).colname eq 'A4' then begin
   		if fn ne '' then begin
       		openw,lun,fn,/get_lun
	   		printf,lun,xrdata
	   		free_lun,lun
		endif
	endif
	if (*pstate).colname eq 'E' then begin
	   if fn ne '' then begin
       		openw,lun,fn,/get_lun
       		printf, lun, 'X	','Y ','Z	', 'Intensity'
       		sz=size((*(*pstate).xdata))
       		for i=0, sz(1)-1 do printf,lun,(*(*pstate).xdata)(i),(*(*pstate).ydata)(i), (*(*pstate).stitchdata)(i),(*(*pstate).intdata)(i)
	   		free_lun,lun
		endif
	endif; for independent variable is E
	if (*pstate).colname eq 'Temp' then begin
   		if fn ne '' then begin
       		openw,lun,fn,/get_lun
	   		printf,lun,xrdata
	   		free_lun,lun
		endif
	endif
end;save_psd_data

;*************************************************
pro Gauss_Init_Guess, x, y, loop, p=p
	p(0)=min(y)
	dis1 = 10000000
	dis2 = 10000000
	for i=0, loop do begin
		if abs(y(i)-max(y)) le dis1 then begin
			dis1 = abs(y(i)-max(y))
			xcen = i
		endif
		if abs(y(i)-(max(y)-min(y))/2) le dis2 then begin
			dis2 = abs(y(i)-(max(y)-min(y))/2)
			xwidth = i
		endif
	endfor
	p(1) = x(Xcen)
	p(2) = 2*abs(x(Xcen)-x(Xwidth))/sqrt(Alog(4))
	p(3) = max(y) * p(2)

end;Gaussian_Init_Guess

pro Gaussian_function,x,p, yfit
	y0=p(0)
	xc=p(1)
	w=p(2)
	A=p(3)
    yfit = y0+A/(w*sqrt(!PI/2))*exp(-2*(x-xc)^2/w^2)

end;fit_gaussian
;********************************************
pro PSD_Gauss_Fit,event, x, y, p=p
	widget_control,event.top,get_uvalue = pstate
	device,get_decomposed = dc
	device,decomposed = 1
	green = ((2L)^8 - 1)*256L
	dy=sqrt(y)
	; Fit it using the MPCURVEFIT
	yfit = mpcurvefit(x,y,dy,p,sigma,/noderivative,  $
   		function_name = 'Gaussian_function',/quiet,status = status)
	fit_result = p

	wset,(*pstate).winvis

	if n_elements(sigma) eq 4 then begin
		oplot,x,yfit,color=green, linestyle = 0,thick = 4.0
		xyouts,(min(x)),(0.95*max(y)),'Background: '+strtrim(string(p(0)),2)+ $
	 		 ' +/- '+strtrim(string(sigma(0)),2), font = 5
		xyouts,(min(x)),(0.9*max(y)),'FWHM: '+strtrim(string(p(2)*sqrt(Alog(4))),2)+  $
	 		 ' +/- '+strtrim(string(sigma(2)),2), font = 5
	 	xyouts,(min(x)),(0.85*max(y)),'Area: '+strtrim(string(p(3)),2)+  $
	 		' +/- '+strtrim(string(sigma(3)),2), font = 5
	 	xyouts,(min(x)),(0.8*max(y)),'Center: '+strtrim(string(p(1)),2)+  $
	 		' +/- '+strtrim(string(sigma(1)),2), font = 5

	endif

end; PSD_Gauss_Fit
;***********************************************************************
;Plot the zoom in funciton in the graphic window
pro Psd_plot_refresh, event
	widget_control, event.top, get_uvalue=pstate
	etype = tag_names(event,/structure_name)
	wset,(*pstate).winvis
	if etype eq 'WIDGET_DRAW' && (*pstate).imageflage eq 1 then begin
		if event.press eq 1 then begin
			(*pstate).press = 1
			(*pstate).zoomed = 0
			(*pstate).pressx = event.x
			(*pstate).pressy = event.y
			(*pstate).currentx = event.x
			(*pstate).currenty = event.y

		endif;press

		if event.release eq 1 then begin
			(*pstate).press = 0

			!x = *(*pstate).xptr
			!y = *(*pstate).yptr

			pressr   = convert_coord((*pstate).pressx,(*pstate).pressy,/device,/to_data)
			releaser = convert_coord((*pstate).currentx,(*pstate).currenty,/device,/to_data)

			minx = min([pressr[0],releaser[0]],max=maxx)
			miny = min([pressr[1],releaser[1]],max=maxy)

			plot,*(*pstate).xdata,*(*pstate).ydata,psym = 4, thick = 4.0, $
				xrange = [minx, maxx], yrange= [miny,maxy], $
				title = (*pstate).plot_title, xtitle = (*pstate).xtitle, $
		 		ytitle= (*pstate).ytitle
		 	errplot, *(*pstate).xdata, (*(*pstate).ydata)-(*(*pstate).errdata),(*(*pstate).ydata) $
		 			+(*(*pstate).errdata)
		 	;image = tvrd(true=1)
			if ptr_valid(zoom_image) then ptr_free,zoom_image
			(*pstate).zoomImage = ptr_new(image)
			(*pstate).zoomed = 0
		endif;release

		;IF THE MOUSE IS HELD PRESSED, THEN PERFORM THE DRAG EVENTS
		if (*pstate).press eq 1 then begin

			;RECORD THE CURRENT MOUSE POSITION
			(*pstate).currentx = event.x
			(*pstate).currenty = event.y
			x1 = [(*pstate).pressx,(*pstate).currentx,(*pstate).currentx,(*pstate).pressx,(*pstate).pressx]
			y1 = [(*pstate).pressy,(*pstate).pressy,(*pstate).currenty,(*pstate).currenty,(*pstate).pressy]

			green = ((2L)^8 - 1)*256L
			tv,*(*pstate).zoomImage,true=1
			plots,x1,y1,color=green,/device
			(*pstate).zoomed = 1
		endif

		;FINALLY CHECK FOR RIGHT MOUSE BUTTON PRESS
		if event.press eq 4 then begin
			(*pstate).zoomed = 0
		   tv,*(*pstate).image,true=1
		endif;event.press
	endif;etype

end; psd_plot_refresh,
;*******************************************************

pro crystalbt7psdbackup,group_leader = group_leader

	if n_elements(group_leader) eq 0 then group_leader = 0L

	if xregistered('PSD_plot') then return
	tlb = widget_base(/row,title = 'BT7 PSD data reduction', group_leader = group_leader)
	col_base = widget_base(tlb,/col, space =2)
	void = widget_button(col_base, value = 'Load Data File', font ='8', uname='LoadDat')
	void = widget_label(col_base, /align_center, font='8', Value = 'View Raw Data')
	RawCh_base = widget_base(col_base, /row, space =15)
	void = widget_label(Rawch_base, value = 'Channel',/align_center)
	Ch_slider = widget_slider(Rawch_base,/drag, maximum = 47, minimum=0, xsize=120, value=10, uname='PSD_CH' )
	RawVar_base = widget_base(col_base, /row, space =10)
	void = widget_label(RawVar_base, value = 'Indep var')
	Var_slider = widget_slider(RawVar_base,/drag, maximum=300, xsize =120, value = 0,uname='VAR')
	void = widget_label(col_base,/align_center, font='8',value = 'Used PSD Channels')
	subrow_base = widget_base(col_base, /row, space =27)
 	void=widget_label(subrow_base,/align_center,value = 'Center Channel')
 	PsdCenter  = widget_text(subrow_base,frame = 6, xsize=10, value='23',/editable,uname='Psd_Center')
 	leftrow_base = widget_base(col_base, /row, space = 40)
 	void=widget_label(leftrow_base,/align_center,value = 'Left Channel')
 	PsdLeft  = widget_text(leftrow_base,value='0',/editable,xsize =10, uname='Psd_Left')
 	rightrow_base = widget_base(col_base, /row, space=32)
	void=widget_label(rightrow_base,/align_center,value = 'Right Channel')
	PsdRight  = widget_text(rightrow_base,value='47',/editable,xsize = 10, uname='Psd_Right')
	void = widget_button(col_base,/align_center, font='8',value = 'Exclude PSD Channels',uname='ExcludeCh')
	BadCh  = widget_text(col_base, value='',/editable,xsize = 10, uname='BadChannel')
	void = widget_button(col_base,value = 'Sum Channel', font='8',uname = 'SUM')
	void = widget_button(col_base,value = 'Stich Pattern', font='8',uname = 'STICH')
	void = widget_button(col_base, value = 'Save Stitch Pattern', font='8',uname = 'SaveStich')
	void = widget_button(col_base,value = 'Quit', font='8',uname = 'QUIT')
	win = widget_draw(tlb,xsize = 700,ysize = 700, $
			/motion_events,/button_events,/tracking_events,	uname = 'WIN')
	widget_control,tlb,/realize
	widget_control,win,get_value = winvis
	window,/free,/pixmap,xsize = 500,ysize = 500
	winpix = !d.window
	stitchdata = ptr_new(/allocate_heap)
	data=ptr_new(/allocate_heap)
	vardata = ptr_new(/allocate_heap)
	state =  {  winvis:		winvis,     $
            winpix:			winpix,     $
            win:			win,		$
            npoints:		0,			$ ;the number of column of data in the file after #Ncolumns
            colname: 		''        , $ ;the independent variable name after #Scan
            unit: 			0L		  , $ ;time or monitor after #reference
            unitname: 		'',			$ ;
            ef:				14.7, 		$ ;
            lattice:		findgen(6), $ ; The lattice parameter and angles
            orient:         indgen(6),	$ ; the orientation plane
            monspacing:     0L,			$ ;The d-spacing of the monochromator
            Ch_axis_Eff:	findgen(48),$ ; instore the channel efficiency
            Ch_Eng_Eff:		findgen(48),$ ; Instore the engery diffence for every channel in 3-axis mode
            plot_title: 	'',         $
            xtitle:			'',			$
            ytitle:			'',			$
            vardata:		vardata,    $ ; store all the other variables
            data: 			data,		$
            var_slider:		var_slider, $
            ch_slider:		ch_slider,  $
            PsdCenter:		PsdCenter,  $
            PsdLeft: 		PsdLeft,  	$
            PsdRight: 		PsdRight,  	$
            BadCh:			BadCh,      $
            stitchdata: 	stitchdata, $
            imageflage:		0,			$   ;image flage =1 if there is a actived image in draw window
            image:			ptr_new(/allocate_heap),$
            zoomImage: 		ptr_new(/allocate_heap), $
            xptr:			ptr_new(!x),$	;KEEP TRACK OF AXES SYSTEM VARIABLES
            yptr:			ptr_new(!y),$
            press:			0,			$	;THIS IS THE WAY TO KEEP TRACK OF DRAGGING
            pressx:			0,			$	;KEEP TRACK OF POSITION OF MOUSE PRESS
            pressy:			0,			$
            currentx:		0,			$	;KEEP TRACK OF POSITION OF CURRENT MOUSE POSITION
            currenty:		0,			$
            zoomed:			0,			$;KEEP TRACK OF WHETHER PLOT IS CURRENTLY ZOOMED.
            xdata:			ptr_new(/allocate_heap), $  ;the x, y data used to draw the current graph
            ydata:			ptr_new(/allocate_heap), $
            errdata:		ptr_new(/allocate_heap), $
            zdata:			ptr_new(/allocate_heap), $
            psd_a4_spacing:  'PSD_A4_Spacing.dat',$
            psd_channeal_Eff: 'PSD_Channeal_Eff.dat',$
            psd_Ef_spacing:'psd_Ef_spacing.dat',$
            psd_directory:'c:\psd\'            
            intdata:		ptr_new(/allocate_heap)		  }
	pstate = ptr_new(state)
	widget_control,tlb,set_uvalue = pstate

	xmanager,'getbt7psd',tlb,/no_block, event_handler = 'PSD_event', cleanup = 'PSD_cleanup'

end; getbt7psd