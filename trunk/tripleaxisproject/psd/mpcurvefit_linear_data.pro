pro linear_function,x,p,yfit
; The model function that defines the
; line.
offset = p[0] & slope = p[1]
yfit = offset+slope*x
end

pro mpcurvefit_linear_data
; Create some fake linear data
offset = 1500.0 & slope = -0.25
npts = 100
x = makepts(xlo = -20.0,xhi = 30.0,npts = npts)
yline = offset+slope*x
; "Color" the data with Poisson noise
y = fltarr(npts)
for i = 0,npts-1 do y[i] = randomn(s,1,poisson = yline[i])
dy = sqrt(y)
plot,x,y,psym = 4
errplot,x,y-dy,y+dy,width = 0.0
; Provide initial guess for the fit parameters
p = [100.0,1.0]
; Fit it using the MPCURVEFIT
yfit = mpcurvefit(x,y,1d/dy^2,p,sigma,/noderivative,  $
   function_name = 'linear_function',/quiet,status = status)
fit_result = p

oplot,x,yfit,psym = 0,thick = 2.0

xyouts, [-10],[100],'This is the fit result'; print the line in the graph
print,"Fit results"
print,"-----------"
print,'Offset: '+strtrim(string(fit_result[0]),2)+ $
   ' +/- '+strtrim(string(sigma[0]),2)
print,'Slope: '+strtrim(string(fit_result[1]),2)+  $
   ' +/- '+strtrim(string(sigma[1]),2)

;; Write the jpeg
;tvlct,r,g,b,/get
;image2d = tvrd()  ; tvrd read the current windows to a 2 by 2 matrix
;help, image2d
;s = size(image2d,/dimensions)
;print, s[0], s[1]
;image24 = bytarr(3,s[0],s[1])
;tvlct,r,g,b,/get
;image24[0,*,*] = 255B - r[image2d]
;image24[1,*,*] = 255B - g[image2d]
;image24[2,*,*] = 255B - b[image2d]
;
;filename = 'c:\idl_training\ornl\documentation\mplinfit.jpg'
;write_jpeg,filename,image24,quality = 100,true = 1
;write_jpeg,fn, tvrd()
;set_plot, ps
;device,
;filename,
;plot;
;close
;set_plot to the windows
;look for the help device for an example of writing ps file


;An example of writing ps file
myDevice =!d.name  ; record the current graphic device
print, myDevice
set_plot, 'PS'     ; set the plot to ps
device, filename='C:\IDLCourse\linearfit.ps', /landscape  ; set the file name
plot, x, yfit, psym=1
;errplot,
;oplot,
device, /close
set_plot, myDevice


end

; fn=diag_pickfile()
; write_