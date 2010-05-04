function bsol,temp,p
    ;Tn,Jt,Nf,Bk=p
    ;print, 'bsol p ', p
    Tn=p(0) & Jt=p(1) & Nf=p(2) & Bk=p(3)
    t=4.0*(Jt/(Jt+1.0))*Tn/temp
    ;print, t
    if (Tn le 0) or (Jt le 0) or  (temp ge Tn) then begin
        xout=0.0
    endif else begin
        xout=zbrentp(0.0,t,FUNC='bfun',temp=temp,myp=p,tol=1e-6)
	endelse
    return, xout
    end

function bfun,x,T,p
;print, x
;print, T
;print, 'bfun p ',p
    Tn=p(0) & Jt=p(1) & Nf=p(2) & Bk=p(3)
    if x eq 0.0 then begin
        B=-1.0;   /* so that it wont find solution at zero */
     endif else begin
        B=(x-3*brill(Jt,x)*(Jt/(Jt+1))*(Tn/T))
    endelse
    return, B
    end

 function brill,j,x
    ;print 'Tn=',Tn,' Jt=',Jt,' Nf=',Nf,' Bk=',Bk
    temp=(2*j+1.0)/2/j
    ;print x
    if x eq 0 then begin
        Br=0.0
    endif else begin
        Br=temp/tanh(temp*x)-1.0/tanh(x/2/j)/2/j
    endelse
    return, Br
	end
 function Intensity,T,p
    Tn=p(0) & Jt=p(1) & Nf=p(2) & Bk=p(3)
    ;print, Tn,Jt,Nf,Bk
    br=brill(Jt,bsol(T,p))
    bout=Bk+Nf*br^2
    return, bout
    end


pro meanfield
p=[50.0,0.5,100.0,0.0]
print, Intensity(20.0,p)
T=findgen(100)+1
Icalc=fltarr(n_elements(T))
for i=0, n_elements(T)-1 do begin
	Icalc(i)=Intensity(T(i),p)
endfor
print, Icalc
print, T
plot, T,Icalc
end