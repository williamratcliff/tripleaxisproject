; $Id: drebin.pro,v 1.4 2005/02/03 18:04:12 azuah Exp $
;###############################################################################
; LICENSE:
;  The software in this file is written by an employee of 
;  National Institute of Standards and Technology 
;  as part of the DAVE software project.
;
;  The DAVE software package is not subject to copyright protection
;  and is in the public domain. It should be considered as an
;  experimental neutron scattering data reduction, visualization, and
;  analysis system. As such, the authors assume no responsibility
;  whatsoever for its use, and make no guarantees, expressed or
;  implied, about its quality, reliability, or any other
;  characteristic. The use of certain trade names or commercial
;  products does not imply any endorsement of a particular product,
;  nor does it imply that the named product is necessarily the best
;  product for the stated purpose. We would appreciate acknowledgment
;  if the DAVE software is used of if the code in this file is
;  included in another product.
;
;###############################################################################
; Written by J.R.D. Copley.
;************************************************************************************************
pro drebin_histo,x_in,z_in,dz_in,x_out,z_out,dz_out
;************************************************************************************************
; This procedure rebins 1-dimensional data. The data to be rebinned are described by the
; arrays x_in, which is dimensioned n_in+1, and by z_in and dz_in, both of which are dimensioned
; [n_in,ng], where n_in is the number of channels and ng is the number of groups, i.e. sets of
; data. The grid onto which the data are to be rebinned is defined by the array x_out and the
; results are placed in z_out and dz_out; x_out is dimensioned n_out+1 and z_out and dz_out are
; dimensioned [n_out,ng] where n_out is the number of output channels. It is up to the user of
; this procedure to ensure that the input arrays are properly dimensioned.
;
compile_opt strictarr
;
;
; Generate number, left limit, right limit and width in x of input channels.
n_in=n_elements(x_in)-1
l_in=x_in[0:n_in-1]
r_in=x_in[1:n_in]
w_in=r_in-l_in
;
; Generate number, left limit, right limit and width in x of output channels.
n_out=n_elements(x_out)-1
l_out=x_out[0:n_out-1]
r_out=x_out[1:n_out]
w_out=r_out-l_out
;
; Determine the number of sets of input data.
ng=n_elements(z_in)/n_in
;
; Create unit vectors for later use.
uvec_nout=intarr(n_out)+1
uvec_nin=intarr(n_in)+1
uvec_ng=intarr(ng)+1
;
;	Calculate integrated counts and associated errors for input channels.
t_in=z_in*(w_in#uvec_ng)
dt_in=dz_in*(w_in#uvec_ng)
;
;	The array element (*f_limPtr)[i,j] gives the fraction of input channel j
;		that contributes to output channel i.
r_limPtr=ptr_new((r_out#uvec_nin) < (uvec_nout#r_in),/no_copy)
l_limPtr=ptr_new((l_out#uvec_nin) > (uvec_nout#l_in),/no_copy)
d_limPtr=ptr_new((*r_limPtr-*l_limPtr) > 0,/no_copy)
f_limPtr=ptr_new(*d_limPtr/(uvec_nout#w_in),/no_copy)
;
;	Do the actual rebinning.
; Integrated counts and associated errors are computed for all output channels.
t_out=(*f_limPtr)#t_in
dt_out=sqrt((*f_limPtr)#dt_in^2)
;
;	Hence count rate and associated error for output channels.
z_out=t_out/(w_out#uvec_ng)
dz_out=dt_out/(w_out#uvec_ng)
;
; Free some pointers.
ptr_free,r_limPtr
ptr_free,l_limPtr
ptr_free,d_limPtr
ptr_free,f_limPtr
;
end


;************************************************************************************************
pro drebin_pts,x_in,z_in,dz_in,x_out,z_out,dz_out
;************************************************************************************************
; This procedure rebins 1-dimensional ***POINTS*** data. The data are described by the
; arrays x_in, which is dimensioned n_in, and by z_in and dz_in, both of which are dimensioned
; [n_in,ng], where n_in is the number of channels and ng is the number of groups, i.e. sets of
; data. The grid onto which the data are to be rebinned is defined by the array x_out and the
; results are placed in z_out and dz_out; x_out is dimensioned n_out+1 and z_out and dz_out are
; dimensioned [n_out,ng] where n_out is the number of output channels. It is up to the user of
; this procedure to ensure that the input arrays are properly dimensioned.
;
compile_opt strictarr
;
;
; Generate number of input channels.
n_in=n_elements(x_in)
;
; Generate number, left limit, right limit and width in x of output channels.
n_out=n_elements(x_out)-1
l_out=x_out[0:n_out-1]
r_out=x_out[1:n_out]
w_out=r_out-l_out
;
; Determine the number of sets of input data.
ng=n_elements(z_in)/n_in
;
; Create unit vectors for later use.
uvec_nout=intarr(n_out)+1
uvec_nin=intarr(n_in)+1
uvec_ng=intarr(ng)+1
;
;	The array element (*f_limPtr)[i,j] is 1 if input value j
; contributes to output channel i, otherwise 0, except that
; if input value j lines up with the boundary
; between output channels i and i+1,
; (*f_limPtr)[i,j]=0.5 and (*f_limPtr)[i+1,j]=0.5.
lPtr=ptr_new(l_out#uvec_nin,/no_copy)
rPtr=ptr_new(r_out#uvec_nin,/no_copy)
xPtr=ptr_new(uvec_nout#x_in,/no_copy)
;
l1_Ptr=ptr_new(*lPtr lt *xPtr,/no_copy)
r1_Ptr=ptr_new(*xPtr lt *rPtr,/no_copy)
f1_Ptr=ptr_new(*l1_Ptr and *r1_Ptr,/no_copy)
;
l2_Ptr=ptr_new(*lPtr eq *xPtr,/no_copy)
r2_Ptr=ptr_new(*xPtr eq *rPtr,/no_copy)
f2_Ptr=ptr_new(*l2_Ptr + *r2_Ptr,/no_copy)
;
ff_Ptr=ptr_new(*f1_Ptr + 0.5 * *f2_Ptr,/no_copy)
;
; Free some pointers.
ptr_free,lPtr
ptr_free,rPtr
ptr_free,xPtr
;
ptr_free,l1_Ptr
ptr_free,r1_Ptr
ptr_free,f1_Ptr
;
ptr_free,l2_Ptr
ptr_free,r2_Ptr
ptr_free,f2_Ptr
;
;	Do the actual rebinning.
; Integrated counts and associated errors are computed for all output channels.
t_out=(*ff_Ptr)#z_in
dt_out=sqrt((*ff_Ptr)#dz_in^2)
;
;	Hence count rate and associated error for output channels.
z_out=t_out/(w_out#uvec_ng)
dz_out=dt_out/(w_out#uvec_ng)
;
; Free the remaining pointer.
ptr_free,ff_Ptr
;
end


;************************************************************************************************
function strincrease,x
;************************************************************************************************
; This function returns 1 if the array x is strictly increasing, otherwise 0.
; "Strictly increasing" means that x(i) < x(i+1) for all i.
;
compile_opt strictarr
;
nel=n_elements(x)
not_decreasing=fix(total(x eq x[sort(x)])+0.001) eq nel
unique=n_elements(x[uniq(x)]) eq nel
return,not_decreasing and unique
end


;************************************************************************************************
pro drebin_xp2h,xpts,xbounds
;************************************************************************************************
; Given nx x-points this routine creates nx+1 x-bounds such that
; for i=1,nx-1 xbounds[i]=0.5*(xpnts[i-1]+xpnts[i]),
; xbounds[0]=xbounds[1]-(xpnts[1]-xpnts[0])
; xbounds[nx]=xbounds[nx-1]+(xpnts[nx-1]-xpnts[nx-2])
;
compile_opt strictarr
;
nx = n_elements(xpts)
mid = 0.5*(xpts[0:nx-2]+xpts[1:nx-1])
wid = xpts[1:nx-1] - xpts[0:nx-2]
xbndi = mid[0] - wid[0]
xbndf = mid[nx-2] + wid[nx-2]
xbounds = [xbndi,mid,xbndf]
end


;************************************************************************************************
pro drebin_zp2h,xbounds,zvals,zheights
;************************************************************************************************
; Given nx x-bounds this routine converts values of z (areas) into z-heights.
;
compile_opt strictarr
;
nx=n_elements(xbounds)
dx=xbounds[1:nx-1]-xbounds[0:nx-2]
ng=n_elements(zvals)/(nx-1)
zheights=zvals/(dx#(bytarr(ng)+1))
end


;************************************************************************************************
pro drebin_zh2p,xbounds,zheights,zvals
;************************************************************************************************
; Given nx x-bounds this routine converts z-heights to values of z, i.e. areas.
compile_opt strictarr
;
nx=n_elements(xbounds)
dx=xbounds[1:nx-1]-xbounds[0:nx-2]
ng=n_elements(zheights)/(nx-1)
zvals=zheights*(dx#(bytarr(ng)+1))
end


;************************************************************************************************
pro drebin,x_in,z_in,dz_in,x_out,z_out,dz_out,err=err,emsg=emsg,$
	histogram=hst,points=pts,to_histogram=to_hst,to_points=to_pts
;************************************************************************************************
;+
; NAME:
;	DREBIN
;
; PURPOSE:
;	Perform rebinning.
; The input data are one or more sets of points or one or more sets of histograms.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;	drebin,x_in,z_in,dz_in,x_out,z_out,dz_out,err=err,emsg=emsg
;
; INPUTS:
;	x_in:   Input x values.
;	z_in:   Input z values.
;	dz_in:  Input dz values.
; If the input data are one or more sets of points, x_in contains the x
;		coordinates of the points.
; If the input data are one or more histograms, x_in contains the boundaries
; 	of the histogram channels.
; In general z_in and dz_in are dimensioned (n,ng) where ng is the number of
;		sets of points and n is the number of points or histogram channels per set.
; If the data are points the dimension of x_in is n.
; If the data are histograms the dimension of x_in is n+1.
;	x_out:  Output x values.
; x_out contains the boundaries of the output histogram channels.
;
;	INPUT KEYWORDS:
; points:	If defined and nonzero, the input data are points data.
;	histogram: If defined and nonzero, the input data are histogrammed data.
; to_points: If defined and nonzero, the output data are points data.
;	to_histogram: If defined and nonzero, the output data are histogrammed data.
; The only permitted combinations are /points,/to_points and /histogram,/to_histogram.
;
; OUTPUTS:
;	z_out:   Output z values.
;	dz_out:  Output dz values.
;
;	OUTPUT KEYWORDS:
; err:   An integer specifying an error code, 0 if no error.
;	emsg:  A string containing an error message, "" if no error.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	See pdf file.
;
; MODIFICATION HISTORY:
;	Written, JRDC, April 2002.
; Modified to accept input keywords and perform one of two types of conversion, from points
; to points or from histogram to histogram, JRDC, November 2002.
;
compile_opt strictarr
;
z_out=0
dz_out=0
;
err=0
emsg=""
;
opcode=8*keyword_set(pts)+4*keyword_set(hst)+2*keyword_set(to_pts)+keyword_set(to_hst)
;
; The only allowed values of opcode are 5 (,/histogram,/to_histogram) and
;		10 (,/points,/to_points).
if (opcode ne 5 and opcode ne 10) then begin
	err=-20
	emsg="DREBIN: Forbidden combination of input keywords."
	return
endif
;
; Keywords err and emsg must both appear in the call to this procedure.
if (not arg_present(err) and not arg_present(emsg)) then return
;
if (not arg_present(err) and arg_present(emsg)) then begin
	; (No point in trying to define err)
	emsg="DREBIN: Keyword err is missing from call."
	return
endif
;
if (arg_present(err) and not arg_present(emsg)) then begin
	err=-2
	; (No point in trying to define emsg)
	return
endif
;
; There must be 6 positional parameters.
if (n_params() ne 6) then begin
	err=-1
	emsg="DREBIN: There must be 6 parameters."
	return
endif
;
; x_in must be a vector.
if (size(x_in,/n_dimensions) ne 1) then begin
	err=1
	emsg="DREBIN: x_in must be a vector."
	return
endif
;
; The vector x_in must be strictly increasing.
if (strincrease(x_in) ne 1) then begin
	err=2
	emsg="DREBIN: The vector x_in must be strictly increasing."
	return
endif
;
; x_in must be of type float or double.
if (size(x_in,/type) ne 4 and size(x_in,/type) ne 5) then begin
	err=3
	emsg="DREBIN: x_in must be of type float or double."
	return
endif
;
xdims=size(x_in,/dimensions)
zdims=size(z_in,/dimensions)
dzdims=size(dz_in,/dimensions)
;
; The dimensions of z_in and dz_in must be identical
if (fix(total(zdims ne dzdims)+0.001)) then begin
	err=11
	emsg="DREBIN: The dimensions of z_in and dz_in must be identical."
	return
endif
;
xdims=xdims[0]
zdims=zdims[0]
;
; If opcode=10, first dimension of z_in must equal dimension of x_in.
if (opcode eq 10 and (xdims ne zdims)) then begin
	err=12
	emsg="DREBIN: Input is points so first dimension of z_in must equal dimension of x_in."
	return
endif
;
; If opcode=5, first dimension of z_in must be 1 less than dimension of x_in.
if (opcode eq 5 and (xdims-1 ne zdims)) then begin
	err=13
	emsg="DREBIN: Input is histogram so first dimension of z_in must be 1 less than dimension of x_in."
	return
endif
;
; z_in and dz_in must be of type float or double.
if ((size(z_in,/type) ne 4 and size(z_in,/type) ne 5) $
	or (size(dz_in,/type) ne 4 and size(dz_in,/type) ne 5)) then begin
	err=14
	emsg="DREBIN: z_in and dz_in must be of type float or double."
	return
endif
;
; x_out must be a vector.
if (size(x_out,/n_dimensions) ne 1) then begin
	err=21
	emsg="DREBIN: x_out must be a vector."
	return
endif
;
; The vector x_out must be strictly increasing.
if (strincrease(x_out) ne 1) then begin
	err=22
	emsg="DREBIN: The vector x_out must be strictly increasing."
	return
endif
;
; x_out must be of type float or double.
if (size(x_out,/type) ne 4 and size(x_out,/type) ne 5) then begin
	err=23
	emsg="DREBIN: x_out must be of type float or double."
	return
endif
;
; Points to points.
redefining_z=0
;
if (opcode eq 10) then begin
	drebin_xp2h,x_in,x_inbounds
	drebin_xp2h,x_out,x_outbounds
	if (redefining_z) then begin
		drebin_zp2h,x_inbounds,z_in,z_inheights
		drebin_zp2h,x_inbounds,dz_in,dz_inheights
		drebin_histo,x_inbounds,z_inheights,dz_inheights,x_outbounds,z_outheights,dz_outheights
		drebin_zh2p,x_outbounds,z_outheights,z_out
		drebin_zh2p,x_outbounds,dz_outheights,dz_out
	endif else begin
		drebin_histo,x_inbounds,z_in,dz_in,x_outbounds,z_out,dz_out
	endelse
	return
endif
;
; Histogram to histogram.
if (opcode eq 5) then begin
	drebin_histo,x_in,z_in,dz_in,x_out,z_out,dz_out
	return
endif
;
return
;
end
