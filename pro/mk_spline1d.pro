function mk_spline1d,x,aa,pder,xloc=xloc,useinit=useinit,$
	tension=tension,verbose=verbose, _extra=e
;+
;function	mk_spline1d
;	returns a spline passing through the points [XLOC,AA]
;
;syntax
;	sy=mk_spline1d(x,aa,pder,xloc=xloc,/useinit,tension=tension,$
;	verbose=verbose)
;
;parameters
;	x	[INPUT array; required] where SY(X) must be computed
;	aa	[INPUT; required] array of y-values through which the
;		spline curve passes
;		* if only one point is given, the array is padded so
;		  that the curve goes to 0 at the endpoints of X
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if supplied in calling sequence.
;
;keywords
;	xloc	[INPUT] x-locations which correspond to AA
;		* if given, size must match AA, or else an error is returned
;		* if not given, a regular grid running from min(X) to max(X)
;		  and containing n_elements(AA) points are used
;	useinit	[INPUT] if set, uses the older IDL functions SPL_INIT()
;		and SPL_INTERP() to construct the spline.  the default
;		is to call SPLINE()
;	tension	[INPUT] the amount of tension to be applied in the call
;		to SPLINE()
;		* identical to the parameter SIGMA of SPINE() -- if close
;		  to 0, effectively cubic spline, otherwise like polynomial
;		* if not given, assumed to be 1
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		* SPLINE: DOUBLE
;		* SPL_INIT: YP0, YPN_1, DOUBLE
;
;usage summary
;	* call as a function
;
;subroutines
;	SPLINE
;	SPL_INIT
;	SPL_INTERP
;
;history
;	vinay kashyap (Sep08)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & na=n_elements(aa)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if na eq 0 then ok='Spline coefficients are not defined'
if ok ne 'ok' then begin
  print,'Usage: sy=mk_spline1d(x,aa,pder,xloc=xloc,/useinit,tension=tension,$'
  print,'       verbose=verbose, /double,yp0=yp0,ypn_1=ypn_1)'
  print,'  returns a spline passing through [XLOC,AA]'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
xmin=min(x,max=xmax)
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
if keyword_set(useinit) then callspline='SPL_INIT()/SPL_INTERP()' else callspline='SPLINE()'
nxloc=n_elements(xloc)
if nxloc eq 0 then begin	;(XLOC is undefined
  xloc=(xmax-xmin)/(na-1)
endif else begin		;XLOC undefined)(XLOC is given
  if nxloc ne na then begin
    message,'XLOC and AA are inconsistent; quitting',/informational
    return,-1L
  endif
endelse
xx=xloc & yy=aa & os=sort(xx) & xx=xx[os] & yy=yy[os]
if na eq 1 then begin
  yy=[0.,yy,0.]
  xx=[xmin,xx,xmax]
endif

if vv gt 20 then message,'Generating spline with '+callspline,/informational

;	compute function
if callspline eq 'SPLINE()' then begin
  tense=1.0 & if keyword_set(tension) then tense=tension[0]
  ix=lindgen(nx) & os=sort(x) & z=x[os] & iz=ix[os]
  sy=spline(xx,yy,z,tense, _extra=e)
  oss=sort(iz) & sy=sy[oss]
endif else begin
  os=sort(xx) & y2=spl_init(xx[os],yy[os], _extra=e)
  sy=spl_interp(xx[os],yy[os],sy,x, _extra=e)
endelse

if vv gt 50 then plot,x,sy,xtitle='X',ytitle='spline'

;	compute partial derivatives
if np ge 3 then begin
  pder=fltarr(nx,na)+1.
  ;	y(x_i) = y_i, so the partials wrt the parameter y_i, (d/dy_i)y_i == 1
endif

return,sy
end
