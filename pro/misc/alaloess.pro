function alaloess,yy,xx,hwidth=hwidth,ysig=ysig,ndeg=ndeg,$
	yerr=yerr,verbose=verbose, _extra=e
;+
;function	alaloess
;	computes a loess/lowess type smoothed function for the
;	given set of input points by fitting polynomials over
;	the neighborhood of each point and using the local value
;	of the polynomial fit as the smoothed curve
;
;	another one in a hoary line of PINTofALE smoothing routines,
;	joining SMOOTHIE, NOISMOOTH, VARSMOOTH, VOORSMOOTH, REGROUP,
;	SPLAC, CONV_RMF, SCRMF, UNKINK, LINEREM, HAARLINE, etc.
;
;syntax
;	yl=alaloess(yy,xx,hwidth=hwidth,ysig=ysig,ndeg=ndeg,$
;	yerr=yerr,verbose=vv)
;
;parameters
;	yy	[INPUT; required] the points that need to be smoothed
;	xx	[INPUT] the abscissae for which YY are defined
;		* if size eq size of YY or +1, taken to be the midpoints
;		  or the grid
;		* if insensible, taken to be array indices, LINDGEN(N(YY))
;
;keywords
;	hwidth	[INPUT] half range over which local polynomial fits are
;		to be performed
;		* if XX are array indices, HWIDTH is forced to be >0
;		* if XX is a predefined grid, the default value is the
;		  minimum of delta(XX)
;	ysig	[INPUT] error on YY
;		* if size matches that of YY, then maps 1-to-1
;		* if scalar, taken to be
;		    <0 : absolute error
;		    >0 and <1 : fractional error
;		    >1 : percentage error
;		    =0 : ignored in call to POLY_FIT()
;		  on YY
;		* otherwise, assumed to be sqrt(abs(YY)+0.75)+1
;	ndeg	[INPUT] degree of the polynomial to be fit
;		* default is 3
;	yerr	[OUTPUT] 1-sigma error at the smooth location from
;		each fit.
;		* WARNING: this is to be taken only as an illustrative
;		  number and not used for anything quantitative, because
;		  the errors in adjacent bins are not independent
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run alaloess	;(requires PEASECOLR)
;
;these are all the smoothing tools in PINTofALE
;	ALALOESS() makes a loess curve
;	CLTSMOOTH() removes -ves from background subtracted spectra or light curves
;	CONV_RMF convolves with specified RMF
;	HAARTRAN() smoothes by filtering on Haar wavelet coefficients
;	LINEREM() removes lines from spectrum by iterative S/N filtering
;	NOISMOOTH() does boxcar accumulation a la Ebeling's asmooth
;	REGROUP() accumulates from one end
;	SCRMF does fast convolution using optimized RMF
;	SMOOTHIE does peak-descent accumulation to build up S/N
;	SPLAC computes a segmented piecewise linear approximations
;	UNKINK() removes sharp kinks from a curve
;	VARSMOOTH() does local smoothing at specified scales
;	VOORSMOOTH() does fast local smoothing at specified scales
;
;history
;	vinay kashyap (MarMMX)
;	some bug fixes, and hardening to poly_fit failures -- where it fails
;	  now simply keep old value (VK; MayMMX)
;	updated to be more robust to NaNs (VK; AprMMV)
;-

;	usage
np=n_params() & ok='ok' & ny=n_elements(yy) & nx=n_elements(xx)
if np eq 0 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y(X) is undefined' else $
  if ny eq 1 then ok='Y(X) must be a vector'
if ok ne 'ok' then begin
  print,'Usage: yl=alaloess(y,x,hwidth=hwidth,ysig=ysig,ndeg=ndeg,$'
  print,'       yerr=yerr,verbose=vv)'
  print,'  returns lowess like smooth curve based on local polynomial fit'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	parameters
yz=0.0+[yy[*]] & xok='units'
case nx of
  ny: xz=xx
  ny+1: xz=0.5*(xx[1:*]+xx)
  else: begin
    xz=lindgen(ny)
    xok='indices'
  endelse
endcase

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
mpoly=3 & if keyword_set(ndeg) then mpoly=abs(fix(ndeg[0]))
;
nye=n_elements(ysig) & sigy=sqrt(abs(yz)+0.75)+1.
if nye eq 1 then begin
  if ysig[0] lt 0 then sigy[*]=abs(ysig[0])
  if ysig[0] gt 0 and ysig[0] lt 1 then sigy=ysig[0]*yz
  if ysig[0] ge 1 then sigy=(ysig[0]/100.)*yz
  if ysig[0] eq 0 then sigy[*]=0
endif
if nye eq ny then sigy=ysig+0.0
;
hrange=min(abs(xz[1:*]-xz))
if keyword_set(hwidth) then hrange=hrange>hwidth[0]

;	step through the array and get the smoothed values
yl=yz & yerr=yl
if 2*hrange le mpoly then begin
  message,'fitting window too small for polynomial',/informational
  message,'returning without any smoothing',/informational
  print,mpoly,hrange
  yerr=sigy
  return,yl
endif
for i=0L,ny-1L do begin
  i0=(i-hrange)>0
  i1=(i+hrange)<(ny-1L)
  if xok eq 'units' then begin
    oo=where(abs(xz-xz[i]) le hrange,moo)
    i0=oo[0] & i1=oo[moo-1L]
  endif
  j=i-i0
  if i1-i0 gt mpoly then begin
    xin=xz[i0:i1] & yin=yz[i0:i1] & yein=sigy[i0:i1]
    if sigy[i] gt 0 then begin
      tmp=poly_fit(xin,yin,mpoly,measure_errors=yein,yfit=yfit,$
      yband=yband,/double,status=ist)
      if ist gt 0 then tmp=poly_fit(xin,yin,mpoly,yfit=yfit,yband=yband,/double,status=ist)
    endif else begin
      tmp=poly_fit(xin,yin,mpoly,yfit=yfit,yband=yband,/double,status=ist)
      if ist gt 0 then tmp=poly_fit(xin,yin,mpoly,/double,status=ist)
    endelse
    if ist eq 0 then begin
      yl[i]=yfit[j] & yerr[i]=yband[j]
    endif
  endif
endfor

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,yl
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example calling sequence

yl=alaloess()

if not keyword_set(ndeg) then ndeg=3
if not keyword_set(hwidth) then hwidth=8
if n_elements(verbose) eq 0 then verbose=1
print,'HWIDTH='+strtrim(hwidth,2)
print,'NDEG='+strtrim(ndeg,2)
print,'VERBOSE='+strtrim(VERBOSE,2)

y=50*abs(randomn(seed))*mk_gauss(findgen(256),128,10,1)+$
  30*abs(randomn(seed))*mk_gauss(findgen(256),165,15,1)
yoff=40*randomu(seed)+10.
yr=long(0*y)
for i=0,256-1 do yr[i]=randomu(seed,poisson=y[i]+yoff)

yl=alaloess(yr,hwidth=hwidth,ndeg=ndeg,yerr=yerr,verbose=verbose)

peasecolr & loadct,3 & peasecolr
plot,yr,psym=1 & oplot,y+yoff,col=1
oplot,yl,col=2 & oplot,yl+yerr,col=2,line=2 & oplot,yl-yerr,col=2,line=2
print,'' & print,'blue: actual; red: smoothed' & print,''

end
