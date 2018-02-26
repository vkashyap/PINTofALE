function unkink,y,n,xloc=xloc,hires=hires,nres=nres,$
	onestop=onestop,verbose=verbose, _extra=e
;+
;function	unkink
;	given a curve with kinks in it, returns a curve by
;	sweeping line segments across that smooths out kinks
;
;syntax
;	ys=unkink(y,n,xloc=xloc,/hires,nres=nres,/onestop,$
;	verbose=verbose)
;
;parameters
;	y	[INPUT; required] the ordinate values of the curve to smooth
;		* cannot handle multi-valued functions.  i.e., curve can be
;		  concave upward or downward, but not leftward or rightward
;	n	[INPUT; default=1] the number of points to cast on either
;		side of a point to compute the smoothing line segment
;
;keywords
;	xloc	[INPUT] the abscissa values at which Y are defined
;		* the default is to use lindgen(n_elements(Y))+1
;		* this is useful only if the grid is non-uniform
;		* if specified, but size does not match that of Y,
;		  then it is ignored
;	hires	[INPUT] if set, resamples the grid into a finer
;		resolution (with NRES additional points in each
;		segment), linearly interpolates Y onto this finer
;		grid, and operates on the high-resolution curve
;	nres	[INPUT] the number of additional points to include
;		in each segment of Y(XLOC)
;		* default is 1
;	onestop	[INPUT] if set, the smoothed curve is returned as
;		the value interpolated on to XLOC by a line segment
;		that is symmetric around XLOC
;		* if not set, an average of the interpolated Ys is
;		  computed from all the line segments that pass through
;		  the given XLOC
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	none
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
;example
;	x=findgen(20) & y=abs(x-10) & plot,x,y,psym=-4
;	oplot,x,unkink(y,10,/hires,nres=5),thick=2
;
;history
;	vinay kashyap (sep08)
;-

;	usage
ok='ok' & np=n_params() & ny=n_elements(y) & nn=n_elements(n)
if np eq 0 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y is undefined' else $
  if ny lt 3 then ok='Y is not kinky'
   if nn gt 1 and nn ne ny then ok='N and Y are incompatible'
if ok ne 'ok' then begin
  print,'Usage: ys=unkink(y,n,xloc=xloc,/hires,nres=nres,$'
  print,'       /onestop,verbose=verbose)'
  print,'  given a curve with kinks in it, generates a smooth curve by'
  print,'  sweeping line segments across it that smooths out the kinks'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
yy=[y[*]]	;one long happy array
mm=lonarr(ny)+1
if nn gt 0 then mm[*]=long(n[0])>1	;constant smoothing
if nn eq ny then mm=long(n)>1	;variable smoothing

;	check keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
xx=lindgen(ny)+1
nx=n_elements(xloc) & if nx eq ny then xx[*]=xloc
if vv gt 10 then plot,xx,yy,psym=-4
;
npt=ny
if keyword_set(hires) then begin
  mres=1 & if keyword_set(nres) then mres=long(nres[0])>1
  npt=ny+(ny-1L)*mres
  ii=lindgen(ny) & iih=interpol(ii,npt)
  xxh=interpol(xx,ii,iih)
  yyh=interpol(yy,xx,xxh)
  mmh=long(interpol(mm,xx,xxh))>1
  xx=xxh & yy=yyh & mm=mmh
  ;	to see that that sequence of interpols does work, try
  ;	x=lindgen(5)^2 & print,x
  ;	ii=lindgen(5) & iii=interpol(ii,5+(5-1)*2) & print,iii
  ;	xx=interpol(x,ii,iii) & print,xx
  ;	plot,ii,x,psym=-1 & oplot,iii,xx,psym=-4,col=2
  ;	oplot,iii,iii^2,col=3 ;<-- that it does not lie on xx is a feature, not a bug
endif

ys=0.*yy & nchk=lonarr(npt)
for i=0L,npt-1L do begin
  nx=mm[i]
  i0=(i-nx)>0L
  i1=(i+nx)<(npt-1L)
  xseg=xx[i0:i1] & yseg=yy[i0:i1]
  ix=(i-i0)
  x0=xx[i0] & x1=xx[i1] & y0=yy[i0] & y1=yy[i1]
  if x1 ne x0 then begin
    slope=(y1-y0)/(x1-x0) & intercept=(x1*y0-x0*y1)/(x1-x0)
    yseg=slope*xseg+intercept
  endif
  if vv gt 100 then oplot,xseg,yseg,color=2
  if keyword_set(onestop) then ys[i]=yseg[ix] else begin
    nchk[i0:i1]=nchk[i0:i1]+1L
    ys[i0:i1]=ys[i0:i1]+yseg
  endelse
endfor
if not keyword_set(onestop) then ys=ys/float(nchk)

;	downsample if required
if keyword_set(hires) then begin
  iil=lindgen(ny)*(mres+1) & ys=ys[iil] & xx=xx[iil]
endif
if vv gt 10 then oplot,xx,ys,thick=3,color=3

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,ys
end
