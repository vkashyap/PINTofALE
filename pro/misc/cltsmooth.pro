function cltsmooth,y,syerr=syerr,swidth=swidth,verbose=verbose, _extra=e
;+
;function	cltsmooth
;	An adaptive smoothing technique that takes advantage
;	of the central limit theorem.  Applicable to background
;	subtracted light curves or spectra that need to have
;	the negative values removed without touching anything
;	that is at large S/N.
;
;	First picks out all the -ves, assumes these are normally
;	distributed, and collects all +ves that are of the same
;	magnitude.  Then computes mean and stderr at all selected
;	points, then at progressively smaller scales until at most
;	32 bins widths, and resets the local mean and stderr iff
;	it is significantly different from that computed at larger
;	scale.
;
;syntax
;	sy=cltsmooth(y,syerr=syerr,swidth=swidth,verbose=verbose)
;
;parameters
;	y	[INPUT; required] the spectrum or light curve
;
;keywords
;	syerr	[OUTPUT] the errors computed at each bin
;		* note that these will not be independent
;	swidth	[OUTPUT] the widths over which the SY and
;		SYERR are computed
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
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
;	vinay kashyap (2015feb)
;	guard against numerical overflow by rescaling by max(Y) (VK; 2016oct)
;-

;	usage
ok='ok' & np=n_params() & ny=n_elements(y)
if np lt 1 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y is not defined' else $
  if ny eq 1 then ok='Y must be a vector' else $
   if ny lt 8 then ok='Y is not large enough'
if ok ne 'ok' then begin
  print,'Usage: sy=cltsmooth(y,syerr=syerr,swidth=swidth,verbose=verbose)'
  print,'  helps to eliminate -ve excursions from background-subtracted spectra or light curves'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
yy=1.0*y	;make sure it is at least float
ynorm=max(yy,/nan) & yy=yy/ynorm

;	first figure out how the -ves are distributed
om=where(yy lt 0,mom)
if mom eq 0 then begin
  if vv ge 1 then message,'all bins are +ve, nothing to do',/informational
  return,y
endif
;
xm=yy[om] & xx=[xm,-xm] & sdx=stddev(xx)

;	now pick out all adjacent bins, and bins adjacent to them
oset=[om-2,om-1,om,om+1,om+2] & oset=oset[uniq(oset,sort(oset))]
yset=yy[oset] & nset=n_elements(oset)

;	ignore any bins where the fluctuation is > 10*sdx
ok=where(abs(yset) le 10.*sdx,mok)
if mok eq 0 then message,'!BUG'
yset=yset[ok] & oset=oset[ok] & nset=mok

;	top-level scale
ymean=mean(yset) & ysdev=stddev(yset)/sqrt(mok)

;	define outputs
sy=yy & sy[oset]=ymean & syerr=0.*yy & syerr[oset]=ysdev & swidth=lonarr(ny)-1L & swidth[oset]=mok
if mok lt 4 then begin
  if vv gt 5 then message,'all done already',/informational
  sy=sy*ynorm & syerr=syerr*ynorm	;put back normalization
  return,sy
endif

;	cascade down the scales
nscale=mok & go_on=1
while go_on do begin
  nscale=nscale/2L & kern=fltarr(nscale)+1./float(nscale)
  ysm=convol(yset,kern) & ysm2=convol(yset^2,kern) & ysd=sqrt(ysm2-ysm^2)
  dsnr=fltarr(nset) & dsnr=(sy[oset]-ysm)/sqrt(ysd^2+syerr[oset]^2)
  oho=where(abs(dsnr) gt 1,moho)
  if moho gt 0 then begin
    ihi=oset[oho]
    sy[ihi]=ysm[oho] & syerr[ihi]=ysm2[oho] & swidth[ihi]=nscale
  endif else go_on=0
  if nscale lt 4 then go_on=0
  !p.multi=[0,2,1]
  plot,yy & oplot,sy,col=2
  plot,oset,dsnr,psym=-1
  !p.multi=0
  stop,nscale,moho,go_on
endwhile

;	put back the normalization
sy=sy*ynorm & syerr=syerr*ynorm

if vv gt 1000 then stop,'HALTing for debugging; type .CON to continue'

return,sy
end
