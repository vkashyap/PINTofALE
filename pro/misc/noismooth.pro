function noismooth,y,ysig=ysig,snrthr=snrthr,snrout=snrout,xsiz=xsiz,$
	sdepth=sdepth,verbose=verbose, _extra=e
;+
;function	noismooth
;	return the input function boxcar smoothed at varying
;	scales such that the value at any point is determined
;	at a sufficiently high S/N.
;
;	beware that this will NOT conserve flux, and the adjacent
;	points (or errors) of the smoothed function will no longer
;	be statistically independent.
;
;syntax
;	ys=noismooth(y,ysig=ysig,snrthr=thresh,snrout=snr,xsiz=xsiz,$
;	sdepth=sdepth,verbose=v)
;
;parameters
;	y	[INPUT; required] array that must be smoothed
;		* note that integer count arrays are returned as
;		  integers, which may not necessarily be what is
;		  optimal.  to return floats instead, simply input
;		  (y*1.0)
;
;keywords
;	ysig	[INPUT] errors on Y -- S/N on smoothed function will
;		be calculated by square-adding
;		* if not given, Gehrel's approximation to Poisson,
;		  sqrt(abs(Y)+0.75)+1 will be assumed
;		* if single-element, then assumed to be
;		  -- additive constant if <0 (i.e., abs(YSIG[0]), or else
;		  -- fractional if <1, or else
;		  -- percentage if <100, and
;		  -- additive constant otherwise
;	snrthr	[INPUT] S/N threshold -- try to make every point lie
;		above this S/N
;		* default is 3
;	snrout	[OUTPUT] final S/N at each point
;	xsiz	[OUTPUT] final smoothing scale at each point
;	sdepth	[INPUT] maximum smoothing scale to consider
;		* default is n_elements(Y)/2-1, which is also the
;		  hardcoded maximum
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	akin to Harald Ebeling's adaptive smoothing of 2-D images
;	(encoded in CIAO CSMOOTH), this program boxcar smooths the
;	input 1-D array over a number of scales and at each bin
;	keeps the result from the smallest scale which puts it
;	over the specified S/N threshold.
;
;examples
;	peasecolr & y=randomu(seed,1000,poisson=1) & y[400:499]=randomu(seed,100,poisson=5)
;	plot,y & oplot,smooth(y,5),col=2,thick=3 & oplot,smooth(y,100),col=29,thick=3
;	ys=noismooth(y,ysig=ysig,snrthr=10,snrout=snr,xsiz=xsiz)
;	plot,y & oplot,ys,col=2,thick=3 & oplot,snr,col=3
;	ys=noismooth(y,ysig=ysig,snrthr=3,snrout=snr,xsiz=xsiz,verbose=10)
;	oplot,ys,col=25,thick=3 & oplot,snr,col=35
;	ys=noismooth(y,ysig=ysig,snrthr=13,snrout=snr,xsiz=xsiz)
;	oplot,ys,col=29,thick=3 & oplot,snr,col=39
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
;	vinay kashyap (OctMM)
;	modified YSIG behavior to handle YSIG[0]>1 (VK; Nov'02)
;-

;	usage
ok='ok' & np=n_params() & ny=n_elements(y)
if np eq 0 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y is undefined' else $
  if ny lt 3 then ok='Y is too small to smooth!'
if ok ne 'ok' then begin
  print,'Usage: ys=noismooth(y,ysig=ysig,snrthr=thresh,snrout=snr,$'
  print,'          xsiz=xsiz,sdepth=sdepth,verbose=v)'
  print,'  return noise-based smoothed function'
  if np ne 0 then message,ok,/info
  if ny gt 0 then return,y else return,-1L
endif

;	inputs
v=0 & if keyword_set(verbose) then v=long(verbose[0]) > 1
;
nye=n_elements(ysig) & sigy=sqrt(abs(Y)+0.75)+1.
if nye eq 0 and v gt 0 then message,$
	'assuming Gehrels approx to Poisson errors',/info
if nye eq 1 then begin
  if ysig[0] le 0 then sigy[*]=abs(ysig[0]) else $
   if ysig[0] lt 1 then sigy[*]=y*ysig[0] else $
    if ysig[0] lt 100 then sigy[*]=y*ysig[0]/100. else $
     sigy[*]=ysig[0]
  if v gt 0 then begin
    ok='ok'
    if ysig[0] le 0 then ok='assuming constant additive errors' else $
     if ysig[0] lt 1 then ok='assuming fractional errors' else $
      if ysig[0] lt 100 then ok='assuming fractional percentage errors' else $
       ok='assuming constant additive errors'
    if ok ne 'ok' and v gt 1 then message,ok,/informational
  endif
endif
;
thresh=3.0 & if keyword_set(snrthr) then thresh=float(snrthr[0])
;
sclmax=ny/2L-1L & if keyword_set(sdepth) then sclmax=long(sdepth[0]) < sclmax

;	outputs
snrout=fltarr(ny)-1. & xsiz=lonarr(ny) & ys=y

;	smooth at increasingly higher bin sizes
scale=0L
os=where(sigy gt 0,mos) & snrout[os]=y[os]/sigy[os]
oo=where(snrout ge 0 and snrout lt thresh,moo)
while moo gt 0 do begin			;{continue smoothing
  scale=scale+1L & ss=2L*scale+1L
  if v gt 0 then kilroy,dot=strtrim(scale,2)+':'
  tmp=smooth(y,ss,/edge_truncate,/nan)
  tmpe=sqrt(smooth((sigy>0)^2,ss,/edge_truncate,/nan)/ss)
  tmps=snrout < thresh
  ok=where(tmpe gt 0 and tmps lt thresh,mok) & tmps[ok]=tmp[ok]/tmpe[ok]
  ys[ok]=tmp[ok] & snrout[ok]=tmps[ok] & xsiz[ok]=ss
  oo=where(snrout ge 0 and snrout lt thresh,moo)
  if v gt 2 then kilroy,dot='['+strtrim(moo,2)+']'
  ;
  ;other stopping rules
  if scale ge sclmax then moo=0L
endwhile				;S/N < THRESH}

return,ys
end
