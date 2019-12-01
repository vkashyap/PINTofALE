function haartran,y,scales,ss,ysig=ysig,thrsig=thrsig,thrloc=thrloc,$
	thrpoi=thrpoi,wy=wy,bwy=bwy,verbose=verbose, _extra=e
;+
;function	haartran
;	Run the input 1D array through a sequence of Haar wavelets,
;	and return the filtered wavelet transform coefficients.
;
;	convolves Y with Haar wavelet at different scales; at each
;	scale, a threshold is applied, either a global one based on
;	the histogram of the wavelet coefficients or a local one based
;	on the propagated error or a local one based on the local
;	background and the probability of fluctuations; the filtered
;	coefficients are returned in a 2D array, the first row of
;	which contains the "reconstructed/denoised" Y.
;
;	a Haar wavelet is simply a boxcar with a moat.
;
;syntax
;	hy=haartran(y,scales,ss,ysig=ysig,thrsig=thrsig,/thrloc,$
;	thrpoi=thrpoi,wy=wy,bwy=bwy,verbose=v)
;
;parameters
;	y	[INPUT; required] the input 1D function, assumed to
;		be on a regular grid
;	scales	[INPUT] scales at which to construct the Haar wavelets
;		* must be integer
;		* if not specified, then taken to be integers going in
;		  powers of 2 from 1..(N(Y)/6)
;		  (6 because the kernel width is actually 3*scale, and
;		  you want to go only so far as to fill half of Y, it
;		  would be useless beyond that)
;		* if array, taken to be the actual wavelet scales given
;		  in bin or pixel units
;		* if scalar, taken to be the _number_ of scales to be
;		  considered, 1..2^(SCALES)
;		* in all cases, maximum scale size must be < N(Y)/6
;		  (scales higher than that are simply ignored)
;	ss	[OUTPUT] the scales at which the wavelet is applied
;
;keywords
;	ysig	[INPUT] error on Y
;		* if not given, and Y are integers and max(Y)>1, then assumed
;		  to be sqrt(abs(Y)+0.75)+1, else 0
;		* if given and is a scalar, then
;		  -- if <0, the abs value is assumed to be a constant error
;		  -- if >0 and <1, assumed to be a constant fractional error
;		  -- if >1 and <100, assumed to be a constant percentage error
;		  -- if >100, assumed to be a constant additive error
;		* if an array and size does not match N(Y), then ignored
;		* used only if THRLOC is set
;	thrsig	[INPUT; default=1] threshold sigma to use for filtering
;		* if THRLOC is set, then filter out all wavelet coefficients
;		  with values < THRSIG*local_error
;		* if THRLOC is not set (default), filter out all wavelet
;		  coefficients with values<(mean_of_dist+THRSIG*stddev_of_dist)
;		* note that this does _not_ imply a probablity of false
;		  detection or anything of that ilk.  this is just a number,
;		  use it as such.
;	thrloc	[INPUT] if set, filters wavelet coefficients based on a
;		locally computed error
;	thrpoi	[INPUT] if set, first estimates a local background at each
;		bin at the given scale, and then sets the local threshold
;		such that at most THRPOI bins are found to be false lines
;		over the range of Y (i.e., sets a probability threshold
;		locally at THRPOI/N(Y))
;		* if set, overrides THRLOC and ignores THRSIG
;	wy	[OUTPUT] _all_ the computed coefficients, without any filtering,
;		for those that may wish to apply a more rigorous filtering,
;		as a 2D array of size (N(SCALES),N(Y))
;	bwy	[OUTPUT] same as WY, but computed using only the moat of the
;		Haar as the kernel.  So this is effectively the local background
;		at each bin.
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
;
;subroutines
;	KILROY
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
;	vinay kashyap (Dec'02)
;	now centers the kernel; added keyword BWY (Nov'14)
;	*really* added keyword BWY (Nov'19)
;-

;	usage
ok='ok' & np=n_params() & ny=n_elements(y)
if np eq 0 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y is undefined' else $
  if ny lt 3 then ok='Y must have at least 3 elements'
if ok ne 'ok' then begin
  print,'Usage: hy=haartran(y,scales,ss,ysig=ysig,thrsig=thrsig,/thrloc,wy=wy,verbose=v)'
  print,'  compute Haar wavelet coefficients, and return the filtered coefficients'
  print,'  as a 2D array of size (N(SCALES)+1,N(Y))'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check SCALES
ns=n_elements(scales)
smax=long(alog(ny/6)/alog(2))	;the "6" because you want the kernel,
				;which is 3*scale, to at best fit in
				;within half the size of the spectrum.
if ns eq 0 then begin			;(NS=0
  ns=smax+1L & ss=2L^(lindgen(ns))
endif else begin			;NS=0)(NS>0
  ssz=size(scales)
  if ssz[0] eq 0 then begin				;(NS=1
    ss=2L^(lindgen(long(scales[0]>1)))
    oi=where(ss lt 2L^(smax),moi)
    if moi gt 0 then ss=ss[oi] else message,'BUG!'
    ns=moi
  endif else begin				;NS=1)(NS>1
    ssz=size(scales,/tname)
    if ssz ne 'INT' and ssz ne 'LONG' then begin
      message,'SCALES must be an integer array. Returning as is.',/informational
      return,y
    endif
    oi=where(scales lt 2L^(smax),moi)
    if moi gt 0 then ss=scales[oi] else begin
      message,'SCALES must be < N(Y)/2.  Returning as is.',/informational
      return,y
    endelse
    ns=moi
  endelse					;NS>1)
endelse					;NS>0)

;	check keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
sy=size(y,/tname)
if sy eq 'INT' or sy eq 'LONG' and max(y) gt 1 then $
	sigy=sqrt(abs(y)+0.75)+1. else sigy=0.*y
nye=n_elements(ysig)
if nye gt 0 then begin
  if nye eq 1 then begin
    if ysig[0] lt 0 then sigy[*]=abs(ysig[0]) else $
     if ysig[0] lt 1 then sigy[*]=y*ysig[0] else $
      if ysig[0] lt 100 then sigy[*]=y*ysig[0]/100. else $
       sigy[*]=ysig[0]
  endif else begin
    if nye eq ny then sigy=ysig else $
      message,'YSIG size does not match Y; ignoring',/informational
  endelse
endif
;
sigthr=1. & if keyword_set(thrsig) then sigthr=float(thrsig[0])

;	define the outputs
wy=fltarr(ns+1,ny)+0*y[0]		;the "0*y[0]" in case output needs to be in double precision
bwy=fltarr(ns+1,ny)+0*y[0]
hy=fltarr(ns+1,ny)+0*y[0]

;	for each scale, get the wavelet coefficients and filter
for i=0L,ns-1L do begin			;{step through the scales
  scale=ss[i]
  kern=[-0.5+fltarr(scale), 1.0+fltarr(scale), -0.5+fltarr(scale)]/(2.*scale)
  bkern=[0.5+fltarr(scale), fltarr(scale), 0.5+fltarr(scale)]/(2.*scale)
  if vv gt 1 then kilroy,dot=strtrim(i,2)+'.'
  if vv gt 3 then kilroy,dot=strtrim(ss[i],2)+':'
  tmp=convol(float(y),bkern,/edge_truncate,center=1) & bwy[i+1,*]=tmp
  tmp=convol(float(y),kern,/edge_truncate,center=1) & wy[i+1,*]=tmp
  if keyword_set(thrpoi) then begin	;(local probability threshold
    message,'sorry, THRPOI not implemented yet',/informational
    return,-1L
  endif else begin			;THRPOI)(THRPOI not set
    if keyword_set(thrloc) then begin	;(local thresholds
      ;tmpe=sqrt(convol(sigy^2,abs(kern),/edge_truncate,center=1))
      tmpe=sqrt(convol(sigy^2,(-kern>0),/edge_truncate,center=1))
      tmpbe=sqrt(convol(sigy^2,(bkern>0),/edge_truncate,center=1))
      oo=where(tmp gt sigthr*tmpe,moo)
    endif else begin			;THRLOC=1)(global threshold
      jnk=moment(tmp)
      oo=where(tmp gt jnk[0]+sigthr*sqrt(jnk[1]),moo)
    endelse				;THRLOC=0)
  endelse				;THRLOC/THRSIG)
  if moo gt 0 then hy[i+1L,oo]=tmp[oo]
endfor					;I=0,NS-1}
for i=0L,ny-1L do hy[0,i]=total(hy[1L:ns,i])
for i=0L,ny-1L do wy[0,i]=total(wy[1L:ns,i])
for i=0L,ny-1L do bwy[0,i]=total(bwy[1L:ns,i])

return,hy
end
