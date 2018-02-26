function voorsmooth,y,s,x,type=type,usekern=usekern,eps=eps,$
	verbose=verbose, _extra=e
;+
;function	voorsmooth
;	returns an array smoothed by convolving with appropriate function
;
;	this is basically a front-end to IDL's CONVOL function, and is
;	a faster alternative to VARSMOOTH (a VARSMOOTH that goes vrooom,
;	if you will), smooths with gaussians, lorentzians, beta-profiles,
;	asymmetric beta-profiles, etc.
;
;syntax
;	sy=voorsmooth(y,s,x,type=type,usekern=usekern,eps=eps,verbose=v,$
;	/fwhm,betap=betap,angle=angle,missing=missing,/perbin)
;	
;parameters
;	y	[INPUT; required] 1D array to be smoothed
;	s	[INPUT; required] smoothing scale.  if an array, only
;		the first element is used.
;		* assumed to be in pixel units unless X is given
;		* the meaning of the scale depends on the function
;		  being used -- sigma for TYPE='Gauss' (unless FWHM
;		  is set), half-width for TYPE='Lorentz', etc.
;		* >>ignored<< if USEKERN is defined, but set to some
;		  non-zero value anyway
;	x	[INPUT; optional] values at which Y=Y(X) is defined
;		* size must match Y
;		  OR must be bin boundaries flanking Y
;		  otherwise X is ignored
;		* if given, S is assumed to be in the same units as X
;		* if dX is non-uniform, Y and X are rebinned to the
;		  finest possible grid and THEN convolved.
;
;keywords
;	type	[INPUT; default='Gaussian'] define the 3-parameter family
;		function to be used to smooth
;		* may be anything accepted by X3MODEL/MK_3MODEL, e.g.,
;		  "Gaussian", "Lorentizan", "Beta=<value>",
;		  "Slant=<angle>,<beta>", etc.
;		* >>ignored<< if USEKERN is defined
;	usekern	[I/O] if defined legally on input, will ignore parameter
;		S and keyword TYPE and use the defined kernel to smooth
;		the input array Y.
;	eps	[INPUT] a small number, 1e-6 by default
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] use this to pass defined keywords to subroutines
;		-- MK_GAUSS: FWHM, MISSING
;		-- MK_LORENTZ: BETAP, MISSING
;		-- MK_SLANT: ANGLE, BETAP, MISSING
;		-- REBINW: PERBIN
;
;restrictions
;	calls subroutines:
;		X3MODEL
;		MK_3MODEL
;		MK_GAUSS
;		MK_LORENTZ
;		MK_SLANT
;		REBINW
;		FINDEX
;
;example
;	make a "spectrum":
;	y=fltarr(500) & y[250]=2.
;	smooth with Gaussian of sigma=10 pix
;	ys=voorsmooth(y,10.,type='Gauss',usekern=usekern)
;	make a new kernel:
;	kern=0.5*(shift(usekern,-30)+shift(usekern,30))
;	smooth again:
;	ys2=voorsmooth(y,10.,type='Gauss',usekern=kern)
;	display:
;	plot,y & peasecolr & oplot,ys,color=2 & oplot,ys2,color=3
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
;	vinay kashyap (2001Jul; derived from varsmooth, but really
;	  a precursor, hence the [Noordwijk influenced] name)
;	bug correction oo v/s o0; now reverse kernel before convol (VK; Oct01)
;	allowed X to be bin boundaries (VK/LL; Sep02)
;	accounted for 1-pix shift due to centering of kernel (VK; May03)
;-

;	usage
ok='ok' & np=n_params() & ny=n_elements(y) & ns=n_elements(s)
nx=n_elements(x)
if np lt 2 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y is undefined' else $
  if ns eq 0 then ok='smoothing scale is undefined' else $
   if ny eq 1 then ok='Y has 1 element; nothing to smooth'
if ok ne 'ok' then begin
  print,'Usage: sy=voorsmooth(y,s,x,type=type,usekern=kern,eps=eps,verbose=v,$'
  print,'       /fwhm,betap=betap,angle=angle,missing=missing,/perbin)'
  print,'  smooth Y(X) at scale S using a defined function'
  if np ne 0 then message,ok,/info
  if ny eq 1 then return,Y else return,-1L
endif

;	keywords
ttyp='Gaussian'	;the default
if keyword_set(type) then begin
  szt=size(type) & nszt=n_elements(szt)
  if szt[nszt-2] eq 7 then ttyp=type[0]
endif
vv=0	;the default
if n_elements(verbose) eq 0 then setsysval,'VERBOSE',vv,/getval else $
 vv=long(verbose[0])
;
nk=n_elements(usekern)
;
if n_elements(eps) eq 0 then eps=1e-6

;	recast variables
yy=[y[*]] & ss=s[0]
if ss eq 0 then begin
  message,'Cannot smooth at scale 0; doing nothing',/info
  return,y
endif

;	if X is given
if nx eq ny or nx eq ny+1L then begin
  if nx eq ny+1L then begin
    xx=0.5*(x[1:*]+x) & dx=x[1:*]-x
  endif else begin
    xx=[x[*]] & dx=xx[1:*]-xx
  endelse
  o0=where(dx ne 0,mo0)
  if mo0 eq nx then begin
    message,'grid seems to be a point; returning w/o doing anything',/info
    return,y
  endif
  udx=dx[uniq(dx,sort(dx))] & nudx=n_elements(udx)
  if nudx gt 1 then begin
    if vv gt 0 then message,'non-uniform gridding; rebinning',/info
    regrid=1
    mindx=min(abs(dx[o0]))
    xrange=[min(x),max(x)] & nbin=long((xrange[1]-xrange[0])/mindx)+1L
    yyr=rebinw(yy,xx,nbin,wrange=xrange, _extra=e)
    oldxx=xx & xx=nbin & yy=yyr
  endif
endif else begin
  if nx gt 0 then message,'X does not match Y; ignoring X',/info
  xx=findgen(ny)
endelse
nxx=n_elements(xx)

;	make the kernel
if nk eq 0 then begin
  ;aa=[xx[ny/2],ss,1.]
  aa=[0.5*(min(xx)+max(xx)),ss,1.]
  call_procedure,'x3model',xx,aa,kern,/norm,type=ttyp, _extra=e

  ;	make the kernel as small as possible by throwing out
  ;	points from either end that are too small to matter.
  ;	(must do this symmetrically, otherwise, as Erica points out,
  ;	there will be unsightly shifts in the smoothed spectrum)
  k=0L & sum_edge=0. & kernorm=total(abs(kern))
  if kernorm eq 0 then begin
    message,'kernel is zero?  returning without smoothing',/informational
    return,y
  endif
  dnorm=sum_edge/kernorm
  while dnorm le eps do begin
    nk=n_elements(kern) & kern=kern[1L:nk-2L]
    sum_edge=sum_edge+abs(kern[0])+abs(kern[nk-2L-1L])
    dnorm=sum_edge/kernorm
    k=k+2L
  endwhile
  if k eq ny or k eq ny-1L then begin
    message,'Kernel is trivial, nothing to smooth with; returning',/info
    return,y
  endif else if vv gt 5 then message,'compressing kernel to '+$
	strtrim(ny-k,2)+' points out of '+strtrim(ny,2),/info
  kernorm=total(abs(kern)) & kern=kern/kernorm	;renormalize
  usekern=kern
endif else begin
  if vv gt 0 then message,'using the input kernel to smooth the function',/info
  kern=[usekern[*]]
endelse

;	convolve
if vv gt 5 then message,'convolving...',/info
sy=convol(yy,reverse(kern),/edge_wrap)

;	rebin to original grid if necessary
if keyword_set(regrid) then begin
  noldxx=n_elements(oldxx) & ndx=n_elements(dx)
  xxr=[oldxx,oldxx[noldxx-1L]+dx[ndx-1L]]
  yyr=rebinw(sy,xx,xxr,wrange=xrange, _extra=e)
  sy=yyr
endif

return,sy
end
