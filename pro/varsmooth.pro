function varsmooth,func,scale,xfunc=xfunc,weight=weight,type=type,$
	steps=steps,gauss=gauss,boxcar=boxcar,nueff=nueff,nutilde=nutilde,$
	_extra=e
;+
;function	varsmooth
;	returns a locally smoothed array
;
;	by default, smooths with a boxcar of specified scale.
;	may also smooth with gaussians, lorentzians, beta-profiles,
;	and asymmetric beta-profiles of varying widths; or obtain
;	stepped averages over varying widths.
;
;syntax
;	sfunc=varsmooth(func,scale,xfunc=xfunc,weight=weight,type=type,$
;	steps=steps,/gauss,/boxcar,nueff=nueff,nutilde=nutilde,$
;	missing=missing,/fwhm,betap=betap,angle=angle)
;
;parameters
;	func	[INPUT; required] 1D array to be smoothed
;		* if multi-D, gets converted to 1D
;	scale	[INPUT] array containing the half-length scales over which
;		FUNC must be smoothed.
;		* if smaller in size than FUNC, fills out with last element
;		* must be in bin units, unless
;		  -- XFUNC and TYPE are both defined
;
;keywords
;	xfunc	[INPUT] x-coordinates for FUNC
;		* size must match FUNC or must be bin boundaries on FUNC,
;		else ignored.
;		* ignored in case of boxcar smoothing
;	weight	[INPUT] allows weighting of adjacent points while boxcar
;		smoothing
;		* if scalar, WEIGHT(*)=(dPIX)^(-WEIGHT(0))
;	type	[INPUT] define the 3-parameter family function to be used
;		to smooth
;		* may be anything accepted by X3MODEL/MK_3MODEL, e.g.,
;		  "Gaussian", "Lorentizan", "Beta=<value>",
;		  "Slant=<angle>,<beta>", etc.
;	steps	[INPUT] if set, ignores the actual values of SCALE
;		and replaces FUNC in bins with contiguous stretches of
;		SCALE with <FUNC> averaged over this contiguous stretch.
;		* only those values of SCALE > STEPS[0] are considered
;		  in the calculation -- those below it are held at one
;		  bin each.
;		* NUEFF and NUTILDE are not calculated for this case
;	gauss	[INPUT] retained for backwards compatibility.  if set,
;		and TYPE is not set, then sets TYPE='Gaussian'
;	boxcar	[INPUT] if TYPE, GAUSS and STEPS are not set, the default
;		is to carry out boxcar type smoothing.  however, the default
;		there is to limit the smoothing scale to the maximum
;		extent of the range, i.e., N(LSCAL)/2.  if BOXCAR is
;		set, this limitation is ignored, and smoothing scales
;		can be as large as possible.  however, note that this
;		does not handle spillage over endpoints gracefully, so
;		_THIS OPTION IS NOT RECOMMENDED_ unless you know what
;		you are doing.
;	nueff	[OUTPUT] the effective number of degrees of freedom
;		in the problem, after smoothing.
;		* the trace of the smoothing matrix
;	nutilde	[OUTPUT] correction to NUEFF
;		* the trace of the product of the transpose of the smoothing
;		  matrix and the matrix
;		* dof ~ 2*NUEFF-NUTILDE
;	_extra	[INPUT] pass defined keywords to
;		MK_GAUSS: MISSING, FWHM
;		MK_LORENTZ: BETAP, MISSING
;		MK_SLANT: ANGLE, BETAP, MISSING
;
;restrictions
;	input array, FUNC, must be evenly sampled.
;	requires subroutines:
;		X3MODEL
;		MK_3MODEL
;		MK_GAUSS
;		MK_LORENTZ
;		MK_SLANT
;		KILROY
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
;	vinay kashyap (Apr97)
;	changed gaussian smoothing to wrap around (VK; May97)
;	corrected bug of SCALE not filling out with last element (VK; Aug98)
;	added keywords XFUNC, STEPS, and TYPE (latter makes GAUSS obsolete)
;	  (VK; MarMM)
;	allowed XFUNC to be bin boundaries (VK/LL; Sep02)
;	now variable rebinning takes effect only if scales are > STEPS[0]
;	  (VK; Feb03)
;	added keywords NUEFF and NUTILDE (VK; Jan06)
;	corrected NUTILDE calculation; added keyword BOXCAR (VK; Mar06)
;-

;	usage
if n_elements(func) eq 0 then begin
  print,'Usage: sfunc=varsmooth(func,scale,xfunc=xfunc,weight=weight,$'
  print,'       type=type,steps=steps,/gauss,nueff=nueff,nutilde=nutilde,$'
  print,'       missing=missing,/fwhm,betap=betap,angle=angle)'
  print,'  returns array smoothed at local scales'
  return,0L
endif

;	cast variables
f=[func(*)] & nf=n_elements(f)			;input function
;
s=0*f & ns=n_elements(scale)			;smoothing half-scales
if ns gt 0 then begin
  s(*)=scale(ns-1)
  if ns le nf then s(0:ns-1)=scale else s=scale(0:nf-1)
endif
;
x=findgen(nf) & usex=0				;x-coordinate of function
nx=n_elements(xfunc)
if nx eq nf or nx eq nf+1L then begin
  if nx eq nf+1L then x=0.5*(xfunc[1:*]+xfunc) else x=[xfunc(*)]
  usex=1
endif else begin
  if nx gt 0 then message,'XFUNC does not match FUNC; ignoring XFUNC',/info
endelse
;
wts=0*f+1. & nw=n_elements(weight)		;weights for measure
if nw eq 1 then wts=(findgen(nf)+1)^(-weight(0))
if nw gt 1 and nw lt nf then begin
  wts(0:nw-1)=weight & wts(nw:*)=weight(nw-1)
endif
if nw ge nf then wts=weights(0:nf-1)
rwts=reverse(wts)			;for time saving

;	initialize
sfunc=0*f					;the output

;	backwards compatibility hack
if keyword_set(gauss) and not keyword_set(type) then type='Gaussian'

if arg_present(nutilde) or arg_present(nueff) then ll=fltarr(nf,nf)

;	smooth with 3-parameter family of functions
if keyword_set(type) then begin
  s0=s(0) & a0=[x(nf/2),s0,1.]
  call_procedure,'x3model',x,a0,fmod,/norm,type=type, _extra=e
  fmod=shift(fmod,-nf/2-1L)/total(fmod)
  for i=0L,nf-1L do begin		;{step through input function
    if i eq 100*long(i/100.) then kilroy; make a mark.
    if usex eq 1 then s0=s(i)-1.	;force recalculation at each point
    if s(i) ne s0 then begin
      a0=[x(i),s(i),1.] & s0=s(i)
      call_procedure,'x3model',x,a0,fmod,/norm,type=type, _extra=e
      fmod=fmod/total(fmod)
    endif else fmod=shift(fmod,1)
    sfunc=sfunc+f(i)*fmod
    if arg_present(nutilde) or arg_present(nueff) then ll(i,*)=fmod
  endfor				;I=0,NF-1}
  if arg_present(nutilde) or arg_present(nueff) then begin
    ii=lindgen(nf) & nueff=total(ll(ii,ii))
    nutilde=total( ((transpose(ll))(ii,ii))*(ll(ii,ii)) )
  endif
  ;
  ;	git
  ;
  return,sfunc
endif

;{OBSOLETE
;;	gaussian smooth
;x=findgen(nf) & z=(x-nf/2)^2/2.
;if keyword_set(gauss) then begin		;{add up gaussians
;  ;DEBUG: plot,f
;  for i=0,nf-1 do begin			;(step through function
;    ;OLD code: g=(x-i)^2/2.
;    g=shift(z,nf/2+i)				;wrap around
;    if s(i) ne 0 then begin
;      g=g/s(i)^2 & norm=1./s(i)/sqrt(2.*!pi)
;    endif else begin
;      g(*)=69. & g(i)=0. & norm=1.
;    endelse
;    g = g < 69.
;    sfunc=sfunc+f(i)*exp(-g)*norm
;    ;DEBUG: oplot,f(i)*exp(-g)*norm;,col=250-i*250/nf
;    ;DEBUG: oplot,sfunc;,col=i*250/nf
;    ;DEBUG: stop,i,s(i)
;  endfor					;I=0,NF-1)
;  ;
;  ;	git
;  ;
;  return,sfunc
;endif						;GAUSS}
;OBSOLETE}

;	make like a step
if keyword_set(steps) then begin
  ;	this part figures out the "bitmap" of the scales
  k=1L & ks=lonarr(nf)+k & smin=float(steps(0))
  for i=1L,nf-1L do begin
    if s(i) eq s(i-1) and s(i) gt smin then ks(i)=ks(i-1) else begin
      k=k+1L & ks(i)=k
    endelse
  endfor
  ;
  ;	and this part gets the averages in each stretch
  for i=1L,k do begin
    oo=where(ks eq i,moo)
    if moo eq 0 then message,'BUG!'
    sfunc(oo)=total(func(oo))/moo
  endfor
  ;
  ;	git
  ;
  return,sfunc
endif

;	boxcar smooth
for i=0,nf-1 do begin				;{step through function
  ;scl=long(float(s(i))/2+0.5)
  ;scl=long(float(s(i))/2)
  scl=s(i)
  im=(i-scl)>0 & ip=(i+scl)<(nf-1)
  dim=i-im & dip=ip-i
  sfunc(i)=f(i)
  if dim gt 0 then sfunc(i)=sfunc(i)+total(rwts(nf-dim:*)*f(im:i-1))
  if dip gt 0 then sfunc(i)=sfunc(i)+total(wts(0:dip-1)*f(i+1:ip))
  if keyword_set(boxcar) then begin
    sfunc(i)=sfunc(i)/(2*s(i)+1.)
    if arg_present(nutilde) or arg_present(nueff) then ll(i,*)=1./(2*s(i)+1.)
  endif else begin
    sfunc(i)=sfunc(i)/(dim+dip+1.)
    if arg_present(nutilde) or arg_present(nueff) then ll(i,*)=1./(dim+dip+1.)
  endelse
endfor						;I=0,NF-1}
if arg_present(nutilde) or arg_present(nueff) then begin
  ii=lindgen(nf) & nueff=total(ll(ii,ii))
  ;nutilde=total((transpose(ll)#ll)(ii,ii))
  nutilde=total( ((transpose(ll))(ii,ii))*(ll(ii,ii)) )
endif

return,sfunc
end
