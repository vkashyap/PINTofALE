pro smoothie,y,outy,outye,ix,yerr=yerr,snrthr=snrthr,verbose=verbose, _extra=e 
;+ 
;procedure       smoothie
;	return the input function smoothed such that a minimum required S/N 
;	is met at each point in the function.  useful for smoothing out the
;	weak wings and continuua in line dominated spectra.  starts from
;	the peak of the lines and steps down the profile, accumulating the
;	bins as necessary.
;
;	note that unlike most adaptive smoothing algorithms, this one
;	conserves flux.
; 
;syntax 
;       smoothie,y,outy,outye,ix,yerr=yerr,snrthr=snrthr,verbose=v
; 
;parameters 
;       y	[INPUT; required] array to be smoothed 
;       outy	[OUTPUT] smoothed array
;       outye	[OUTPUT] errors on smoothed array
;		* WARNING: these do not include correlation errors,
;		  and so do not represent a rigorous estimate.  i.e.,
;		  even though correctly propagated for each bin, they
;		  cannot be used as independent errors on the values in
;		  adjacent bins.
;       ix	[OUTPUT] bin indices -- tells which bins are grouped
;		* all bins which have the same value of IX are in the same group
;
;keywords 
;       yerr	[INPUT] error on Y
;		* default is 1.+sqrt(Y+0.75)
;		* the default is set that way for consistency with
;		  standard PINTofALE behavior, but SQRT(Y) is a better
;		  option for this algorithm.
;       snrthr	[INPUT] S/N threshold below which to smooth bins
;		* default is 3.0
;	verbose	[INPUT] controls chatter
;	
;	_extra	[JUNK] here only to prevent crashing
; 
;example 
;	y=50*(mk_gauss(findgen(256),128,10,1)+0.01) & yr=long(0*y)
;	for i=0,256-1 do yr[i]=randomu(seed,poisson=y[i])
;	smoothie,yr,outy,outye,ix,yerr=1+sqrt(0.75+yr),snrthr=4.0
;	plot,y & oplot,yr,psym=10,col=2 & oplot,outy,col=3
;	** notice the wings **
;
;	y=randomu(seed,1000,poisson=1)
;	y[400:499]=randomu(seed,100,poisson=5)
;	smoothie,y,outy,outye,snrthr=5
;	plot,y & oplot,outy,col=2,thick=3
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
;	VK/LL (Nov 04)
;	cleaned up I/O (VK; Jan05)
;	now returns a flat average if total(y)/sqrt(total(ysig^2))<SNRTHR
;	  (VK; Oct06)
;	bugfix: unnecessarily crashing because of a variable name mismatch
;	  when too few counts to meet SNR threshold (VK; Oct09)
;- 

;       usage
ok='ok' & np=n_params(0) & ny=n_elements(y)
if np eq 0 then ok='Insufficient inputs' else $
 if ny eq 0 then ok='Y is undefined' else $
  if ny lt 3 then ok='Y is too small'
if ok ne 'ok' then begin
  print, 'Usage: smoothie,y,outy,outye,ix,yerr=yerr,snrthr=snrthr,$'
  print,'        verbose=verbose'
  if np ne 0 then message,ok,/informational
  return
endif

;	check inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
nye=n_elements(yerr) & ysig=1.+sqrt(y+0.75)
if nye ne 0 then begin
  if nye eq ny then ysig=yerr else begin
    if nye eq 1 then begin
      if yerr[0] gt 0 and yerr[0] lt 1 then ysig=yerr[0]*y else $
       if yerr[0] ge 1 and yerr[0] lt 100 then ysig=yerr[0]*y/100. else $
	if yerr[0] ge 100 then ysig=0.*y+yerr[0] else $
	 if yerr[0] lt 0 then ysig=0.*y+abs(yerr[0])
    endif else message,$
	'YERR incompatible with Y; assuming poisson',/informational
  endelse
endif
snrthresh=3.0 & if keyword_set(snrthr) then snrthresh=float(snrthr[0])>0

;	outputs and other variables
kx=findgen(ny)
outy=1.*y & outye=ysig & ix=lonarr(ny)	;outputs
ok=where(ix eq 0,mok)		;look only at these bins
ty=0. & tye2=0. & ibin=1L	;initialize local variables

;       feassibility test  
if total(y)/sqrt(total(ysig^2)) lt snrthresh then begin 
   message, 'The required S/N is unrealistic with the input y',/informational
   outy[*]=mean(y) & outye[*]=sqrt(total(ysig^2)/ny)
   return
endif

;	smooth
while mok gt 0 do begin		;{bins remain to be smoothed
  if ty eq 0 then begin
    ;	start new group
    ty=max(outy[ok],j) & i=ok[j]
    ty=outy[i] & tye2=outye[i]^2
    ix[i]=ibin
    in1 = i & ip1 = i
  endif else begin
  while in1 ge 0 and (ix[[in1]])[0] ne 0 do in1=in1-1L
  while ip1 lt ny and (ix[[ip1]])[0] ne 0 do ip1=ip1+1L
  if vv gt 10 then print, in1, ip1
     if in1 ge 0 then begin
         ty=ty+outy[in1] & tye2=tye2+outye[in1]^2
        ix[in1]=ibin
    endif
    if ip1 lt ny then begin
        ty=ty+outy[ip1] & tye2=tye2+outye[ip1]^2
        ix[ip1]=ibin
    endif 
  endelse
  ok=where(ix eq 0,mok)
  tstthr=ty & if tye2 gt 0 then tstthr=ty/sqrt(tye2)
  if tstthr ge snrthresh or mok eq 0 then begin 
    oo=where(ix eq ibin,moo)
    if moo gt 0 then begin
      outy[oo]=ty/moo & outye[oo]=sqrt(tye2)/moo
    endif
    ty=0. & tye=0. & ibin=ibin+1L
  endif
endwhile			;MOK>0}

end
