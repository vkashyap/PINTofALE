function hawaliner,counts,wave,wrange=wrange,lsffwhm=lsffwhm,maxkern=maxkern,$
	clev=clev,maxiter=maxiter,fullerr=fullerr,verbose=verbose,$
	hawastr=hawastr,$
	kkcont=kkcont,ctline=ctline,wscale=wscale,$
	wcont=wcont,ewcont=ewcont,$
	lpos=lpos,lerrpos=lerrpos,lflx=lflx,lerrflx=lerrflx,$
	lwdt=lwdt,lerrwdt=lerrwdt, _extra=e
;+
;function	hawaliner
;	returns the wavelengths at which spectral lines are found
;	to be located in the high-resolution counts spectrum
;
;syntax
;	linewvl=hawaliner(counts,wave,wrange=wrange,lsffwhm=lsffwhm,$
;	maxkern=maxkern,clev=clev,maxiter=maxiter,/fullerr,verbose=verbose,$
;	hawastr=hawastr,kkcont=kkcont,ctline=ctline,wscale=wscale,$
;	wcont=wcont,ewcont=ewcont,lpos=lpos,lerrpos=lerrpos,$
;	lflx=lflx,lerrflx=lerrflx,$
;	lwdt=lwdt,lerrwdt=lerrwdt,type=type)
;
;parameters
;	counts	[INPUT; required] counts spectrum
;		* this routine was developed for integer counts spectra, but
;		  currently no explicit assumption of Poisson distributions
;		  is made, and it _ought_ to work for spectral intensities
;		  in any units -- but no guarantees!
;	wave	[INPUT; required] wavelength grid for the spectrum
;		* if size matches that of COUNTS, taken to be mid-bin values
;		* if size exceeds COUNTS by 1, taken to be grid boundaries
;		* if size is double that of COUNTS, taken to be two
;		  concatenated arrays of lower and upper boundaries of grid
;		* it is assumed that WAVE is on a regular grid, i.e., that
;		  WAVE[1:*]-WAVE is constant
;
;keywords
;	wrange	[INPUT] range of wavelengths to consider in the analysis
;		* if not given, or is not a 2D array, a range that includes
;		  99% of the total counts that includes the peak of the
;		  spectrum is used
;	lsffwhm	[I/O] fwhm of the spectral lines, in same units as WAVE
;		* forced to be greater than the FWHM of the strongest line
;		  in the spectrum
;		* if given, ignores all kernels that are smaller than this
;		* it is assumed that LSFFWHM does not vary significantly
;		  across the spectrum
;		* if undefined, then the program will calculate it internally
;		  (by finding the highest point in the spectrum and walking
;		  down the slopes) and this calculated will be returned on
;		  output
;		* TO REITERATE: if defined on input, will be unchanged on output
;	maxkern	[INPUT] maximum kernel size to use, in bins
;		* by default, we start with a kernel of size 3 and increase
;		  the scale such that the central +ve part of the next scale
;		  is as large as the full extent of the previous scale
;		* if not given, set to a size corresponding to 30*LSFFWHM
;	clev	[INPUT] the confidence level at which to filter the
;		correlation coefficients
;		* default is 0.95
;		* if < 0, abs(CLEV) is used
;		* if > 1 and < 100, then assumed to be given as a percentage
;		* if > 100, then 1-1/CLEV is used
;	maxiter	[INPUT] maximum number of iterations
;		* default is 10
;		* hardcoded minimum is 2
;	fullerr	[INPUT] if set, calls MCERROR to derive better error bars
;		* WARNING: this can slow the program to a crawl
;	hawastr	[OUTPUT] a structure that contains all sorts of useful arrays:
;		{filtered WAVE, filtered COUNTS, wavelength range,
;		scales, LSF width, estimated continuum, filtered continuum,
;		line pos, line pos error, line flux, line flux error,
;		line width, line width error}
;	kkcont	[OUTPUT] the continuua calculated at each scale
;		* this is an array of size [WVL,SCALES], where WVL
;		  is a filtered version of WAVE (see HAWASTR.X)
;	ctline	[OUTPUT] the line "spectrum" devoid of continuum
;		* on same grid as WAVE
;	wcont	[OUTPUT] the continuum "spectrum" devoid of lines
;		* on same grid as WAVE
;	ewcont	[OUTPUT] error on WCONT
;		* on same grid as WAVE
;		* note that the errors on each bin are highly correlated
;		  with those at adjacent bins and should not be assumed to
;		  be independent
;	wscale	[INPUT] the kernel sizes to be used in the calculation
;		* the actual width of the wavelets will be 3*WSCALE
;		* if not given, will be calculated as
;		  2*LSFFWHM/dWAVE*[1,3,...,NBIN/3]/3
;	lpos	[OUTPUT] the positions of the detected lines
;		* (same as primary output)
;	lerrpos	[OUTPUT] errors on LPOS
;	lflx	[OUTPUT] the fluxes of the detected lines
;	lerrflx	[OUTPUT] errors on LFLX
;	lwdt	[OUTPUT] the widths of the detected lines
;	lerrwdt	[OUTPUT] errors on LWDT
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		FIT_LEVMAR : TYPE
;
;description
;	carries out a wavelet-like analysis, by computing the correlations
;	of the counts spectrum with Haar-like wavelets, and identifies bins
;	that are likely to have lines, finds the nearest local maximum
;	in the data corresponding to each identified line, and computes
;	their mean location
;
;subroutines
;	HIPD_INTERVAL()
;	MID2BOUND()
;	GETLOCMAX()
;	FIT_LEVMAR
;	MCERROR
;	ADJUSTIE()
;	LEVMARQ
;	LMCOEFF
;	LNPOISSON()
;
;history
;	Vinay Kashyap (Apr2007)
;	added call to MCERROR, included SIGY in call to FIT_LEVMAR,
;	  changed CTCONT to KKCONT, added keyword EWCONT (VK; Jun2007)
;	various tweaks and bug fixes, changed behavior of WSCALE
;	  and LSFFWHM (VK; Jul2008)
;
;about the name
;	It is a Haar-like Wavelet-based line finding program, so works
;	for the purpose of taking up a unique spot in namespace.  Well,
;	unique at least until the Indians decide to make a luxury airliner.
;-

;	usage
ok='ok' & np=n_params() & nc=n_elements(counts) & nw=n_elements(wave)
if np lt 2 then ok='Insufficient parameters' else $
 if nc eq 0 then ok='COUNTS array undefined' else $
  if nw eq 0 then ok='WAVE array undefined' else $
   if nw ne nc and nw ne nc+1L and nw ne 2L*nc then $
   	ok='COUNTS and WAVE arrays incompatible' else $
    if nc lt 9 then ok='arrays too small to bother'
if ok ne 'ok' then begin
  print,'Usage: linewvl=hawaliner(counts,wave,wrange=wrange,lsffwhm=lsffwhm,$'
  print,'       maxkern=maxkern,clev=clev,maxiter=maxiter,/fullerr,verbose=verbose,$'
  print,'       hawastr=hawastr,kkcont=kkcont,ctline=ctline,wscale=wscale,$'
  print,'       wcont=wcont,lpos=lpos,lerrpos=lerrpos,lflx=lflx,lerrflx=lerrflx,$
  print,'       lwdt=lwdt,lerrwdt=lerrwdt,type=type)
  print,'  returns the wavelengths at which spectral lines are found'
  print,'  to be located in the high-resolution counts spectrum'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
ct=float(counts) & wvl=findgen(nc) & wgrid=findgen(nc+1L)-0.5
if nw eq nc then begin	;(WAVE are mid-bin values
  wvl=float(wave)
  wgrid=mid2bound(wvl, _extra=e)
endif			;NW==NC)
if nw eq nc+1L then begin	;(WAVE are grid boundaries
  wgrid=float(wave)
  wvl=0.5*(wgrid[1:*]+wgrid)
endif				;NW==NC+1)
if nw eq 2L*nc then begin	;(WAVE are BIN_LO,BIN_HI arrays
  w1=wave[0L:nc-1L] & w2=wave[nc:*]
  if w1[0] lt w2[0] then wgrid=[w1[0],w2] else wgrid=[w1,w2[0]]
  wvl=0.5*(w1[0L:nc-1L]+w2[nc:*])
endif				;NW==2*NC)

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
nwr=n_elements(wrange)
if nwr ne 2 then $	; choose an interval that includes 99% of the counts
  wvlrange=hipd_interval(ct,wvl,clev=0.99) else wvlrange=wrange
wmin=min(wvlrange,max=wmax)
ow=where(wgrid ge wmin and wgrid le wmax,mow)
if mow eq 0 then begin
  message,'WRANGE selects nothing',/informational
  return,-1L
endif
wvl=wvl[ow] & ct=ct[ow] & iwoffset=ow[0]
dw=abs(median(wgrid[1:*]-wgrid))	;bin width, assumed constant
;
lsfwdt=3.
;	if LSFFWHM
;	figure out the width of the strongest line and
;	force this to be the minimum FWHM, regardless of
;	what is specified
ctmax=max(ct,imx) & go_on=1 & kplus=1L & kminus=1L
while go_on do begin
  i0=imx-kminus > 0
  i1=imx+kplus < (mow-1L)
  c0=ct[i0] & c1=ct[i1]
  if c0 gt 0.5*ctmax then kminus=kminus+1L
  if c1 gt 0.5*ctmax then kplus=kplus+1L
  if c0 le 0.5*ctmax and c1 le 0.5*ctmax then go_on=0
  if i0 eq 0 and i1 eq mow-1L then go_on=0
endwhile
wminus=wvl[i0] ;+(wvl[imx]-wvl[i0])*ctmax/2./(ctmax-ct[i0])
wplus=wvl[i1] ;+(wvl[imx]-wvl[i1])*ctmax/2./(ctmax-ct[i1])
;	eh, the interpolation doesn't work because linear
;	approximation causes errors of ~2x.  this way, errors
;	are at most ~2bins, and conservative.  if you have a
;	better value, specify it in the call.
lsfwdt=abs(wplus-wminus)/dw
if keyword_set(lsffwhm) then begin
  ;lsfwdt=abs(float(lsffwhm[0]))/dw > lsfwdt
  lsfwdt=abs(float(lsffwhm[0]))/dw > 1
endif else lsffwhm=lsfwdt*dw
lsfsigma=lsfwdt*dw/2.355
;
;kernmax=lsfwdt*100.
kernmax=mow/3 > (lsfwdt*100.)
if keyword_set(maxkern) then kernmax=abs(float(maxkern[0]))
kernmax = kernmax < mow
;
itermax=10 & if keyword_set(maxiter) then itermax=long(maxiter[0])>2
;
if not keyword_set(clev) then clevel=0.95 else begin
  clevel=0.0+clev[0]
  if clevel lt 0 then clevel=abs(clevel)
  if clevel ge 1 and clevel lt 100 then clevel=clevel/100.
  if clevel ge 100 then clevel=1.D - 1.D/clevel
endelse

;	so, how many kernels can we fit in here?
kk=[long(2*lsfwdt)>3] & while max(kk) lt kernmax/3 do kk=[kk,3L*max(kk)]
if n_elements(wscale) ne 0 then begin
  if keyword_set(wscale) then kk=3*(long(wscale[*])) > 3
endif
nkern=n_elements(kk)

;	arrays to store the intermediate products
ww=fltarr(mow,nkern)	;correlation coefficients at first pass
ol=lonarr(mow,nkern)-1L	;indices to keep track of where the lines are
oc=lonarr(mow,nkern)-1L	;indices to keep track of where the lines ain't
kkcont=fltarr(mow,nkern)	;the "background"
kkconterr=fltarr(mow,nkern)	;error on kkcont
hpd=fltarr(2,nkern)	;highest density intervals at first pass
ctc=fix(ct) & ctl=ctc	;continuum and lines

;	get the correlations with the Haar wavelets
for i=0L,nkern-1L do begin		;{for each scale
  ;	the wavelet kernel is [-1/2,+1,-1/2]
  nk=kk[i]/3
  kern=[-1./2./nk+fltarr(nk),$
	 1./nk+fltarr(nk),$
	-1./2./nk+fltarr(nk)]
  negkern=(-kern)>0
  if vv gt 5 then kilroy,dot=strtrim(nk,2)
  ;
  go_on=1 & kiter=0L & ctj=ct
  if vv gt 10 then plot,wvl,ctj,/xs
  while go_on eq 1 do begin		;{iterate on lopping off lines
    if vv gt 5 then kilroy
    oldctj=ctj & kiter=kiter+1L
    wk=convol(ctj,kern,/edge_wrap)
    hpdint=hipd_interval(wk,/fsample,clev=clevel)
    o1=where(wk gt hpdint[1],mo1,complement=xo1)
    if kiter eq 1 then ww[*,i]=wk
    if kiter eq 1 then hpd[*,i]=hpdint
    if mo1 gt 0 then ol[o1,i]=o1
    if mo1 lt mow then oc[xo1,i]=xo1
    bgk=convol(ctj,negkern,/edge_wrap)
    ebgk=sqrt(convol(ctj,negkern^2,/edge_wrap))
    if mo1 gt 0 then ctj[o1]=bgk[o1]
    ;
    if vv gt 15 then oplot,wvl,ctj,col=(kiter mod 8)+1
    if total(abs(ctj-oldctj))/total(ctj) lt 1e-3 then go_on=0
    if total(abs(ctj-oldctj)) lt 1 then go_on=0
    if kiter ge itermax then go_on=0
    ;
    if mo1 gt 0 then begin
      ctc[o1]=!values.F_NAN
      if not keyword_set(oll) then oll=o1 else oll=[oll,o1]
    endif
    if mo1 lt mow then begin
      if not keyword_set(occ) then occ=xo1 else occ=[occ,xo1]
    endif
  endwhile				;GO ON}
  ;
  wk=convol(ctj,kern,/edge_wrap)
  o1=where(abs(wk) gt bgk+3*(sqrt(bgk+0.75)+1.),mo1,complement=xo1)
  if mo1 gt 0 then begin
    ol[o1,i]=o1
    oc[xo1,i]=xo1
    ctj[o1]=bgk[o1]
    ctc[o1]=!values.F_NAN
    if not keyword_set(oll) then oll=o1 else oll=[oll,o1]
    if mo1 lt mow then begin
      if not keyword_set(occ) then occ=xo1 else occ=[occ,xo1]
    endif
    if vv gt 15 then oplot,wvl[o1],ctj[o1],psym=1
  endif
  kkcont[*,i]=bgk
  kkconterr[*,i]=ebgk
  if vv gt 100 then begin
    jll=oll[uniq(oll,sort(oll))]
    tmp=fix(ct) & tmp[jll]=!values.F_NAN
    plot,wvl,ct,/xs & oplot,wvl,ct-tmp,col=2
    if vv gt 2000 then stop,'HALTing; type .CON to continue'
  endif
endfor					;I=0,NKERN-1}
ctline=fix(ct)-ctc
;tmp=0*ct & tmp[iwoffset:iwoffset+mow-1L]=ctline[*] & ctline=tmp	;what does this line do???
if keyword_set(oll) then oll=oll[uniq(oll,sort(oll))]
if keyword_set(occ) then occ=occ[uniq(occ,sort(occ))]

if vv gt 50 then begin
  plot,wvl,ct,/xs & oplot,wvl,ctl,col=2
  wait,(0.01*vv)<2
endif

;	now find all the line locations
;sct=smooth(smooth(ct,lsfwdt/2,/edge_truncate),lsfwdt/2,/edge_truncate)
sct=smooth(smooth(ct,lsfwdt,/edge_truncate),lsfwdt,/edge_truncate)
oline=where(ctl gt 0,moline) & iline=lonarr(moline)
;ilm=getlocmax(wvl,sct,width=lsfwdt,sigmay=sqrt(sct),nsigma=0.1,/flattop)
ilm=getlocmax(wvl,sct,width=lsfwdt,sigmay=sqrt(sct),nsigma=1,/flattop)
;
for i=0L,moline-1L do begin
  tmp=min(abs(oline[i]-ilm),imn) & iline[i]=ilm[imn]
endfor
uniqiline=iline[uniq(iline,sort(iline))]
lpos=wvl[uniqiline] & nline=n_elements(uniqiline)
if vv gt 0 then message,strtrim(nline,2)+' lines found',/informational

;	compute the continuum
wcont=0.*ct & for i=0,nkern-1L do wcont=wcont+kkcont[*,i] & wcont=wcont/nkern
;tmp=0*ct & tmp[iwoffset:iwoffset+mow-1L]=wcont[*] & wcont=tmp		;what does this line do???
ewcont=0.*ct & for i=0,nkern-1L do ewcont=ewcont+kkconterr[*,i]^2 & ewcont=sqrt(ewcont)/nkern
;tmp=0*ct & tmp[iwoffset:iwoffset+mow-1L]=ewcont[*] & ewcont=tmp	;what does this line do???

if vv gt 1500 then stop,'HALTing; type .CON to continue'

;	now do some post processing
if arg_present(hawastr) then begin
  lflx=0.*lpos & lwdt=lflx
  lerrpos=lflx & lerrflx=lflx & lerrwdt=lflx
  for i=0L,nline-1L do begin
    j=uniqiline[i]
    i0=(j-2*lsfwdt) > 0
    i1=(j+2*lsfwdt) < (mow-1L)
    xx=wvl[i0:i1] & yy=ct[i0:i1] & bg=wcont[i0:i1] & sigy=sqrt(yy+0.75)+1.
    ymx=max(yy,imx)
    ties=[$
	'A0=(A0>('+strtrim(lpos[i]-lsfsigma,2)+'))<('+strtrim(lpos[i]+lsfsigma,2)+')',$
	'A1=(A1>('+strtrim(lsfsigma/2,2)+'))<('+strtrim(5*lsfsigma,2)+')',$
	'A2=(A2>0)']
    aa=[lpos[i],lsfsigma,ymx-bg[imx]]
    fit_levmar,xx,yy-bg,aa,yfunc,sig=sigy,erra=erra,chisq=chisq,/dumb,ties=ties,$
	funcs='libmodel',verbose=(vv<1), _extra=e
    fit_levmar,xx,yy-bg,aa,yfunc,sig=sigy,erra=erra,chisq=chisq,/dumb,ties=ties,$
	funcs='libmodel',verbose=(vv<1), _extra=e
    lpos[i]=aa[0] & lwdt[i]=aa[1] & lflx[i]=aa[2]
    lerrpos[i]=erra[0] & lerrwdt[i]=erra[1] & lerrflx[i]=erra[2]
    if vv gt 5 then begin
      plot,xx,yy-bg,psym=10,title='line # '+strtrim(i+1,2)
      oplot,xx,yfunc,color=2
      if vv gt 200 then wait,(0.001*vv)<0.5
    endif
    aa0=aa
    if keyword_set(fullerr) then begin
      mcerror,xx,yy-bg,aa0,errav,ysig=sigy,/dumb,ties=ties,$
	funcs='libmodel',type='gauss',algo='LevMarq+SVD',verbose=0
      lerrpos[i]=errav[0] & lerrwdt[i]=errav[1] & lerrflx[i]=errav[2]
    endif
  endfor
  hawastr=create_struct('X',wvl,'Y',fix(ct),'XMIN',wmin,'XMAX',wmax,$
	'SCALES',kk,'LSFWIDTH',lsfwdt,'CONTINUUM',wcont,$
	'YCONT',wcont[ow],'YCONTERR',ewcont[ow],$
	'LINEPOS',lpos,'LINEPOSERR',lerrpos,$
	'LINEFLX',lflx,'LINEFLXERR',lerrflx,$
	'LINEWDT',lwdt,'LINEWDTERR',lerrwdt)
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,lpos
end
