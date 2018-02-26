function linerem,lamda,spec,sig=sig,cell=cell,nsigma=nsigma,$
	bkgval=bkgval,bkgerr=bkgerr,quiet=quiet,posve=posve,negve=negve,$
	_extra=e
;+
;function	linerem
;	remove lines from input spectrum and return the "cleaned" spectrum
;
;syntax
;	cleanspec=linerem(lamda,spec,sig=sig,cell=cell,nsigma=nsigma,$
;	bkgval=bkgval,bkgerr=bkgerr,/quiet,/posve,/negve)
;
;parameters
;	lamda	[INPUT; required] wavelengths at which spectrum
;		is defined.
;	spec	[INPUT; optional] the spectrum.
;		NOTE: if not given, LAMDA is taken to be SPEC and
;		the array indices are taken to be LAMDA
;
;keywords
;	sig	[INPUT] error at each point; if not given, the
;		errors are taken to be 1+sqrt(abs(SPEC)+0.75).
;	nsigma	[INPUT; default: 4] multiple of SIG to consider
;		as a threshold for detection of features
;	cell	[INPUT] 1D filter to use in computing the background
;		* default is [1,1,0,0,0,1,1]
;		* if scalar, then [IC+1,ICC,IC+1], where
;		  IC=intarr(2*CELL), ICC=intarr(2*CELL+1)
;	bkgval	[OUTPUT] final background values at each bin
;		NOTE: This is essentially a smoothed version of
;		the cleaned spectrum!
;	bkgerr	[OUTPUT] error estimates on BKGVAL
;	quiet	[INPUT] if set, doesn't show, doesn't tell
;	posve	[INPUT] if set, removes only +ve deviations
;	negve	[INPUT] if set, removes only -ve deviations
;	_extra	[JUNK] here only to prevent crashes
;
;description
;	1. make sure that the spectrum is defined on a regular grid.
;	   (if not, rebin [NOT IMPLEMENTED!])
;	2. convolve the spectrum with background cell to determine
;	   local background.
;	3. also propagate errors (assume gaussian; if poisson, use
;	   gaussian approximation)
;	4. compare local value with local background
;	5. flag those bins which are significantly different from local
;	   background (use +-NSIGMA*SIG as a threshold value)
;	6. reset the flagged bin values to local background values
;	7. repeat 2-6 until no new bins are flagged
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
;	vinay kashyap (Oct96)
;	changed keyword SIGMA to SIG (VK; Jun98)
;	CELL was being altered if input as scalar -- corrected; altered
;	  behavior of output when no deviations are found (VK; Oct98)
;	now if QUIET>1, won't ask to quit if many bins are found (VK; 99Aug)
;	changed POS and NEG to POSVE and NEGVE (VK; MMJul)
;-

np=n_params(0)
if np eq 0 then begin
  print, 'Usage: clean_spectrum=linerem(lamda,spectrum,sig=sig,nsigma=nsigma,$'
  print, '       cell=cell,bkgval=bkgval,bkgerr=bkgerr,/quiet,/posve,/negve)'
  print, '   removes lines from spectrum and returns continuum'
  return,-1L
endif

;	relabel
x=lamda & nx=n_elements(x) & ny=n_elements(spec)
if ny ne nx then begin
  y=x & x=lindgen(nx)
endif else y=spec

;	check keywords
if not keyword_set(sig) then sig=sqrt(abs(y)+0.75)+1.
if not keyword_set(nsigma) then nsigma=4.
if keyword_set(cell) then begin
  csize=cell
  nc=n_elements(cell)
  if nc eq 1 then begin
    cc=abs(cell(0))
    ic=intarr((2*cc)>1) & icc=intarr(2*cc+1) & csize=[ic+1,icc,ic+1]
    nc=n_elements(csize)
  endif
  if nc gt nx then begin
    message,'Background cell too large',/info & return,y
  endif
endif else csize=[1,1,0,0,0,1,1]
if not keyword_set(quiet) then quiet=0
if keyword_set(posve) and keyword_set(negve) then begin
  posve=0 & negve=0
endif

;	initialize
norm=total(csize)
if norm le 0. then begin
  c1='background cell has zero area! ignoring normalization'
  message,c1,/info
endif

;	check to see that the spectrum is defined on a uniform grid
dx=x(1:*)-x & ddx=dx(uniq(dx,sort(dx))) & nddx=n_elements(ddx)
odx=where(abs(dx-median(dx))/median(dx) gt 1e-3,modx)
if modx gt 1 then begin
  c1='spectrum not on uniform grid... convolutions may result in nonsense'
  message,c1,/info
  if not quiet then print, 'there are '+strtrim(nddx,2)+' bin sizes:',ddx
  if quiet lt 2 then begin
    c1='hit any key to continue, q to return, x to stop' & print,c1
    c1=get_kbrd(1)
  endif else c1=''
  if strlowcase(c1) eq 'q' then return,y
  if strlowcase(c1) eq 'x' then begin
    print,'there are '+strtrim(nx,2)+' bins and '+strtrim(nddx,2)+$
	' unique bin widths'
    help,x,y,dx,ddx
    stop,'type RETURN,Y to return sans change'
  endif
endif

;	convolve spectrum with cell to get local background
bkg=convol(y,csize,/edge_truncate) & if norm ne 0 then bkg=bkg/norm

;	propagate errors
bge=convol(sig^2,abs(csize),/edge_truncate) & if norm ne 0 then bge=bge/norm
bge=sqrt(bge)

;	find significant deviations
dely=(y-bkg) & nok=where(abs(dely) gt nsigma*abs(bge))
if keyword_set(posve) then nok=where(dely gt nsigma*abs(bge))
if keyword_set(negve) then nok=where(dely lt -nsigma*abs(bge))
if nok(0) eq -1 then begin
  message,'no significant deviations.. spectrum is the continuum?',/info
  bkgval=bkg & bkgerr=bge
  return,y					;hey.
endif
y(nok)=bkg(nok)
inew=lonarr(nx) & inew(nok)=1

;	show
if not quiet then begin
  print,strtrim(long(total(inew)),2)+' bins reset'
  print,'plotting lambda v/s spectrum, line-removed spectrum, bkg w. errors'
  plot,x,spec,psym=10,/xs & oplot,x,y,col=100 & oplot,x,bkg,col=150
  oplot,x,bkg+nsigma*bge,col=150,line=2
  oplot,x,bkg-nsigma*bge,col=150,line=2
endif

;	iterate
y1=y & bkg1=bkg & sigg=sig & sig(nok)=bge(nok)
while total(inew) gt 0 do begin
  ;	all as above...
  bkg1=convol(y1,csize,/edge_truncate) & if norm ne 0 then bkg1=bkg1/norm
  bge1=convol(sigg^2,abs(csize),/edge_truncate)
  if norm ne 0 then bge1=bge1/norm
  bge1=sqrt(bge1)
  dely=(y1-bkg1) & nok=where(abs(dely) gt nsigma*abs(bge1))
  if keyword_set(posve) then nok=where(dely gt nsigma*abs(bge))
  if keyword_set(negve) then nok=where(dely lt -nsigma*abs(bge))
  inew=lonarr(nx)
  if nok(0) ne -1 then begin
    inew(nok)=1 & y1(nok)=bkg1(nok) & sigg(nok)=bge(nok)
  endif
  if not quiet then begin
    print,strtrim(long(total(inew)),2)+' bins reset in this iteration'
    plot,x,spec,psym=10,/xs & oplot,x,y1,col=100 & oplot,x,bkg1,col=150
    oplot,x,bkg1+nsigma*bge1,col=150,line=2
    oplot,x,bkg1-nsigma*bge1,col=150,line=2
  endif
endwhile

bkgval=bkg1 & bkgerr=sigg

return,y1
end
