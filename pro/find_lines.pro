function find_lines,cts,wvl,hwhm=hwhm,effar=effar,wvlar=wvlar,ctb=ctb,bkgscal=bkgscal,$
	exptime=exptime,verbose=verbose, _extra=e
;+
;function	find_lines
;		detect lines in an X-ray counts spectrum using a combination of wavelet
;		transforms, statistics, and heuristics, and return the list of found
;		lines in a structure of the form
;			{WLOC, EWLOC, OBSCT, ESTBG, pvalBG, ESTCONT, NETCT, NETCTERR, SNR,
;			FLUX, FLUXERR, WIDTH, SCALE, QUALITY}
;
;syntax
;	lines=find_lines(cts,wvl,hwhm=hwhm,effar=effar,wvlar=wvlar,ctb=ctb,bkgscal=bkgscal,exptime=exptime,verbose=verbose)
;
;parameters
;	CTS	[INPUT; required] the counts of the spectrum, must be integers
;	WVL	[INPUT] the wavelength grid over which CTS is defined
;		* if not given, assumed to be array indices
;		* if size is N(CTS), assumed to be mid-bin points, else
;		  if size is N(CTS)+1, assumed to be bin boundaries, else
;		  ignored
;
;keywords
;	hwhm	[INPUT] the approximate half-width at half-max of an isolated line
;		* this is used only to filter out narrower spikes, as well as to
;		  combine the counts around the nominal location to get line intensity
;		  estimates
;		* default is 3 times median WVL bin widths
;	effar	[INPUT] effective area, used to compute flux
;	wvlar	[INPUT] wavelength grid over which EFFAR is defined
;		* size must match EFFAR, otherwise it is as though it was not given
;		* if WVLAR is not given, then EFFAR is assumed to be defined over
;		  the range minmax(WVL), and is then interpolated to have the same
;		  size as CTS
;	ctb	[INPUT] counts spectrum of the background
;		* can be scalar or array
;		* if scalar, assumed to be constant in all bins
;		* if array, size must match CTS, or else is ignored
;		  -- can be integer or float
;	bkgscal	[INPUT; default=1] the ratio of the area in which CTB is
;		collected to the area in which CTS is collected
;	exptime	[INPUT; default=1] the exposure time of the observation, used only
;		to divide into FLUX
;	verbose	[INPUT; controls chatter]
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;
;history
;	Vinay Kashyap (MMXIX.XI)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(cts) & nw=n_elements(wvl)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='CTS undefined' else $
  if ns lt 10 then ok='too few bins in this "spectrum"'
if ok ne 'ok' then begin
  print,'Usage: lines=find_lines(cts,wvl,hwhm=hwhm,effar=effar,wvlar=wvlar,ctb=ctb,bkgscal=bkgscal,exptime=exptime,verbose=verbose)'
  print,'  detect lines and return their properties as a structure with fields'
  print,'  {WLOC, EWLOC, OBSCT, ESTBG, pvalBG, ESTCONT, NETCT, NETCTERR, SNR, FLUX, FLUXERR, WIDTH, SCALE, QUALITY}
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords and parameters
vv=0L & if keyword_set(verbose) then vv=(long(verbose[0]))>1L
;
wmid=findgen(ns)
if nw eq ns then wmid=wvl
if nw eq ns+1L then wmid=(wvl[1:*]+wvl)/2.
delw=abs(median(wmid[1:*]-wmid))
;
lsfwdt=1
if keyword_set(hwhm) then lsfwdt=ceil(abs(hwhm[0])/delw)
;
nea=n_elements(effar)
if nea eq ns then areff=effar else begin
  nwea=n_elements(wvlar)
  if nwea ne nea then begin
    areff=(interpol(effar,findgen(nea)*float(ns)/float(nea),findgen(ns))>(min(effar)))<(max(effar))
  endif else begin
    areff=(interpol(effar,wvlar,wmid)>(min(effar)))<(max(effar))
  endelse
endelse
;
nb=n_elements(ctb)
bct=0.*cts
if nb eq 1 then bct=bct+ctb[0]
if nb eq ns then bct=ctb
;
backscal=1.
if keyword_set(bkgscal) then backscal=bkgscal[0]
;
duration=1.D
if keyword_set(exptime) then duration=exptime[0]
;
inicon,fundae=ff

;this is the heuristic:
;	carry out a moat wavelet analysis (haartran) for a series of scales going from 1 bin to roofn(10*HWHM,2)
;	select a non-overlapping subset of the output (e.g., scales 1,2,4,8,16,etc), filter out negative values,
;	  and add them up to get cumulative wavelet signal
;
;	find all the local maxima in this cumulative signal (this is the seed set from which lines will be selected)
;	  these all have Quality=1
;
;	all local maxima for which adjacent bins have non-zero global threshold wavelet signals have Quality+=2
;
;	find the scaled background under the source in each bin
;	compute the Poisson likelihood of obtaining CTS counts or higher, given the scaled background
;	all local maxima for which p(\ge{CTS}|bg) < 0.1 have Quality+=4
;
;	compute SCALE, the scales at which the wavelet response is maximum, for each bin
;	compute minSCALE, the first scale at which wavelet response reaches a maximum, for each bin
;
;	find all local minima of SCALE along WVL
;	all local maxima that match these SCALE local minima have Quality+=8
;
;	find all local minima of minSCALE along WVL
;	all local maxima that match these minSCALE local minima have Quality+=16

;	set the scales
;maxscal=roofn(long(10.*lsfwdt),2) & scales=lindgen(maxscal)+1L & isc=2L^lindgen(long(alog(maxscal)/alog(2))+1L)
minscal=(lsfwdt/2)>1 & maxscal=roofn(long(10.*lsfwdt),2)
scales=lindgen(maxscal-minscal+1L)+minscal
isc=2L^lindgen((long(alog(maxscal)/alog(2)+0.5)-long(alog(minscal)/alog(2)+0.5))+1)*2L^long(alog(minscal)/alog(2)+0.5)

;	compute wavelet coefficients
hy=haartran(cts,scales,ss,ysig=fxe,wy=wy,bwy=bwy,verbose=vv)

;	compute cumulative +ve wavelet signal
iwy=wy[isc,*]>0 & ciwy=total(iwy,1)
ibwy=bwy[isc,*] & cibwy=total(ibwy,1)

;	compute cumulative global threshold filtered wavelet signal
ihy=hy[isc,*]>0 & cihy=total(ihy,1)

;	find local maxima
tmp=ciwy
if lsfwdt gt 1 then begin & smoothcts=smooth(float(cts),lsfwdt,/edge_wrap) & tmp=smoothcts+ciwy & endif
;omx=where(ciwy gt shift(ciwy,1) and ciwy gt shift(ciwy,-1),momx)
omx=where(tmp gt shift(tmp,1) and tmp gt shift(tmp,-1),momx)
qual=lonarr(ns) & qual[omx]=1L & ilocmax=qual

;	which of these are "good maxima"?
;omxp=omx+1L & omxn=omx-1L & thy=cihy[omxp]+cihy[omxn] & omxg=where(thy gt 0 and cihy gt 0 and ilocmax eq 1,momxg)
thy=shift(cihy,1)+shift(cihy,-1) & omxg=where(thy gt 0 and cihy gt 0 and ilocmax eq 1,momxg)
if momxg gt 0 then qual[omxg]=qual[omxg]+2L

;	estimate p(\ge{CTS}|BCT)
pBG=dblarr(ns)+1.
bct=bct/backscal
for i=0L,ns-1L do begin
  if bct[i] gt 0 then begin
    if cts[i] gt 0 then begin
      kk=lindgen(cts[i])
      lnpoi=kk*alog(bct[i])-bct[i]-lngamma(kk+1.D)
      mxlnpoi=max(lnpoi) & lnpoi=lnpoi-mxlnpoi & pBG[i]=1.D - total(exp(lnpoi))*exp(mxlnpoi)
    endif else pBG[i]=1.0
  endif
endfor
omxgp=where(pBG lt 0.1 and ilocmax ge 1,momxgp) & if momxgp gt 0 then qual[omxgp]=qual[omxgp]+4L

;	find BestSCALE and FirstSCALE for each bin
sclbest=lonarr(ns)+maxscal & sclfirst=lonarr(ns)+maxscal & wmxbest=dblarr(ns) & wmxfirst=dblarr(ns)
for i=0L,ns-1L do begin
  hc=hy[1:*,i] & wmxbest[i]=max(hc,imx) & if wmxbest[i] gt 0 then sclbest[i]=ss[imx]
  wmxfirst[i]=wmxbest[i] & sclfirst[i]=sclbest[i]	;defaults
  dhc=hc[1:*]-hc & oo=where(dhc lt 0,moo)
  if moo gt 2 then begin & if oo[1]-oo[0] eq 1 and oo[2]-oo[1] eq 1 then sclfirst[i]=ss[oo[0]+1] & endif
  o0=where(hc eq 0,mo0)
  if mo0 gt 0 then begin
    if o0[0] gt 0 then begin
      hc2=hc[0:o0[0]] & wmxfirst[i]=max(hc2,imx) & sclfirst[i]=ss[imx]
    endif else if hy[0,i] eq 0 then sclfirst[i]=max(ss)
  endif
endfor

;	flag bins which have low BestSCALE and FirstSCALE
;omxgpb=where((sclbest eq 1 or sclbest eq 2) and ilocmax ge 1,momxgpb) & if momxgpb gt 0 then qual[omxgpb]=qual[omxgpb]+8L
;omxgpb1=where((sclfirst eq 1 or sclfirst eq 2) and ilocmax ge 1,momxgpb1) & if momxgpb1 gt 0 then qual[omxgpb1]=qual[omxgpb1]+16L
omxgpb=where((sclbest le 2*lsfwdt) and ilocmax ge 1,momxgpb) & if momxgpb gt 0 then qual[omxgpb]=qual[omxgpb]+8L
omxgpb1=where((sclfirst le 2*lsfwdt) and ilocmax ge 1,momxgpb1) & if momxgpb1 gt 0 then qual[omxgpb1]=qual[omxgpb1]+16L

;to compute flux:
;	at each line location, collect the counts within ±2*hwhm
;	estimate a continuum based on BWY[maxscale,*]
;	compute net counts in the line, divide by EA and exposure time to get flux

ok=where(qual eq 31,mok)
if mok eq 0 then begin
  message,'No lines found',/informational
  return,create_struct('WLOC',-1,'eWLOC',-1,'OBSCT',-1,'ESTBG',-1,'pvalBG',-1,'ESTCONT',-1,'NETCT',-1,'NETCTERR',-1,'SNR',-1,'FLUX',-1,'FLUXERR',-1,'WIDTH',-1,'SCALE',-1,'QUALITY',-1)
endif

;	catch surrounding bins and merge any instances of 
iok=lonarr(ns)
for i=0L,mok-1L do begin
  j=ok[i] & j0=(j-2*lsfwdt)>0 & j1=(j+2*lsfwdt)<(ns-1L)
  iok[j0:j1]=1
endfor
o1=where(iok eq 1,mo1) & kok=lonarr(ns) & k=1L
if o1[0] eq 0 then kok[o1[0]]=1
for i=1L,mo1-1L do begin
  j=o1[i] & if iok[j-1L] eq 0 then k=k+1L & kok[j]=k
endfor
nlin=max(kok) & srcct=lonarr(nlin) & 
wloc=fltarr(nlin) & ewloc=wloc
obsct=lonarr(nlin) & estbg=fltarr(nlin) & pvalBG=dblarr(nlin)
estcont=fltarr(nlin)
netct=fltarr(nlin) & enetct=netct & snr=netct & flux=netct & eflux=netct
nbins=lonarr(nlin) & detqual=bytarr(nlin) & detscal=lonarr(nlin)
for il=0L,nlin-1L do begin
  ol=where(kok eq il+1L,mo1)
  wloc[il]=total(wmid[ol]*cts[ol])/total(cts[ol])
  ewloc[il]=sqrt(total(wmid[ol]^2*cts[ol])/total(cts[ol])-wloc[il]^2)
  obsct[il]=total(cts[ol])
  estbg[il]=total(bct[ol])
    kk=lindgen(obsct[il]) & lnpoi=kk*alog(estbg[il])-estbg[il]-lngamma(kk+1.D)
    mxlnpoi=max(lnpoi) & lnpoi=lnpoi-mxlnpoi
  pvalBG[il]=1.D - total(exp(lnpoi))*exp(mxlnpoi)
  estcont[il]=total(bwy[maxscal-minscal+1,ol])
  netct[il]=obsct[il]-estcont[il]
  enetct[il]=sqrt(obsct[il]+estcont[il])
  snr[il]=netct[il]/enetct[il]
  flux[il]=(ff.h*ff.c*1e8/wloc[il])*netct[il]/mean(areff[ol])/duration
  eflux[il]=(ff.h*ff.c*1e8/wloc[il])*enetct[il]/mean(areff[ol])/duration
  nbins[il]=mo1
  detscal[il]=min(sclfirst[ol])
  detqual[il]=max(qual[ol])
endfor
lines=create_struct('WLOC',wloc,'eWLOC',ewloc,'OBSCT',obsct,'ESTBG',estbg,'pvalBG',pvalBG,'ESTCONT',estcont,'NETCT',netct,'NETCTERR',enetct,'SNR',snr,'FLUX',flux,'FLUXERR',eflux,'WIDTH',nbins,'SCALE',detscal,'QUALITY',detqual)

if vv gt 1000 then begin
  help,name='momx*'
  ;DEBUG plot,wmid,cts & oplot,wmid,ciwy,col=2,thick=2 & oplot,wmid,cihy,col=3,thick=2 & oplot,wmid[ok],cts[ok],psym=-7,thick=2 & oplot,wmid,cibwy,col=4,thick=2
  ;DEBUG oo=where(wmid ge 11.5 and wmid le 11.6,moo) & for i=0,moo-1 do if qual[oo[i]] eq 31 then print,oo[i],wmid[oo[i]],cts[oo[i]],qual[oo[i]],thy[oo[i]],ciwy[oo[i]],cihy[oo[i]],sclbest[oo[i]],sclfirst[oo[i]]
  stop,'halting; type .CON to continue'
endif

return,lines
end
