function colourcurve,NH,arffil,rmffil,logT,Sband=Sband,Hband=Hband,$
	n_e=n_e,logP=logP,abund=abund,$
	ldbdir=ldbdir,ceroot=ceroot,cdbdir=cdbdir,ioneqf=ioneqf,$
	lstr=lstr,cstr=cstr,verbose=verbose, _extra=e
;+
;function	colourcurve
;	computes color (log(S)-log(H)) for given NH, for various temperatures,
;	for a given ARF/RMF combo
;
;syntax
;	colcurv=colourcurve(NH,arffil,rmffil,logT,Sband=Sband,Hband=Hband,$
;	n_e=n_e,logP=logP,abund=abund,$
;	ldbdir=ldbdir,ceroot=ceroot,cdbdir=cdbdir,ioneqf=ioneqf,$
;	lstr=lstr,cstr=cstr,verbose=verbose,$
;	chidir=chidir, cocofil=cocofil,chdwvl=chdwvl,/twoph,/noff,/nofb,$
;	nhne=nhne, fH2=fH2,He1=He1,HeII=HeII,/Fano,/wam,/bam,/mam,/vion,/noHeH,$
;	icrstr=icrstr,ionfrac=ionfrac, vfkydir=vfkydir,,nconsec=nconsec,$
;	EBV=EBV,R_V=R_V,LMC2=LMC2,AVGLMC=AVGLMC,ExtCurve=ExtCurve,$
;	fmgamma=fmgamma,fmx0=fmx0,fmc1=fmc1,fmc2=fmc2,fmc3=fmc3,fmc4=fmc4)
;
;parameters
;	NH	[INPUT; required] H I column density in [cm^-2]
;	arffil	[INPUT; required] name of ARF file
;		* if blank, ignores RMFFIL and proceeds to compute
;		  flux ratios rather than counts ratios
;		  - hardcoded range of [15 kev, 200 Ang], edit the code, or
;		    keep SBAND and HBAND in this range
;	rmffil	[INPUT; required] name of RMF file
;		* quits if ARFFIL and RMFFIL are set but files cannot be read
;	logT	[INPUT] temperatures at which to compute colors [log10([K])]
;		* if not given, or if input is not understandable, assumed
;		  to be 5:8:0.1
;
;keywords
;	Sband	[INPUT] 2-element float array of energies for the soft band [keV]
;		* default is [0.5,2]
;		* if 1-element, then alters the nearest bound.  e.g., if
;		  Sband=0.3, changes to [0.3,2]
;		  Sband=1.24, changes to [1.24,2]
;		  Sband=1.26, changes to [0.5,1.26]
;		  Sband=4, changes to [0.5,4]
;		* if integer, assumed to be channel number (cf. RMFFIL)
;	Hband	[INPUT] 2-element float array of energies for the hard band [keV]
;		* default is [2.0,7.0]
;		* if 1-element, then alters the nearest bound (cf. SBAND)
;		* if integer, assumed to be channel number (cf. RMFFIL)
;	n_e	[INPUT] passed w/o comment to RD_LINE() and RD_CONT()		  
;	logP	[INPUT] passed w/o comment to RD_LINE() and RD_CONT()		  
;	abund	[INPUT] abundances
;		* default is to use solar coronal
;	ldbdir	[INPUT] line emissivity database
;		* default is '$APED'
;	ceroot	[INPUT] what type of continuum to use
;		* default is 'apec'
;		* options are 'cie' and 'chianti'
;		  (but latter has problems; still being investigated)
;	cdbdir	[INPUT] continuum emissivity database
;		* default is '$CONT'
;	ioneqf	[INPUT] ion balance file
;		* default is 'ioneq/chianti.ioneq'
;		* not used for APED/ATOMDB
;	lstr	[I/O] structure that contains line emissivity database
;	cstr	[I/O] structure that contains continuum emissivity database
;		* if LSTR and CSTR are given on input, then the emissivity
;		  databases are not read in
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		[RD_IONEQ] CHIDIR
;		[RD_CONT] COCOFIL, CHDWVL, TWOPH, NOFF, NOFB
;		[ISMTAU] FH2, HE1, HEII, FANO, VION, BAM, MAM, NOHEH, ICRSTR
;		[IONABS] ICRSTR, IONFRAC, VFKYDIR, NCONSEC
;		[FM_UNRED] EBV, R_V, LMC2, AVGLMC, EXTCURVE, FMGAMMA, FMX0,
;		           FMC1, FMC2, FMC3, FMC4
;		[LINEFLX] NHNE
;-

;	usage
ok='ok' & np=n_params() & nnh=n_elements(NH) & nT=n_elements(logT)
nea=n_elements(arffil) & sza=size(arffil,/type)
nrm=n_elements(rmffil) & szr=size(rmffil,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if nnh eq 0 then ok='NH is not defined' else $
  if nnh gt 1 then ok='NH must be scalar' else $
   if nea eq 0 then ok='ARFFIL is not defined' else $
    if nrm eq 0 then ok='RMFFIL is not defined' else $
     if sza ne 7 then ok='ARFFIL must be a filename' else $
      if szr ne 7 then ok='RMFFIL must be a filename' else $
       if nea gt 1 then ok='ARFFIL must be a scalar' else $
        if nrm gt 1 then ok='RMFFIL must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: colcurv=colourcurve(NH,arffil,rmffil,logT,Sband=Sband,Hband=Hband,$'
  print,'       n_e=n_e,logP=logP,abund=abund,$'
  print,'       ldbdir=ldbdir,ceroot=ceroot,cdbdir=cdbdir,ioneqf=ioneqf,$'
  print,'       lstr=lstr,cstr=cstr,verbose=verbose,$'
  print,'       chidir=chidir, cocofil=cocofil,chdwvl=chdwvl,/twoph,/noff,/nofb,$'
  print,'       nhne=nhne, fH2=fH2,He1=He1,HeII=HeII,/Fano,/wam,/bam,/mam,/vion,/noHeH,$'
  print,'       icrstr=icrstr,ionfrac=ionfrac, vfkydir=vfkydir,,nconsec=nconsec,$'
  print,'       EBV=EBV,R_V=R_V,LMC2=LMC2,AVGLMC=AVGLMC,ExtCurve=ExtCurve,$'
  print,'       fmgamma=fmgamma,fmx0=fmx0,fmc1=fmc1,fmc2=fmc2,fmc3=fmc3,fmc4=fmc4)'
  print,'  computes color (log(S)-log(H)) for given NH, for various temperatures'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	initialize
inicon,fundae=ff
EM0=1.D

;	inputs
tlog=findgen(31)*0.1+5
if nt ne 0 then begin	;(NT>0
  if nt eq 1 then begin	;(NT=1
    if size(logt,/type) lt 4 then begin	;(LOGT is integer
      nn=logT[0] & delT=3./float(nn) & tlog=findgen(nn)*delT+5.
    endif else begin			;integer)(float
      if size(logt,/type) ne 7 then tlog=[logT[0]]
    endelse				;LOGT is float)
  endif else tlog=logT	;NT=1)
endif			;NT>0)
nT=n_elements(tlog)
;
if strtrim(arffil[0],2) ne '' then begin	;(arffil is not ''
fil=file_search(arffil[0],count=nfil)
  if nfil eq 0 then begin
    message,arffil+': file not found; quitting',/informational
    return,-1L
  endif
  arf=rdarf(arffil[0],arstr)
  elo=arstr.ELO & ehi=arstr.EHI & whi=ff.kevang/elo & wlo=ff.kevang/ehi & ww=0.5*(wlo+whi) & os=sort(ww)
  ww=ww[os] & wlo=wlo[os] & whi=whi[os] & elo=elo[os] & ehi=ehi[os] & arf=arf[os]
  wgrid=[wlo,max(whi)] & wgridmid=0.5*(wgrid[1:*]+wgrid) & egrid=ff.kevang/wgrid
  ;
  fil=file_search(rmffil[0],count=nfil)
  if nfil eq 0 then begin
    message,rmffil+': file not found; quitting',/informational
    return,-1L
  endif
  rmstr=rd_ogip_rmf(rmffil[0],effar=rspea,rmfimg=rmfimg)
  nchan=n_elements(rmstr.EMN)
  emn=rmstr.EMN & emx=rmstr.EMX & wmx=ff.kevang/emn & wmn=ff.kevang/emx
endif else begin				;arffil is not '')(is blank
  ;	compute flux-ratios rather than counts-ratios
  ;	make up numbers for spectral grid, etc.
  ;	-- not very important what it is as long as it covers range of interest
  ;	   and if it doesn't edit this line (someday will become a keyword)
  wmin=ff.kevang/15. & wmax=200. & nwbin=16384L & dwbin=(wmax-wmin)/float(nwbin-1L) & wgrid=findgen(nwbin)*dwbin+wmin
  wgridmid=(wgrid[1:*]+wgrid)/2.
  arf=0.*wgridmid+1.
  egrid=ff.kevang/wgrid
  ee=ff.kevang/wgrid & emn=ee[0L:nwbin-2L] & emx=ee[1:*]
  nchan=nwbin
endelse						;arffil is '')

;	keywords
;
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if n_elements(abund) ne 30 then abund=getabund('coronal')
if not keyword_set(ldbdir) then ldbdir='$APED'
if not keyword_set(ceroot) then ceroot='apec'
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(ioneqf) then ioneqf='ioneq/chianti.ioneq'
rdlem=1 & if n_tags(lstr) gt 0 then rdlem=0
rdcem=1 & if n_tags(cstr) gt 0 then rdcem=0
;
bandS=[0.5,2.0] & ibandS=[0,nchan/2-1L]
tmp=min(abs(bandS[0]-emn),imn) & ibandS[0]=imn
tmp=min(abs(bandS[1]-emx),imx) & ibandS[1]=imx
mS=n_elements(Sband) & szS=size(Sband,/type)
;
if mS eq 1 then begin	;(Sband is a scalar
  if szS lt 4 then begin
    xx=Sband[0] & dx0=abs(xx-ibandS[0]) & dx1=abs(xx-ibandS[1])
    if dx0 le dx1 then ibandS[0]=xx else ibandS[1]=xx	;update defaults
  endif else begin
    xx=Sband[0] & dx0=abs(xx-bandS[0]) & dx1=abs(xx-bandS[1])
    if dx0 le dx1 then bandS[0]=xx else bandS[1]=xx	;update defaults
  endelse
endif			;mS=1)
if mS eq 2 then begin
  if szS lt 4 then begin
    ibandS=Sband
    smin=min(ibandS,max=smax) & ibandS=[smin,smax]
  endif else begin
    bandS=Sband
    smin=min(bandS,max=smax) & bandS=[smin,smax]
  endelse
endif
if mS gt 2 then begin
  message,'Sband cannot be of size > 2; quitting',/informational
  return,-1L
endif
if szS ge 4 then begin	;(convert to channels if in energy
  tmp=min(abs(bandS[0]-emn),imn)
  omn=where(bandS[0] ge emn,momn) & if momn gt 0 then ibandS[0]=omn[0] else ibandS[0]=imn
  tmp=min(abs(bandS[1]-emx),imx)
  omx=where(bandS[1] le emx,momx) & if momx gt 0 then ibandS[1]=omx[0] else ibandS[1]=imx
endif			;szS.ge.4)
;
bandH=[2.0,7.0] & ibandH=[nchan/2,nchan-1L]
tmp=min(abs(bandH[0]-emn),imn) & ibandH[0]=imn
tmp=min(abs(bandH[1]-emx),imx) & ibandH[1]=imx
mH=n_elements(Hband) & szH=size(Hband,/type)
;
if mH eq 1 then begin	;(Hband is a scalar
  if szH lt 4 then begin
    xx=Hband[0] & dx0=abs(xx-ibandH[0]) & dx1=abs(xx-ibandH[1])
    if dx0 le dx1 then ibandH[0]=xx else ibandH[1]=xx	;update defaults
  endif else begin
    xx=Hband[0] & dx0=abs(xx-bandH[0]) & dx1=abs(xx-bandH[1])
    if dx0 le dx1 then bandH[0]=xx else bandH[1]=xx	;update defaults
  endelse
endif			;mH=1)
if mH eq 2 then begin
  if szH lt 4 then begin
    ibandH=Hband
    Hmin=min(ibandH,max=Hmax) & ibandH=[Hmin,Hmax]
  endif else begin
    bandH=Hband
    Hmin=min(bandH,max=Hmax) & bandH=[Hmin,Hmax]
  endelse
endif
if mH gt 2 then begin
  message,'Hband cannot be of size > 2; quitting',/informational
  return,-1L
endif
if szH ge 4 then begin	;(convert to channels if in energy
  tmp=min(abs(bandH[0]-emn),imn)
  omn=where(bandH[0] ge emn,momn) & if momn gt 0 then ibandH[0]=omn[0] else ibandH[0]=imn
  tmp=min(abs(bandH[1]-emx),imx)
  omx=where(bandS[1] le emx,momx) & if momx gt 0 then ibandH[1]=omx[0] else ibandH[1]=imx
endif			;szH.ge.4)
;
emin=emn[ibandS[0]] & emax=emx[ibandH[1]]

;	compute optical depths
tau=ismtau(wgridmid,NH=NH[0], _extra=e) & trans=exp(-tau)

;	read in the line emissivity database
if keyword_set(rdlem) then begin
  lemis=rd_line(atom,logP=logP,n_e=n_e,wrange=ff.kevang/[emax,emin],dbdir=ldbdir,fstr=lstr,verbose=vv)
  cc=strupcase(ldbdir)
  ;	add ion balance if not ATOMDB, and remove abundances if ATOMDB
  if strpos(cc,'APED',0) lt 0 and strpos(cc,'ATOMDB',0) lt 0 and strpos(cc,'APEC',0) lt 0 then $
    lemis=fold_ioneq(ff,lstr.Z,lstr.JON,logt=lstr.LOGT,verbose=vv) else $
    apedance,lemis,lstr.Z
  lstr.LINE_INT=lemis
endif
lwvl=abs(lstr.WVL)
nlT=n_elements(lstr.LOGT)
lnrg=ff.h*ff.c*1e8/lwvl

;	read in the continuum emissivity database
if keyword_set(rdcem) then begin
  cemis=rd_cont(ceroot,logP=logP,n_e=n_e,wrange=ff.kevang/[emax,emin],dbdir=cdbdir,fcstr=cstr,verbose=vv,abund=abund, _extra=e)
endif
cwvl=cstr.midWVL
;this can be numerically unstable, use hack of next line instead --> cww=cstr.WVL & cdw=abs(cww[1:*]-cww)
cwb=mid2bound(cwvl) & cdw=abs(cwb[1:*]-cwb)
ncT=n_elements(cstr.LOGT)
cnrg=ff.h*ff.c*1e8/cwvl

;	correct line emissivities for abundances
labund=abund[lstr.Z-1]

;	at each temperature,
;	compute fluxes, make spectrum, convolve with RMF, get counts in S and H bands, compute ratio
colcurv=fltarr(nT)
for iT=0L,nT-1L do begin

  dT=min(abs(tlog[iT]-lstr.LOGT),ioT)
  if dT lt 0.01 then begin	;(dT<0.01
    lemis=reform((lstr.LINE_INT)[ioT,*])
  endif else begin		;dT<0.01)(dt>0.01
    message,'this feature has not been implemented; '+$
    	'can only use input LOGT that matches the logT grid in databases',$
	/informational
    return,-1L
    ;if ioT gt 0 and ioT lt nlT-1L then begin
    ;  o0=where(lstr.logT le tlog[iT],mo0) & o1=where(lstr.LOGT gt tlog[iT],mo1)
    ;endif else begin
    ;endelse
  endelse			;dT>0.01)

  dT=min(abs(tlog[iT]-cstr.LOGT),ioT)
  if dT lt 0.01 then begin	;(dT<0.01
    cemis=reform((cstr.CONT_INT)[ioT,*])
    cemis=cemis*cdw	;[1e23 ergs cm^3/s/A]*[A]
  endif else begin		;dT<0.01)(dt>0.01
    message,'this feature has not been implemented; '+$
    	'can only use input LOGT that matches the logT grid in databases',$
	/informational
    return,-1L
    ;if ioT gt 0 and ioT lt nlT-1L then begin
    ;  o0=where(cstr.logT le tlog[iT],mo0) & o1=where(cstr.LOGT gt tlog[iT],mo1)
    ;endif else begin
    ;endelse
  endelse			;dT>0.01)

  ;	line and continuum fluxes
  lflx=lemis*EM0*labund/lnrg/1d23 	;[1e23 ergs cm^3/s]*[cm^-5]/[erg/ph]/1e23 = [1e23 ph/s/cm^2]
  cflx=cemis*EM0/cnrg/1d23 	;[1e23 ergs cm^3/s]*[cm^-5]/[erg/ph]/1e23 = [1e23 ph/s/cm^2]

  ;	make spectra
  lspec=hastogram(lwvl,wgrid,wts=lflx)	;[ph/s/cm^2/bin]
  cspec=rebinw(cflx,cwvl,wgrid,/perbin)	;[ph/s/cm^2/bin]

  ;	combine and include effective area and optical depth transmission factor
  spec=(lspec+cspec)*arf*trans	;[ct/s/bin]

  ;	push through RMF
  if n_tags(rmstr) gt 0 then begin
    conv_rmf,egrid,spec,chan,ctspec,rmstr
  endif else begin
    chan=egrid & ctspec=spec
  endelse

  ;	counts in S and H bands
  ctS=total(ctspec[ibandS[0]:ibandS[1]])
  ctH=total(ctspec[ibandH[0]:ibandH[1]])

  ;	hardness ratio
  colcurv[iT]=alog10(ctS)-alog10(ctH)

  ;	debug
  if vv gt 5 then print,iT,tlog[iT],ctS,ctH,colcurv[iT]

endfor

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,colcurv
end
