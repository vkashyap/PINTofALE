function starflux,lstr,cstr,tlog=tlog,DEM=DEM,abund=abund,wrange=wrange,$
	NH=NH,effar=effar,wvlar=wvlar,dist=dist,radius=radius,noph=noph,$
	lflux=lflux,cflux=cflux, _extra=e, ikeV=ikeV
;+
;function	*flux
;	returns total flux in the passband of interest at the specified
;	temperatures OR for the given DEM
;
;syntax
;	flx=starflux(lstr,cstr,tlog=tlog,DEM=DEM,abund=abund,wrange=wrange,$
;	    NH=NH,effar=effar,wvlar=wvlar,dist=dist,radius=radius,/noph,$
;	    lflux=flux,cflux=cflux,fH2=fH2,He1=He1,HeII=HeII,/Fano)
;
;a word on the units
;	the line emissivities are in [ergs cm^3/s], and the continuum
;	emissivities are also converted to this.  the [ergs] are converted
;	to [ph] before output unless specified otherwise.
;	thus, output will have units:
;		[(ergs|ph) cm^3/s][EM][RADIUS^2][DIST^-2][AREA]
;	e.g., if EM is [cm^-3] at source, and DIST is given, then output
;	is [ph/cm^2/s], the flux at the telescope, which is fine.
;	on the other hand, if DEM is [cm^-5] (n^e*V/(4*pi*D^2) at the
;	telescope) and DIST are given, we'd get [ph/cm^4/s] which would
;	be meaningless.  instead, RADIUS may be given to get [ph/s].
;	*******************************************************************
;	IT IS UP TO THE USER TO MAKE SURE THAT THE OUTPUT UNITS MAKE SENSE!
;	*******************************************************************
;
;parameters
;	lstr	[INPUT; required] structure containing line emissivity data
;		(see RD_LINE or RD_LIST)
;		* BEWARE: ion balance is assumed to be included!
;	cstr	[INPUT; required] structure containing continuum emissivities
;		(see RD_CONT)
;
;keywords
;	tlog	[INPUT] EITHER [A] the temperatures at which to compute the
;		output fluxes XOR [B] the temperatures at which the DEM is
;		defined
;		* if n(DEM) < 2, then [A]
;		* if n(DEM) > 1 .AND. n(DEM) = n(TLOG), then [B]
;	DEM	[INPUT] EITHER [A] volume emission measure [cm^-3 or cm^-5]
;		XOR [B] differential emission measure [cm^-3/logK or cm^-5/logK]
;		* if [B] then n(DEM) .MUST=. {n(TLOG) or n(LSTR.LOGT)}
;		* default is 1e14
;	abund	[INPUT] abundances
;		* default: Grevesse et al. (1992)
;	wrange	[INPUT] passband in which to compute flux
;		* if not given, assumes larger of range(LSTR.WVL,CSTR.WVL)
;	NH	[INPUT] ISM column density [cm^-2]
;	effar	[INPUT] effective area [cm^2]
;	wvlar	[INPUT] wavelengths at which effective area is defined [Ang]
;		* n(EFFAR) .MUST=. n(WVLAR)
;		* if not given, or illegal, output will be [ph/s/cm^2]
;	dist	[INPUT] distance to star, in [pc]
;		* BUT, if > 1e6, taken to be in [cm]
;		* if given, divides computed flux by 4*pi*DIST^2 before output
;		* BEWARE: no cosmological corrections are made!
;	radius	[INPUT] radius of star, in [cm]
;		* BUT, if < 1e6, taken to be in [RSol]
;		* if given, multiplies flux by 4*pi*RADIUS^2 before output
;		* BEWARE: no scale-height/volume effects included!
;	noph	[INPUT] if set, does not convert [ergs] to [ph]
;	lflux	[OUTPUT] flux due to lines only
;	cflux	[OUTPUT] flux due to continuum only
;		* LFLUX and CFLUX will >NOT< have DIST and RADIUS included!!
;	_extra	[INPUT] pass defined keywords to
;		ISMTAU [fH2,He1,HeII,Fano]
;	ikeV	[JUNK] does nothing.  only to avoid passing it on to ISMTAU
;
;restrictions
;	requires subroutines ISMTAU, GETABUND, MK_DEM
;	does not call LINEFLX or GRATFLX, so perhaps a wee bit faster,
;	  but on the other hand, cannot handle grating effective areas either.
;
;history
;	vinay kashyap (1999May)
;	changed KEV to IKEV (VK; DecMM)
;-

;	usage
ok='ok'
ml=n_elements(lstr) & mc=n_elements(cstr)
nl=n_tags(lstr)     & nc=n_tags(cstr)
ltg=[''] & if nl gt 1 then ltg=strupcase(tag_names(lstr))
ctg=[''] & if nc gt 1 then ctg=strupcase(tag_names(cstr))
il1=where(ltg eq 'LINE_INT',ml1) & ic1=where(ctg eq 'CONT_INT',mc1)
il2=where(ltg eq 'LOGT',ml2)     & ic2=where(ctg eq 'LOGT',mc2)
il3=where(ltg eq 'WVL',ml3)      & ic3=where(ctg eq 'WVL',mc3)
il4=where(ltg eq 'Z',ml4)        & ic4=where(ctg eq 'MIDWVL',mc4)
if ml eq 0 then ok='Line emissivity structure must be supplied' else $
 if mc eq 0 then ok='Continuum emissivity structure must be supplied' else $
  if nl eq 0 then ok='LSTR not a structure' else $
   if nc eq 0 then ok='CSTR not a structure' else $
    if nl lt 5 then ok='LSTR not in correct format' else $
     if nc lt 5 then ok='CSTR not in correct format' else $
      if ml1 eq 0 then ok='LSTR.LINE_INT missing' else $
       if ml2 eq 0 then ok='LSTR.LOGT missing' else $
	if ml3 eq 0 then ok='LSTR.WVL missing' else $
	 if ml4 eq 0 then ok='LSTR.Z missing' else $
	  if mc1 eq 0 then ok='CSTR.CONT_INT missing' else $
	   if mc2 eq 0 then ok='CSTR.LOGT missing' else $
	    if mc3 eq 0 then ok='CSTR.WVL missing' else $
	     if mc4 eq 0 then ok='CSTR.MIDWVL missing'
if ok ne 'ok' then begin
  print,'Usage: flx=starflux(lstr,cstr,tlog=tlog,DEM=DEM,abund=abund,$'
  print,'       wrange=wrange,NH=NH,effar=effar,wvlar=wvlar,dist=dist,$'
  print,'       radius=radius,/noph,fH2=fH2,He1=He1,HeII=HeII,/Fano)'
  print,'  return fluxes at given temperature'
  if n_params() gt 0 then message,ok,/info
  return,0.
endif

;	check keywords
;
mt=n_elements(Tlog)
if mt eq 0 then logT=lstr.LOGT else logT=tlog & mt=n_elements(logT)
;
VEM=1			;begin by assuming VEM is given
md=n_elements(DEM)
if md eq 0 then DiffEM=1d14		;VEM
if md eq 1 then DiffEM=DEM(0)		;still VEM
if md gt 1 then begin
  if md eq mt or md eq n_elements(lstr.LOGT) then begin
    DiffEM=DEM & VEM=0			;not VEM
    if md eq mt then begin
      if mt ne n_elements(lstr.LOGT) then $
	DiffEM=mk_dem('interpolate',lstr.LOGT,pardem=logT,indem=DEM)
    endif
  endif else begin
    message,'Input DEM is not understandable; summing up and using as VEM',/info
    DiffEM=total(DEM)			;VEM
  endelse
endif
;
if n_elements(abund) lt 30 then abund=getabund('grevesse et al.')
;
if n_elements(wrange) ne 2 then begin
  wlmin=min(lstr.WVL,max=wlmax) & wcmin=min(cstr.WVL,max=wcmax)
  wmin = wlmin < wcmin
  wmax = wlmax < wcmax
  wr=[wmin,wmax]
endif else wr=wrange
if wr(1) lt wr(0) then wr=reverse(wr)
;
etaul=0.*lstr.WVL+1.
etauc=0.*cstr.midWVL+1.
;
areffl=0.*lstr.WVL+1. & arwvll=abs(lstr.WVL)
areffc=0.*cstr.midWVL+1. & arwvlc=abs(cstr.midWVL)
nea=n_elements(effar) & nwa=n_elements(wvlar) & ow=sort(wvlar)
;
if not keyword_set(noph) then begin
  h=6.626176e-27 & c=2.9979e10
endif

;	now figure out array size of output
if VEM eq 0 then nout=1L else nout=mt
;	figure out nearest values of the temperatures at which you
;	can calculate the fluxes
iT=lonarr(nout)-1L & tt=lstr.LOGT
nt=n_elements(tt) & mdt=median((tt(1:*)-tt))
for i=0L,nout-1L do begin
  tmp=min(abs(tt-logT(i)),imn) & iT(i)=imn
  if abs(tt(imn)-logT(i)) gt mdt/2. then begin
    ok='approximating logT='+strtrim(logT(i),2)+' by '+strtrim(tt(imn),2)
    message,ok,/info
  endif
endfor

;	filter on wavelength
owl=where(abs(lstr.WVL) ge wr(0) and abs(lstr.WVL) le wr(1),mowl)
owc=where(cstr.WVL ge wr(0) and cstr.WVL le wr(1),mowc)
if mowl gt 0 then begin
  wvll=abs(lstr.WVL) & zz=lstr.Z
  wvll=wvll(owl) & zz=zz(owl) & zab=abund(zz-1)
  if keyword_set(NH) then etaul=exp(-ismtau(wvll,NH=NH, _extra=e))
  if nea gt 0 and nea eq nwa then $
	areffl=interpol(effar(ow),wvlar(ow),wvll) > 0
  if not keyword_set(noph) then nrgl=h*c*1e8/abs(wvll)
endif
if mowc gt 0 then begin
  wvlc=cstr.midWVL & ww=cstr.WVL
  wvlc=wvlc(owc) & ww=ww(owc) & dw=ww(1:*)-ww
  if keyword_set(NH) then etauc=exp(-ismtau(wvlc,NH=NH, _extra=e))
  if nea gt 0 and nea eq nwa then $
	areffc=interpol(effar(ow),wvlar(ow),wvlc) > 0
  if not keyword_set(noph) then nrgc=h*c*1e8/abs(wvlc)
endif

;	outputs
flx=fltarr(nout) & lflux=flx & cflux=flx

;	populate the outputs
ot=lindgen(nt)
emisl=reform((lstr.LINE_INT)(*,owl))		;[1e-23 ergs cm^3/s]
emisc=reform((cstr.CONT_INT)(*,owc))		;[1e-23 ergs cm^3/s/A]
for i=0L,nout-1L do begin
  kilroy; was here.	;make a mark
  if VEM eq 1 then ot=iT(i)
  if mowl gt 0 then begin
    flxl=fltarr(mowl)
    for j=0L,mowl-1L do begin			;{for each line
      tmp=(1d-23*DiffEM)*emisl(ot,j)*zab(j)
      tmp([0,nt-1])=0.5*tmp([0,nt-1])
      flxl(j)=mdt*alog(10.)*total(tmp)		;[ergs cm^3/s][EM]
    endfor					;J=0,MOWL-1}
    flxl=flxl*areffl*etaul			;[ergs cm^3/s][EM][AREA]
    if not keyword_set(noph) then flxl=flxl/nrgl	;[ergs..]->[ph..]
    lflux(i)=total(flxl)
  endif
  if mowc gt 1 then begin
    flxc=fltarr(mowc)
    for j=0L,mowc-2L do begin			;{for each bin
      tmp=(1d-23*DiffEM)*emisc(ot,j)*dw(j)	;[ergs cm^3/s]
      tmp([0,nt-1])=0.5*tmp([0,nt-1])
      flxc(j)=mdt*alog(10.)*total(tmp)		;[ergs cm^3/s][EM]
    endfor					;J=0,MOWC-1}
    flxc=flxc*areffc*etauc			;[ergs cm^3/s][EM][AREA]
    if not keyword_set(noph) then flxc=flxc/nrgc	;[ergs..]->[ph..]
    cflux(i)=total(flxc)
  endif
endfor
;
flx=lflux+cflux			;[(ergs|ph) cm^3/s][EM][AREA]

;	distance corrections
if keyword_set(dist) then begin
  pc2cm=3.084d18 & if dist gt 1e6 then pc2cm=1.
  oo=where(flx gt 0,moo)
  if moo gt 0 then begin
    tmp=alog10(flx(oo))-alog10(4.*!dpi)-2.*alog10(dist)-2.*alog10(pc2cm)
    flx(oo)=10.^(tmp)		;[(ergs|ph) cm^3/s][EM][AREA][DIST^-2]
  endif
endif
if keyword_set(radius) then begin
  cm2rsun=1.D & if dist lt 1e6 then cm2rsun=6.969d10
  tmp=alog10(flx(oo))+alog10(4.*!dpi)+2.*alog10(radius*cm2rsun)
  flx=10.^(tmp)		;[(ergs|ph) cm^3/s][EM][AREA][RADIUS^2][DIST^-2]
endif

return,flx
end
