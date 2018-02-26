function thermal_emis,wvlar,effar,tgrid,ldbdir=ldbdir,cdbdir=cdbdir,$
	ioneqf=ioneqf,abund=abund,cieroot=cieroot,n_e=n_e,NH=NH,$
	lstr=lstr,cstr=cstr,EM0=EM0,$
	powl=powl,powc=powc,nrgt=nrgt,pht=pht,verbose=verbose,$
	_extra=e
;+
;function	thermal_emis
;	compute and return the thermal response over the passband
;	for various temperatures given the instrument effective areas,
;	in units of [ct/s]
;
;syntax
;	pow=thermal_emis(wvlar,effar,tgrid,ldbdir=ldbdir,$
;	cdbdir=cdbdir,ioneqf=ioneqf,abund=abund,cieroot=cieroot,$
;	n_e=n_e,NH=NH,lstr=lstr,cstr=cstr,EM0=EM0,verbose=verbose,$
;	chidir=chidir,fH2=fH2,He1=He1,HeII=HeII)
;
;parameters
;	wvlar	[INPUT; required] wavelength grid for effective areas
;	effar	[INPUT] effective area at each WVLAR in units of [cm^2 ct/ph]
;		* if not given, assumed to be unity
;		* if given
;		  -- and is scalar or single element, is expanded at same
;		     value to match the size of WVLAR
;		  -- otherwise must match size of WVLAR
;	tgrid	[I/O] the log(T[K]) grid over which the output
;		is to be defined
;		* if not specified on input, will be set to be
;		  the same as that in the emissivity tables,
;		  the PoA default, findgen(81)*0.05+4
;
;keywords
;	ldbdir	[INPUT; '$CHIANTI'] line emissivity database directory
;	cdbdir	[INPUT; '$CONT'] continuum emissivity database directory
;	ioneqf	[INPUT; 'ioneq/chianti.ioneq'] location of the
;		ion fraction tables
;	abund	[INPUT; getabund('Grevesse et al')] abundances
;	cieroot	[INPUT; 'cie'] root name for the files in $CONT
;	n_e	[INPUT; 1e9] electron number density in the plasma
;	NH	[INPUT; 0] H column density to apply
;	lstr	[OUTPUT] line emissivities structure, as read in
;		from RD_LINE()
;	cstr	[OUTPUT] continuum emissivities structure, as read in
;		from RD_CONT()
;	EM0	[INPUT; default=1d10] multiples the emissivities by
;		this value to derive [erg/s] or [erg/cm^2/s] depending
;		on whether the input has units [cm^-3] or [cm^-5]
;		* must be a scalar, otherwise the total of an array is used
;	powl	[OUTPUT] the response due solely to lines
;	powc	[OUTPUT] the response due solely to continuum
;		* both POWL and POWC are identical in size and units
;		  to the primary output
;	nrgt	[OUTPUT] the total energy in [erg/s/cm^2] over the same
;		energy grid (i.e., excluding conversion to [ph] and
;		multiplication by EFFAR)
;	pht	[OUTPUT] the total photon counts in [ph/s/cm^2] over the
;		same energy grid (i.e., excluding multiplication by EFFAR)
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		FOLD_IONEQ : CHIDIR
;		ISMTAU : FH2, HE1, HEII
;
;history
;	vinay kashyap (May2017; based on hinode_xrt_emis.pro)
;-

;	usage
ok='ok' & np=n_params() & nw=n_elements(wvlar) & nea=n_elements(effar) & nt=n_elements(tgrid)
if np eq 0 then ok='Insufficient parameters' else $
 if nw eq 0 then ok='WVLAR is undefined'
if ok ne 'ok' then begin
  print,'Usage: pow=thermal_emis(wvlar,effar,tgrid,ldbdir=ldbdir,$
  print,'	cdbdir=cdbdir,ioneqf=ioneqf,abund=abund,cieroot=cieroot,$
  print,'	n_e=n_e,NH=NH,lstr=lstr,cstr=cstr,EM0=EM0,verbose=verbose,$
  print,'	chidir=chidir,fH2=fH2,He1=He1,HeII=HeII)
  print,'  compute and return power loss at different temperatures functions given
  if np gt 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
WAVE=wvlar
wmin=min(WAVE,max=wmax)
if nea eq 0 then EFF_AREA=0.*wvlar+1.0 else begin
  if nea eq 1 then EFF_AREA=0.*wvlar+effar[0] else begin
    if nea eq nw then EFF_AREA=effar else begin	;(could also have a "or nea eq nw-1" here, but later
      message,'EFFAR and WVLAR are incompatible; exiting',/informational
      return,-1L
    endelse					;NEA .ne. NW)
  endelse
endelse

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
if not keyword_set(ldbdir) then ldbdir='$CHIANTI'
if not keyword_set(ioneqf) then ioneqf='ioneq/chianti.ioneq'
if n_elements(abund) lt 30 then abund=getabund('grevesse et al.')
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(cieroot) then cieroot='cie'
if not keyword_set(n_e) then n_e=1e9
if not keyword_set(NH) then NH=0
if not keyword_set(EM0) then EM_0=1d10 else begin
  nEM=n_elements(EM0)
  if nEM eq 1 then EM_0=EM0[0] else EM_0=total(EM0,/nan)
endelse

;	some other useful stuff
if vv ge 5 then peasecolr
inicon,fundae=fundae
xtitle='log T' & ytitle='[ct s!u-1!n]'

;	read in line emissivities
if n_tags(lstr) eq 0 then begin
  lconf=rd_line(atom,n_e=n_e,wrange=[wmin,wmax],dbdir=ldbdir,fstr=lstr,verbose=verbose)
  ;	if APED, then remove the Anders & Grevesse abundances
  if strpos(strlowcase(ldbdir),'ape') ge 0 then begin
    apedance,lconf,lstr.Z & lstr.LINE_INT=lconf
  endif
  ;	if not APED, apply ion balances
  if strpos(strlowcase(ldbdir),'ape') lt 0 and strpos(strlowcase(ioneqf),'none') lt 0 then begin
    lconf=fold_ioneq(lstr.LINE_INT,lstr.Z,lstr.JON,eqfile=ioneqf,verbose=verbose, _extra=e)
    lstr.LINE_INT=lconf
  endif
endif
lwvl=abs(lstr.WVL) & Z=lstr.Z & nlw=n_elements(Z)
nrgl=(fundae.h*fundae.c*1e8)/lwvl	;[erg/ph]

;	read in continuum emissivities
if n_tags(cstr) eq 0 then begin
  cconf=rd_cont(cieroot,n_e=n_e,wrange=[wmin,wmax],dbdir=cdbdir,abund=abund,fcstr=cstr,verbose=verbose,_extra=e)
endif
cww=mid2bound(cstr.midWVL) & cwvl=0.5*(cww[1:*]+cww) & ncw=n_elements(cwvl)
nrgc=(fundae.h*fundae.c*1e8)/cwvl	;[erg/ph]
cdw=abs(cww[1:*]-cww)
cconf=cstr.CONT_INT & for i=0L,n_elements(cstr.LOGT)-1L do cconf[i,*]=cconf[i,*]*cdw
;	recast the continuum on the same temperature grid as the line emissivities
if n_elements(cstr.LOGT) ne n_elements(lstr.LOGT) then begin
  tmp=rebinx(cconf,cstr.LOGT,lstr.LOGT,verbose=vv)
  if vv gt 9 then stop,'HALTing; type .CON to continue'
  cconf=tmp
endif

;	any intrinsic absorption?
if keyword_set(NH) then begin
  lopt=exp(-ismtau(lwvl,NH=NH,/bam,abund=abund,verbose=(vv<1),_extra=e))
  copt=exp(-ismtau(cwvl,NH=NH,/bam,abund=abund,verbose=(vv<1),_extra=e))
endif else begin
  lopt=fltarr(nlw)+1.
  copt=fltarr(ncw)+1.
endelse
;lopt=lopt*(1d-23/nrgl)
;copt=copt*(1d-23/nrgc)

;	output
if n_elements(tgrid) eq 0 then tgrid=lstr.LOGT
nT=n_elements(tgrid) & nmaxT=n_elements(lstr.LOGT) & jT=findex(lstr.LOGT,tgrid)
powt=dblarr(nT) & powl=dblarr(nT) & powc=dblarr(nT)
nrgt=dblarr(nT) & pht=dblarr(nT)
warf=WAVE[0:nW-1L,0] & larf=fltarr(nlw) & carf=fltarr(ncw)
larf=(interpol(EFF_AREA,warf,lwvl)>0)<(max(EFF_AREA))
carf=(interpol(EFF_AREA,warf,cwvl)>0)<(max(EFF_AREA))

;	compute the response
zab=abund[Z-1]
for it=0L,nT-1L do begin
  if vv gt 0 then kilroy
  jj=round(jT[iT])
  jj0=(jj>0)<(nmaxT-1L) & jj1=(jj+1)<(nmaxT-1L)	;no extrapolations!
  djT=abs(jT[iT]-jj)
  ;
  ;	intensities
  linint_erg=1d-23*lopt*reform((lstr.LINE_INT)[jj0,*])*EM_0*zab	;[erg/s/cm^2]
  conint_erg=1d-23*copt*reform((cconf)[jj0,*])*EM_0		;[erg/s/cm^2]
  ;	interpolate if necessary
  if djT gt 0.001 then begin
    if jj ge 0 and jj lt nmaxT-2 then begin
      linint_erg = djT*linint_erg[jj0,*] + (1.-djT)*linint_erg[jj1,*]
      conint_erg = djT*conint_erg[jj0,*] + (1.-djT)*conint_erg[jj1,*]
    endif
    kilroy,dot='*'
  endif

  ;	energy fluxes
  nrgt[it]=total(linint_erg,/nan)+total(conint_erg,/nan)

  ;	photon fluxes
  linint_ph=linint_erg/nrgl & conint_ph=conint_erg/nrgc	;[ph/s/cm^2]
  pht[it]=total(linint_ph,/nan)+total(conint_ph,/nan)

  ;	count rates
  linint=larf*linint_ph & conint=carf*conint_ph	;[ct/s]
  powl[it]=total(linint,/nan)
  powc[it]=total(conint,/nan)
  powt[it]=powl[it]+powc[it]

  ;linint=lopt*larf*reform((lstr.LINE_INT)[jj0,*])*EM_0*zab
  ;conint=copt*carf*reform((cconf)[jj0,*])*EM_0
  ;if djT gt 0.001 and $
  ;   jj ge 0 and jj lt nmaxT-2 then begin	;(interpolate between slices
  ;  linint=djT*linint + (1.-djT)*lopt*reform((lstr.LINE_INT)[jj1,*])*EM_0*zab
  ;  conint=djT*conint + (1.-djT)*copt*reform((cconf)[jj1,*])*EM_0
  ;endif						;JT>0.001)

  if it gt 0 and vv ge 5 then begin
    plot,tgrid[0:it],powt[0:it],/yl,yr=max(powt[0:it])*[1e-5,2],$
	xtitle=xtitle,ytitle=ytitle,_extra=e
    oplot,tgrid[0:it],powc[0:it],line=1
    oplot,tgrid[0:it],powl[0:it],line=2
  endif
endfor

if vv ge 5 then begin
  plot,tgrid,powt,/yl,yr=max(powt)*[1e-5,2],$
    xtit=xtitle,ytit=ytitle,_extra=e
  oplot,tgrid,powc,line=1
  oplot,tgrid,powl,line=2
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,powt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example

print,getpoadef('ARDB')
print,getpoadef('LDBDIR')
print,getpoadef('CDBDIR')
print,getpoadef('CEROOT')
print,getpoadef('ABREF')
print,getpoadef('CHIDIR')
print,getpoadef('IONEQF')

tgrid=findgen(31)*0.05+6.

rd_pimms_file,get_pimms_file('chandra','acis-s'),effar,wvlar,/wave
powt=thermal_emis(wvlar,effar,tgrid,ldbdir=!LDBDIR,cdbdir=!CDBDIR,ioneqf=!IONEQF,abund=getabund(!ABREF),cieroot=!CEROOT,n_e=1e9,NH=NH,lstr=lstr,cstr=cstr,EM0=EM0,verbose=!verbose,chidir=!CHIDIR)

end
