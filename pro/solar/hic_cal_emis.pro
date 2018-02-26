function hic_cal_emis,filter,tgrid,ldbdir=ldbdir,cdbdir=cdbdir,$
	ioneqf=ioneqf,abund=abund,cieroot=cieroot,n_e=n_e,NH=NH,$
	hiceff=hiceff,lstr=lstr,cstr=cstr,toel=toel,toph=toph,toDN=toDN,$
	EM0=EM0,pixsize=pixsize,hicl=hicl,hicc=hicc,verbose=verbose,$
	_extra=e
;+
;function	hic_cal_emis
;	compute and return the combined emissivities of hic 
;	filter in an array of form [Ntemp,Nfilt] in units of
;	[1d-23 ergs cm^5/s] for each pixel
;
;syntax
;	hicemis=hic_cal_emis(filter,tgrid,ldbdir=ldbdir,$
;	cdbdir=cdbdir,ioneqf=ioneqf,abund=abund,cieroot=cieroot,$
;	n_e=n_e,NH=NH,hiceff=hiceff,lstr=lstr,cstr=cstr,$
;	/toel,/toph,/toDN,EM0=EM0,pixsize=pixsize,verbose=verbose,$
;	chidir=chidir,fH2=fH2,He1=He1,HeII=HeII)
;
;parameters
;	filter	[JUNK] scalar or array of filter names
;		* since Hi-C only used the 193 filter, this variable
;		  is completely ignored, and 193 is assumed regardless
;		  of input.  It is present here merely for compatibility
;		  with hinode_xrt_emis and sdo_aia_emis
;	tgrid	[I/O] the log(T[K]) grid over which the output
;		is to be defined`
;		* if not specified on input, will be set to be
;		  the same as that in the emissivity tables,
;		  the PoA default, findgen(81)*0.05+4
;
;keywords
;	ldbdir	[INPUT; '$CHIANTI'] line emissivity database directory
;	cdbdir	[INPUT; '$CONT'] continuum emissivity database directory
;	ioneqf	[INPUT; 'ioneq/mazzotta_etal.ioneq'] location of the
;		ion fraction tables
;	abund	[INPUT; getabund('Grevesse et al')] abundances
;	cieroot	[INPUT; 'cie'] root name for the files in $CONT
;	n_e	[INPUT; 1e9] electron number density in the plasma
;	NH	[INPUT; 0] H column density to apply
;	hiceff	[I/O] the HIC effective areas; input can be any of the
;		following:-
;		- a structure that contains the following fields --
;		  {NAME, CHANNELS, Awvl}
;		  where Awvl is a structure that contains the fields
;		  {NAME, WAVE, EA, PLATESCALE, UNITS}
;		  which can be generated using the SSW routine
;		  hiceff = HIC_GET_RESPONSE(/area)
;		- the name of a save file that contains the previously
;		  generated structure, named either HICEFF or EFF
;		- a flag that indicates that this must be calculated
;		  in situ by calling HIC_GET_RESPONSE() (requires SSW)
;		  - this can be accomplished by setting /HICEFF, but
;		    that will prevent it from being returned up.  instead,
;		    it is better to do something like
;		    hiceff=1 & hicemis=hic_emis(...,hiceff=hiceff,...)
;		* on output, will always contain the structure that
;		  was read in or calculated
;	lstr	[OUTPUT] line emissivities structure, as read in
;		from RD_LINE()
;	cstr	[OUTPUT] continuum emissivities structure, as read in
;		from RD_CONT()
;	toel	[INPUT; default=0] if set, returns the output in units
;		of [el cm^5/s]
;	toph	[INPUT; default=0] if set, returns the output in units
;		of [ph cm^5/s]
;	toDN	[INPUT; default=0] if set, returns the output in units
;		of [DN cm^5/pix]
;		* NOTE: TODN overrides TOPH overrides TOEL
;	EM0	[INPUT; default=1] if set, multiples the output by
;		this value to derive [ph/s] or [ph/cm^2/s] depending
;		on whether the input has units [cm^-3] or [cm^-5]
;	pixsize	[INPUT; default=0.6 arcsec] the size of a pixel on
;		the detector
;	hicl	[OUTPUT] the response due solely to lines
;	hicc	[OUTPUT] the response due solely to continuum
;		* both HICL and HICC are identical in size and units
;		  to the primary output, HICEMIS
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		FOLD_IONEQ : CHIDIR
;		ISMTAU : FH2, HE1, HEII
;		HIC_GET_RESPONSE : DN, PHOT
;
;history
;	Srividya Subramanian (Jul2013; based on sdo_aia_emis.pro)
;	made compatible with other mission_instrument_emis routines
;	  (VK; Jul2013)
;-

forward_function hic_get_response

;	usage
ok='ok' & np=n_params() ;& nf=n_elements(filter) & szf=size(filter,/type)
if np eq 0 then ok='Insufficient parameters' ;else $
 ;if nf eq 0 then ok='FILTER is undefined' else $
  ;if szf ne 7 then ok='FILTER must be a char string'
if ok ne 'ok' then begin
  print,'Usage: hicemis=hic_cal_emis(whatever,tgrid,ldbdir=ldbdir,$'
  print,'       cdbdir=cdbdir,ioneqf=ioneqf,abund=abund,cieroot=cieroot,$'
  print,'       n_e=n_e,NH=NH,hiceff=hiceff,lstr=lstr,cstr=cstr,$'
  print,'	/toel,/toph,/toDN,EM0=EM0,pixsize=pixsize,verbose=verbose,$'
  print,'       chidir=chidir,fH2=fH2,He1=He1,HeII=HeII)'
  print,'  compute and return emissivity functions for Hi-C'
  if np gt 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
if not keyword_set(ldbdir) then ldbdir='$CHIANTI'
if not keyword_set(ioneqf) then ioneqf='ioneq/chianti.ioneq'
if n_elements(abund) lt 30 then abund=getabund('grevesse et al.')
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(cieroot) then cieroot='cie'
if not keyword_set(n_e) then n_e=1e9
if not keyword_set(NH) then NH=0
if not keyword_set(EM0) then EM0=1.D
pixsizx=0.103 & pixsizy=pixsizx
if keyword_set(pixsize) then begin
  npix=n_elements(pixsize) & szp=size(pixsize,/type)
  if szp[0] le 5 then begin
    if pixsize[0] ne 0 then pixsizx=abs(pixsize[0])
    if npix ge 2 then $
     if pixsize[1] ne 0 then pixsizy=abs(pixsize[1]) else pixsizy=pixsizx
  endif else $
    message,'ignoring input PIXSIZE; using default',/informational
endif

if not keyword_set(hiceff) then begin	;(hiceff is not set
  ; look for save file in !ARDB
  defsysv,'!ARDB',ardbdir,exists=ivar	;!ARDB is defined if ivar=1
  if ivar eq 0 then begin	;(!ARDB undefined
    help,/source_files,output=cc
    ii=strpos(cc,'hic_cal_emis.pro') & oo=where(ii ge 0,moo)
    if moo eq 0 then begin	;(this file cannot be found?
      ardbdir='/data/fubar/SCAR/ardb'
    endif else begin		;file not found)(file found
      ccc=cc[oo[0]]
      iii=strpos(ccc,'pro') & i0=strpos(ccc,' ',/reverse_search)
      ardbdir=filepath('',root_dir=strmid(ccc,i0+1,iii-i0-1),subdir='ardb')
    endelse			;file found)
  endif else setsysval,'!ARDB',ardbdir,/getval	;!ARDB undefined)
  message,'restoring HICEFF from '+filepath('hic_eff.sav',root_dir=ardbdir),$
  	/informational
  restore,filepath('hic_eff.sav',root_dir=ardbdir),verbose=(vv<1)
  hiceff=hicea_dn
  if n_tags(hiceff) eq 0 then begin	;(hiceff doesn't exist
    if n_tags(eff) eq 0 then begin	;(eff also doesn't exist
      message,'Hi-C effective areas not found; quitting',/informational
      return,-1L
    endif else hiceff=eff		;eff)
  endif					;hiceff)
endif else begin			;hiceff is not set)(is set
  if n_tags(hiceff) eq 0 then begin	;(hiceff not a structure
    if size(hiceff,/type) eq 7 then begin	;(hiceff is a filename
      message,'restoring HICEFF from '+hiceff[0],/informational
      restore,hiceff[0],/verbose
      if n_tags(hicea_dn) ne 0 then hicaeff=hicea_dn
      if n_tags(hiceff) eq 0 then begin	;(hiceff doesn't exist
        if n_tags(eff) eq 0 then begin	;(eff also doesn't exist
          message,'HIC effective areas not found; quitting',/informational
          return,-1L
        endif else hiceff=eff		;eff)
      endif				;hiceff)
    endif else begin			;hiceff not filename)(calculate
      message,'calling hic_get_response() to compute HICEFF',$
      	/informational
      hiceff=hic_get_response(/area,_extra=e)	;REQUIRES SSW!
      message,'calling hic_get_response() to compute HICEFF',$
      	/informational
    endelse				;calculate hiceff in situ)
  endif 				;hiceff not a structure)
endelse					;hiceff)

intoel=0 & if keyword_set(toel) then intoel=1
intoph=0 & if keyword_set(toph) then intoph=1
if keyword_set(intoph) then intoel=0	;override TOEL
intoDN=0 & if keyword_set(toDN) then intoDN=1
if keyword_set(intoDN) then intoph=0	;override TOPH
if keyword_set(intoDN) then intoel=0	;override TOEL

;	extract some useful arrays from HICEFF
;
;<-- someday: here put in a check to make sure HICEFF is kosher -->
;

wave=hiceff.WAVE
eff_area=hiceff.EFFAREA

wmin=min(wave,max=wmax)

;	some other useful stuff
if vv ge 5 then peasecolr
inicon,fundae=fundae
fov=pixsizx*pixsizy
fov=fov*((!pi/180.)*(1./3600.))^2/(4.*!dpi)	;FOV of 1" pixels covers this many steradians/pix
ergperDN=fundae.ergev*3.65*15.	;[erg/eV]*[eV/e]*[e/DN]
if not keyword_set(toDN) then ergperDN=1.	;NO conversion to DN, even if
						;"ergperDN" is used downstream!
if keyword_set(intoel) then ergperDN=fundae.ergev*3.65	;[erg/eV]*[eV/e]
ergperDN=1d23*ergperDN		;to account for the 1e-23 in the emissivities
ergperDN=ergperDN/fov
xtitle='log T' & ytitle='[10!u-23!n ergs cm!u5!n s!u-1!n]'
if keyword_set(intoel) then ytitle='[electrons cm!u5!n s!u-1!n]'
if keyword_set(intoph) then ytitle='[ph cm!u5!n s!u-1!n]'
if keyword_set(intoDN) then ytitle='[DN cm!u5!n s!u-1!n]'

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
if keyword_set(toph) then nrgl=(fundae.h*fundae.c*1e8)/lwvl	;[erg/ph]

;	read in continuum emissivities
if n_tags(cstr) eq 0 then begin
  cconf=rd_cont(cieroot,n_e=n_e,wrange=[wmin,wmax],dbdir=cdbdir,abund=abund,fcstr=cstr,verbose=verbose,_extra=e)
endif
cww=mid2bound(cstr.midWVL) & cwvl=0.5*(cww[1:*]+cww) & ncw=n_elements(cwvl)
if keyword_set(toph) then nrgc=(fundae.h*fundae.c*1e8)/cwvl	;[erg/ph]
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
if keyword_set(toph) then lopt=lopt*(1d-23/nrgl)
if keyword_set(toph) then copt=copt*(1d-23/nrgc)

CHANNEL_NAME=['193'] & nchan=n_elements(channel_name)
ifilt=intarr(nchan)
ofilt=where(ifilt ge 0,mofilt)
; if mofilt eq 0 then begin
;   message,'No filters selected; quitting',/informational
;   return,-1L
; endif

;	output
if n_elements(tgrid) eq 0 then tgrid=lstr.LOGT
nT=n_elements(tgrid) & nmaxT=n_elements(lstr.LOGT) & jT=findex(lstr.LOGT,tgrid)
mf=1	;mofilt
hicemis=dblarr(nT,mf) & hicl=dblarr(nT,mf) & hicc=dblarr(nT,mf)
totemis=dblarr(nT)
totlemis=dblarr(nT)
totcemis=dblarr(nT)
warf=WAVE & larf=fltarr(nlw) & carf=fltarr(ncw)
larf=(interpol(EFF_AREA,warf,lwvl)>0)<(max(EFF_AREA))
carf=(interpol(EFF_AREA,warf,cwvl)>0)<(max(EFF_AREA))

;	compute the response
zab=abund[Z-1]
for it=0L,nT-1L do begin
  if vv gt 0 then kilroy
  jj=fix(jT[iT])
  jj0=(jj>0)<(nmaxT-1L) & jj1=(jj+1)<(nmaxT-1L)	;no extrapolations!
  djT=abs(jT[iT]-jj)
  ;
  ;	intensities
  linint=lopt*reform((lstr.LINE_INT)[jj0,*])*EM0*zab
  conint=copt*reform((cconf)[jj0,*])*EM0
  if djT gt 0.001 and $
     jj ge 0 and jj lt nmaxT-2 then begin totlemis	;(interpolate between slices
    linint=djT*linint + (1.-djT)*lopt*reform((lstr.LINE_INT)[jj1,*])*EM0*zab
    conint=djT*conint + (1.-djT)*copt*reform((cconf)[jj1,*])*EM0
  endif						;JT>0.001)
  ;	fluxes
  totemis[it]=(total(linint)+total(conint))/ergperDN
  totlemis[it]=(total(linint)+0*total(conint))/ergperDN
  totcemis[it]=(0*total(linint)+total(conint))/ergperDN
  for jf=0,mf-1 do hicemis[it,jf]=$
	(total(linint*larf[*,jf])+total(conint*carf[*,jf]))/ergperDN
  for jf=0,mf-1 do hicl[it,jf]=(total(linint*larf[*,jf]))/ergperDN
  for jf=0,mf-1 do hicc[it,jf]=(total(conint*carf[*,jf]))/ergperDN
  if it gt 0 and vv ge 5 then begin
    plot,tgrid[0:it],totemis[0:it,0],/yl,yr=max(hicemis[0:it,*])*[1e-5,2],$
	xtitle=xtitle,ytitle=ytitle,_extra=e
    for i=0,mf-1 do oplot,tgrid,hicemis[*,i],col=1+(i mod 8)
    if vv ge 10 then begin
      dy=(!y.crange[1]-!y.crange[0])/(n_elements(CHANNEL_NAME)+5)
      for i=0,mf-1 do xyouts,!x.crange[0]+0.1,10.^(!y.crange[1]-(i+1)*dy),CHANNEL_NAME[ofilt[i]],color=1+(i mod 8)
    endif
  endif
endfor

if vv ge 5 then begin
  plot,tgrid,totemis,/yl,yr=max(hicemis)*[1e-5,2],$
    xtit=xtitle,ytit=ytitle,_extra=e
  for i=0,mf-1 do oplot,tgrid,hicemis[*,i],col=1+(i mod 8)
  dy=(!y.crange[1]-!y.crange[0])/(n_elements(CHANNEL_NAME)+5)
  for i=0,mf-1 do xyouts,!x.crange[0]+0.1,10.^(!y.crange[1]-(i+1)*dy),$
	CHANNEL_NAME[ofilt[i]],color=1+(i mod 8)
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,hicemis
end
