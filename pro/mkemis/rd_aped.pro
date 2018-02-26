function rd_aped,apedfil,wvls,temps,edens,Zs,ions,levs,$
	apeddir=apeddir,xtnames=xtnames, _extra=e
;+
;function	rd_aped
;	read in a line list from APED and return the 3D
;	emissivities [WVLS,TEMPS,EDENS]
;
;syntax
;	emiss=rd_aped(apedfil,wvls,temps,edens,Zs,ions,levs,$
;	apeddir=apeddir,xtnames=xtnames)
;
;parameters
;	apedfil	[INPUT; required] FITS file, output of APEC
;		* assumed to be in following format --
;		  [PARAMETERS][cols kT,EDensity,Nelement,Nline]
;		  	where kT is in [keV], EDensity is in [cm^-3]
;		  [L{logT}_{logNe}][cols Lambda,Lambda_Err,Epsilon,Epsilon_Err,Element,Ion,UpperLev,LowerLev]
;		  	where {logT} and {logNe} are in '(f3.1)' or '(f4.1)',
;		  	Lambda* are in [\AA], and Epsilon* are in [ph cm^3/s]
;	wvls	[OUTPUT; required] all the wavelengths listed in APEDFIL [\AA]
;	temps	[OUTPUT; required] all the temperatures in APEDFIL [log[K]]
;	edens	[OUTPUT; required] all the e-densities in APEDFIL [cm^-3]
;	Zs	[OUTPUT; required] atomic numbers corresponding to each WVLS
;	ions	[OUTPUT; required] ionic state of each of ZS
;	levs	[OUTPUT; required] level designations for each line [2,N(WVLS)]
;
;keywords
;	apeddir	[INPUT] directory containing APEDFIL
;		* if not set, looks for APEDFIL in $cwd
;	xtnames	[OUTPUT] all the extension names
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	for APED format as of Apr 2000 only
;	subroutines: KILROY
;
;history
;	vinay kashyap (AprMM)
;	bug correction (thanks to Antonio Maggio) -- fundae.e and
;	  fundae.kB were not defined in lines 90 and 92 (VK; Aug02)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(apedfil)
szf=size(apedfil) & nszf=n_elements(szf)
if np lt 7 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='APED file undefined' else $
  if szf[nszf-2] ne 7 then ok='APED FITS file needed'
if ok ne 'ok' then begin
  print,'Usage: emiss=rd_aped(apedfil,wvls,temps,edens,Zs,ions,levs,$'
  print,'             apeddir=apeddir,xtnames=xtnames)'
  print,'  read in line emissivities from APED line list'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	keywords
daped='' & nd=n_elements(apeddir) & szd=size(apeddir) & nszd=n_elements(szd)
if nd ne 0 then begin
  if szd[nszd-2] eq 7 then begin
    daped=strtrim(apeddir[0],2) & dalen=strlen(daped)
    if strmid(daped,dalen-1,1) ne '/' then daped=daped+'/'
  endif
endif

;	initialize
faped=daped+apedfil[0]
xtparam=1

;	read the [PARAMETERS] extension
tmp=readfits(faped,hpar,exten=xtparam)
xtnam=strupcase(strtrim(sxpar(hpar,'EXTNAME'),2))
if xtnam ne 'PARAMETERS' then begin
  message,faped+': in unknown format',/info & return,-1L
endif
;
nflds=sxpar(hpar,'TFIELDS')
for i=0,nflds-1 do begin
  colnam=strupcase(strtrim(sxpar(hpar,'TTYPE'+strtrim(i+1,2)),2))
  if colnam eq 'KT' then allkT=fits_get(hpar,tmp,i+1)
  if colnam eq 'EDENSITY' then allN_e=fits_get(hpar,tmp,i+1)
  if colnam eq 'NLINE' then allNline=fits_get(hpar,tmp,i+1)
endfor

;	how many points in the grid?
kTs=allkT(uniq(allkT,sort(allkT))) & nkTs=n_elements(kTs)
edens=allN_e(uniq(allN_e,sort(allN_e))) & neds=n_elements(edens)
;	convert the kTs to log(Temperature)s
;temps=alog10(erg(erg(kTs*1e3,/eV),/degK,/inv))
;	convert [keV] -> [K]
t=kTs*1e3	;[keV] -> [eV]
t=t*1.6021892e-19	;[eV] -> [J]
t=t*1e7		;[J] -> [erg]
t=t/1.3806620e-16	;[erg] -> [K]
temps=alog10(t)

;	initialize outputs
numext=nkTs*neds
xtnames=strarr(nkTs,neds)

;	step high and wide over APEDFIL
for iext=2L,numext+1L do begin		;{read in from each extension
  tmp=readfits(faped,htmp,exten=iext,pointlun=pointlun,/silent)
  ntmp=n_elements(tmp) & szt=size(tmp) & nszt=n_elements(szt)
  if abs(iext) gt 0 and ntmp eq 1 and szt(nszt-2) eq 2 and $
	tmp(0) eq -1 then begin
    message,faped+': bad file?',/info & return,-1L
  endif

  ;	get the extension name and figure out which (T,N_e) point it refers to
  xtnam=strupcase(strtrim(sxpar(htmp,'HDUNAME'),2))
  tt=float(strmid(xtnam,1,3)) & dd=float(strmid(xtnam,5))
  mint=min(abs(temps-tt),itmin) & mind=min(abs(alog10(edens)-dd),idmin)
  xtnames[itmin,idmin]=xtnam
  print,xtnam,itmin,idmin

  ;	extract the wavelengths, emissivities, elements, ions, and levels
  icol_w=0 & icol_em=2 & icol_Z=4 & icol_i=5 & icol_ul=6 & icol_ll=7
  ww=fits_get(htmp,tmp,icol_w+1)
  ee=fits_get(htmp,tmp,icol_em+1)
  zz=fits_get(htmp,tmp,icol_z+1)
  ii=fits_get(htmp,tmp,icol_i+1)
  ul=fits_get(htmp,tmp,icol_ul+1)
  ll=fits_get(htmp,tmp,icol_ul+1)

  ;	store as necessary
  nww=n_elements(ww)
  if keyword_set(oink) then begin	;(all variables have been defined
    for j=0L,nww-1L do begin	;{check each wavelength for duplication
      dw=abs(wvls-ww[j])
      ow=where(dw lt 1e-3 and Zs eq zz[j] and ions eq ii[j] and $
	reform(levs[0,*]) eq ll[j] and reform(levs[1,*]) eq ul[j],mow)
      if mow eq 0 then begin	;(new entry
	wvls=[wvls,ww[j]] & Zs=[Zs,zz[j]] & ions=[ions,ii[j]]
	mw=n_elements(wvls)
	if mw eq 100L*long(mw/100.) then kilroy,dot=' '+strtrim(mw,2)
	;emiss=reform([emiss[*],dblarr(nkTs*neds)],mw,nkTs,neds)
	;levs=reform([levs[*],ll[j],ul[j]],2,mw)
	tmp=dblarr(mw,nkTs,neds) & for k=0L,mw-2L do tmp[k,*,*]=emiss[k,*,*]
	tmp[mw-1L,*,*]=dblarr(nkTs,neds) & emiss=tmp
	emiss[mw-1L,itmin,idmin]=ee[j]
	tmp=intarr(2,mw) & for k=0,1 do tmp[k,0:mw-2L]=levs[k,*]
	tmp[0,mw-1L]=ll[j] & tmp[1,mw-1L]=ul[j] & levs=tmp
      endif else begin		;)(old entry
	tmp=min(dw,iwmn)
	emiss[iwmn,itmin,idmin]=emiss[iwmn,itmin,idmin]+ee[j]
	;if mow gt 1 then message,'DR line being subsumed',/info
	if mow gt 1 then kilroy,dot='*'
      endelse			;MOW)
    endfor			;J=0,NWW-1}
  endif else begin			;)(not
    wvls=ww & Zs=zz & ions=ii
    emiss=dblarr(nww,nkTs,neds) & emiss[*,itmin,idmin]=ee[*]
    levs=intarr(2,nww) & levs[0,*]=ll[*] & levs[1,*]=ul[*]
    oink=1
  endelse				;first time initialize)
endfor					;IEXT=2,NUMEXT+1}

return,emiss
end
