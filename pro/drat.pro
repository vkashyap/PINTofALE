function drat,line1,line2,edens,DEM,logT,fx1=fx1,fx2=fx2, _extra=e
;+
;function	drat
;	return the ratio of fluxes in one line relative to another line
;	for a variety of densities
;
;syntax
;	frat=drat(line1,line2,edens,DEM,logT,fx1=fx1,fx2=fx2,$
;	mapping=mapping,chifil=chifil,chidir=chidir,eqfile=eqfile,$
;	abund=abund,/noph,effar=effar,wvlar=wvlar,NH=NH,fH2=fH2,$
;	He1=He1,HeII=HeII,/Fano,dbdir=dbdir,sep=sep,prefix=prefix)
;
;parameters
;	line1	[INPUT; required] string describing the line/wavelength-range
;		that defines the numerator
;	line2	[INPUT; required] string describing the line/wavelength-range
;		that defines the denominator
;		* LINE1 and LINE2 must be in the format used by RD_LIST, i.e.,
;		  	Z ION <sep> WAVE <sep> DBDIR <sep> DESCRIPTION
;		  and WAVE = "WVL"/"WVL +- dW"/"WMIN,WMAX"/"WMIN-WMAX"
;	edens	[INPUT; required] electron densities at which to compute
;		the flux ratios
;	DEM	[INPUT; required] differential emission measure to use in
;		computing the fluxes
;	logT	[INPUT; required] temperatures at which DEM is defined
;		* size must match that of DEM
;
;keywords
;	fx1	[OUTPUT] the fluxes due to line 1 at each EDENS
;	fx2	[OUTPUT] the fluxes due to line 2 at each EDENS
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		RD_LIST: MAPPING,DBDIR,SEP,PREFIX
;		FOLD_IONEQ: CHIFIL,EQFILE
;		RD_IONEQ: CHIDIR
;		LINEFLX: ABUND,NOPH,EFFAR,WVLAR
;		ISMTAU: NH,fH2,He1,HeII,FANO,BAM
;		BAMABS: ABUND
;
;subroutines
;	RD_LIST
;	RD_LINE
;	FOLD_IONEQ
;	LINEFLX
;	ISMTAU
;	BAMABS
;	CAT_LN
;
;history
;	vinay kashyap (MarMM)
;-

;	usage
ok='ok' & np=n_params() & nden=n_elements(edens)
n1=n_elements(line1) & szl1=size(line1) & nszl1=n_elements(szl1)
n2=n_elements(line2) & szl2=size(line2) & nszl2=n_elements(szl2)
ndem=n_elements(DEM) & nT=n_elements(logT)
if np lt 5 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='Line 1 undefined' else $
  if n2 eq 0 then ok='Line 2 undefined' else $
   if szl1(nszl1-2) ne 7 then ok='Line 1 in unknown format' else $
    if szl2(nszl2-2) ne 7 then ok='Line 2 in unknown format' else $
     if nden eq 0 then ok='Densities undefined' else $
      if nden lt 2 then ok='Insufficient densities specified' else $
       if ndem eq 0 then ok='DEM undefined' else $
        if nT eq 0 then ok='logT in DEM(logT) undefined' else $
         if ndem ne nT then ok='DEM(logT) and logT incompatible'
if ok ne 'ok' then begin
  print,'Usage: frat=drat(line1,line2,edens,DEM,logT,fx1=fx1,fx2=fx2,$'
  print,'       mapping=mapping,chifil=chifil,chidir=chidir,eqfile=eqfile,$'
  print,'       abund=abund,/noph,effar=effar,wvlar=wvlar,NH=NH,fH2=fH2,$'
  print,'       He1=He1,HeII=HeII,/Fano,dbdir=dbdir,sep=sep,prefix=prefix)'
  print,'  return the ratio of fluxes in two different lines for various'
  print,'  electron densities'
  if np ne 0 then message,ok,/info
  return,0.
endif

;	define outputs
frat=fltarr(nden)-1. & fx1=fltarr(nden) & fx2=fltarr(nden)

;	compute fluxes
for i=0L,nden-1L do begin			;{for each density
  n_e=edens(i)
  lstr1=rd_list(line1,/incieq,n_e=n_e, _extra=e)
  lstr2=rd_list(line2,/incieq,n_e=n_e, _extra=e)
  if n_elements(lstr1.LOGT) ne nT then begin
    D_E_M=mk_dem('interpolate',logT=lstr1.LOGT,indem=DEM,pardem=logT)
    log_T=lstr1.LOGT
  endif else begin
    D_E_M=DEM & log_T=logT
  endelse
  flux1=lineflx(lstr1.LINE_INT,log_T,lstr1.WVL,lstr1.Z,DEM=D_E_M, _extra=e)
  flux2=lineflx(lstr2.LINE_INT,log_T,lstr2.WVL,lstr2.Z,DEM=D_E_M, _extra=e)
  rtau1=exp(-ismtau(lstr1.WVL, _extra=e))
  rtau2=exp(-ismtau(lstr2.WVL, _extra=e))
  fx1(i)=total(flux1*rtau1) & fx2(i)=total(flux2*rtau2)
endfor						;I=0,NDEN-1}

;	make the flux ratios
oo=where(fx2 gt 0,moo) & if moo gt 0 then frat(oo)=fx1(oo)/fx2(oo)

return,frat
end
