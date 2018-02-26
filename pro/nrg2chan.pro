pro nrg2chan,nrg,flx,chan,spec,rmstr,effar=effar,nrgar=nrgar, _extra=e
;+
;procedure	nrg2chan
;	convolves input energy spectrum with an OGIP-style response matrix
;
;syntax
;	nrg2chan,nrg,flx,chan,spec,rmstr,effar=effar,nrgar=nrgar
;
;parameters
;	nrg	[INPUT; required] bin boundary values at which input FLX is
;		defined.
;		* usually [keV], depends on RMF file
;	flx	[INPUT; required] fluxes in [ph/bin/(s/cm^2/sr/...)]
;	chan	[OUTPUT; required] bin boundary values at which output SPEC
;		is defined
;		* same units as NRG
;	spec	[OUTPUT; required] output, convolved spectrum
;	rmstr	[INPUT; required] structure containing all the relevant
;		data (e.g., see output of RD_OGIP_RMF):
;		{NNRG,ELO,EHI,NCHAN,EMN,EMX,N_GRP,F_CHAN,N_CHAN,MATRIX}
;		* ELO and EHI refer to range of photon energies at which
;		  response is valid
;		* EMN and EMX refer to bin boundaries for each channel
;		* N_GRP, F_CHAN, and N_CHAN refers to number of groups
;		  of non-zero data, beginning channel, and number of
;		  channels in each group
;
;keywords
;	effar	[INPUT] effective areas [cm^2]
;	nrgar	[INPUT] energies at which EFFAR are defined
;		* same units as NRG
;		* sizes of EFFAR and NRGAR must match
;		* range of NRGAR must overlap NRG
;		* if EFFAR and NRGAR are legal, FLX is multiplied by
;		  EFFAR (over intersection of NRGAR and NRG; zeroed
;		  elsewhere) before convolving with response matrix
;		* of course, this multiplication makes sense only if
;		  units on FLX are [ph/bin/cm^2/...]
;	_extra	[JUNK] ignore - here only to prevent crashing
;
;restrictions
;	requires subroutine BINERSP
;	requires output of RD_OGIP_RMF
;
;history
;	vinay kashyap (Apr99; based on CONV_RMF)
;-

;	usage
ok='ok'
np=n_params() & mnrg=n_elements(nrg) & mflx=n_elements(flx)
mrm=n_tags(rmstr)
if np lt 5 then ok='missing parameters' else $
 if mnrg eq 0 then ok='requires input spectrum' else $
  if mflx eq 0 then ok='requires spectrum' else $
   if mnrg ne mflx and mnrg ne (mflx+1L) then ok='input spectrum garbled' else $
    if mrm lt 7 then ok='requires output of RD_OGIP_RMF'
if ok ne 'ok' then begin
  print,'Usage: nrg2chan,nrg,flx,chan,spec,rmstr,effar=effar,nrgar=nrgar'
  print,'  convolve spectrum with response'
  if np gt 0 then message,ok,/info
  return
endif

;	check inputs
ee=[nrg(*)] & fe=[flx(*)]
if mnrg eq mflx then begin
  ;	assume input NRG are mid-bin values and decipher the bin boundaries
  de=nrg(1:*)-nrg
  ee=[nrg-0.5*de,nrg(mnrg-1)+0.5*de(mnrg-2)*[-1,1]]
  mnrg=mnrg+1L
  if mnrg ne mflx+1L then message,'BUG!'
endif
;
trm=tag_names(rmstr) & i_MATRIX=-1L
i_NNRG=-1L & i_ELO=-1L & i_EHI=-1L
i_NCHAN=-1L & i_EMN=-1L & i_EMX=-1L
for i=0L,n_elements(trm)-1L do begin
  if trm(i) eq 'NNRG' then i_NNRG=i
  if trm(i) eq 'ELO' then i_ELO=i
  if trm(i) eq 'EHI' then i_EHI=i
  if trm(i) eq 'NCHAN' then i_NCHAN=i
  if trm(i) eq 'EMN' then i_EMN=i
  if trm(i) eq 'EMX' then i_EMX=i
  if trm(i) eq 'MATRIX' then i_MATRIX=i
endfor
ok='ok'
if i_NNRG lt 0 then ok='Input structure missing NNRG; returning'
if i_ELO lt 0 then ok='Input structure missing ELO; returning'
if i_EHI lt 0 then ok='Input structure missing EHI; returning'
if i_NCHAN lt 0 then ok='Input structure missing NCHAN; returning'
if i_EMN lt 0 then ok='Input structure missing EMN; returning'
if i_EMX lt 0 then ok='Input structure missing EMX; returning'
if i_MATRIX lt 0 then ok='Input structure missing MATRIX; returning'
if ok ne 'ok' then begin
  message,ok,/info & return
endif

;	include effective areas?
mnar=n_elements(nrgar) & mear=n_elements(effar)
if mnar eq mear and mear gt 0 then begin	;(effective area defined
  arnrg=ee & areff=0.*arnrg+1.
  if mnar eq 1 then areff(*)=effar(0) else begin
    minea=min(effar,max=maxea) & minna=min(nrgar,max=maxna)
    areff=(interpol(effar,wvlar,ee) > (minea)) < (maxea)
    oo=where(arnrg lt minna or arnrg gt maxna,moo)
    if moo gt 0 then areff(oo)=0.
  endelse
  fe=fe*areff	;assumes input FLX is in units of [ph/bin/cm^2/...]
endif						;EFFAR(NRGAR))

;	define outputs
chan=[rmstr.EMN,rmstr.EMX(rmstr.NCHAN-1L)]
spec=fltarr(rmstr.NCHAN)

;	find binning for input NRG
ie=binersp(ee,bstr=rmstr)

;	distribute each input FLX into PHAs
for inrg=0L,rmstr.NNRG-1L do begin
  rsp=fltarr(rmstr.NCHAN) ;& whee,spoke=spoke,/moveit
  ;
  ngrp=rmstr.N_GRP(inrg) & jbeg=0
  szn=size(rmstr.N_CHAN)
  for ig=0,ngrp-1 do begin
    if szn(0) gt 1 then begin
      ibeg=rmstr.F_CHAN(ig,inrg) & iw=rmstr.N_CHAN(ig,inrg)
    endif else begin
      ibeg=rmstr.F_CHAN(inrg) & iw=rmstr.N_CHAN(inrg)
    endelse
    ;someone please tell me if this is needed or not!!
    ;needed for ASCA, not for XTE??
    ibeg=ibeg-1			;IDL index correction
    if iw gt 0 then rsp(ibeg:ibeg+iw-1)=rmstr.MATRIX(jbeg:jbeg+iw-1,inrg)
    jbeg=jbeg+iw
  endfor
  ;
  oi=where(ie eq inrg,moi)
  if moi gt 0 then spec=spec+rsp*total(fe(oi))
endfor

return
end
