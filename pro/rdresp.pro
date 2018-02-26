function rdresp,rspfil,ostr,effar=effar, _extra=e
;+
;function	rdresp
;	reads in an OGIP-compliant response matrix file and returns the
;	full response matrix
;
;syntax
;	rsp=rdresp(rspfil,ostr,effar=effar,fchcol=fchcol,/shift1)
;
;parameters
;	rspfil	[INPUT; required] name of response file
;	ostr	[OUTPUT] structure containing all the information for the
;		response matrix: {NNRG, ELO, EHI, NCHAN, EMN, EMX, N_GRP,
;		F_CHAN, N_CHAN, MATRIX, FIRSTCHAN}
;		(see RD_OGIP_RMF() for description)
;		* output response matrix is in the form RSP(NNRG,NCHAN)
;
;keywords
;	effar	[OUTPUT] returns the effective area as a function of ELO
;		NOTE: if an RMF file is read instead of an RSP file,
;		EFFAR are actually the efficiencies
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		RD_OGIP_RMF: FCHCOL,SHIFT1
;
;restrictions
;	requires IDLASTRO library routine MRDFITS
;	requires RD_OGIP_RMF
;
;history
;	vinay kashyap (Jul97)
;	added keyword SHIFT1 (VK; JanMMI)
;	removed SHIFT1, added call to RD_OGIP_RMF, cleaned up (VK; Nov2001)
;-

;	usage
ok='ok' & nrsp=n_elements(rspfil) & szr=size(rspfil) & nszr=n_elements(szr)
if nrsp eq 0 then ok='Name of Response Matrix File required' else $
 if nrsp gt 1 then ok='Cannot handle multiple RSPs!' else $
  if szr[nszr-2] ne 7 then ok='input must be filename'
if ok ne 'ok' then begin
  print,'Usage: rsp=rdresp(rspfil,ostr,effar=effar,fchcol=fchcol,/shift1)'
  print,'  read OGIP-compliant RMF and return response matrix'
  if nrsp ne 0 then message,'	'+ok,/info
  return,0.
endif

ostr=rd_ogip_rmf(rspfil[0],effar=effar, _extra=e)
if n_tags(ostr) eq 0 then return,0.

nnrg=ostr.NNRG & elo=ostr.ELO & ehi=ostr.EHI
n_grp=ostr.N_GRP & f_chan=ostr.F_CHAN & n_chan=ostr.N_CHAN
matrix=ostr.MATRIX & firstchan=ostr.FIRSTCHAN 
nchan=ostr.NCHAN & emn=ostr.EMN & emx=ostr.EMX

if keyword_set(FIRSTCHAN) then message,$
	'assuming that input RMF is 1-based, not 0-based',/info

;	read into full response matrix and figure out effective area
rsp=fltarr(nnrg,nchan) & effar=fltarr(nnrg) & szn=size(n_chan)
for inrg=0,nnrg-1 do begin
  ngrp=n_grp[inrg] & jbeg=0
  for ig=0,ngrp-1 do begin
    if szn[0] gt 1 then begin
      ibeg=f_chan[ig,inrg] & iw=n_chan[ig,inrg]
    endif else begin
      ibeg=f_chan[inrg] & iw=n_chan[inrg]
    endelse
    if keyword_set(firstchan) then ibeg=ibeg-1	;IDL index correction
    if iw gt 0 then rsp[inrg,ibeg:ibeg+iw-1]=matrix[jbeg:jbeg+iw-1,inrg]
    if iw gt 0 then effar[inrg]=effar[inrg]+total(matrix[jbeg:jbeg+iw-1,inrg])
    jbeg=jbeg+iw
  endfor
endfor

return,rsp
end
