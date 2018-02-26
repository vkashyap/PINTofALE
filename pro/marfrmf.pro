function marfrmf,arffil,rmffil,ostr,effar=effar,noext=noext,nofull=nofull,$
	_extra=e
;+
;function	marfrmf
;	returns the instrument response matrix computed as the product of
;	the redistribution matrix with the ancillary response.
;
;syntax
;	rsp=marfrmf(arffil,rmffil,ostr,effar=effar,/noext,/nofull,$
;	fchcol=fchcol,/shift1)
;
;parameters
;	arffil	[INPUT; required] OGIP-compliant Ancillary Response File
;		* appends a ".arf" if missing (overridden by NOEXT)
;	rmffil	[INPUT; required] OGIP-compliant Redistribution Matrix File
;		* appends a ".rmf" if missing (overridden by NOEXT)
;	ostr	[OUTPUT] structure containing the axes information for the
;		response matrix -- see RD_OGIP_RMF()
;
;keywords
;	effar	[OUTPUT] returns the effective area as a function of ELO
;	noext	[INPUT] normally, if an extension is missing from ARFFIL
;		or RMFFIL, the program will append the appropriate extension.
;		if this keyword is set, no such action will occur.
;	nofull	[INPUT] if set, bypasses RDRESP() and the accompanying
;		computation of the full response matrix (i.e., with the OGIP
;		compression removed).  this is useful for those cases where
;		the uncompressed matrix is too large for memory.
;		* if this is set, the output will be simply the same as
;		  OSTR.MATRIX, which is the compressed version of the
;		  effective area weighted redistribution matrix
;	_extra	[INPUT] pass defined keywords to subroutines
;		RD_OGIP_RMF: FCHCOL, SHIFT1
;
;restrictions
;	* requires subroutines RD_OGIP_RMF, RDRESP and RDARF
;	* requires IDLASTRO library of routines
;
;history
;	vinay kashyap (Jul97)
;	allowed for keyword SHIFT1 (VK; JanMMI)
;	added keyword NOFULL; now multiplication also works on OSTR.MATRIX
;	  (VK; Nov2001)
;	changed output in case of /NOFULL from 1.0 to OSTR.MATRIX (VK; Oct2008)
;-

;	usage
narf=n_elements(arffil) & nrmf=n_elements(rmffil)
if narf ne 1 or nrmf ne 1 then begin
  print,'Usage: rsp=marfrmf(arffil,rmffil,ostr,effar=effar,/noext,/shift1)'
  print,'  returns RSP = ARF*RMF'
  return,-1L
endif

;	need to add extensions?
farf=arffil[0] & frmf=rmffil[0]
if not keyword_set(noext) then begin
  iextarf=strpos(farf,'.',0) & iextrmf=strpos(frmf,'.',0)
  if iextarf lt 0 then farf=strtrim(farf,2)+'.arf'
  if iextrmf lt 0 then frmf=strtrim(frmf,2)+'.rmf'
endif

;	read ARF
arf=rdarf(farf,astr) & a_elo=astr.elo & naelo=n_elements(a_elo)

;	read RMF
if keyword_set(nofull) then begin
  ostr=rd_ogip_rmf(frmf,effar=effar,_extra=e)
  rmf=1.0 & rsp=rmf
endif else begin
  rmf=rdresp(frmf,ostr,effar=effar,_extra=e)
endelse
r_elo=ostr.elo & nrelo=n_elements(r_elo)

;	sanity check
if naelo ne nrelo then begin
  message,farf+' and '+frmf+' are not compatible.  returning RMF',/info
  return,rmf
endif

;	do the multiplication
matrix=ostr.MATRIX
for i=0L,nrelo-1L do begin
  tmp=(ostr.MATRIX)[*,i] * arf[i] & matrix[*,i]=tmp[*]
endfor
ostr.MATRIX=matrix
if not keyword_set(nofull) then begin
  rsp=rmf & for i=0L,nrelo-1 do rsp[i,*]=rsp[i,*]*arf[i]
endif else rsp=matrix

;	get the effective areas
if not keyword_set(nofull) then begin
  effar=fltarr(nrelo)+1. & effar=reform(effar#rsp)
endif else begin
  effar=fltarr(nrelo)+1.
  for i=0L,nrelo-1L do effar[i]=total(matrix[*,i])
endelse

return,rsp
end
