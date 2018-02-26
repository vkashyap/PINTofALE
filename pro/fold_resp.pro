function fold_resp,resp,iE,flx,binE=binE, _extra=e
;+
;function	fold_resp
;	returns the PH spectrum as observed with an instrument defined
;	by the specified response
;
;syntax
;	ph=fold_resp(resp,iE,flx,/binE,fchcol=fchcol,/shift1)
;
;parameters
;	resp	[INPUT; required] either of:
;		A: 2D array (N_NRG,N_PH) describing the response to a photon
;		B: name of OGIP-compliant FITS file containing the response
;	iE	[INPUT; required] either of:
;		A': bin indices of energy array corresponding to rows of the
;		    response matrix
;		B': energy of photons which must be binned according to
;		    the response matrix resolution
;	flx	[INPUT] flux at given iE -- a weighting function
;		* if scalar or 1-element vector, essentially acts as
;		  "normalization"
;
;keywords
;	binE	[INPUT] if set, *and* RSP is a filename, IE is assumed to
;		be actual energies and is binned appropriately.
;	_extra	[INPUT ONLY] pass defined keywords to subroutines:
;		RD_OGIP_RMF: FCHCOL, SHIFT1
;
;restrictions
;	* requires subroutines RD_OGIP_RMF, RDRESP, BINERSP,
;	  and the IDLASTRO library
;
;history
;	vinay kashyap (Jul97)
;	corrected bug with bad binning in NRG, now handles large RESP
;	  (VK; Mar99)
;-

;	usage
szr=size(resp) & szE=size(iE) & n_r=n_elements(szr) & n_E=n_elements(szE)
nfx=n_elements(flx)
ok=''		;assume all OK
if szr(1) eq 0 then ok='Response not defined' else $
 if szE(1) eq 0 then ok='Photon energies not defined' else $
  if szr(0) ne 2 then begin
    if szr(n_r-1) ne 1 then ok='cannot understand Response' else $
     if szr(n_r-2) ne 7 then ok='Need name of response matrix file'
  endif ;else if szr(1) ne szE(1) then ok='RESP and E not compatible'
if ok ne '' then begin
  message,ok,/info
  print,'Usage: ph=fold_resp(resp,iE,flx,/binE)'
  print,'  returns PH spectrum'
  return,-1L
endif

;	if RESP=filename, call RDRESP
jE=[iE(*)]
resp_was_file=0
if szr(n_r-1) eq 1 and szr(n_r-2) eq 7 then begin
  rsp=rdresp(resp,bstr,_extra=e) & szr=size(rsp)
  ;	if IE are actual energies, bin them to position indices
  if keyword_set(binE) then jE=binersp(iE,bstr=bstr)
  resp_was_file=1
endif else rsp=temporary(resp)
szr=size(rsp)

;	get weighting function
;nnrg=szE(n_E-1) & fx=fltarr(nchan)+1.
nnrg=n_elements(jE) & nchan=szr(2)
fx=fltarr(nnrg)+1.
if nfx eq 1 then fx(*)=flx(0)
if nfx gt 1 and nfx le nnrg then fx(0)=flx
if nfx gt nnrg then fx(*)=flx(0:nfx-1)

;	do the additions
ph=fltarr(nchan)
;for i=0,nchan-1 do if jE(i) ge 0 then ph=ph+rsp(jE(i),*)*fx(i)
if nchan le nnrg then begin
  for i=0,nchan-1 do begin
    ;oo=where(jE eq i,moo)
    ;if moo gt 0 then ph=ph+rsp(jE(i),*)*total(fx(oo))
    ph(i)=total(reform(rsp(jE,i))*fx)
  endfor
endif else begin
  for i=0,nnrg-1 do begin
    if jE(i) ge 0 then ph=ph+rsp(jE(i),*)*fx(i)
  endfor
endelse

if resp_was_file eq 0 then resp=temporary(rsp)

return,ph
end
