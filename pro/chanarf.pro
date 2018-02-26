function chanarf,rmfstr,effar,spec=spec,chnrg=chnrg,verbose=verbose, _extra=e
;+
;function	chanarf
;	compute the average effective area over each channel, as a
;	weighted average of the effective area from all the energies
;	that contribute to that channel.
;
;	makes a difference mainly for low and medium spectral
;	resolution data.  for high-res grating data, not so much.
;
;syntax
;	carf=chanarf(rmfstr,effar,spec=spec)
;
;parameters
;	rmfstr	[INPUT; required] response matrix structure, in the
;		same format as returned from RD_OGIP_RMF()
;	effar	[INPUT] effective area as a function of energy [cm^2]
;		* assumed to be on the same grid as the input energies
;		  of RMFSTR, i.e., RMFSTR.ELO and RMFSTR.EHI
;		* ignored if size doesn't match RMFSTR
;
;keywords
;	spec	[INPUT] a spectrum to further weight the contribution
;		of EFFAR
;		* assumed to be on the same grid as RMFSTR.ELO and RMFSTR.EHI
;		* ignored if size doesn't match EFFAR
;		* ideally should be in units of [ph/keV/...]
;	chnrg	[OUTPUT] the average energy of a photon in this channel
;		* it should be similar to 0.5*(RMFSTR.EMN+RMFSTR.EMX), and
;		  serves as a check on the gain calibration
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing
;
;history
;	vinay kashyap (MMVI.IX)
;-

;	usage
ok='ok' & np=n_params() & nr=n_elements(rmfstr) & nrt=n_tags(rmfstr)
if np eq 0 then ok='Insufficient parameters' else $
 if nr eq 0 then ok='RMFSTR is undefined' else $
  if nrt eq 0 then ok='RMFSTR must be a structure; see RD_OGIP_RMF()'
if ok ne 'ok' then begin
  print,'Usage: carf=chanarf(rmfstr,effar,spec=spec,chnrg=chnrg,$'
  print,'       verbose=verbose)'
  print,'  compute the average effective area over each channel'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
ok='ok' & rtags=tag_names(rmfstr)
iNNRG=where(strpos(rtags,'NNRG') ge 0)
if iNNRG[0] ne -1 then NNRG=rmfstr.NNRG else ok='RMFSTR.NNRG not found'
iELO=where(strpos(rtags,'ELO') ge 0)
if iELO[0] ne -1 then ELO=rmfstr.ELO else ok='RMFSTR.ELO not found'
iEHI=where(strpos(rtags,'EHI') ge 0)
if iEHI[0] ne -1 then EHI=rmfstr.EHI else ok='RMFSTR.EHI not found'
iNCHAN=where(strpos(rtags,'NCHAN') ge 0)
if iNCHAN[0] ne -1 then NCHAN=rmfstr.NCHAN else ok='RMFSTR.NCHAN not found'
iEMN=where(strpos(rtags,'EMN') ge 0)
if iEMN[0] ne -1 then EMN=rmfstr.EMN else ok='RMFSTR.EMN not found'
iEMX=where(strpos(rtags,'EMX') ge 0)
if iEMX[0] ne -1 then EMX=rmfstr.EMX else ok='RMFSTR.EMX not found'
iN_GRP=where(strpos(rtags,'N_GRP') ge 0)
if iN_GRP[0] ne -1 then N_GRP=rmfstr.N_GRP else ok='RMFSTR.N_GRP not found'
iF_CHAN=where(strpos(rtags,'F_CHAN') ge 0)
if iF_CHAN[0] ne -1 then F_CHAN=rmfstr.F_CHAN else ok='RMFSTR.F_CHAN not found'
iN_CHAN=where(strpos(rtags,'N_CHAN') ge 0)
if iN_CHAN[0] ne -1 then N_CHAN=rmfstr.N_CHAN else ok='RMFSTR.N_CHAN not found'
iMATRIX=where(strpos(rtags,'MATRIX') ge 0)
if iMATRIX[0] ne -1 then MATRIX=rmfstr.MATRIX else ok='RMFSTR.MATRIX not found'
iFIRSTCHAN=where(strpos(rtags,'FIRSTCHAN') ge 0)
if iFIRSTCHAN[0] ne -1 then FIRSTCHAN=rmfstr.FIRSTCHAN else ok='RMFSTR.FIRSTCHAN not found'
if ok ne 'ok' then begin
  message,ok,/informational
  if vv gt 0 then message,'Returning without further ado',/informational
  return,-1L
endif
;
areff=total(matrix,1)
nea=n_elements(effar) & if nea eq nnrg then areff=areff*effar
;
wts=ehi-elo
nsp=n_elements(spec) & if nsp eq nnrg then wts=wts*spec

;	make a full 2D image of the RMF
rmfimg=fltarr(nchan,nnrg) & szn=size(n_chan)
for inrg=0L,nnrg-1L do begin
  ngrp=n_grp[inrg] & jbeg=0
  for ig=0,ngrp-1 do begin
    if szn[0] gt 1 then begin
      ibeg=f_chan[ig,inrg] & iw=n_chan[ig,inrg]
    endif else begin
      ibeg=f_chan[inrg] & iw=n_chan[inrg]
    endelse
    if keyword_set(firstchan) then ibeg=ibeg-1        ;IDL index correction
    if iw gt 0 then rmfimg[ibeg:ibeg+iw-1,inrg]=matrix[jbeg:jbeg+iw-1,inrg]
    jbeg=jbeg+iw
  endfor
endfor

;	now compute the average weighted ARF in each detector channel
carf=fltarr(nchan)
for ic=0L,nchan-1L do carf[ic]=total(rmfimg[ic,*]*wts*areff)/total(rmfimg[ic,*]*wts)
chnrg=fltarr(nchan) & nrg=0.5*(elo+ehi)
for ic=0L,nchan-1L do chnrg[ic]=total(rmfimg[ic,*]*wts*nrg)/total(rmfimg[ic,*]*wts)

;	show
if vv gt 10 then begin
  plot,elo,effar,thick=3,line=1 & oplot,emn,carf,thick=3,line=2
  if vv gt 1000 then stop
endif

return,carf
end
