function extreffcorr,rmf,arffile,rmfstr=rmfstr,nrgar=nrgar,effar=effar,$
	verbose=verbose, _extra=e
;+
;function	extreffcorr
;	Chandra grating effective areas do not include the effect of the finite
;	size of the spectrum extraction region.  These are included within the
;	grating RMFs.  This function computes and returns a correction factor
;	array that can be multipled with the gARF (or does the correction and
;	returns the corrected effective area if the gARF is supplied).
;
;syntax
;	gcor=extreffcorr(rmf,arffile,rmfstr=rmfstr,nrgar=nrgar,effar=effar,$
;	verbose=verbose)
;
;parameters
;	rmf	[INPUT; required] RMF from which to figure out the
;		extraction efficiency correction
;		* can be either a filename, or an actual RMF as read in with rd_ogip_rmf()
;	arfile	[INPUT] name of gARF
;		* if given, reads in the ARF and modifies the output by the ARF
;		  such that the returned value is the gARF with extraction efficiency
;		  correction applied
;		* the energy grid in ARFFILE is what is used to define the spectrum
;		  array
;		* if not given, and also if NRGAR is not given, uses the ELO and EHI
;		  grid of the RMF
;
;keywords
;	rmfstr	[OUTPUT] if RMF is a file, the rmf that is read in is returned in
;		this variable
;		* if RMF is a structure, it, too, is returned in this variable
;	nrgar	[INPUT] if given, and ARFFILE is not given, defines the input flat
;		spectrum over this energy grid
;	effar	[INPUT] if given, and ARFFILE is not given, assumed to hold the
;		uncorrected ARF that must be modified with gCOR
;	verbose	[INPUT] controls chatter
;
;	_extra	[JUNK] here only to prevent crashing
;
;subroutines
;	CONV_RMF
;
;example
;	If you run chandra_repro on a grating observation, it will generate .arf and .rmf
;	files in repro/tg/.  Pick a pair and call them arffile and rmffile.  Then, the
;	plain extraction efficiencies are
;		gcor0=extreffcorr(rmffile,rmfstr=rmfstr) ; first run
;		gcor0alt=extreffcorr(rmfstr)	;subsequent runs
;		arfcor=extreffcorr(rmfstr,arffile)	;to automatically return corrected gARF
;	You can check that it worked by plotting them
;		arf=rdarf(arffile,arstr)	;to read in the gARF separately
;		arfarr=interpol(arf,arstr.ELO,rmfstr.ELO)	;just in case the grids differ
;		plot,12.398521/arstr.ELO,arf,xtitle='Wvl',ytitle='EA [cm!u2!n]'	;the untouched gARF
;		oplot,12.398521/rmfstr.ELO,arfarr*gcor0 	;corrected for extraction efficiency
;		oplot,12.398521/rmfstr.ELO,arfarr*gcor0alt	;should be identical to above
;		oplot,12.398521/rmfstr.ELO,arfcor       	;should be identical to above
;
;history
;	Vinay Kashyap (2015sep)
;-

;	usage
ok='ok' & np=n_params() & nr=n_elements(rmf) & szr=size(rmf,/type) & nrt=n_tags(rmf)
if np lt 1 then ok='Insufficient parameters' else $
 if nr eq 0 then ok='RMF is undefined' else $
  if nr gt 1 then ok='RMF must be scalar' else $
   if szr ne 7 and nrt eq 0 then ok='RMF must be filename or output of rd_ogip_rmf()'
if ok ne 'ok' then begin
  print,'Usage: gcor=extreffcorr(rmf,arffile,rmfstr=rmfstr,verbose=verbose)'
  print,'  compute extraction region correction for Chandra grating effective areas'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	grok the inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if szr eq 7 then rmfstr=rd_ogip_rmf(rmf[0]) else begin
  tnam=tag_names(rmf)
  o0=where(strpos(tnam,'ELO') ge 0 or $
     strpos(tnam,'EHI') ge 0 or $
     strpos(tnam,'MATRIX') ge 0 or $
     strpos(tnam,'EMN') ge 0 or $
     strpos(tnam,'EMX') ge 0 or $
     strpos(tnam,'F_CHAN') ge 0 or $
     strpos(tnam,'N_CHAN') ge 0,mo0)
  if mo0 ne 7 then begin
    message,'Input RMF does not appear to be an RMF',/informational
    return,-1L
  endif
  rmfstr=rmf
endelse
nrg=(rmfstr.ELO+rmfstr.EHI)/2.			;default grid
;
na=n_elements(arffile)
if na gt 0 then begin
  if na gt 1 then message,'only one ARF filename at a time, please',/informational
  sza=size(arffile[0],/type)
  if sza ne 7 then begin
    message,'ARFFILE is not a filename',/informational
    message,'taking energy grid from RMF',/informational
  endif else begin
    effar=rdarf(arffile[0],arfstr)
    nrgar=(arfstr.ELO+arfstr.EHI)/2.
  endelse
endif else begin
  ;	check whether nrgar and effar are passed in as keywords
  nnar=n_elements(nrgar) & nnea=n_elements(effar)
  if nnar gt 1 then begin
    if nnar ne nnea then effar=0.*nrgar+1.0
  endif else begin
    nrgar=nrg & effar=0.*nrgar+1.			;default flat effective area
  endelse
endelse
flx=0.*nrg+1.
if not keyword_set(nrgar) then begin			;flat EFFAR by default
  nrgar=nrg
  effar=0.*nrgar+1.
endif

;	convolve flat spectrum with RMF
conv_rmf,nrg,flx,chan,gcor,rmfstr,verbose=vv

;	correct EFFAR
tmp=(interpol(effar,nrgar,nrg)>0)<(max(effar))
gcor=gcor*tmp

if vv gt 1000 then stop,'halting; type .CON to continue'

return,gcor
end
