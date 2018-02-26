pro conv_rmf,nrg,flx,chan,spec,rmf,effar=effar,nrgar=nrgar,rmfstr=rmfstr,$
	verbose=verbose, _extra=e
;+
;procedure	conv_rmf
;	convolves input energy spectrum with the response matrix
;
;	unlike the combination of RDRESP and FOLD_RESP, this one goes
;	easy on memory at the expense of speed.
;
;syntax
;	conv_rmf,nrg,flx,chan,spec,rmf,effar=effar,nrgar=nrgar,$
;	rmfstr=rmfstr,fchcol=fchcol,/shift1,verbose=verbose
;
;parameters
;	nrg	[INPUT; required] mid-bin values at which input FLX is
;		defined.
;		* usually [keV], depends on RMF file
;		* if size = N(FLX)+1, assumed to be bin-boundaries
;	flx	[INPUT; required] fluxes in [ph/bin/(s/cm^2/sr/...)]
;	chan	[OUTPUT; required] bin boundary values at which output SPEC
;		is defined
;		* same units as NRG
;		* WARNING: the values that populate this array are taken
;		  from the EMN and EMX values in the RMF.  These are not
;		  exact values, and simply represent approximate bounds
;		  on the detector channels that the input spectrum gets
;		  redistributed into.  They are to be used mainly for
;		  plotting purposes, and should not be used for scientific
;		  purposes unless the user recognizes their limitations.
;		  In other words, the array indices of SPEC are the
;		  detector channels, and the CHAN array is a mapping from
;		  channels to energy, as defined in the EBOUNDS extension
;		  of the RMF
;	spec	[OUTPUT; required] output, convolved spectrum
;	rmf	[INPUT; required] name of OGIP-compliant FITS file containing
;		the response matrix
;		* may also be a structure, such as the output of RD_OGIP_RMF()
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
;	rmfstr	[OUTPUT] structure containing info of response matrix,
;		the output of RD_OGIP_RMF()
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines:
;		RD_OGIP_RMF: FCHCOL, SHIFT1
;		REBINW : SLOWOK
;
;restrictions
;	requires IDLAstro library
;	requires RD_OGIP_RMF()
;	requires REBINW()
;	requires BINERSP()
;
;these are all the smoothing tools in PINTofALE
;	ALALOESS() makes a loess curve
;	CLTSMOOTH() removes -ves from background subtracted spectra or light curves
;	CONV_RMF convolves with specified RMF
;	HAARTRAN() smoothes by filtering on Haar wavelet coefficients
;	LINEREM() removes lines from spectrum by iterative S/N filtering
;	NOISMOOTH() does boxcar accumulation a la Ebeling's asmooth
;	REGROUP() accumulates from one end
;	SCRMF does fast convolution using optimized RMF
;	SMOOTHIE does peak-descent accumulation to build up S/N
;	SPLAC computes a segmented piecewise linear approximations
;	UNKINK() removes sharp kinks from a curve
;	VARSMOOTH() does local smoothing at specified scales
;	VOORSMOOTH() does fast local smoothing at specified scales
;
;history
;	vinay kashyap (Apr99)
;	added keyword SHIFT1 (VK; JanMMI)
;	bug: incorrect usage of WVLAR instead of NRGAR in effar block
;	  (Erica R; Aug01)
;	deleted keyword SHIFT1, cleaned up and speeded up (VK; Nov2001)
;	bug: when spectrum size < RMF size and RMF has multiple groups,
;	  was using only the first group (VK; Nov'02)
;	bug: was crashing when NRG was outside the bounds of the RMF
;	  in cases where the RMF was at a higher resolution (LL/VK; Apr'03)
;	added keyword VERBOSE (VK; MarMMV)
;	bug: when input spectrum energy range was smaller than RMF range,
;	  was assuming RMF was at higher resolution even if it wasn't
;	  (VK; Aug'07)
;	bug: wasn't checking for frequency beating between input spectrum
;	  grid and input RMF grid in above case; now if it does, and calls
;	  REBINW() and forces input spectrum to the right grid (VK; Aug'07)
;	bug: EFFAR multiplication was being ignored when input spectrum had
;	  to be rebinned to match RMF grid (VK; Oct'07)
;-

;	usage
ok='ok'
np=n_params() & mnrg=n_elements(nrg) & mflx=n_elements(flx)
szr=size(rmf) & nszr=n_elements(szr)
nrmf=n_elements(rmf) & mrmf=n_tags(rmf)
if np lt 5 then ok='missing parameters' else $
 if mnrg eq 0 then ok='NRG is undefined: require spectrum grid' else $
  if mflx eq 0 then ok='FLX is undefined: require spectrum' else $
   if mnrg ne mflx and mnrg ne (mflx+1L) then ok='input spectrum garbled' else $
    if nrmf eq 0 then ok='RMF is undefined' else $
     if szr[nszr-2] ne 7 and mrmf eq 0 then $
     ok='RMF must be either a filename or a response matrix structure '
if ok ne 'ok' then begin
  print,'Usage: conv_rmf,nrg,flx,chan,spec,rmf,effar=effar,nrgar=nrgar,$'
  print,'       rmfstr=rmfstr,fchcol=fchcol,/shift1,verbose=verbose'
  print,'  convolve spectrum with response'
  if np gt 0 then message,ok,/info
  return
endif

;	verbosity
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

;	check inputs
ee=[nrg[*]] & fe=[flx[*]]
if mnrg eq mflx+1 then begin
  ;	assume that NRG are bin boundaries and obtain the mid-bin values
  ee=[0.5*(nrg[1:*]+nrg)]
  mnrg=mnrg-1L
  if mnrg ne mflx then message,'BUG!'
endif
;if mnrg eq mflx then begin
;  ;	assume input NRG are mid-bin values and decipher the bin boundaries
;  de=nrg[1:*]-nrg
;  ee=[nrg-0.5*de,nrg[mnrg-1]+0.5*de[mnrg-2]*[-1,1]]
;  mnrg=mnrg+1L
;  if mnrg ne mflx+1L then message,'BUG!'
;endif

;	include effective areas?
mnar=n_elements(nrgar) & mear=n_elements(effar)
if mnar eq 0 and mear ne 0 then begin
  message,'WARNING: EFFAR is given, but NRGAR is not; ignoring EFFAR',$
	/informational
endif
if mnar eq mear and mear gt 0 then begin	;(effective area defined
  arnrg=ee & areff=0.*arnrg+1.
  if mnar eq 1 then areff[*]=effar[0] else begin
    minea=min(effar,max=maxea) & minna=min(nrgar,max=maxna)
    areff=(interpol(effar,nrgar,ee) > (minea)) < (maxea)
    oo=where(arnrg lt minna or arnrg gt maxna,moo)
    if moo gt 0 then areff[oo]=0.
  endelse
  fe=fe*areff	;assumes input FLX is in units of [ph/bin/cm^2/...]
endif						;EFFAR(NRGAR))

;	read in response matrix
if mrmf eq 0 then rmfstr=rd_ogip_rmf(rmf[0]) else rmfstr=rmf
if n_tags(rmfstr) eq 0 then return
;
NNRG=rmfstr.NNRG & ELO=rmfstr.ELO & EHI=rmfstr.EHI
N_GRP=rmfstr.N_GRP & F_CHAN=rmfstr.F_CHAN & N_CHAN=rmfstr.N_CHAN & FIRSTCHAN=rmfstr.FIRSTCHAN
MATRIX=rmfstr.MATRIX
NCHAN=rmfstr.NCHAN & EMN=rmfstr.EMN & EMX=rmfstr.EMX
szn=size(N_CHAN)

;	define outputs
chan=[EMN,EMX[nchan-1L]] & spec=fltarr(nchan)

;	find binning for input NRG
;ie=binersp(ee,bstr=rmfstr)
ie=binersp(ee,elo,ehi)

if keyword_set(FIRSTCHAN) and vv gt 5 then message,$
	'assuming that input RMF is 1-based, not 0-based',/info

;	which is at the higher resolution? the spectrum or the RMF grid?
spechires=1	;the default
if min(nrg) lt min(nnrg) and max(nrg) gt max(nnrg) and mflx le nnrg then spechires=0
if min(abs(nrg[1:*]-nrg)) gt max(abs(rmfstr.EHI-rmfstr.ELO)) then spechires=0

;	distribute each input FLX into PHAs
;rspimg=fltarr(nchan,nnrg)
nnoi=intarr(nnrg) & ffe=fltarr(nnrg)
if spechires gt 0 then begin	;(spectrum at higher resolution than RMF
  ;	check if there will be problems with bin sizes
  ;	and frequency beating between spectrum grid and
  ;	RMF grid, and if there is, rebin the input spectrum
  ;	to the RMF grid
  deltaspec=median(abs(nrg[1:*]-nrg)) & deltaRMF=median(EHI-ELO)
  specrmf=deltaspec/deltaRMF
  if (specrmf-fix(specrmf)) gt 1./float(nnrg) then begin
    if vv gt 0 then message,'input spectrum will be rebinned to match the RMF grid',/informational
    newnrg=[ELO,max(EHI)] & newnrg=newnrg[sort(newnrg)]
    ;newflx=rebinw(flx,nrg,newnrg,/perbin, _extra=e)
    newflx=rebinw(fe,ee,newnrg,/perbin, _extra=e)
    ee=newnrg & if n_elements(newnrg) eq n_elements(newflx)+1 then ee=0.5*(newnrg[1:*]+newnrg)
    fe=newflx
    ie=binersp(ee,ELO,EHI)
  endif
  for inrg=0L,nnrg-1L do begin
    rsp=fltarr(nchan)
    ngrp=n_grp[inrg] & jbeg=0
    for ig=0L,ngrp-1L do begin
      if szn[0] gt 1 then begin
	ibeg=f_chan[ig,inrg] & iw=n_chan[ig,inrg]
      endif else begin
	ibeg=f_chan[inrg] & iw=n_chan[inrg]
      endelse
      if keyword_set(firstchan) then ibeg=ibeg-1	;IDL index correction
      if iw gt 0 then rsp[ibeg:ibeg+iw-1L]=matrix[jbeg:jbeg+iw-1L,inrg]
      jbeg=jbeg+iw
    endfor
    ;
    oi=where(ie eq inrg,moi)
    nnoi[inrg]=moi	;diagnostic array
    if moi gt 0 then spec=spec+rsp*total(fe[oi])
    if moi gt 0 then ffe[inrg]=total(fe[oi])	;diagnostic array
    ;if moi gt 0 then begin
    ;  if vv gt 500 then begin
    ;	 opl=where(spec gt 0.01*max(spec),mopl)
    ;	 if mopl gt 0 then begin
    ;	   plot,spec[opl],title=inrg,/xs,psym=10
    ;	   oplot,rsp[opl]*max(spec[opl])/max(rsp[opl]),line=1,psym=10
    ;	   wait,0.2
    ;	   cc=get_kbrd(0)
    ;	   if cc eq 'x' then stop,'HALTing; type .CON to continue'
    ;	 endif
    ;  endif
    ;endif
    ;rspimg[*,inrg]=rsp
  endfor
endif else begin		;SPECHIRES>0)(RMF at higher resolution than spectrum
  ;the following is not necessary for the case where the input
  ;spectrum is at a lower resolution, because each bin of the
  ;input _will_ be accounted for anyway.  problem will arise
  ;if the RMF varies radically across a bin, but if that is so,
  ;the user should input a spectrum at better resolution, and
  ;it is not our job to guess.
  ;{
  ;;	check if there will be problems with bin sizes
  ;;	and frequency beating between spectrum grid and
  ;;	RMF grid, and if there is, rebin the input spectrum
  ;;	to the RMF grid
  ;deltaspec=median(abs(nrg[1:*]-nrg)) & deltaRMF=median(EHI-ELO)
  ;specrmf=deltaspec/deltaRMF
  ;if (specrmf-fix(specrmf)) gt 1./float(nnrg) then begin
  ;  if vv gt 0 then message,'input spectrum will be rebinned to match the RMF grid',/informational
  ;  newnrg=[ELO,max(EHI)] & newnrg=newnrg[sort(newnrg)]
  ;  newflx=rebinw(flx,nrg,newnrg,/perbin, _extra=e)
  ;  ee=newnrg & if n_elements(newnrg) eq n_elements(newflx)+1 then ee=0.5*(newnrg[1:*]+newnrg)
  ;  fe=newflx
  ;  ie=binersp(ee,ELO,EHI)
  ;endif
  ;}
  oi=where(ie ge 0,moi)
  for i=0L,moi-1L do begin
    ii=oi[i]
    inrg=ie[ii] & rsp=fltarr(nchan)
    ngrp=n_grp[inrg] & jbeg=0
    for ig=0L,ngrp-1L do begin
      if szn[0] gt 1 then begin
	ibeg=f_chan[ig,inrg] & iw=n_chan[ig,inrg]
      endif else begin
	ibeg=f_chan[inrg] & iw=n_chan[inrg]
      endelse
      if keyword_set(firstchan) then ibeg=ibeg-1	;IDL index correction
      if iw gt 0 then rsp[ibeg:ibeg+iw-1L]=matrix[jbeg:jbeg+iw-1L,inrg]
      jbeg=jbeg+iw
    endfor
    ;
    spec=spec+fe[ii]*rsp
  endfor
endelse				;SPECHIRES=0)

;if vv gt 1000 then stop,'HALTing; type .CON to continue'

return
end
