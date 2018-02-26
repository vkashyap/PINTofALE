pro licospec,wgrid,lspec,cspec,at_logT,lrf=lrf,lstr=lstr,cstr=cstr,$
	lemis=lemis,cemis=cemis,abund=abund,lcthr=lcthr,constr=constr,$
	ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,verbose=verbose,$
	allah=allah,_extra=e
;+
;procedure	licospec
;	compute the contributions of line and continuum
;	emissivities over a specified wavelength grid.
;
;syntax
;	licospec,wgrid,lspec,cspec,at_logT,lrf=lrf,lstr=lstr,cstr=cstr,$
;	lemis=lemis,cemis=cemis,abund=abund,lcthr=lcthr,constr=constr,$
;	ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,verbose=verbose,$
;	pres=pres,logP=logP,n_e=n_e,chifil=chifil,chidir=chidir,$
;	eqfile=eqfile,logT=logT,DEM=DEM,effar=effar,wvlar=wvlar,$
;	type=type,/norm,/fwhm,betap=betap,angle=angle,vrot=vrot
;
;parameters
;	wgrid	[INPUT; required] wavelength bin boundaries over
;		which to compute line and continuum contributions
;		* must be in same units as the wavelengths stored
;		  in LDBDIR, i.e., Angstroms
;	lspec	[OUTPUT; required] line emission at each wavelength bin
;	cspec	[OUTPUT; required] continuum emission at each wavelength bin
;		* note that these are essentially spectra, and include
;		  the effects of specified abundance, DEM, and any input
;	 	  RMFs, LRFs, or effective areas.
;	at_logT	[INPUT] if set to a scalar, assumes a delta-function
;		emission measure to compute the appropriate fluxes;
;		the default is to assume a DEM that covers the range
;		4..8 in logT (but note that any arbitrary DEM can be
;		defined using keywords LOGT and DEM).
;		* NOT TESTED
;
;keywords
;	lrf	[INPUT] if scalar, the width of the line response function
;		NOTE: the function used to define the LRF can be set
;		using the keyword TYPE (e.g., TYPE='gauss', TYPE='beta=2.5')
;		* if float, the width in the same units as WGRID
;		* if integer, the width in number of bins
;		* if RMF structure (e.g., output of RD_OGIP_RMF()), use
;		  it to convolve the output and place in channels defined
;		  by WGRID indices (i.e., if LRF is defined as an RMF,
;		  it _must_ match the wavelength grid)
;	lstr	[I/O] line emissivity structure out of RD_LINE()*FOLD_IONEQ()
;		* if defined on input, RD_LINE() and FOLD_IONEQ() will not
;		  be called unless bounds of WGRID spill out, in which case
;		  the user will be asked what to do
;		* unlike LEMIS below, abundances are _not_ included
;		  (unless in case of APED)
;	cstr	[I/O] continuum emissivity structure out of RD_CONT() --
;		* if defined on input, RD_CONT() will not be called unless
;		  bounds of WGRID spill out, in which case the user will be
;		  asked what to do
;	ldbdir	[INPUT; '$CHIANTI'] line emissivity database directory
;	cdbdir	[INPUT; '$CONT'] continuum emissivity database directory
;	ceroot	[INPUT; 'cie'] prefix of files containing the continuum
;		emissivities
;	abund	[INPUT] element abundances to use while computing the
;		continuum emissivities, the line and continuum fluxes,
;		combining line emissivities into wavelenth bins, etc.
;		* if not defined, calls GETABUND(!ABREF)
;		* if !ABREF is not defined, adopts Grevesse et al.
;	lemis	[OUTPUT] line emissivities recast to be on the input
;		wavelength grid, i.e., an array of size N(LOGT)xN(WGRID)
;		and with units [1e-23 ergs cm^3/s]
;		* the output line emissivities _include_ the specified
;		  ion balances and abundances
;	cemis	[OUTPUT] continuum emissivities recast to be on the input
;		wavelength grid, i.e., an array of size N(LOGT)xN(WGRID)
;		and with units [1e-23 ergs cm^3/s]
;		* NOTE: the /Ang usually present in continuum emissivities
;		  is missing here
;	lcthr	[INPUT] if set, combines the continuum emissivity from all
;		the bins with LSPEC/CSPEC less than LCTHR and returns it
;		in CONSTR
;		* if negative, then combines all bins with CSPEC < LCTHR
;	constr	[OUTPUT] continuum emissivity structure with the emissivities
;		marginalized over wavelength
;       allah   [INPUT] if set, no bounds of WGRID are not checked against 
;               LSTR and the structure is used as is 
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		RD_LINE: PRES, LOGP, N_E
;		FOLD_IONEQ: CHIFIL, CHIDIR, EQFILE
;		RD_CONT: PRES, LOGP, N_E, ABUND
;		LINEFLX: LOGT, DEM, ABUND, EFFAR, WVLAR
;		LIBMODEL: TYPE, NORM, FWHM, BETAP, ANGLE, VROT
;
;example
;	;	make a wavelength grid
;	wgrid=findgen(2001)*0.001+18. & lrf=0.02 & type='beta=2.5'
;	;	call LICOSPEC
;	licospec,wgrid,lspec,cspec,verbose=10,lstr=lstr,cstr=cstr,$
;	lcthr=0.1,constr=constr,lemis=lemis,cemis=cemis,n_e=1e10,$
;	lrf=lrf,type=type
;	;	plot line to continuum flux ratio
;	plot,wgrid,lspec/cspec,/ylog,$
;	xtitle='wavelength [Ang]',ytitle='line-to-continuum'
;	;	plot line and continuum integrated emissivities
;	tmpl=dblarr(81) & for i=0,80 do tmpl[i]=total(lemis[i,*])
;	tmpc=dblarr(81) & for i=0,80 do tmpc[i]=total(cemis[i,*])
;	plot,lstr.LOGT,tmpl,/yl & oplot,cstr.LOGT,tmpc,thick=3,line=2
;	;	integrated emissivities over a sub range
;	oo=where(lspec/cspec lt 0.1)
;	tmpl=dblarr(81) & for i=0,80 do tmpl[i]=total(lemis[i,oo])
;	tmpc=dblarr(81) & for i=0,80 do tmpc[i]=total(cemis[i,oo])
;	plot,lstr.LOGT,tmpl,/yl & oplot,cstr.LOGT,tmpc,thick=3,line=2
;
;subroutines
;	CONV_RMF [RD_OGIP_RMF]
;	FOLD_IONEQ [RD_IONEQ [READ_IONEQ], WHEE]
;	GETABUND
;	HASTOGRAM [KILROY]
;	LIBMODEL [MK_GAUSS, MK_LORENTZ, MK_ROGAUSS [MK_GAUSS], MK_SLANT]
;	LINEFLX [GETABUND, WHEE]
;	MID2BOUND
;	RD_CONT [SETSYSVAL [LEGALVAR], SYMB2ZION [LAT2ARAB]]
;	RD_LINE [INICON, SETSYSVAL [LEGALVAR], SYMB2ZION [LAT2ARAB]]
;	REBINW [FINDEX]
;	SETSYSVAL [LEGALVAR]
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Dec2003)
;	added keywords LEMIS,CEMIS,LCTHR,CONSTR,ABUND,LDBDIR,CDBDIR,CEROOT;
;	  corrected bug with CSTR.CONT_INT being degraded in repeated calls;
;	  now -asks- if input wavelength range has changed and LSTR/CSTR do
;	  not span this range (VK; Jan2004)
;	cleaned up behavior of CSTR.CONT_INT for compatibility with RD_CONT();
;	  now output LEMIS and CEMIS add up in photon space (VK; Apr2004)
;       added ALLAH keyword to facilitate use of ZRESP keyword in
;         SOLAR_TRESP() (LL; Sep2005) 
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;-

;	usage
ok='ok' & np=n_params() & nw=n_elements(wgrid)
if np eq 0 then ok='Insufficient parameters' else $
 if np lt 1 then ok='Insufficient parameters: WGRID is missing' else $
  if np lt 2 then ok='Insufficient parameters: LSPEC is missing' else $
   if np lt 3 then ok='Insufficient parameters: CSPEC is missing' else $
    if nw eq 0 then ok='WGRID is undefined' else $
     if nw eq 1 then ok='WGRID must have bin _boundaries_'
if ok ne 'ok' then begin
  print,'Usage: licospec,wgrid,lspec,cspec,at_logT,lrf=lrf,lstr=lstr,cstr=cstr,$'
  print,'       lemis=lemis,cemis=cemis,abund=abund,lcthr=lcthr,constr=constr,$'
  print,'       verbose=verbose,$'
  print,'       pres=pres,logP=logP,n_e=n_e,chifil=chifil,chidir=chidir,$'
  print,'       eqfile=eqfile,logT=logT,DEM=DEM,effar=effar,wvlar=wvlar,$'
  print,'       type=type,/norm,/fwhm,betap=betap,angle=angle,vrot=vrot'
  print,'  compute the continuum and line contributions over a specified'
  print,'  wavelength grid'
  if np ne 0 then message,ok,/informational
  return
endif

;	figure out keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;	wavelength range
wmin=min(wgrid,max=wmax) & wmid=0.5*(wgrid[1:*]+wgrid)
;
if not keyword_set(ldbdir) then setsysval,'LDBDIR',ldbdir,/getval
if not keyword_set(cdbdir) then setsysval,'CDBDIR',cdbdir,/getval
if not keyword_set(ceroot) then setsysval,'CEROOT',ceroot,/getval
if not keyword_set(ldbdir) then ldbdir='$CHIANTI'
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(ceroot) then ceroot='cie'
if n_elements(abund) lt 30 then begin
  if not keyword_set(abref) then setsysval,'ABREF',abref,/getval
  if not keyword_set(abref) then abref='grevesse et al.'
  abund=getabund(abref)
endif
;	LSTR
if n_tags(lstr) ne 0 then begin
  lwmin=min(abs(lstr.WVL),max=lwmax)
  if not keyword_set(allah) then begin 
    if wmin lt lwmin or wmax gt lwmax then begin
      message,'Wvl range = '+strtrim(wmin,2)+'-'+strtrim(wmax,2),/informational
      message,'LSTR '+strtrim(lwmin,2)+'-'+strtrim(lwmax,2)+' does not cover input range',/informational
      c1='Reread line emissivities? [Y/n] > ' & read,prompt=c1,c1
      if c1 eq '' or strmid(strlowcase(strtrim(c1,2)),0,1) eq 'y' then lstr=0
    endif
  endif
endif
if n_tags(lstr) eq 0 then begin
  tmp=rd_line(atom,wrange=[wmin,wmax],dbdir=ldbdir,fstr=lstr,verbose=vv,$
	_extra=e)
  if ldbdir ne '$APED' then tmp=fold_ioneq(lstr.LINE_INT,lstr.Z,lstr.JON,logT=lstr.LOGT,verbose=vv,$
	_extra=e)
  lstr.LINE_INT=tmp
endif
LINE_INT=lstr.LINE_INT
;	CSTR
if n_tags(cstr) ne 0 then begin
  cwmin=min(abs(cstr.WVL),max=cwmax)
  if wmin lt cwmin or wmax gt cwmax then begin
    message,'CSTR '+strtrim(cwmin,2)+'-'+strtrim(cwmax,2)+' does not cover input range',/informational
    c1='Reread continuum emissivities? [Y/n] > ' & read,prompt=c1,c1
    if c1 eq '' or strmid(strlowcase(strtrim(c1,2)),0,1) eq 'y' then cstr=0
  endif
endif
if n_tags(cstr) eq 0 then begin
  tmp=rd_cont(ceroot,wrange=[wmin,wmax],dbdir=cdbdir,fcstr=cstr,verbose=vv,$
	abund=abund, _extra=e)
endif else begin
  ;	if CSTR is not the right type, quit.
  ;cww=mid2bound(cstr.midWVL) & dcww=cww[1:*]-cww
endelse
CONT_INT=cstr.CONT_INT
;	convert continuum emissivities from [.../A] to [.../bin]
cww=mid2bound(cstr.midWVL) & dcww=cww[1:*]-cww
tmp=CONT_INT
for i=0L,81L-1L do tmp[i,*]=tmp[i,*]*dcww & CONT_INT=tmp
;
contags=tag_names(cstr)
;	AT_LOGT
if n_elements(at_logT) gt 0 then begin
  logT=findgen(81)*0.05+4. & tmp=min(abs(logT-at_logT),iT)
  DEM=dblarr(81) & DEM[iT]=1d14
endif
;	LCTHR
lct=1d20
if is_keyword_set(lcthr) then lct=lcthr[0]

;	convolving function, if necessary
if n_tags(lrf) eq 0 then begin
  if keyword_set(lrf) then begin
    if size(lrf,/type) eq 1 then w=findgen(n_elements(wgrid)-2) else $
	w=wmid
    a=[0.5*(min(w)+max(w)),lrf[0],1] & libmodel,w[1:*],a,convf, _extra=e
    convf=float(convf)/total(convf)
  endif
endif

;	output emissivities
nT=n_elements(lstr.LOGT) & lemis=dblarr(nT,nw-1L) & cemis=lemis
abZ=abund[lstr.Z-1L]

;	get the fluxes	[ph/s/bin/...]
lfx=lineflx(LINE_INT,lstr.LOGT,lstr.WVL,lstr.Z,$
	logT=logT,DEM=DEM,abund=abund, _extra=e)
cfx=lineflx(CONT_INT,cstr.LOGT,cstr.midWVL,$
	logT=logT,DEM=DEM, _extra=e)

lin_erg2ph=6.626176d-27*2.9979d10*1d8/abs(lstr.WVL)
con_erg2ph=6.626176d-27*2.9979d10*1d8/abs(cstr.WVL)
;	rebin to the right grid
if n_tags(lrf) gt 0 then begin		;(RMF is available
  xgrid=lrf.elo
  if xgrid[0] gt xgrid[1] then xgrid=[max(lrf.ehi),xgrid] else $
	xgrid=[xgrid,max(lrf.ehi)]
  lsp=hastogram(12.3985/abs(lstr.WVL),xgrid,wts=lfx,verbose=vv)
  conv_rmf,xgrid,lsp,chan,ltmp,lrf
  lsp=rebinw(ltmp,12.3985/xgrid,wgrid,/perbin)
  lspec=lsp
  csp=rebinw(cfx,12.3985/cww,xgrid,/perbin,verbose=vv, _extra=e)
  conv_rmf,xgrid,csp,chan,ctmp,lrf
  csp=rebinw(ctmp,12.3985/xgrid,wgrid,/perbin)
  cspec=csp
  ;
  xmid=0.5*(xgrid[1:*]+xgrid)
  grid_ph2erg=6.626176d-27*2.9979d10*1d8/(12.3985/xmid)
  for iT=0L,nT-1L do begin
    if vv gt 0 then print,'logT = ',lstr.logT[iT]
    lsp=hastogram(12.3985/abs(lstr.WVL),xgrid,$
	wts=abZ*reform((LINE_INT)[iT,*])/lin_erg2ph,verbose=vv)
    conv_rmf,xgrid,lsp,chan,ltmp,lrf
    lsp=rebinw(ltmp,12.3985/xgrid,wgrid,/perbin)
    lemis[iT,*]=lsp*grid_ph2erg
    csp=rebinw(reform((CONT_INT)[iT,*])/con_erg2ph,12.3985/cww,xgrid,/perbin,$
	verbose=vv, _extra=e)
    conv_rmf,xgrid,csp,chan,ctmp,lrf
    csp=rebinw(ctmp,12.3985/xgrid,wgrid,/perbin)
    cemis[iT,*]=csp*grid_ph2erg
  endfor
endif else begin			;RMF)(no RMF
  lsp=hastogram(abs(lstr.WVL),wgrid,wts=lfx,verbose=vv)
  csp=rebinw(cfx,cww,wgrid,/perbin,verbose=vv, _extra=e)
  if keyword_set(lrf) then begin
    ltmp=convol(lsp,convf,/edge_truncate)
    ctmp=convol(csp,convf,/edge_truncate)
    lsp=ltmp & csp=ctmp
  endif
  lspec=lsp & cspec=csp
  ;
  grid_ph2erg=6.626176d-27*2.9979d10*1d8/wmid
  for iT=0L,nT-1L do begin
    if vv gt 0 then print,'logT = ',lstr.logT[iT]
    lsp=hastogram(abs(lstr.WVL),wgrid,$
	wts=abZ*reform((LINE_INT)[iT,*])/lin_erg2ph,verbose=vv)
    csp=rebinw(reform((CONT_INT)[iT,*])/con_erg2ph,cww,wgrid,/perbin,$
	verbose=vv, _extra=e)
    if keyword_set(lrf) then begin
      lsp=convol(lsp,convf,/edge_truncate)
      csp=convol(csp,convf,/edge_truncate)
    endif
    lemis[iT,*]=lsp*grid_ph2erg & cemis[iT,*]=csp*grid_ph2erg
  endfor
endelse					;no RMF)

;	filter over wavelength acc. to LCTHR
if is_keyword_set(lcthr) then begin
  if lcthr[0] gt 0 then oo=where(lspec/cspec le lct,moo) else $
	oo=where(cspec le abs(lct),moo)
  if moo eq 0 then begin
    message,'Line-to-Continuum threshold too stringent; including all',/informational
    moo=n_elements(cspec) & oo=lindgen(moo)
  endif
  ;	collapse the emissivity array
  temis=cemis[*,oo] & twvl=wmid[oo] & conemis=dblarr(nT)
  conwvl=total(twvl)/moo
  for iT=0L,nT-1L do conemis[iT]=total(temis[iT,*]*twvl)/conwvl
  wvl_range=minmax(wmid[oo]) & delta_wvl=wvl_range[1]-wvl_range[0]
  conemis=conemis/delta_wvl	;[1e-23 ergs cm^3/s/Ang]
  ;	construct the emissivity structure
  ;	NOTE: must be in exact same format as in RD_CONT()
  constr=create_struct('CONT_INT',conemis,'LOGT',cstr.LOGT,'WVL',[wmin,wmax],$
	  'TkeV',cstr.TkeV,'EkeV',12.3985/[wmin,wmax],'midWVL',conwvl)
endif

return
end
