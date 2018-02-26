function potem_tool,idstr,flux,logT=logT,abund=abund,ldir=ldir,dwvl=dwvl,$
	okmult=okmult,kettle=kettle,outid=outid,multid=multid,$
	wdth=wdth,verbose=verbose, _extra=e
;+
;function	potem_tool
;	return emission measures EM(LOGT,WVL) for all lines that have been ID'd
;	(but see keyword OKMULT)
;
;	for each temperature, assume a delta-function emission measure,
;	compute line fluxes seen through some instrument, account for
;	interstellar absorption, and scale to match the observed flux.
;	
;	in case of multiple IDs for a single observed line, the emissivities
;	of the IDs are added together, and assumed for all practical purposes
;	to be the strongest line modified by the other IDs.
;	in case of IDs spanning different elements, the emissivities are
;	added as before, but weighted with their relative abundances.
;
;syntax
;	EM=potem_tool(idstr,flux,logT=logT,abund=abund,ldir=ldir,$
;	dwvl=dwvl,/okmult,/kettle,outid=outid,wdth=wdth, pres=pres,$
;	logP=logP,n_e=n_e, chifil=chifil,chidir=chidir,eqfile=eqfile,$
;	NH=NH,/noph,defEM=defEM, thresh=thresh,/temp,/ikev,effar=effar,$
;	wvlar=wvlar,fh2=fh2,He1=He1,HeII=HeII,/fano, level=level,$
;	verbose=verbose)
;
;parameters
;	idstr	[INPUT; required] the ID structure that comes out of LINEID
;	flux	[INPUT] observed fluxes
;		* if not given, taken from IDSTR.(*).FLUX
;		  (if zero, forcibly set to 1 [(ph|erg)/s/cm^2])
;		* size is forced to match IDSTR.WVL
;		  if <, filled out with 1s
;		  if >, extra elements are ignored
;
;keywords
;	logT	[INPUT] temperatures at which to compute the emission measures
;		* default is findgen(81)*0.05+4.
;	abund	[INPUT] abundances
;		* default is to use Anders & Grevesse
;	ldir	[INPUT] string array containing line database directory
;		names, one for each ID.
;		* if given, will NOT use the embedded emissivities in
;		  IDSTR.EMIS (default is to use them)
;		* if size < N(ALLMATCHES), the last element is replicated
;		  as often as necessary
;		* to selectively use IDSTR.(#).EMIS, set to '0' within a
;		  string array (e.g., ['$CHIANTI','$SPEX','0','0','$SCAR']
;		* if set to a number, assumes "$SCAR"
;	dwvl	[INPUT; default=5e-3] search the line database for the
;		wavelength of interest with this assumed roundoff error
;	okmult	[INPUT] if set, then multiple IDs (if any) are >not<
;		collapsed, and the output will be of the form EM(LOGT,ALLWVL)
;		* automatically set if size of FLUX matches the total number
;		  of matches
;	kettle	[INPUT] if set, does NOT do Pottasch corrections
;	outid	[OUTPUT] if defined on input, returns the IDs with new
;		emissivities in an ID structure
;		* multiple IDs are "collapsed" unless OKMULT is set
;	multid	[OUTPUT] number of IDs for each observed line (this is
;		of use in, e.g., PLOTEM)
;	wdth	[OUTPUT] the full width in logT of the emissivity curve,
;		as calculated by POTTASCH()
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] use to pass defined keywords to subroutines:
;		RD_LINE: PRES,LOGP,N_E,VERBOSE
;		FOLD_IONEQ: CHIFIL,CHIDIR,EQFILE,VERBOSE
;		FLUX_TO_EM: NH,NOPH,DEFEM,THRESH
;		LINEFLX: TEMP,IKEV,EFFAR,WVLAR
;		ISMTAU: FH2,HE1,HEII,FANO
;		POTTASCH: LEVEL
;
;subroutines
;	INICON
;	RD_LIST
;	    RD_LINE
;	    FOLD_IONEQ
;	        READ_IONEQ
;	            RD_IONEQ (CHIANTI routine)
;	FLUX_TO_EM
;	    LINEFLX
;	    ISMTAU
;	POTTASCH
;
;history
;	vinay kashyap (Nov98; "concatenated" id_to_emis and flux_to_em)
;	changed default behavior of FLUX; added keyword MULTID; changed
;	  OUTID.().Z and OUTID.().ION to floats in case multiple IDs have
;	  been collapsed to one emissivity (VK; Dec98)
;	added field FLUXERR to output IDSTR (VK; Mar99)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	no need to read in ion balance info in POTTASCH again (VK; JanMMI)
;	added keyword WDTH to carry the emissivity width info outside
;	  (VK; SepMMVII)
;-

;	usage
nid=n_tags(idstr)
if nid eq 0 then begin
  print,'Usage: EM=potem_tool(idstr,flux,logT=logT,abund=abund,ldir=ldir,$'
  print,'       dwvl=dwvl,/okmult,/kettle,outid=outid,wdth=wdth,
  print,'       RD_LINE:pres,logp,n_e; FOLD_IONEQ:chifil,chidir,eqfile;'
  print,'       FLUX_TO_EM:NH,noph,defEM,thresh; LINEFLX:temp,ikev,effar,wvlar;'
  print,'       ISMTAU:fh2,He1,HeII,fano; POTTASCH:level,'
  print,'       verbose=verbose)
  print,'  constructs emission measure predictions for each observed line'
  return,[-1L]
endif

;	trivial errors check
nf=n_elements(flux) & nT=n_elements(logT)
tid=tag_names(idstr) & owvl=idstr.(0) & nwvl=n_elements(owvl)
ok='ok'
if tid(0) ne 'WVL' then ok='input structure of unknown form' else $
 if nwvl ne nid-1 then ok='ID structure not in standard format' else $
  if tid(1) eq 'WVL_COMMENT' then ok='ID structure contains no data'
if ok ne 'ok' then begin
  em=[-1.D] & if nT gt 0 then em=dblarr(nT)-1
  message,ok,/info & return,em
endif

;	how many matched wavelengths, really?
multid=intarr(nwvl) & mwvl=0L
for i=0,nwvl-1 do multid(i)=n_elements(IDSTR.(i+1).WVL)
mwvl=long(total(multid))

;	are the fluxes given?
fx=fltarr(nwvl)
for i=0,nwvl-1 do begin		;{default is to use fluxes from IDSTR
  tmpid=idstr.(i+1) & fx(i)=total(tmpid.FLUX)
  if fx(i) eq 0 then fx(i)=1.	;if zero, set to 1
endfor				;I=0,NWVL-1}
if nf gt 0 and nf lt nwvl then fx(0:nf-1)=flux
if nf ge nwvl then begin
  if nf eq mwvl then begin
    okmult=1
    k=0L
    for i=0,nwvl-1 do begin
      tmpid=idstr.(i+1) & mw=multid(i)
      fx(i)=total(flux(k:k+mw-1L))
      k=k+mw
    endfor
  endif else fx=flux(0:nwvl-1)
endif

;	keywords
if nT eq 0 then begin			;define temperature grid
  nT=81L & logT=findgen(nT)*0.05+4.
endif
;
if n_elements(abund) lt 30 then abund=getabund('anders & grevesse')
;
if not keyword_set(dwvl) then dwvl=5e-3
;
;	define line database directories
ldbdir=strarr(mwvl)+'0' & ndir=n_elements(ldir)
if ndir ne 0 then begin
  szd=size(ldir) & nszd=n_elements(szd) & ndir=szd(nszd-1)
  if szd(nszd-2) ne 0 then ldbdir(*)='$SCAR'	;e.g., if /LDIR is set
  if szd(nszd-2) eq 7 then begin
    kdir=0L
    for i=0,nwvl-1 do begin
      mw=n_elements(idstr.(i+1).wvl)
      if ndir eq nwvl then ldbdir(kdir:kdir+mw-1)=ldir(i) else $
	ldbdir(kdir)=ldir([kdir])
      kdir=kdir+mw
    endfor
  endif
endif
;
v=0 & if keyword_set(verbose) then v=long(verbose(0))>1

;	initialize
atom=1 & rom=1 & inicon,atom=atom,roman=rom & atom=['X',atom] & rom=['',rom]
		;atomic symbols and ionization states, allow for "Unknown"
sep='	'				;separator, here a <tab>
if arg_present(outid) then outid=create_struct('WVL',owvl)

;	define the output
if keyword_set(okmult) then em=dblarr(nT,mwvl) else em=dblarr(nT,nwvl)

;	now step through the IDs and compute the emissivities
;basically, for each observed wavelength, get the IDs, and get the
;line emissivity for each line, then fold in ion balance, then
;interpolate to the required grid, then call FLUX->EM to get the emission
;measures, then apply the pottasch correction.
kwvl=0L
for i=0,nwvl-1 do begin
  idi=idstr.(i+1)
  ww=idi.WVL & zz=idi.Z & jon=idi.ION & labl=idi.LABL & mw=n_elements(ww)
  fjmx=dblarr(mw)
  if v gt 0 then print,'ID '+strtrim(i+1,2)+':'
  if v gt 5 then junk=cat_id(idstr,pick=[i])

  if zz(0) eq 0 or jon(0) eq 0 or labl(0) eq 'Unknown' then begin
    ;	ID is unknown -- nothing to calculate!
    em(*,kwvl)=dblarr(nT)-1.
    if keyword_set(outid) then outid=create_struct(outid,'ID'+strtrim(i,2),idi)
  endif else begin			;(compute EMs
    if keyword_set(okmult) then emis=dblarr(nT,mw) else emis=dblarr(nT)
    for j=0,mw-1 do begin		;{for each match
      ;
      if ldbdir(kwvl+j) eq '0' then begin	;(use IDSTR.(#).EMIS
	nam_idi=tag_names(idi)
	oi=where(nam_idi eq 'EMIS',moi)
	if moi eq 0 then ldbdir(kwvl+j)='$SCAR' else begin
	  ff=idi.EMIS & tlog=idi.LOGT & ff=reform(ff(*,j))
	endelse
      endif					;use IDSTR.(#).EMIS)
      if ldbdir(kwvl+j) ne '0' then begin	;(go read databse
	lnlst=atom(zz(j))+rom(jon(j))+sep+strtrim(ww(j),2)+' ? '+$
		strtrim(dwvl,2)+sep+ldbdir(kwvl+j)
	lstr=rd_list(lnlst,sep=sep,/incieq, _extra=e)
	ff=lstr.LINE_INT & tlog=lstr.LOGT
        if n_elements(lstr.WVL) gt 1 then begin
	  ;		too many?  collapse them
	  ff=0*logT
	  for k=0,n_elements(lstr.WVL)-1 do ff=ff+reform((lstr.LINE_INT(*,k)))
        endif
      endif					;read databse)
      ;
      if ff(0) le -1 then begin
	;		if line is missing, no help for it
	tlog=logT & ff=0*logT
      endif
      ;	interpolate emissivities from TLOG to LOGT
      fmin=min(ff,max=fmax) & fjmx(j)=fmax
      tmp=((interpol(ff,tlog,logT)>(fmin))<(fmax))
      if keyword_set(okmult) then emis(*,j)=tmp(*) else $
	emis(*)=emis(*)+abund(zz(j)-1)*tmp(*)
    endfor				;J=0,MW-1}

    fjmax=max(fjmx,jmx) & tfjmx=total(fjmx)
    if keyword_set(okmult) then begin
      emis=emis;/abund(zz(jmx)-1)
      if tfjmx gt 0 then fxj=fx(i)*fjmx/total(fjmx) else $
	fxj=fltarr(mw)+fx(i)/mw
    endif else begin
      emis=emis/abund(zz(jmx)-1)
      fxj=fx(i)
      if tfjmx gt 0 then ww=total(fjmx*ww)/tfjmx else ww=total(ww)/mw
      if mw gt 1 then zz=[float(zz(jmx))] else zz=[zz(jmx)]
      if mw gt 1 then jon=[float(jon(jmx))] else jon=[jon(jmx)]
    endelse
    tmp=flux_to_em(emis,flux=fxj,logT=logT,wvl=ww,Z=zz, abund=abund, _extra=e)
    if keyword_set(okmult) then em(*,lindgen(mw)+kwvl)=tmp else $
	em(*,i)=tmp
    if keyword_set(kettle) then pot=0*zz+1. else $
    	pot=pottasch(emis,logT,wdth=wdth, _extra=e)
    	;pot=pottasch(emis,logT,Z=zz,ion=jon, _extra=e)
    if v ge 10 then print,'Pottasch correction: ',pot
    npot=n_elements(pot)
    for k=0,npot-1 do begin
      if keyword_set(okmult) then begin
        oo=where(em(*,kwvl+k) ge 0,moo)
        if moo gt 0 and pot(k) gt 0 then em(oo,kwvl+k)=em(oo,kwvl+k)/pot(k)
      endif else begin
        oo=where(em(*,i) ge 0,moo)
        if moo gt 0 and pot(k) gt 0 then em(oo,i)=em(oo,i)/pot(k)
      endelse
    endfor
  endelse				;ID known)
  if arg_present(outid) then begin
    idwvl=create_struct('WVL',ww,'Z',zz,'ION',jon,'LABL',labl,'FLUX',fxj,$
	'FLUXERR',0*fxj,'LOGT',logt,'EMIS',emis)
    outid=create_struct(outid,'ID'+strtrim(i+1,2),idwvl)
  endif
  kwvl=kwvl+mw
endfor

return,em
end
