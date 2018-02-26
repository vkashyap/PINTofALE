pro munge_list,list,idstr,lnstr,wvl,flux,fluxerr,dbdir=dbdir,sep=sep,$
	prefix=prefix,comment=comment, _extra=e
;+
;procedure	munge_list
;	chomps through ascii file containing observed wavelengths/fluxes
;	and possible matches to them, then generates the appropriate
;	line database and ID structures (see RD_LINE and LINEID)
;
;syntax
;	munge_list,list,idstr,lnstr,wvl,flux,fluxerr,dbdir=dbdir,$
;	sep=sep,prefix=prefix,comment=comment, mapping=mapping,pres=pres,$
;	logp=logp,n_e=n_e,desig=desig,econf=econf,chifil=chifil,chidir=chidir,$
;	eqfile=eqfile,DEM=DEM,abund=abund,/noph,effar=effar,wvlar=wvlar
;
;parameters
;	list	[INPUT; required] file containing wavelengths, fluxes, and
;		matches.
;		* see format below
;	idstr	[OUTPUT] LINEID type structure
;	lnstr	[OUTPUT] RD_LINE type structure
;	wvl	[OUTPUT] all the observed wavelengths
;	flux	[OUTPUT] all the measured wavelengths
;	fluxerr	[OUTPUT] errors on FLUX
;
;keywords
;	dbdir	[I/O] directory in which to look for files
;		* will be used ONLY if not given in LNLIST
;		* if not specified, *default* will be set at run-time to
;		  whatever was used last
;		* on output, will contain complete list of line database
;		  directories used
;	sep	[INPUT; default=<tab>] separator used to delineate fields
;		in LNLIST
;	prefix	[INPUT; default='*'] any line beginning with this letter
;		is considered a comment
;		* this is in addition to "%", "#", ";", "/"
;	comment	[OUTPUT] any commented lines are returned in this variable
;		as a string array
;	_extra	[INPUT ONLY] use this to pass defined keywords to
;		-- RD_LIST (MAPPING)
;		-- RD_LINE (PRES,LOGP,N_E,DESIG,ECONF,VERBOSE)
;		-- FOLD_IONEQ (CHIFIL,CHIDIR,EQFILE,VERBOSE)
;		-- LINEFLX (DEM,ABUND,NOPH,EFFAR,WVLAR)
;
;input file format
;	* one entry per line
;	* each observed wavelength and related IDs must be in a block
;	  of the form
;	  	/{ <sep> WVL <sep> FLUX <sep> FLUXERR
;	  	ELEM1 ION1 <sep> WAVE1 <sep> SOURCE1 <sep> DESCRIPTION1
;	  	ELEM2 ION2 <sep> WAVE2 <sep> SOURCE2 <sep> DESCRIPTION2
;	  	...
;	  	ELEMn IONn <sep> WAVEn <sep> SOURCEn <sep> DESCRIPTIONn
;		// human readable comments
;	  	/}
;	* the beginning and ending "/{" and "/}" are REQUIRED
;	* WVL, FLUX, and FLUXERR are optional, and will be guessed at from
;	  the IDs:
;	  WVL is set to 0 or the WAVE of the strongest line;
;	  FLUX is the total of the peak emissivities of the IDs;
;	  FLUXERR is 1% of FLUX (if non-zero), or 1 otherwise.
;	* obviously, WVL is required if FLUX is given, and WVL and FLUX
;	  are required if FLUXERR is given.
;	* of the IDs, ELEM is required, ION is useful, WAVE and SOURCE are
;	  optional (obviously, WAVE is required if SOURCE is given, and
;	  WAVE and SOURCE are required if DESCRIPTION is given.)
;	* due to a shortcoming in ZION2SYMB, it is best if "ELEM ION" is
;	  given as "ELEMENT ION", e.g., "Fe XII" or "FeXII" or "Fe 12"
;	  and not "26 12" or "26_12"
;	* <SEP> is a separator, usually a <tab>, but had better not be any
;	  of "," "+-" "-" and "?"
;	* WAVE may be (any combination of):-
;	  exact match		-- "WVL"
;	  inexact match		-- "WVL +- dW" or "WVL ? dW"
;	  range			-- "WMIN,WMAX" or "WMIN-WMAX"
;	  enclosed in brackets	-- "()" "[]" "{}"
;	* lines prefixed by PREFIX, "%", "#", ";" and "/" are ignored
;	  (except for "/{", "/}", and "//" -- the last is treated as
;	  scribbles to go into the NOTES field of the ID)
;	* if no IDs are found, WVL is flagged as "UNKNOWN"
;	* emissivities of IDs not attached to any WVL are returned in LNSTR
;
;subroutines
;	RD_LIST
;	  RD_LINE
;	    SYMB2ZION
;	      LAT2ARAB
;	  FOLD_IONEQ
;	    RD_IONEQ
;	      READ_IONEQ
;	  CAT_LN
;	LINEFLX
;	  GETABUND
;	  WHEE
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Nov98)
;	updated to include FLUXERR in IDSTR (VK; Dec98?)
;	updated to include NOTES in IDSTR (VK; Mar99)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	corrected bug where LNSTR was being concatenated if defined
;	  on input; now input is ignored (VK; Jul2008)
;	corrected bug where LNSTR was not being output (VK; Jan2015)
;-

;	usage
np=n_params()
if np eq 0 then begin
  print,'Usage: munge_list,list,idstr,lnstr,wvl,flux,fluxerr,dbdir=dbdir,$'
  print,'       sep=sep,prefix=prefix,comment=comment; RD_LIST:mapping;'
  print,'       RD_LINE:pres,logp,n_e,desig,econf; FOLD_IONEQ:chifil;'
  print,'       RD_IONEQ:chidir,eqfile; LINEFLX:DEM,abund,noph,effar,wvlar'
  print,'  return IDs and emissivities for given observed list of lines and IDs'
  return
endif

ss='	'	;default separator
if keyword_set(sep) then ss=sep

pfx='*'		;default comment prefix
if keyword_set(prefix) then pfx=string(prefix)

;	stupid user tricks
sz=size(list) & nsz=n_elements(sz)
if sz(nsz-2) ne 7 then begin
  ;	LIST is not a string, and most likely not the filename
  message,'Require name of file as input',/info & return
endif
;
if ss eq ',' or ss eq '-' or ss eq '+-' or ss eq '?' then begin
  message,'harrumph.  illegal field separator: <'+sep+'>',/info
  return
endif

;	initialize
comment=['']		;store all comments
id=0			;id=1 ==> we're handling an ID
kid=-1L			;which'th ID?
lnlist=0 & idlist=0	;ascii lists
lnstr=0

;	read each line in LIST and munge
openr,ulst,list,/get_lun
while not eof(ulst) do begin			;{read each line
  line='' & readf,ulst,line ;& stop,line
  c0=strmid(strtrim(line,2),0,1) & com=0
  if c0 eq '%' then com=1		;yes, a comment
  if c0 eq '#' then com=1		;yes, a comment
  if c0 eq ';' then com=1		;yes, a comment
  if c0 eq pfx then com=1		;yes, a comment
  if com eq 1 then begin
    comment=[comment,line]
    goto,go001				;yes, a goto
  endif

  ;	ok, not a comment
  c1=strmid(strtrim(line,2),1,1)

  if c0 ne '/' and np eq 3 then begin	;(add to the kitty
    if keyword_set(dbdir) then odbdir=dbdir
    tmp=rd_list(line,dbdir=odbdir,sep=sep,prefix=prefix, _extra=e)
    if n_tags(lnstr) eq 0 then lnstr=tmp else lnstr=cat_ln(lnstr,tmp)
  endif							;loose match)

  if c0 ne '/' and id eq 1 then begin		;(just accumulate for now
    if is_keyword_set(idlist) then idlist=[idlist,line] else idlist=[line]
  endif						;ID)

  ;	but wait a while -- we still have to deal with the "special comment"
  if c0 eq '/' then begin

    if c1 ne '{' and c1 ne '}' then begin
      comment=[comment,line]
      if c1 eq '/' then begin		;line starts with "//"
	if scratch eq '' then scratch=strmid(line,2,strlen(line)-2) else $
	  scratch=scratch+string("12b)+string("15b)+$
		strmid(line,2,strlen(line)-2)
      endif
      goto,go001			;skip the rest of the agonizing
    endif

    if id eq 1 and c1 eq '{' then begin		;(nested IDs?  ignore
      message,"IDs appear to be nested?  Reading user's mind",/info
      msg='abracadabra' & bmsg=byte(msg) & nb=n_elements(bmsg)
      for i=0,nb-1 do begin
        kilroy,dot=bmsg(i) & j=fix(randomu(seed)*5+1)
	for k=0,j do kilroy
        ;wait,randomu(seed)
      endfor
      kilroy,dot='!'
      message,'	OK, calling this a >new< ID set',/info
      lookagain=1 & goto,go002		;(go and close off the ID
      go004: lookagain=0		;and come right back here)
    endif					;ID=1 and C1="{")

    ;	but if an initial "{" is missing..
    if id eq 0 and c1 eq '}' then begin		;(huh, where'd this come from?
      message,'Missing observation? ..Ignoring',/info
    endif					;ID=0 and C1="}")

    if id eq 0 and c1 eq '{' then begin		;(read observed WVL and FLUX
      id=1
      scratch=''	;any comments to add?
      cc=str_sep(line,ss) & ncc=n_elements(cc)	;extract the fields
      if ncc ge 2 then wvl=float(cc(1)) else wvl=-1.
      if ncc ge 3 then flx=float(cc(2)) else flx=-1.
      if ncc ge 4 then ferr=float(cc(3)) else ferr=-1.
      kid=kid+1L
      if kid eq 0 then owvl=wvl else owvl=[owvl,wvl]
      if kid eq 0 then oflx=flx else oflx=[oflx,flx]
      if kid eq 0 then ofer=ferr else ofer=[ofer,ferr]
      flstr=0		;reset the line emissivity structure
    endif					;ID=0 and C1="{")

    if id eq 1 and c1 eq '}' then begin		;(end this ID
      go002:		;{this block must be touched no matter what
      id=0
      ;
      ;	here construct the KIDth component of IDSTR
      if keyword_set(idlist) then begin		;(IDs have been read in
        ;	read in the emissivities
	if keyword_set(dbdir) then odbdir=dbdir
	flstr=rd_list(idlist,/incieq,dbdir=odbdir,sep=sep,prefix=prefix,$
	  _extra=e)
        if n_tags(lnstr) eq 0 then lnstr=flstr else lnstr=cat_ln(lnstr,flstr)
	if not keyword_set(odir) then odir=[odbdir] else odir=[odir,odbdir]
	;
	;	get the label for the line
	desig=flstr.DESIG & econf=flstr.CONFIG & nw=n_elements(flstr.WVL)
	labl=strarr(2,nw)
	if n_elements(desig) eq 2*nw then labl=desig
	if n_elements(econf) eq 2*nw then labl='('+econf+') '+labl
	;
	;	get the line fluxes
	flux=lineflx(flstr.LINE_INT,flstr.logT,flstr.WVL,flstr.Z, _extra=e)
	tfx=total(flux) & if tfx le 0 then flux(*)=1.
	;
	;	reset the observed wavelength and/or flux
	if owvl(kid) lt 0 then owvl(kid)=total(flux*flstr.WVL)/total(flux)
	if oflx(kid) le 0 then oflx(kid)=total(flux) else $
		flux=flux*oflx(kid)/total(flux)
	flxerr=0*flux
	if ofer(kid) lt 0 then begin
	  ;default error is 1% of flux
	  if oflx(kid) ne 0 then ofer(kid)=0.01*abs(oflx(kid)) else ofer(kid)=0.
	endif else begin
	  ;flxerr=flux*ofer(kid)/total(flux)
	  flxerr=0*flux+sqrt(ofer(kid)^2/n_elements(flux))
	endelse
					;force theoretical to match observed
	;	make LINEID-style structure
	logT=flstr.logT
	idkid=create_struct('WVL',flstr.wvl,'Z',flstr.Z,'ION',flstr.ION,$
	  'LABL',labl,'FLUX',flux,'FLUXERR',flxerr,$
	  'LOGT',flstr.logT,'EMIS',flstr.LINE_INT)
	if scratch ne '' then idkid=create_struct(idkid,'NOTES',scratch)
	idlist=0
      endif else begin				;)(unknown IDs
	labl=strarr(2,1) & labl(0)='Unknown'
	if not keyword_set(logT) then logT=findgen(81)*0.05+4.
	emis=0*logT
	idkid=create_struct('WVL',[-abs(wvl)],'Z',[0],'ION',[0],$
	  'LABL',labl,'FLUX',[0.],'FLUXERR',[0.],'LOGT',logT,'EMIS',emis)
	if scratch ne '' then idkid=create_struct(idkid,'NOTES',scratch)
      endelse					;FLSTR)
      if kid eq 0 then idwvl=create_struct('ID0',idkid) else $
        idwvl=create_struct(idwvl,'ID'+strtrim(kid,2),idkid)

      if keyword_set(cleanup) then goto,go003		;get out..
      if keyword_set(lookagain) then goto,go004		;.. NOW!}
    endif					;ID=1 and C1="}")
  endif

  go001: ;get here directly if LINE was commented out
endwhile					;!EOF(ULST)}
close,ulst,/all & free_lun,ulst

;	clean up
if id eq 1 then begin	;(file ended before closing off the ID
  cleanup=1 & goto,go002	;(just repeat the closing block above
  go003: cleanup=0		;and come right back here)
endif			;if ID is still set)

;	the ID structure
if n_tags(idwvl) gt 0 then idstr=create_struct('WVL',owvl,idwvl) else begin
  message,'no IDs read in',/info & return
endelse
if np eq 1 then help,cat_id(idstr)

;	the line database directories
if keyword_set(odir) then dbdir=odir

;	sundry output
wvl=owvl
flux=oflx
fluxerr=ofer

return
end
