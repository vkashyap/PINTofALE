function rd_list,lnlist,incieq=incieq,mapping=mapping,dbdir=dbdir,sep=sep,$
	prefix=prefix,comment=comment,brklst=brklst,eps=eps,verbose=verbose,$
	_extra=e
;+
;function	rd_list
;	read in line emissivities for a list of specified lines and returns
;	the database in a standard structure format (RD_LINE for details)
;
;syntax
;	lstr=rd_list(lnlist,/incieq,desig=desig,econf=econf,mapping=mapping,$
;	dbdir=dbdir,brklst=brklst,sep=sep,prefix=prefix,comment=comment,$
;	pres=pres,logP=logP,n_e=n_e,chifil=chifil,chidir=chidir,eqfile=eqfile,$
;	eps=eps,verbose=verbose)
;
;parameters
;	lnlist	[INPUT; required] line list
;		* must be a string scalar or array
;		* if a line contains a "/" or "./" or ".." in the first field,
;		  assumed to be a filename unless the "/" is followed by
;		  another "/" ("//") or a "*" ("/*"), in which case it is
;		  assumed to be a comment.
;		* any line prefixed by "%", "#", ";", "/*", "*/", "//",
;		  and PREFIX are assumed to be comments
;		* either the file or the direct input must describe the
;		  line as:-
;		 	"Z ION <sep> WAVE <sep> SOURCE <sep> DESCRIPTION"
;		  where Z is the atomic symbol, ION is the ionic state,
;		  WAVE is in [Ang], and SOURCE is DBDIR (see RD_LINE)
;		* Z is required
;		  ION is useful
;		  WAVE, SOURCE, and DESCRIPTION are optional
;		  but WAVE is required if SOURCE is given
;		  and WAVE and SOURCE are required if DESCRIPTION is given
;		* if multiple candidates are available, DESCRIPTION may be
;		  used to choose from among them (see MAPPING)
;		* <SEP> is a separator, usually a <tab>, but had better
;		  not be any of "," "+-" "-" and "?"
;		* WAVE may be (any combination of):-
;		  exact match		-- "WVL"
;		  inexact match		-- "WVL +- dW" or "WVL ? dW"
;		  range			-- "WMIN,WMAX" or "WMIN-WMAX"
;		  enclosed in brackets	-- "()" "[]" "{}"
;		* DESCRIPTION is basically treated as a set of <SP>-separated
;		  keywords, and if a given field is preceded by a "!" any
;		  matches based on that field would be _excluded_ (e.g.,
;		  "!DR" would not match any DR line)
;
;keywords
;	incieq	[INPUT] if set, output emissivities will include ion-balance
;	mapping	[INPUT] this keyword controls what happens when an entry in
;		LNLIST matches multiple instances of (Z,ION,WAVE) in the
;		input database.
;		* if not set, or set to 0, returns all the emissivities,
;		  exactly as found, in separate rows,
;		* but if not set and yet DESCRIPTION is set, then assumed
;		  to be +ve
;		* if set, and is +ve, then
;		  -- looks through the DESCRIPTION field of LNLIST, and
;		     finds the closest match in the database and returns
;		     only that.
;		  -- if "closest match" is ambiguous, then default is set
;		     to the strongest of the closest matches, and the user
;		     is asked to verify.
;		  -- if DESCRIPTION is not set, or the given DESCRIPTION
;		     succeeds in eliminating >all< the matches, then asks
;		     which one (the default being the strongest line in
;		     the set).
;		* if set, and is -ve, then adds the emissivities of all
;		  the extra lines to that of the strongest line in the
;		  set.
;		* NOTE: if "WVL?dW" is set, and multiple lines are found
;		  in the database, MAPPING is overridden and the user is
;		  asked to choose anyway.
;	dbdir	[I/O] directory in which to look for files
;		* will be used ONLY if not given in LNLIST
;		* if not specified, *default* will be set at run-time to
;		  whatever was used last
;		* on output, will be a string array of the DBDIRs used
;		  NOTE: this last feature will soon be eliminated -- use
;		  keyword BRKLST instead
;	sep	[INPUT] separator used to delineate fields in LNLIST
;		* the default is to separate the fields with several blanks
;		  (<tab>)
;	prefix	[INPUT; default='*'] any line beginning with this letter
;		is considered a comment
;		* this is in addition to "%", "#", ";", "/"
;	comment	[OUTPUT] any commented lines are returned in this variable
;		as a string array
;	brklst	[OUTPUT] structure containing a breakdown of all the
;		fields in LNLIST, in the form
;			{SPECIES, WMIN, WMAX, DBDIR, DESCRIPTION}
;	eps	[INPUT; default=1e-5] a small number
;		* used as a "round-off" error for wavelength precision
;		  if the WAVE field does not specify a range
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[INPUT ONLY] use this to pass defined keywords to
;		-- RD_LINE:
;		   N_E		electron density [cm^-3], default is 1e9
;				_explicitly set to 0_ to use PRES/LOGP
;		   PRES		pressure [cm^-3 K]
;		   LOGP		alog10(pressure [cm^-3 K])
;		   DESIG	to read in level designations, set by default
;		   ECONF	to read in e configurations, set by default
;		-- FOLD_IONEQ:
;		   CHIFIL	to read in from CHIANTI style database
;		   EQFILE	name of ion balance file
;		-- RD_IONEQ:
;		   CHIDIR	path to CHIANTI installation
;
;restrictions
;	syntax incompatible for IDL <5
;	function calls incompatible for idl <5.3
;	requires subroutines:
;	  RD_LINE, SYMB2ZION, LAT2ARAB, FOLD_IONEQ, RD_IONEQ, READ_IONEQ,
;	  CAT_LN, LINEFLX, GETABUND, RDABUND, SYZE, INICON, IS_KEYWORD_SET
;
;history
;	vinay kashyap (Jun98)
;	added keyword COMMENT; changed behavior of SETASK (VK; Oct98)
;	DBDIR now contains all directories on output (VK; Nov98)
;	changed input to FOLD_IONEQ from ION to JON (VK; 99May)
;	allowed unfound elements to be read in also, up to a point (VK; 99Jul)
;	added keyword MAPPING, enhanced LNLIST, deleted SETASK, converted
;	  to IDL5.3 (VK; MarMM)
;	now quits if IDL v <5.3 (VK; DecMM)
;	force MAPPING=1 if DESCRIPTION is set; added keyword BRKLST as
;	  first step towards eliminating DBDIR as output; added keyword
;	  VERBOSE and EPS; changed DESIG/ECONF behavior to be set by default
;	  (VK; JanMMI)
;	changed default behavior to read in from constant-density database
;	  at N_E=1e9 if nothing is specified (VK; Jan'03)
;	added ability to exclude fields via DESCRIPTION (VK; Apr'03)
;	changed default SEP to "any series of one or more space or tab
;	  characters" (VK; Feb'04)
;	changed SEP back to <tab> (VK; Sep'05)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;-

;	usage
nl=n_elements(lnlist)
if nl eq 0 then begin
  print,'Usage: lstr=rd_list(lnlist,/incieq,desig=desig,econf=econf,$'
  print,'       mapping=mapping,dbdir=dbdir,brklst=brklst,sep=sep,$'
  print,'       prefix=prefix,comment=comment,logP=logP,n_e=n_e,$'
  print,'       chifil=chifil,chidir=chidir,eqfile=eqfile,eps=eps,$'
  print,'       verbose=verbose)'
  print,'  return line emissivities for given line list'
  print,'  LNLIST is of the form'
  print,'  	"Z ION <sep> WAVE <sep> SOURCE <sep> DESCRIPTION"'
  return,-1L
endif

;	quit if IDL v < 5.3
if float(strmid(!version.RELEASE,0,3)) lt 5.3 then begin
  message,'Requires IDL v5.3 or higher.  Returning',/informational
  return,-1L
endif

;ss='[ ' + STRING(9B) + ']+'	;default separator (one or more <sp> or <tab>)
ss='	'	;default separator
if keyword_set(sep) then ss=sep

pfx='*'		;default comment prefix
if keyword_set(prefix) then pfx=string(prefix)

mm=0		;default mapping of multiple matches
nm=n_elements(mapping) & szm=size(mapping) & nszm=n_elements(szm)
if nm gt 0 then begin
  mm=1		;MAPPING has been set
  if szm(nszm-2) le 5 then mm=mapping[0]	;as given
endif

v=0		;default verbosity
if keyword_set(verbose) then v=long(verbose[0]) > 1

epsilon=1e-5	;default precision
if keyword_set(eps) then epsilon=eps[0] > 0

;	stupid user tricks
sz=size(lnlist) & nsz=n_elements(sz)
if sz(nsz-2) ne 7 then begin		;not a string
  c1='Input must be in the form "ATOM ION <tab> WVL <tab> DBDIR"'
  message,c1,/informational
  return,-1L
endif
;
if ss eq ',' or ss eq '-' or ss eq '+-' or ss eq '?' then begin
  message,'harrumph.  illegal field separator: <'+sep+'>',/informational
  return,-1L
endif

;	build up the compleat line list
llst=0
for i=0L,nl-1L do begin
  cc=strsplit(lnlist[i],ss,/extract)
  c0=strmid(strtrim(cc[0],2),0,1) & c01=strmid(strtrim(cc[0],2),0,2)
  if c01 eq './' or c01 eq '..' then c0='/'			;path to file
  if c01 eq '/*' or c01 eq '//' or c01 eq '*/' then c0=''	;not a file
  if c0 eq '/' then begin			;(lnlist[i] is a file
    openr,ulst,lnlist[i],/get_lun	;{read from file
    while not eof(ulst) do begin
      cc='' & readf,ulst,cc & cc=strtrim(cc,2)
      if not is_keyword_set(llst) then llst=cc else llst=[llst,cc]
    endwhile
    close,ulst & free_lun,ulst		;file I/O}
  endif	else begin				;c0='/')
    if not is_keyword_set(llst) then llst=lnlist[i] else llst=[llst,lnlist[i]]
  endelse
endfor
;
;	get rid of commented lines
nl=n_elements(llst) & icom=intarr(nl)
for i=0L,nl-1L do begin
  c0=strmid(strtrim(llst[i],2),0,1) & c01=strmid(strtrim(cc[0],2),0,2)
  if c0 eq '%' or c0 eq '#' or c0 eq ';' or c0 eq '/' then icom[i]=1
  if c01 eq '//' or c01 eq '/*' or c01 eq '*/' then icom[i]=1
  if icom[i] eq 0 and keyword_set(pfx) then begin
    npfx=strlen(pfx) & icom[i]=strcmp(llst[i],pfx,npfx)
    ;icom[i]=strmatch(llst[i],pfx+'*')
  endif
endfor
oc=where(icom gt 0,moc) & if moc gt 0 then comment=llst[oc]
ok=where(icom eq 0,mok) & if mok gt 0 then llst=llst[ok] & nl=mok
if nl eq 0 then return,-1L

;	and now go and get em
for i=0L,nl-1L do begin			;{for each listed line
  if v gt 2 then message,'Working on: '+llst[i],/informational
  fstr=1 & imap=mm
  if not keyword_set(dbdir) then dbdir=0
  ;
  ;	parse LNLIST
  cc=strsplit(llst[i],ss,/extract) & ncc=n_elements(cc)
  ip=4 & if ncc lt ip then desc='' else desc=strtrim(cc[ip-1],2)
  ip=3 & if ncc lt ip then ldb=dbdir[0] else ldb=strtrim(cc[ip-1],2)
  ip=2 & if ncc lt ip then ww='1-1000' else ww=strtrim(cc[ip-1],2)
  ip=1 & if ncc lt ip then atom=0 else atom=strtrim(cc[ip-1],2)
  if atom eq '*' or atom eq '0' then atom='all'	;short-hand for "all atoms"
  ;
  ;	parse WW
  c=strtrim(ww,2)
  j=strpos(c,'[',0) & if j ge 0 then c=strmid(c,j+1,strlen(c)-j-1)
  j=strpos(c,']',0) & if j ge 0 then c=strmid(c,0,j)
  j=strpos(c,'(',0) & if j ge 0 then c=strmid(c,j+1,strlen(c)-j-1)
  j=strpos(c,')',0) & if j ge 0 then c=strmid(c,0,j)
  j=strpos(c,'{',0) & if j ge 0 then c=strmid(c,j+1,strlen(c)-j-1)
  j=strpos(c,'}',0) & if j ge 0 then c=strmid(c,0,j)
  j=0							; W (exact match)
  if strpos(c,'+-',0) ge 0 then j=1			; W +- dW
  if strpos(c,'-',0) ge 0 and j ne 1 then j=2		; W1 - W2
  if strpos(c,',',0) ge 0 then j=3			; W1,W2
  if strpos(c,'?',0) ge 0 then j=4			; W ? dW (verify)
  case j of					;{for each case
    1: begin
      cc=strsplit(c,'+-',/extract) & ncc=n_elements(cc)
      w=float(cc[0]) & if ncc eq 1 then dw=1e-3*w else dw=float(cc[1])
      w1=w-dw & w2=w+dw
    end
    2: begin
      cc=strsplit(strtrim(c,2),'-',/extract) & ncc=n_elements(cc)
      w1=float(cc[0]) & if ncc eq 1 then w2=w1+1. else w2=float(cc[1])
    end
    3: begin
      cc=strsplit(c,',',/extract) & ncc=n_elements(cc)
      w1=float(cc[0]) & if ncc eq 1 then w2=w1+1. else w2=float(cc[1])
    end
    4: begin
      cc=strsplit(c,'?',/extract) & ncc=n_elements(cc)
      w=float(cc[0]) & if ncc eq 1 then dw=1e-3*w else dw=float(cc[1])
      w1=w-dw & w2=w+dw
      imap=1 & desc=''	;force verification and ignore description
    end
    else: begin
      w1=float(c)-epsilon & w2=w1+2*epsilon
    end
  endcase					;J}
  ;
  ;	output
  if n_elements(bATOM) eq 0 then begin
    bATOM=atom & bWMIN=w1 & bWMAX=w2 & bDBDIR=ldb & bDESC=desc
  endif else begin
    bATOM=[bATOM,atom] & bWMIN=[bWMIN,w1] & bWMAX=[bWMAX,w2]
    bDBDIR=[bDBDIR,ldb] & bDESC=[bDESC,desc]
  endelse
  ;
  ;	call RD_LINE
  ff=rd_line(atom,wrange=[w1,w2],fstr=fstr,dbdir=ldb,verbose=v,$
	/desig,/econf,n_e=1e9, _extra=e)
  dbdir=ldb
  if not keyword_set(outdir) then outdir=[ldb] else outdir=[outdir,ldb]
  logT=fstr.logT & wvl=fstr.logT & Z=fstr.Z & ion=fstr.ION & jon=fstr.JON
  desig=fstr.DESIG & config=fstr.CONFIG & src=fstr.SRC
  ;	call FOLD_IONEQ
  if keyword_set(incieq) then begin
    ff=fold_ioneq(ff,Z,JON,logt=LOGT,verbose=v,_extra=e) & fstr.LINE_INT=ff
  endif
  ;
  ;	if needed, choose among the selections
  nw=n_elements(fstr.WVL)
  if nw eq 1 and v ge 10 then begin
    message,'Replacing line : '+llst[i],/informational
    junk=cat_ln(fstr,/comm)
  endif
  if nw gt 1 then begin

    ;	if a description is available, and MAPPING is not set,
    ;	that would be a mistake.
    if keyword_set(desc) and not keyword_set(imap) then imap=1

    ;	first find the strongest line
    if keyword_set(imap) then begin
      eemx=-1. & iemx=0L
      for j=0L,nw-1L do begin
        tmp=reform((fstr.LINE_INT)[*,j])
        if max(tmp) gt eemx then begin
	  eemx=max(tmp) & iemx=j
        endif
      endfor
    endif

    ;	add extra emissivities to strongest line
    if imap lt 0 then begin	;(IMAP<0
      ee=0.*reform((fstr.LINE_INT)[*,iemx])
      for j=0L,nw-1L do ee=ee+reform((fstr.LINE_INT)[*,j])
      fstr=cat_ln(fstr,pick=[iemx]) & fstr.LINE_INT=ee
    endif			;IMAP<0)

    ;	find the best match
    if imap gt 0 then begin	;(IMAP>0
      jm=-1L		;default is to ask
      if keyword_set(desc) then begin	;(anything to match to?
        ;	first make a string out of the ECONF and DESIG
        dbdesc=strarr(nw)
        cd=fstr.DESIG & ce=fstr.CONFIG & ncd=n_elements(cd) & nce=n_elements(ce)
        for j=0L,nw-1L do begin
	  if ncd gt 2 then c1=cd[0,j]+' '+cd[1,j] else c1=''
	  if nce gt 2 then c2=ce[0,j]+' '+ce[1,j] else c2=''
	  dbdesc[j]=strtrim(c1+' '+c2,2)
        endfor
        ;	now find the closest match
        lndesc=strsplit(desc,/extract) & nld=n_elements(lndesc)
	km=lonarr(nw)
	for j=0L,nld-1L do begin	;{for each keyword descriptor
	  c2=strtrim(lndesc[j],2) & i2=1L
	  if strpos(c2,'!',0) eq 0 then begin	;(this field must be EXCLUDED
	    c2=strmid(c2,1,strlen(c2)-1) & i2=0L
	  endif					;!C2)
	  for k=0L,nw-1L do begin	;{match with each line
	    c1=strtrim(dbdesc[k],2)
	    if i2 eq 1 then begin	;(regular matching
	      if keyword_set(c1) and keyword_set(c2) then $
	    	jj=strmatch(c1,'*'+c2+'*',/fold_case) else jj=0
	    endif else begin		;I2=1)(exclusionary matching
	      if keyword_set(c1) and keyword_set(c2) then $
	    	jj=1-strmatch(c1,'*'+c2+'*',/fold_case)
	    endelse			;I2=0)
	    if jj gt 0 then km[k]=km[k]+1
	  endfor			;K=0,NW-1}
	endfor				;J=0,NLD-1}
	ok=where(km gt 0,mok)
	if mok gt 0 then begin	;(yes, virginia, there is a "best" match
	  jj=max(km,jm)
	  oo=where(ok eq jm,moo)
	  if moo eq 1 then fstr=cat_ln(fstr,pick=[jm],/comm) else begin
	    message,'a surfeit of matches',/informational
	    ffstr=cat_ln(fstr,pick=jm)
	    ;	reset the default selection (in next section) to
	    ;	strongest of the matched lines
	    eemx=0.
	    for j=0L,moo-1L do begin
	      tmp=reform((ffstr.LINE_INT)[*,j])
	      if max(tmp) gt eemx then begin
		eemx=max(tmp) & iemx=jm[j]
	      endif
	    endfor
	    ;	back to square one
	    jm=-1		;too much of a good thing..
	  endelse
	endif			;MOK>0)
      endif				;DESC)
      ;
      if jm lt 0 then begin	;(no "best" match exists, so ask user
	print,'*************************************************'
	print,llst[i],' :-'
	print,'*************************************************'
	help,cat_ln(fstr,/comm)
	print,'type <sp>-separated indices of lines to include in output:'
	print,'<CR> to choose default, - to keep all'
	c='' & read,prompt='['+strtrim(iemx,2)+']> ',c & c=strtrim(c,2)
	if c eq '' then jm=[iemx] else $
	 if strmid(c,0,1) eq '-' then jm=lindgen(nw) else $
	  jm=fix(strsplit(c,/extract))
	fstr=cat_ln(fstr,pick=jm)
      endif			;JM<0)
    endif			;IMAP>0)

  endif
  ;
  ;	accumulate
  if (fstr.line_int)[0] gt -1 then begin
    if n_tags(linstr) eq 0 then linstr=fstr else $
	linstr=cat_ln(linstr,fstr, _extra=e)
  endif else begin
    ;	special handling of unknown, unfound, or unsound elements
    ;	no help if first elements are the unks, but afterwards just
    ;	interpolate onto existing grid
    if n_tags(linstr) gt 0 then begin
      linstr=cat_ln(linstr,fstr, _extra=e)
    endif else begin
      ;	what to do, what to do?
    endelse
  endelse
endfor					;I=0,NL-1}

;	outputs
if arg_present(outdir) then dbdir=outdir	;list of directories
if n_tags(linstr) eq 0 then linstr=fstr		;line structure
if arg_present(brklst) then brklst=create_struct('SPECIES',bATOM,$
	'WMIN',bWMIN,'WMAX',bWMAX,'DBDIR',bDBDIR,'DESCRIPTION',bDESC)

return,linstr
end
