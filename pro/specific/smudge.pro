pro smudge,loc,label,keys=keys,only=only,nisi=nisi,range=range,$
	infile=infile,sep=sep,prefix=prefix,okeV=okeV,help=help,$
	keV=keV,ikeV=ikeV, _extra=e
;+
;procedure	smudge
;	returns labels for features at given locations
;
;	the input file is expected to be in the following format:-
;		LOCATION <sep> UNIT <sep> LABEL <sep> KEYWORDS
;	where
;	  <sep> is a separator, usually a <tab>
;	  LOCATION is a number, usually the wavelength in [Ang]
;	  UNIT specifies the units for LOCATION, e.g., "Ang", "keV", etc.
;	  LABEL is a string label describing the feature at LOCATION
;	  KEYWORDS are a set of keywords that facilitate pinpoint extraction,
;	    usually the instrument name, the origin of the feature, etc.
;	Note that:
;	  LOCATION is required, and everything else is optional
;	  Lines prefixed by "%", "#", ";", "/" are ignored
;	  the following UNITs are understood: Ang, eV, keV, MeV, Hz, MHz, GHz,
;	    wave-#, /cm, cm^-1, cm**-1, mu, micron, mm, cm, meter, km
;
;syntax
;	smudge,x,label,keys=keys,only=only,nisi=nisi,range=range,$
;	infile=infile,sep=sep,prefix=prefix,/okeV,/help
;
;parameters
;	loc	[OUTPUT] float array of location of features
;	label	[OUTPUT] string array of labels for features
;
;keywords
;	keys	[INPUT] string array of keywords.  if given, all features
;		answering to these keywords are returned
;	only	[INPUT] string array of keywords.  while KEYS does a logical
;		OR, ONLY does a logical AND on the keywords
;	nisi	[INPUT] return all >except< those features matching the
;		specified keywords
;		* KEYS is checked first, ONLY next, NISI last, and the
;		  precedence is therefore in reverse order
;	range	[INPUT] if set, returns only those LABEL(LOC)s that lie
;		within RANGE
;		* if 1-element vector or scalar, looks for exact match
;		* WARNING: range check is carried out >after< conversion
;		  to keV if OKEV (see below) is set, so RANGE must then be
;		  given in units of [keV] (otherwise it should be in [Ang])
;	infile	[INPUT] ascii file in which the features are listed
;		* default: $SCARDIR/ardb/smudge.lst
;	sep	[INPUT] separator used to delineate fields in INFILE
;		* default is <tab>
;	prefix	[INPUT] any line beginning with this is considered a comment
;		* default: "*"
;		* this is in addition to "%", "#", ";", "/"
;		* beware: "<sp>" or "<tab>" won't work!
;	okeV	[INPUT] if set, LOC are returned in units of keV
;		* default is to return LOC in units of Ang
;		* KEV is accepted as a synonym for historical reasons
;	help	[INPUT] if set, prints out usage and quits
;	ikeV	[INPUT; DEFUNCT] a synonym for OKEV, still kept
;		(temporarily) only for backward compatibility
;
;	_extra	[JUNK] here only to prevent crashing the program
;
;example usage
;	smudge,/help
;	smudge,wvl,lbl & xr=[min(wvl),max(wvl)] & plot,[0],xr=xr,/xlog
;		xyouts,wvl,randomu(seed,n_elements(wvl)),lbl,orient=90
;	smudge,keys=['edge','gap'],only=['M']
;	smudge,keys='gap',only=['S0','S1'],/okeV
;	smudge,keys='edge gap',range=[1,6]
;	smudge,keys='edge gap',nisi='edge',range=[1,6]
;	smudge,range=0.69,/okev
;
;history
;	vinay kashyap (Nov98)
;	added warning in case of bad SEP (VK; Dec98)
;	changed "caldb" to "ardb" (VK; DecMM)
;	added keyword IKEV in preparation to getting rid of KEV (VK; JanMMI)
;	changed keyword IKEV to OKEV (VK; Jul02)
;-

;	initialize
loc=[-1.] & label=['']
h=6.626176d-27 & c=2.99792458d10

;	usage
if keyword_set(help) then begin
  print,'Usage: smudge,x,label,keys=keys,only=only,nisi=nisi,range=range,$'
  print,'       infile=infile,sep=sep,prefix=prefix,/okeV,/help'
  print,'  extract labels for features from file'
  return
endif

;	check input file
szf=size(infile) & nszf=n_elements(szf)
if szf(nszf-2) eq 7 then begin		;(INFILE is a string
  filnam=strtrim(infile(0),2)
endif else begin
  scardir=getenv('SCAR')
  if scardir(0) eq '' then scardir=getenv('SCARDIR')
  if scardir(0) eq '' then scardir='/data/fubar/SCAR'
  filnam=scardir+'/ardb/smudge.lst'
endelse					;INFILE)

;	check keywords
;	separator between fields in input ascii file
ss='	' & if keyword_set(sep) then ss=sep
;	anything prefixed by this is a comment
prfx='*' & if keyword_set(prefix) then prfx=strtrim(prefix,2)
;	keywords to combine with OR
if keyword_set(keys) then begin		;(OR keys
  keyOR=[''] & nk=n_elements(keys)
  for i=0,nk-1 do begin
    c1=str_sep(keys(i),',') & nc1=n_elements(c1)	;break apart at ","
    for j=0,nc1-1 do begin
      c2=str_sep(c1(j),' ') & nc2=n_elements(c2)	;break at "<sp>"
      for k=0,nc2-1 do begin
	c3=str_sep(c2(k),'	') & nc3=n_elements(c3)	;break at "<tab>"
	for l=0,nc3-1 do begin
	  c4=strtrim(c3(l),2) & if c4 ne '' then keyOR=[keyOR,c4]
	endfor
      endfor
    endfor
  endfor
  if n_elements(keyOR) eq 1 then keyOR=0 else keyOR=keyOR(1:*)
endif else keyOR=0			;ORs)
;	keywords to combine with AND
if keyword_set(only) then begin		;(AND keys
  keyAND=[''] & nk=n_elements(only)
  for i=0,nk-1 do begin
    c1=str_sep(only(i),',') & nc1=n_elements(c1)	;break apart at ","
    for j=0,nc1-1 do begin
      c2=str_sep(c1(j),' ') & nc2=n_elements(c2)	;break at "<sp>"
      for k=0,nc2-1 do begin
	c3=str_sep(c2(k),'	') & nc3=n_elements(c3)	;break at "<tab>"
	for l=0,nc3-1 do begin
	  c4=strtrim(c3(l),2) & if c4 ne '' then keyAND=[keyAND,c4]
	endfor
      endfor
    endfor
  endfor
  if n_elements(keyAND) eq 1 then keyAND=0 else keyAND=keyAND(1:*)
endif else keyAND=0			;ANDs)
;	keywords to explicitly exclude
if keyword_set(nisi) then begin		;(NOT keys
  keyNOT=[''] & nk=n_elements(nisi)
  for i=0,nk-1 do begin
    c1=str_sep(nisi(i),',') & nc1=n_elements(c1)	;break apart at ","
    for j=0,nc1-1 do begin
      c2=str_sep(c1(j),' ') & nc2=n_elements(c2)	;break at "<sp>"
      for k=0,nc2-1 do begin
	c3=str_sep(c2(k),'	') & nc3=n_elements(c3)	;break at "<tab>"
	for l=0,nc3-1 do begin
	  c4=strtrim(c3(l),2) & if c4 ne '' then keyNOT=[keyNOT,c4]
	endfor
      endfor
    endfor
  endfor
  if n_elements(keyNOT) eq 1 then keyNOT=0 else keyNOT=keyNOT(1:*)
endif else keyNOT=0			;NOTs)
;	check for range
rng=[0.,1e7]		;the default
if keyword_set(range) then begin
  nr=n_elements(range)
  if nr eq 1 then rng=range(0)+[-1e-4,1e-4]
  if nr ge 2 then rng=[min(range),max(range)]
endif

openr,uin,filnam,/get_lun			;open INFILE for input
while not eof(uin) do begin			;{step through INFILE

  line='' & readf,uin,line & line=strtrim(line,2)	;read line
  c=strmid(line,0,1) & comment=0
  if line eq '' then comment=1
  if c eq '%' or c eq '#' or c eq ';' or c eq '/' then comment=1
  if strpos(line,prfx,0) eq 0 then comment=1
  if comment eq 1 then goto,go001	;{skip rest of agonizing

  flds=str_sep(line,ss) & nf=n_elements(flds)		;break into fields
  pos=float(flds(0))		;LOCATION

  if nf eq 1 then begin				;(no labels??
    message,'Missing UNIT and LABEL for: ',/info
    print,flds(0)
    if n_elements(str_sep(line,'	')) gt 1 or $
     n_elements(str_sep(line,' ')) gt 1 or $
      n_elements(str_sep(line,',')) gt 1 then message,$
       '	Possible that the separator is incorrectly set',/info
    c1='type any key to continue, x/z to halt, q to quit'
    print,c1 & c1=get_kbrd(1)
    if strlowcase(c1) eq 'x' or strlowcase(c1) eq 'z' then stop,$
	'HALTING!  type .CON to continue'
    if strlowcase(c1) eq 'q' then return
    lbl=''
  endif						;no labels)

  if nf ge 2 then begin				;(get the units right!
    units=flds(1) & oku=-1 & lbl=0

    ;	find out what the units are..
    unit='Angstrom'			;the default
    iun=strpos(strlowcase(units),'ev',0)
    if iun ge 0 then begin
      unit='eV'
      iun=strpos(strlowcase(units),'kev',0) & if iun ge 0 then unit='keV'
      iun=strpos(strlowcase(units),'mev',0) & if iun ge 0 then unit='MeV'
    endif
    iun=strpos(strlowcase(units),'hz',0)
    if iun ge 0 then begin
      unit='Hz'
      iun=strpos(strlowcase(units),'mhz',0) & if iun ge 0 then unit='MHz'
      iun=strpos(strlowcase(units),'ghz',0) & if iun ge 0 then unit='GHz'
    endif
    iun=strpos(strlowcase(units),'mu',0) & if iun ge 0 then unit='micron'
    iun=strpos(strlowcase(units),'micron',0) & if iun ge 0 then unit='micron'
    iun=strpos(strlowcase(units),'cm',0)
    if iun ge 0 then begin
      unit='cm'
      iun=strpos(strlowcase(units),'/',0) & if iun ge 0 then unit='wave'
      iun=strpos(strlowcase(units),'^',0) & if iun ge 0 then unit='wave'
      iun=strpos(strlowcase(units),'-1',0) & if iun ge 0 then unit='wave'
    endif
    iun=strpos(strlowcase(units),'wave',0) & if iun ge 0 then unit='wave'
    iun=strpos(strlowcase(units),'meter',0) & if iun ge 0 then unit='meter'
    if strtrim(units,2) eq 'mm' then unit='mm
    iun=strpos(strlowcase(units),'km',0) & if iun ge 0 then unit='meter'

    ;	.. and translate said units to [Ang]
    case unit of			;{translate location to [Ang]
      'eV': if pos ne 0 then pos=12.3985/(pos/1e3)	;[eV] -> [Ang]
      'keV': if pos ne 0 then pos=12.3985/pos		;[keV] -> [Ang]
      'MeV': if pos ne 0 then pos=12.3985/(pos*1e3)	;[MeV] -> [Ang]
      'Hz': if pos ne 0 then pos=(1e8*c)/pos		;[Hz] -> [Ang]
      'MHz': if pos ne 0 then pos=(1e8*c)/(pos*1e6)	;[MHz] -> [Ang]
      'GHz': if pos ne 0 then pos=(1e8*c)/(pos*1e9)	;[GHz] -> [Ang]
      'wave': if pos ne 0 then pos=1e8/pos		;[wave#] -> [Ang]
      'micron': pos=1e4*pos				;[micron] -> [Ang]
      'mm': pos=1e9*pos					;[mm] -> [Ang]
      'cm': pos=1e8*pos					;[cm] -> [Ang]
      'meter': pos=1e10*pos				;[m] -> [Ang]
      'km': pos=1e13*pos				;[km] -> [Ang]
      else: begin					;[Ang]
	iun=strpos(strlowcase(units),'ang',0)
	if iun lt 0 then begin
	  message,'['+units+'] not understood; assuming [Ang]',/info
	  message,'	acceptable units are:',/info
	  message,'	eV, keV, MeV, Hz, MHz, GHz, wave#, micron, mm, cm, meter, km',/info
	  print,string("7b)
	endif
      end
    endcase				;UNIT}

  endif					;UNITS)

  if nf ge 3 then begin				;(LABEL
    if not keyword_set(lbl) then lbl=flds(2) else lbl=lbl+ss+flds(2)
    kw=flds(2)		;if no keywords field is present, use this itself
  endif						;LABEL)

  if nf ge 4 then begin				;(KEYWORDS
    for i=3,nf-1 do kw=kw+flds(i)
  endif						;KEYWORDS)

  if keyword_set(keyOR) then begin		;(OR search
    ik=n_elements(keyOR) & kk=strlowcase(keyOR) & incthis=0
    for i=0,ik-1 do begin
      j=strpos(strlowcase(kw),kk(i),0) & if j ge 0 then incthis=1
    endfor
  endif else incthis=1				;KEYS)

  if keyword_set(keyAND) then begin		;(AND search
    ik=n_elements(keyAND) & kk=strlowcase(keyAND)
    for i=0,ik-1 do begin
      j=strpos(strlowcase(kw),kk(i),0) & if j lt 0 then incthis=0
    endfor
  endif						;ONLY)

  if keyword_set(keyNOT) then begin		;(NOT search
    ik=n_elements(keyNOT) & kk=strlowcase(keyNOT)
    for i=0,ik-1 do begin
      j=strpos(strlowcase(kw),kk(i),0) & if j ge 0 then incthis=0
    endfor
  endif						;NISI)

  ;	place in output
  if incthis ne 0 then begin
    loc=[loc,pos] & label=[label,lbl]
  endif

  go001:				;get here direct if commented out}
endwhile					;!EOF(UIN)}
close,uin & free_lun,uin			;close INFILE

;	reformat output
n=n_elements(loc) & if n eq 1 then return
ii=lindgen(n-1L)+1L & loc=loc(ii) & label=label(ii)
if keyword_set(ikeV) then okeV=ikeV
if keyword_set(okeV) or keyword_set(keV) then begin
  oo=where(loc ne 0,moo)
  if moo gt 0 then loc(oo)=12.3985/loc(oo)	;[Ang] -> [keV]
endif

;	check range
if keyword_set(range) then begin
  oo=where(loc ge rng(0) and loc le rng(1),moo)
  if moo eq 0 then return
  loc=loc(oo) & label=label(oo) & n=moo+1L
endif

np=n_params()
if keyword_set(okeV) then unit='keV' else unit='Ang'
if np lt 2 then begin
  for i=0L,n-2L do print,'@'+string(loc(i),'(g11.4)')+' '+unit+': '+label(i)
endif

return
end
