pro zion2symb,z,ion,symbol,ziform=ziform,oform=oform, _extra=e
;+
;procedure	zion2symb
;		convert atomic number and ionic state to spectroscopic symbol
;
;syntax
;	zion2symb,z,ion,symbol,ziform=ziform
;
;parameters
;	z	[INPUT; required] atomic number
;	ion	[INPUT] ionic state (0=undefined, 1=neutral, 2=singly
;		ionized, etc.)
;	symbol	[OUTPUT] string containing element designation
;		(e.g., "he_2", "He II", "HeII", "he 2", "he", "He",
;		depending on ZIFORM)
;
;keywords
;	ziform	[INPUT] string containing "Z" and "ION" to describe the
;		output format:-
;		  if "Z" is capitalized, symbol as is (e.g., "He", not "he").
;		  if "z", symbol in small letters (e.g., "fe", "ca", etc).
;		  if "ion" is in lower case, ion in arabic numerals.
;		  if "ION" or "Ion", ion in roman numerals.
;		  symbol and ion separated by whatever "Z" and "Ion" are
;		    separated by.
;		* default is "ZION"
;	oform	[INPUT] kept here only for backwards-compatibility
;		* WARNING -- will be phased out >>very<< soon
;	_extra	[JUNK] here only to prevent crashing the system
;
;restrictions
;	cannot return the "+-" format (e.g., "O^+5")
;
;requires
;	INICON
;
;history
;	vinay kashyap (Nov98)
;	OFORM bug correction (VK; Apr99)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	allowed OFORM to be "Z<sep>number" and "Z<sep>ion" (VK; MMJul)
;	changed keyword OFORM to ZIFORM; converted to IDL5 notation (VK; AugMM)
;-

;	usage
nz=n_elements(Z) & nion=n_elements(ion)
if nz eq 0 then begin
  print,'Usage: zion2symb,z,ion,symbol,ziform=ziform'
  print,'  return spectroscopic symbol for given atomic number and ionic state'
  return
endif

;	initialize
atom=1 & rom=1 & inicon,atom=atom,roman=rom
natom=n_elements(atom) & nrom=n_elements(rom)
atom=['X',atom]		;always allow for the unknown
rom=['',rom]		;clearly the romans didn't use IDL

;	consistency checks
zz=[Z] & jon=0*zz & if nion eq 0 then ion=jon
if nion gt 0 and nion le nz then jon[0:nion-1]=ion[*] else jon[*]=ion[0:nz-1]
oo=where(zz gt natom or zz lt 1,moo) & if moo gt 0 then zz[oo]=0
oo=where(jon gt nrom or jon lt 1,moo) & if moo gt 0 then jon[oo]=0

;	output array
symbol=strarr(nz)

;	components
cz=atom[zz] & ci=rom[jon]

;	now decode the output format
fmtZ='Z' & fmtC='' & fmtI='ION'
if keyword_set(oform) and keyword_set(ziform) then begin
  message,'CONFUSION: keywords OFORM and ZIFORM both present',/info
  message,'ignoring OFORM',/info
endif
if keyword_set(oform) and not keyword_set(ziform) then begin
  message,'WARNING: keyword OFORM is now obsolete.  use keyword ZIFORM',/info
  ZIFORM=OFORM
endif
if keyword_set(ziform) then begin		;(ZIFORM
  fmt=ziform[0]
  kz=strpos(strupcase(fmt),'Z',0) & ki=strpos(strupcase(fmt),'I',0)
  if kz ge 0 then begin				;(Z format
    cc=strmid(fmt,kz,1)
    if cc eq 'z' then fmtZ='z'
  endif else begin				;)(Z not given
    fmtZ='' & cz[*]=''
  endelse					;Z FORM)
  if ki ge 0 then begin				;(I format
    cc=strmid(fmt,ki,1)
    if cc eq 'i' then fmtI='ion'
    if cc eq 'I' then fmtI='ION'
  endif else begin				;Ion)(I=number?
    ki=0
    kii=-1 & ii=0 & fmtI='num'
    while kii lt 0 and ii lt 10 do begin	;{what number
      kii=strpos(strupcase(fmt),strtrim(ii,2),kz+1) & ii=ii+1
    endwhile					;0-9}
    if kii ge 0 then ki=kii else begin		;(I=??
      fmtI='' & ci[*]=''
    endelse					;I form)
  endelse					;I FORM)
  if ki gt 0 then fmtC=strmid(fmt,kz+1,ki-kz-1)
endif						;ZIFORM)

;	figure out what ZIFORM was saying..
if fmtZ eq 'z' then cz=strlowcase(cz)
if fmtI eq 'ion' then ci=strlowcase(ci)
if fmtI eq 'num' then ci=strtrim(jon,2)

;	output
symbol=cz+fmtC+ci

;	print?
cc='Z:'+strtrim(zz,2)+' ION:'+strtrim(jon,2)+':= '+symbol
if n_params() lt 3 then for i=0,nz-1 do message,cc[i],/info

return
end
