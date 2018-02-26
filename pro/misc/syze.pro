function syze,c,oc=oc, _extra=e
;+
;function	syze
;	returns a vector that contains the size and type information a la
;	IDL information routine SIZE.  The diffrence is that this one
;	operates on strings only
;
;	output is always of longward type.
;	the first element contains the number of elements in C,
;	the next elements contain the string-length of each element, one
;	  element per element of C,
;	the last but one element contains the type code:
;		0 : undefined,
;		1 : I*1
;		2 : I*2
;		3 : I*4
;		4 : R*4
;		5 : R*8
;		-1 : expression
;		-2 : equation
;		7 : everything else
;	the last element contains the total number of characters in C
;
;syntax
;	sc=syze(c,oc=oc)
;
;parameters
;	c	[INPUT] string expression to decode
;
;keywords
;	oc	[OUTPUT] converted output, where applicable
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Dec98)
;	now catches sexagesimal (VK; Aug03)
;-

;	initialize
sc=lonarr(3)
cxp=byte('([{/*^}])')
ceq=byte('=')
cpt=byte('.')
cbt=byte('bB')
clw=byte('lL')
cfl=byte('eE')
cdb=byte('dD')
cpm=byte('+-')
csx=byte(':')

;	usage
np=n_params()
if np eq 0 then begin
  print,'Usage: sc=syze(c,oc=oc)'
  print,'  returns decoded size information for string C'
  return,sc
endif

;	check input
nc=n_elements(c)
if nc eq 0 then return,[nc,0L,1L] else sc=[nc,lonarr(nc),7L,0L]

;	populate the output array
for i=0L,nc-1L do begin			;{for each element in C

  ;	how long is each element?
  cc=strtrim(c(i),2)		;discard leading and trailing spaces
  nl=strlen(cc) & sc(i+1)=nl

  ;	of what type is each element?
  b=byte(cc)		;break each character into bytes

  ;	check for each type
  ctyp=7	;default
  if b(0) ge 47 and b(0) le 57 or b(0) eq cpt(0) or $
     b(0) eq cpm(0) or b(0) eq cpm(1) then begin
    ;	begins with a number or a point or a +-
    ctyp=2	;I*2, by default
    if b(0) ne cpm(0) and b(0) ne cpm(1) then begin
      if nl gt 5 then ctyp=3		;change the default to I*4
    endif else if nl gt 6 then ctyp=3	;change the default to I*4
    if b(nl-1L) eq cbt(0) or b(nl-1L) eq cbt(1) then ctyp=1	;byte
    if b(nl-1L) eq clw(0) or b(nl-1L) eq clw(1) then ctyp=3	;longword
    oo=where(b eq cpt(0),moo)
      if moo eq 1 then ctyp=4					;float
      if moo gt 1 then ctyp=7
    oo=where(b eq cfl(0) or b eq cfl(1),moo)			;float
      if moo eq 1 then ctyp=4
      if moo gt 1 then ctyp=7
    oo=where(b eq cdb(0) or b eq cdb(1),moo)
      if moo eq 1 then ctyp=5					;double
      if moo gt 1 then ctyp=7
  endif
  oxp=where(b eq cxp(0) or b eq cxp(1) or b eq cxp(2) or $
	b eq cxp(3) or b eq cxp(4) or b eq cxp(5) or $
	b eq cxp(6) or b eq cxp(7) or b eq cxp(8),moxp)		;expression?
  oeq=where(b eq ceq(0),moeq)		;equation?
  osx=where(b eq csx(0),mosx)		;sexagesimal
  if moxp gt 0 then ctyp=-1		;expression
  if moeq gt 0 then ctyp=-2		;equation
  if mosx gt 0 then ctyp=7		;sexagesimal
  ;
  if i eq 0 then sc(nc+1)=ctyp else begin
    if ctyp gt sc(nc+1) then begin
      message,'input array is of mixed type',/info
      sc(nc+1)=ctyp
    endif
    if ctyp lt 0 then begin
      message,'input is an array of expressions',/info
      sc(nc+1)=-1				;once an expression..
    endif
  endelse
endfor					;I=0,NC-1}

;	and the final element of the output..
sc(nc+2)=total(sc(1:nc))

;	convert if needed
if sc(nc+1) gt 0 and sc(nc+1) lt 7 then begin
  case sc(nc+1) of
    1: oc=byte(fix(c))
    2: oc=fix(c)
    3: oc=long(c)
    4: oc=float(c)
    5: oc=double(c)
    else: oc=c
  endcase
endif

return,sc
end
