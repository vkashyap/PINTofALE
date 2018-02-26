function aria,p1,p2,p3,p4,pickar=pickar,comment=comment, _extra=e
;+
;function	aria
;	puts together the effective areas and corresponding wavelengths
;	into a single structure and returns that structure
;
;syntax
;	A=aria(p1,p2,p3,p4,pickar=pickar,comment=comment)
;
;parameters
;	P1..P4	[INPUT] see USAGE below for description
;
;keywords
;	pickar	[INPUT] if vector, AND only one structure (and nothing else)
;		is input, then returns all the components specified
;		* if scalar, on the other hand, returns ALL BUT the
;		  specified component
;	comment	[INPUT] tacks on comment if given
;	_extra	[JUNK] here only to prevent crashing the program
;
;usage
;	the effective area is defined by the wavelengths, the area at said
;	wavelengths, and the grating order for which this is valid.  so if
;	AREA(WVL;ORDER) are input, out should come a new structure of the
;	form {AREA:AREA, WVL:WVL, ORDER:ORDER, COMMENT:COMMENT}.
;		OUTPUT.(0).AREA		<== AREA
;		OUTPUT.(0).WVL		<== WVL
;		OUTPUT.(0).ORDER	<== ORDER
;		OUTPUT.(0).COMMENT	<== COMMENT
;
;	If such a structure already exists, the supplied AREA(WVL;ORDER) must
;	be concatenated to it.  hence:
;
;	to create a brand new structure
;		arstr=aria(AREA,WVL,ORDER)
;
;	to append AREA(WVL;ORDER) to pre-existing area structure (note that
;	it doesn't matter just where "OLDARSTR" is in the calling sequence)
;		arstr=aria(OLDARSTR, AREA,WVL,ORDER)
;		arstr=aria(AREA, OLDARSTR, WVL,ORDER)
;		arstr=aria(AREA,WVL, OLDARSTR, ORDER)
;		arstr=aria(AREA,WVL,ORDER, OLDARSTR)
;
;	to append ARSTR2 to ARSTR1 by calling this program recursively, once
;	for each segment of ARSTR2
;		arstr=aria(ARSTR1,ARSTR2)
;
;	to selectively choose (with vector PICKAR)
;	or delete (with scalar PICKAR)
;	compoenents (works only when just one parameter is supplied)
;		arstr=aria(OLDARSTR,PICKAR=PICKAR)
;
;subroutines
;	ARIA (recursive calls while adding two structures)
;
;history
;	vinay kashyap (Nov98)
;	changed keyword PICK to PICKAR (VK; JunMM)
;-

;	usage
np=n_params()
n1=n_elements(p1) & n2=n_elements(p2) & n3=n_elements(p3) & n4=n_elements(p4)
t1=n_tags(p1) & t2=n_tags(p2) & t3=n_tags(p3) & t4=n_tags(p4)
ok='ok'
if np lt 1 then ok='Not enough parameters supplied' else $
 if n1 eq 0 and n2 eq 0 then ok='undefined parameters' else $
  if n1+n2+n3+n4 eq 0 then ok='all parameters appear to be empty'
astr={AREA:0., WVL:0., ORDER:-1L, COMMENT:'help'}
astr=create_struct('ARIA_',astr)
if ok ne 'ok' then begin
  print,'Usage: A=aria(area,wvl,order,comment=comment)          OR'	;type=0
  print,'       A=aria(oldA, area,wvl,order,comment=comment)    OR'	;type=1
  print,'       A=aria(area,wvl,order, oldA,comment=comment)    OR'	;type=1
  print,'       A=aria(area, oldA, wvl,order,comment=comment)   OR'	;type=1
  print,'       A=aria(area,wvl, oldA, order,comment=comment)   OR'	;type=1
  print,'       A=aria(area,wvl,oldA, order,comment=comment)    OR'	;type=1
  print,'       A=aria(oldA, pickar=pickar,comment=comment)     OR'	;type=1
  print,'       A=aria(A1,A2,comment=comment)                   OR'	;type=2
  print,'  concatenate effective areas'
  message,ok,/info
  return,astr
endif

;	figure out which case
ti=0 & type=0
if t1 gt 0 then ti=ti+1
if t2 gt 0 then ti=ti+1
if t3 gt 0 then ti=ti+1
if t4 gt 0 then ti=ti+1
if keyword_set(comment) then comment=strtrim(comment,2) else comment=''
case ti of				;{depends on how many structures
  0: begin				;{TYPE=0
    type=0
    if n1 eq 0 then ok='AREAs must be supplied' else $
     if n2 eq 0 then ok='Need Area(WVL)' else $
      if n1 ne n2 then ok='AREA and WVL incompatible' else $
       if n3 gt 1 then ok='do not understand how ORDER can be an array'
    if n3 eq 0 then p3=1
  end					;AREA(WVL;ORDER)}
  1: begin				;{TYPE=1
    type=1
    if t1 gt 0 then begin
      A0=p1
      subtype=1
      if n2+n3+n4 eq 0 then begin
	subtype=-1
      endif else begin
        if n2 eq 0 or n3 eq 0 then ok='AREA(WVL) needs to be given' else $
         if n2 ne n3 then ok='AREA and WVL are incompatible' else $
	  if n4 gt 1 then ok='ORDER must be a scalar'
        if n4 eq 0 then p4=1
        astr=A0
      endelse
    endif
    if t2 gt 0 then begin
      A0=p2
      if n1 eq 0 or n3 eq 0 then ok='AREA(WVL) needs to be given' else $
       if n1 ne n3 then ok='AREA and WVL are incompatible' else $
	if n4 gt 1 then ok='ORDER must be a scalar'
      if n4 eq 0 then p4=1
      subtype=2
      astr=A0
    endif
    if t3 gt 0 then begin
      A0=p3
      if n1 eq 0 or n2 eq 0 then ok='AREA(WVL) needs to be given' else $
       if n1 ne n2 then ok='AREA and WVL are incompatible' else $
	if n4 gt 1 then ok='ORDER must be a scalar'
      if n4 eq 0 then p4=1
      subtype=3
      astr=A0
    endif
    if t4 gt 0 then begin
      A0=p4
      if n1 eq 0 or n2 eq 0 then ok='AREA(WVL) needs to be given' else $
       if n1 ne n2 then ok='AREA and WVL are incompatible' else $
	if n3 gt 1 then ok='ORDER must be a scalar'
      if n3 eq 0 then p3=1
      subtype=4
      astr=A0
    endif
  end					;oldA+AREA(WVL;ORDER)}
  2: begin				;{TYPE=2
    type=2
    if t1 gt 0 then A0=p1 else $
     if t2 gt 0 then A0=p2 else $
      if t3 gt 0 then A0=p3
    if t2 gt 0 then A1=p2
    if t3 gt 0 then A1=p3
    if t4 gt 0 then A1=p4
    astr=A0
  end					;A1+A2}
  else: ok='at most two structures, please'	;can't handle
endcase					;TI}
if ok ne 'ok' then begin
  message,ok,/info & return,astr
endif

case type of				;{assign variables and make structure
  0: begin
    astr={AREA:p1, WVL:p2, ORDER:p3, COMMENT:comment}
    return,create_struct('ARIA0',astr)
  end
  1: begin
    if subtype eq -1 then begin		;(pick selectively..
      nA=n_tags(A0) & j=0L
      szp=size(pickar) & nszp=n_elements(szp)
      setdele=0				;by default, delete none
      if szp(nszp-2) ne 0 and szp(0) eq 0 then setdele=1	;defined scalar
      if szp(nszp-2) eq 0 then pickar=lindgen(nA)
      for i=0,nA-1 do begin
	incthis=1			;include this component
	if setdele eq 1 then begin
	  if pickar(0) eq i then incthis=0	;exclude this component
	endif else begin
	  ok=where([pickar] eq i,mok)
	  if mok eq 0 then incthis=0	;exclude this component
	endelse
	tmp=A0.(i) & nt=n_tags(tmp)
	if nt lt 3 then incthis=0	;exclude this component
	if incthis then begin
	  area=tmp.(0)
	  wvl=tmp.(1)
	  order=tmp.(2)
	  if nt gt 3 then comm=tmp.(3) else comm=comment
	  aa={AREA:area, WVL:wvl, ORDER:order, COMMENT:comm}
	  tname='ARIA'+strtrim(j,2)
	  if j eq 0 then astr=create_struct(tname,aa) else $
	   astr=create_struct(astr,tname,aa)
	  j=j+1
	endif
      endfor
      return,astr
    endif				;SUBTYPE=-1)

    if subtype eq 1 then begin
      area=p2 & wvl=p3 & order=p4
    endif
    if subtype eq 2 then begin
      area=p1 & wvl=p3 & order=p4
    endif
    if subtype eq 3 then begin
      area=p1 & wvl=p2 & order=p4
    endif
    if subtype eq 4 then begin
      area=p1 & wvl=p2 & order=p3
    endif

    tnam=tag_names(A0)
    j=-1L & moj=1
    while moj gt 0 do begin
      j=j+1 & arnam='ARIA'+strtrim(j,2) & oj=where(arnam eq tnam,moj)
    endwhile
    astr={AREA:area, WVL:wvl, ORDER:order, COMMENT:comment}

    return,create_struct(A0,arnam,astr)
  end
  2: begin
    astr=A0 & nt=n_tags(A1)
    for i=0,nt-1 do begin
      tmp=A1.(i) & mt=n_tags(tmp)
      if mt eq 0 then begin
	message,'(SECOND).('+strtrim(i,2)+') not a structure',/info
      endif else begin
	if mt ge 3 then begin
	  area=tmp.(0)
	  wvl=tmp.(1)
	  order=tmp.(2)
	  if mt gt 3 then comm=tmp.(3) else comm=comment
	  astr=aria(astr,area,wvl,order,comment=comm)
	endif
      endelse
    endfor
    return,astr
  end
  else: message,'bug!'
endcase					;by type}

return,astr
end
