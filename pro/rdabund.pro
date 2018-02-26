function rdabund,filnam,comment=comment,sep=sep,inorm=inorm,norm=norm,$
	_extra=e
;+
;function	rdabund
;	return abundance values read from file
;
;	disk file must have one entry per line, with either atomic symbol
;	or atomic number denoting the element (or both), and the actual
;	abundance also.  neither the order of the columns nor the order
;	of the rows matter.  the fields must be separated by "SEP", usually
;	a <tab>.  missing elements are filled out with values from Anders &
;	Grevesse (see GETABUND).
;
;syntax
;	abund=rdabund(filnam,comment=comment,sep=sep,inorm=inorm,norm=norm)
;
;parameters
;	filnam	[INPUT; required] name of file containing the abundances
;
;keywords
;	comment	[OUTPUT] normally, any line prefixed by ";" "%" "#" or "/"
;		are ignored.  these lines are stored as a string array and
;		output via this variable
;	sep	[INPUT] single-character field separator
;		* if set to "." "+" "-" or any of the alphanumerics, is ignored
;		* if not set, tries <sp>, <tab>, "," and ":" in sequence
;	inorm	[INPUT] format of the abundances in the input file
;		* see GETABUND for description
;		* if INORM is defined, first converted to n(Z)/n(H), then
;		  converted back to whatever NORM is set to
;	norm	[INPUT] exactly as in GETABUND
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	GETABUND
;	SYZE
;
;history
;	vinay kashyap (Dec98)
;-

;	initialize
defabu=getabund('anders & grevesse') & nabu=n_elements(defabu)
abund=defabu

;	usage
nf=n_elements(filnam) & szf=size(filnam) & nszf=n_elements(szf)
if nf ne 1 or szf(nszf-2) ne 7 then begin
  print,'Usage: abund=rdabund(filename,comment=comment,sep=sep)'
  print,'  returns abundances read in from file'
  return,defabu
endif

;	check separator
cno='+-.0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
ss=''
if keyword_set(sep) then begin
  ;	check SEP is not one of the verboten
  szs=size(sep) & nszs=n_elements(szs)
  if szs(nszs-2) eq 7 then begin
    bno=byte(cno) & nno=n_elements(bno)
    bss=byte(sep) & nss=n_elements(bss) & iss=intarr(nss)+1
    for i=0,nss-1 do for j=0,nno-1 do if bss(i) eq bno(j) then iss(i)=0
    oi=where(iss gt 0,moi)
    if moi gt 0 then begin
      iss=iss(oi) & bss=bss(oi) & nss=moi
      for i=0,moi-1 do if iss(i) eq 1 then ss=ss+string(bss(i))
    endif
  endif
endif
if not keyword_set(ss) then begin
  ss=' 	,:'					;<sp><tab>,:
  bss=byte(ss) & nss=n_elements(bss) & iss=intarr(nss)+1
endif

;	read from file
openr,udsk,filnam,/get_lun			;{open for input
ilin=0L
while not eof(udsk) do begin			;{read file
  nextline: line=''
  readf,udsk,line & line=strtrim(line,2)
  c1=strmid(line,0,1)
  if c1 eq '#' or c1 eq ';' or c1 eq '%' or c1 eq '/' then begin
    if keyword_set(comment) then comment=[comment,line] else comment=line
    goto,nextline
  endif
  ilin=ilin+1L
 
  k=0L
  while k lt nss do begin	;{check the separators in sequence
    csep=strmid(ss,k,1)
    cc=str_sep(line,csep) & ncc=n_elements(cc)
    if ncc gt 1 then k=nss else k=k+1L
  endwhile			;K<NSS}
  if ncc lt 2 then begin
    print,filnam+':	'+line
    message,'incorrect separator or insufficient columns',/info
    goto,nextline
  endif

  zz=0 & zab=0. & iz=0
  for i=0,ncc-1 do begin
    if cc(i) ne '' then begin
      szc=syze(cc(i),oc=oc)
      if szc(2) eq 7 and szc(3) le 3 then begin
	symb2zion,cc(i),zz,jon & iz=1		;atomic symbol
      endif
      if szc(2) eq 2 then begin
	if keyword_set(zz) then begin
	  if i gt 0 and oc ne ilin and iz eq 0 and oc le nabu then $
		zz=oc				;atomic number
	endif else zz=oc
      endif
      if szc(2) eq 4 or szc(2) eq 5 then zab=oc	;abundance
    endif
  endfor

  if zz gt 0 then begin
    if zz gt nabu then begin
      tmp=fltarr(zz) & tmp(0)=defabu & defabu=tmp
      defabu = defabu > (min(defabu(where(defabu gt 0))))
      tmp=fltarr(zz) & tmp(0)=abund & abund=tmp & nabu=zz
    endif
    if zz le nabu then abund(zz-1)=zab
  endif else begin
    if keyword_set(comment) then comment=[comment,line] else comment=line
  endelse

  if keyword_set(inorm) then begin
    ;	convert to "true" abundances
    if inorm(0) gt 0 and inorm(0) lt 1 then zab=zab/inorm(0)
    if inorm(0) gt 1 then zab=10.^(zab-inorm(0))
    if inorm(0) lt 0 and inorm(0) ge -1 then zab=zab*defabu(zz-1)
    if inorm(0) lt -1 then zab=10.^(zab)*defabu(zz-1)
	;WARNING: in the last 2 cases, if zz > n_elements(defabu)@initial,
	;the output will be garbage!
    abund(zz-1)=zab	;overwrite
  endif

endwhile					;EOF(udsk)}
close,udsk & free_lun,udsk			;close file}

;	renorm
if keyword_set(norm) then begin
  if norm(0) gt 0 and norm(0) le 1 then abund=norm(0)*abund
  if norm(0) gt 1 then begin
    ok=where(abund gt 0,mok)
    oo=where(abund eq 0,moo)
    if mok gt 0 then abund(ok)=alog10(abund(ok))+norm(0)
    if mok gt 0 then minab=min(abund(ok)) else minab=0.
    if moo gt 0 then abund(oo)=minab-10.
  endif
  if norm(0) lt 0 and norm(0) ge -1 then begin
    ok=where(abund gt 0 and defabu gt 0,mok)
    oo=where(abund eq 0 or defabu eq 0,moo)
    if mok gt 0 then abund(ok)=abund(ok)/defabu(ok)
    if moo gt 0 then abund(oo)=0.
  endif
  if norm(0) lt -1 then begin
    ok=where(abund gt 0 and defabu gt 0,mok)
    oo=where(abund eq 0 or defabu eq 0,moo)
    if mok gt 0 then abund(ok)=alog10(abund(ok)/defabu(ok))
    if mok gt 0 then minab=min(abund(ok)) else minab=0.
    if moo gt 0 then abund(oo)=minab-10.
  endif
endif

return,abund
end
