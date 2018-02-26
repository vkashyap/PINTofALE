function rdpha,infile,ext=ext, _extra=e
;+
;function	rdpha
;	return a structure containing the variables in an XSPEC style
;	PHA file.  (doesn't really have to be a PHA file.  any FITS
;	table file will do.)
;
;parameters
;	infile	[INPUT; required] name of PHA file
;
;keywords
;	ext	[INPUT] if set, only reads in the specified extension
;	_extra	[JUNK] here only to avoid crashing the program
;
;restrictions
;	* requires the IDLASTRO library
;	* empty extensions are treated as EOFs
;	* can only decipher tables, not images
;	* cannot handle Type II PHA files
;
;history
;	vinay kashyap (Jan MM)
;-

;	usage
ok='ok' & nfil=n_elements(infile) & szf=size(infile) & nszf=n_elements(szf)
if nfil eq 0 then ok='Input PHA file undefined' else $
 if nfil gt 1 then ok='Cannot handle array of input files' else $
  if szf(nszf-2) ne 7 then ok='filename illegible'
if ok ne 'ok' then begin
  print,'Usage: pha=rdpha(infile)'
  print,'  read in XSPEC-style PHA file into IDL structure'
  if n_params() ne 0 then message,ok,/info
  return,-1L
endif

;	keywords
iext=0L & if keyword_set(ext) then iext=-long(ext)

;	initialize
pha=0 & cols=''

;	read in the file
while 1 do begin		;{infinite loop, break using GOTO
  tmp=readfits(infile,htmp,exten=abs(iext))

  ntmp=n_elements(tmp) & szt=size(tmp) & nszt=n_elements(szt)
  if abs(iext) gt 0 and ntmp eq 1 and szt(nszt-2) eq 2 and tmp(0) eq -1 then $
	goto,fileend	;(nothing to see here, move along
  
  xtnam=strupcase(strtrim(sxpar(htmp,'EXTNAME'),2))
  if not keyword_set(pha) then pha=create_struct(xtnam+'_',htmp) else $
    pha=create_struct(pha,xtnam+'_',htmp)

  xttyp=strupcase(strtrim(sxpar(htmp,'XTENSION'),2))
  if xttyp eq 'TABLE' or xttyp eq 'BINTABLE' or xttyp eq 'A3DTABLE' then begin
  	;(this extension contains a binary table, so extract the columns
    nflds=sxpar(htmp,'TFIELDS')
    for i=0,nflds-1 do begin		;{extract each column separately
      colnam=strupcase(strtrim(sxpar(htmp,'TTYPE'+strtrim(i+1,2)),2))
      ;	check if this column name is a legal IDL variable
      j=legalvar(colnam) & if j eq 0 then colnam='COL_'+strtrim(i+1,2)
      if not keyword_set(cols) then cols=[colnam] else begin
	oo=where(colnam eq cols,moo)
	if moo gt 0 then colnam=colnam+'_'+xtnam	;make unique columns
	cols=[cols,colnam]
      endelse
      pha=create_struct(pha,colnam,fits_get(htmp,tmp,i+1))
    endfor				;I=0,NFLDS-1}
  endif else begin			;)(not a table, just extract blind
    pha=create_struct(pha,'DATA_'+xtnam,tmp)
  endelse				;XTTYP)

  if iext lt 0 then goto,fileend	;(skip if needed!
  iext = iext + 1L	;look in next extension
endwhile			;infinite read-in loop}
fileend:	;the GOTO comes out here))

return,pha
end
