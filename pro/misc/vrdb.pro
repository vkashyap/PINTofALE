pro vrdb,dbstr,pick,cols=cols,row=row,wfield=wfield, _extra=e
;+
;procedure	vrdb
;	print selected rows and columns from RDB table (see RDB.PRO)
;
;syntax
;	vrdb,dbstr,pick,cols=cols,/row,wfield=wfield
;
;parameters
;	dbstr	[INPUT; required] structure as read in by RDB.PRO
;	pick	[INPUT] array of position indices to print out
;		* default is [0]
;
;keywords
;	cols	[INPUT] array of column indices to print out
;		* default is to print all except DBSTR.HEAD
;		* if string array, then assumed to contain column names
;		  -- names do not have to be exact matches
;	row	[INPUT] if set, prints out the fields in a row
;		* default is to print them out in a column.
;		* -ve values will force all possible rows to be output
;	wfield	[INPUT] array of integers specifying the width at which
;		each column must be printed
;		* if not specified, tacks on +2 to the lengths deduced
;		  from the very first row to be printed
;	_exra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (FebMM)
;-

;	usage
ok='ok' & np=n_params() & mdb=n_elements(dbstr) & ndb=n_tags(dbstr)
if np eq 0 then ok='Insufficient parameters' else $
 if mdb eq 0 then ok='Input is undefined' else $
  if ndb eq 0 then ok='Input is not a structure'
if ok ne 'ok' then begin
  print,'Usage: vrdb,dbstr,pick,cols=cols,/row'
  print,'  print selected rows and columns of RDB table'
  if np ne 0 then message,ok,/info
  return
endif

;	inputs
dbnam=tag_names(dbstr)
irow=[0L] & if n_elements(pick) ge 1 then irow=[long(pick)]

;	error check and initialize
oc=where(dbnam eq 'HEAD',moc)
if moc eq 0 then message,'Input does not appear to be an RDB table',/info

;	keywords
;	COLS
ncols=n_elements(cols) & icols=lindgen(ndb)
if moc ne 0 then begin		;(if "HEAD" is a column name
  ok=where(dbnam ne 'HEAD',mok)
  if mok ne 0 then icols=icols(ok) else begin
    message,'Input table appears to be empty',/info
    return	;why bother with rest of stuff?
  endelse
endif				;MOC.NE.0)
if ncols ne 0 then begin		;(COLS has been defined
  szc=size(cols) & nszc=n_elements(szc)
  if szc(nszc-2) eq 7 then begin	;(column names are given
    icols=-1L
    for i=0L,ncols-1L do begin		;{for all specified columns
      oj=lonarr(ndb)
      for j=0L,ndb-1L do oj(j)=strpos(dbnam(j),strupcase(cols(i)),0)
      oo=where(oj ge 0,moo)
      if moo ne 0 then icols=[icols,oo]
    endfor				;I=0,NCOLS-1}
    mcols=n_elements(icols)
    if mcols eq 1 then begin
      message,'None of the column names matched',/info
      return	;nothing to print
    endif else begin
      icols=icols(1:*)
      icols=icols(uniq(icols,sort(icols)))
    endelse
  endif else begin			;)(column indices are given
    icols=[long(cols)]
    ok=where(icols ge 0 and icols lt ndb,mok)
    if mok eq 0 then begin
      message,'No columns in specified index range',/info
      return	;nothing to print
    endif else icols=icols(ok)
    icols=icols(uniq(icols,sort(icols)))
  endelse				;SZC)
endif					;NCOLS.NE.0)
ncol=n_elements(icols)
tabhead=strarr(ncol) & for j=0L,ncol-1L do tabhead(j)=dbnam(icols(j))
;
;	WFIELD
ww=intarr(ncol) & mw=n_elements(wfield)
if mw gt 0 then begin
  if mw eq 1 then ww(*)=fix(wfield(0)) else ww=fix(wfield(icols))
endif

;	check on number of rows
mrow=lonarr(ncol)
for i=0L,ncol-1L do mrow(i)=n_elements(dbstr.(icols(i)))
maxrow=max(mrow,min=minrow)
oo=where(irow lt 0,moo)
if moo ne 0 then begin
  message,'printing out all possible rows',/info
  irow=lindgen(minrow)
endif
oo=where(irow ge maxrow,moo)
if moo ne 0 then begin
  message,'Ignoring rows beyond legal range',/info
  oo=where(irow lt maxrow,moo)
  if moo gt 0 then irow=irow(oo) else begin
    message,'None of specified rows are legal',/info & return
  endelse
endif

;	print them out
nrow=n_elements(irow)
;
for i=0L,nrow-1L do begin		;{for each picked row
  ir=irow(i) & cc=strarr(ncol)
  for j=0L,ncol-1L do begin
    ic=icols(j)		;column index
    val=(dbstr.(ic))([ir])	;value in current column for selected row
    if ir ge mrow(j) then cval='*' else $
	cval=strtrim(val,2)	;convert to string
    cc(j)=cval			;store in array
    if i eq 0 then begin	;to define the formatting width, etc
      lcval=strlen(cval) > (strlen(tabhead(j)))
      if ww(j) le 2 then ww(j)=lcval+2
    endif
  endfor
  if i eq 0 then begin
    wmax0=max(ww) & wmax1=strlen(strtrim(max(irow),2))
    cw='(a'+strtrim(ww,2)+')' & cwmax='(a'+strtrim(wmax0+wmax1+2,2)+')'
  endif
  if keyword_set(row) then begin
    if i eq 0 then begin
      chead=''
      for j=0L,ncol-1L do chead=chead+string(tabhead(j),cw(j))
      print,chead
    endif
    cline=''
    for j=0L,ncol-1L do cline=cline+string(cc(j),cw(j))
    print,cline
  endif else begin
    for j=0L,ncol-1L do print,string(tabhead(j)+'['+strtrim(ir,2)+']',cwmax)+$
	' = '+strtrim(cc(j),2)	;string(cc(j),cw(j))
    print,''
  endelse
endfor					;I=0,NROW-1}

return
end
