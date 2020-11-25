function rdb,rdbfil,cols=cols,comm=comm,fldsep=fldsep,rdbdir=rdbdir,$
	forced=forced,verbose=verbose,_extra=e
;+
;function	rdb
;	return a structure containing the output of ascii table
;	in RDB format
;
;	RDB format is assumed to be:
;	a bunch of comments preceded by "#"s in first column
;	a header stating the names of each column
;	a format description specifying the type of the column
;	(S for strings, N for numeric)
;	and each column is separated by a "<tab>"
;
;	some special things understood by this program are:
;	-- blank lines are treated as comments
;	-- different symbols (e.g., "%") may be used as comment characters
;	-- output comes out as double or float or long depending on what is
;	   in the first row
;	-- double precision may be forced by prepending a ">" to "N"
;
;syntax
;	tab=rdb(filename,cols=cols,comm=comm,fldsep=fldsep,rdbdir=rdbdir,$
;	/forced,verbose=verbose)
;
;parameters
;	rdbfil	[INPUT; required] name of rdb file
;		* the ".rdb" extension is not necessary
;		* a trailing "." is removed, so if there's a file without
;		  an extension, set RDBFIL to "filename."
;		* and if there's a file named "filename.", set RDBFIL to
;		  "filename.." and so on.
;
;keywords
;	cols	[INPUT] string or integer array of columns to extract
;		* if string, column names must match *exactly*, but
;		  case is not important
;		* if integers, remember that this is IDL, so first column
;		  is indexed as 0
;		* if /COLS is set, returns all columns, so if you want
;		  just the second column, remember to say COLS=[1] and
;		  not COLS=1
;	comm	[INPUT] comment character
;		* default is a "#"
;	fldsep	[INPUT] field separator
;		* default is a <TAB>
;	rdbdir	[INPUT] if given, prepends RDBFIL with this pathname
;	forceD	[INPUT] if set, forces all Numeric columns to be double precision
;	verbose	[INPUT] if set, blabbers incoherently along the way
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	STR_2_ARR
;	CREATE_STRUCT
;	KILROY
;	LEGALVAR
;
;restrictions
;	works only on UNIX
;
;history
;	vinay kashyap (MIM.XII)
;	illegal column names are now checked for (VK; MM.I)
;	bug: preceding spaces turned values into zeros; converted all
;	  array subscripting to IDL5 (VK; JulMM)
;	changed call to STR_2_ARR from STR2ARR (VK; MMV.IV)
;	now looks at only the first character of format column field
;	  to decide type (VK; MMXIII.VIII)
;	now basically returns only long integer and doubles (VK; MMXIII.XII)
;	added keyword FORCED (VK; MMXX.VII)
;	now handles NaNs and Infs more gracefully (VK; MMXX.XI)
;-

;	usage
ok='ok' & szf=size(rdbfil) & nszf=n_elements(szf)
if szf(nszf-2) ne 7 then ok='illegible filename' else $
 if szf(0) gt 0 then ok='cannot handle array of RDB files'
if ok ne 'ok' then begin
  print,'Usage: t=rdb(filename,cols=cols,fldsep=fldsep,rdbdir=rdbdir,/forced,verbose=verbose)'
  print,'  return structure containing tabulated values'
  if n_params() gt 0 then message,ok,/info
  return,{HEAD: ok}
endif

;	input
ixtn=rstrpos(rdbfil,'.')
if ixtn lt 0 then rdbfil=rdbfil+'.rdb' else begin
  xtn=strmid(rdbfil,ixtn+1)
  if strlowcase(xtn) ne 'rdb' then message,$
	rdbfil+': input file may not be RDB compatible',/info
  if xtn eq '' then rdbfil=strmid(rdbfil,0,ixtn)	;remove trailing "."s
endelse

;	keywords
;
allcols=1	;default is to return all the columns
mcols=n_elements(cols)
if mcols gt 0 then begin	;(select out specific columns
  allcols=0
  szc=size(cols) & nszc=n_elements(szc) & ctyp=szc(nszc-2)
  case ctyp of			;{string, or integer?
    7: ccols=[strtrim(cols[*],2)]	;string
    1: jcols=[long(cols[*])]		;byte
    2: begin				;integer
      if szc[0] eq 0 and cols[0] eq 1 then begin
	;special case: what if /COLS were set?
	allcols=1	;then return all columns
      endif else jcols=[long(cols[*])]
    end
    3: jcols=[long(cols[*])]		;long integer
    4: jcols=[long(cols[*])]		;float
    5: jcols=[long(cols[*])]		;double
    else: allcols=1			;extract all
  endcase			;CTYP}
endif				;MCOLS>0)
;
comment='#'	;default is a "#"
ok='ok' & szc=size(comm) & nszc=n_elements(szc)
if szc(nszc-2) eq 7 then begin
  lc=strlen(comm[0])
  if lc gt 1 then ok='comment must be single character; ignoring input' else $
   if szc[0] ge 1 then ok='comment characters input as array; keeping first'
  comment=comm[0]
endif else if szc(nszc-2) ne 0 then ok=$
	'comment character must be of type string; ignoring input'
if ok ne 'ok' then message,ok,/info
;
;
fsep='	'	;default is a <TAB>
ok='ok' & szs=size(fldsep) & nszs=n_elements(szs)
if szs(nszs-2) eq 7 then begin
  ls=strlen(fldsep[0])
  if ls gt 1 then ok='separator must be single character; ignoring input' else $
   if szs[0] ge 1 then ok='field separators input as array; keeping only first'
  fsep=fldsep[0]
endif else if szs(nszs-2) ne 0 then ok=$
	'field separator must be of type string; ignoring input'
if ok ne 'ok' then message,ok,/info
;
rdir=''		;default is nothing
szd=size(rdbdir) & nszd=n_elements(szd)
if szd(nszd-2) eq 7 then begin
  if szd[0] ge 1 then message,$
	'RDBDIR input as array; keeping only first element',/info
  rdir=rdbdir[0]
  if strmid(rdir,0,1) eq '$' then rdir=getenv(strmid(rdbdir[0],1,strlen(rdbdir[0])-1))
  if strmid(rdir,strlen(rdir)-1,1) ne '/' then rdir=rdir+'/'
endif
;
v=1		;at least report a wee bit
if keyword_set(verbose) then v=fix(verbose) > 1

;	read from:
frdb=rdir+rdbfil & header='# HISTORY: Filename = '+frdb
if v ge 1 then message,'Reading from file: '+frdb,/info

;	this block currently commented out because of a FINDFILE bug:
;	FINDFILE doesn't understand "~"
;;	check for existence
;tmp=findfile(frdb,count=nfil)
;if nfil ne 1 then begin
;  message,frdb+': file not found',/info & return,{HEAD: frdb}
;endif

;	figure out how many rows of data
cmd='grep -v "^'+comment+'" '+frdb+' | grep -v ^$ | wc -l'
spawn,cmd,cline & nline=long(cline[0])-2L	;-2 because of colname & fmt
if v ge 2 then print,'Number of rows in datafile:',nline

;	all the above was frippery.  this is where the action occurs.
openr,urdb,frdb,/get_lun	;(open for input
while not eof(urdb) do begin	;{read from file

  line='' & readf,urdb,line
  if v ge 5 then print,line

  c0=strmid(strtrim(line,2),0,1)	;does first character denote comment?
  if c0 eq '' then c0=comment		;blank lines are always "comments"
  if c0 eq comment then begin
    if not keyword_set(header) then header=line else header=[header,line]
    goto,rdnextln			;{(oh humor me.
  endif

  ;	first uncommented line is list of columns
  if not keyword_set(colhd) then begin
    colhd=1	;no more column heads
    colnam=strtrim(str_sep(line,fsep),2) & ncols=n_elements(colnam)
    colnamdef='col_'+strtrim(1+lindgen(ncols),2)
    icols=lonarr(ncols)+1	;1: return column; 0: ignore column
    if not keyword_set(allcols) then begin
      icols[*]=0
      for i=0,mcols-1 do begin
	if keyword_set(jcols) then icols(jcols[i])=1 else begin
	 oo=where(strlowcase(ccols[i]) eq strlowcase(colnam),moo)
	 if moo eq 0 then message,acols[i]+': Column not found' else icols[oo[0]]=1
	endelse
      endfor
    endif
    oo=where(icols gt 0,mcols)
    header=[header,'# COLUMNS: '+line]

    ;	if any of the column names are illegal, use generic names
    ll=legalvar(colnam)
    for i=0,mcols-1 do begin
      j=ll[i]
      if j eq 0 then begin
	message,'Replacing name of column '+strtrim(i+1,2)+' [ '+colnam[i]+$
		' ] by [ '+colnamdef[i]+' ] ',/info
	colnam[i]=colnamdef[i]
      endif
    endfor
    goto, rdnextln			;{(same goto as above.
  endif

  ;	second uncommented line is column format
  if not keyword_set(colfm) then begin
    colfm=1	;no more column formats
    colfrm=str_sep(line,fsep) & nfrm=n_elements(colfrm)
    mfrm=nfrm < ncols
    if mfrm ne ncols then message,$
	'Mismatch between number of columns and column formats',/info
    fmtcode=strarr(ncols)+'S'	;default is to extract as string
    for i=0,mfrm-1 do if strpos(strmid(strtrim(strupcase(colfrm[i]),2),0,1),'N',0) ge 0 then $
	fmtcode[i]='N'
    ;if nfrm gt ncols then message,$
    ;	'more format codes than column names; ignoring',/info else $
    ;	for i=0,nfrm-1 do $
    ;	if strpos(strupcase(colfrm[i]),'N',0) ge 0 then fmtcode[i]='N'
    header=[header,'# FORMATS: '+line]
    goto, rdnextln
  endif

  ;	data begin at third uncommented line
  ;	figure out details of column format and define output structure here
  if not keyword_set(kline) then begin
    kline=1L	;we're earnestly reading the data from now on
    data=strarr(ncols)
    cold=str_sep(line,fsep) & ndata=n_elements(cold)
    if ndata le ncols then data[0]=cold else data=cold[0:ncols-1L]

    for i=0L,ncols-1L do begin		;{go through each column

      if icols[i] gt 0 then begin	;(read this column?
	val=data[i]

	if fmtcode[i] eq 'N' then begin	;(figure out the true format
	  ;	if there's a "." or "e", it's R4 or R8, else I2 or I4
	  ;		(well actually the default is now I4, so never I2)
	  ;	if "e" and abs(exponent) > 20, then R8
	  ;	if I and:
	  ;		longer than 4 chars, or
	  ;		format code says longer than 5 chars, or
	  ;		format code has a ">" sign,
	  ;	then I4
	  ;	if R and:
	  ;		longer than 8 chars, or
	  ;		format code says longer than 8 chars, or
	  ;		format code has a ">" sign,
	  ;	then R8

	  if not keyword_set(forced) then fmt='I4' else fmt='R8'	;default is I4 unless FORCED is set
	  idot=strpos(data[i],'.',0) & iexp=strpos(strlowcase(data[i]),'e',0)
	  if idot ge 0 or iexp ge 0 then fmt='R8'	;default double
	  if iexp ge 0 then begin
	    xpo=abs(float(strmid(data[i],iexp+1)))
	    if xpo lt 20 then fmt='R4'
	  endif
	  ifmt=strpos(strupcase(colfrm[i]),'N',0) & igt=strpos(colfrm[i],'>',0)
	  tmp=byte(strmid(colfrm[i],0,ifmt)) & flen=4
	  if tmp[0] ge 48 and tmp[0] le 57 then $
		flen=fix(strmid(colfrm[i],0,ifmt))
	  dlen=strlen(strtrim(data[i],2))
	  if (dlen gt 4 or flen gt 4 or igt ge 0) and fmt eq 'I2' then fmt='I4'	;it is never I2
	  if (dlen gt 8 or flen gt 8 or igt ge 0) and fmt eq 'R4' then fmt='R8'
	  ;	check for NaNs and Infs
	  iNaN=strpos(strlowcase(data[i]),'nan',0)
	  iInf=strpos(strlowcase(data[i]),'inf',0)
	  iInfty=strpos(strlowcase(data[i]),'infinity',0)
	  if iNan ge 0 or iInf ge 0 then fmt='R4'
	  if iInfty ge 0 then fmt='R8'

	  if keyword_set(forced) then begin
	    fmt='R8'
	    if forced[0] eq 4 then fmt='R4'	;undocumented 
	    if forced[0] eq 2 then fmt='I4'	;undocumented 
	  endif

	  cc=strlowcase(strtrim(val,2))
	  case cc of
	    'nan': if fmt eq 'R8' then val=!values.D_NAN else val=!values.F_NAN
	    'inf': if fmt eq 'R8' then val=!values.D_INFINITY else val=!values.F_INFINITY
	    'infinity': if fmt eq 'R8' then val=!values.D_INFINITY else val=!values.F_INFINITY
	    else: begin
	      if fmt eq 'I2' then val=(str_2_arr(strtrim(data[i],2)))[0]
	      if fmt eq 'I4' then val=(str_2_arr(strtrim(data[i],2),/i4))[0]
	      if fmt eq 'R4' then val=(str_2_arr(strtrim(data[i],2),/r4))[0]
	      if fmt eq 'R8' then val=(str_2_arr(strtrim(data[i],2),/r8))[0]
	    end
	  endcase

	  fmtcode[i]=fmt
          if v ge 5 and i eq ncols-1L then begin
	    print,'COLUMN / FORMAT:'
            for k=0,mcols-1 do print,colnam[k],' : ',fmtcode[k]
	  endif
	endif				;FMTCODE)

	;	make array of right size
	var=make_array(nline,size=size([val]),value=val)
	;	define the output structure
	if n_tags(rdbtab) eq 0 then rdbtab=create_struct(colnam[i],var) else $
	  rdbtab=create_struct(rdbtab,colnam[i],var)

      endif				;ICOLS[I]>0)

    endfor				;I=0,NCOLS-1}
    goto, rdnextln			;{(again the same goto
  endif

  ;	and now chug through the rest of the table
  kline=kline+1L		;line number
  rdbcol=1L			;column number
  if kline eq 100*fix(kline/100.) then kilroy; was here.
  data=strarr(ncols)
  cold=str_sep(line,fsep) & ndata=n_elements(cold)
  if ndata le ncols then data[0]=cold else data=cold[0:ncols-1L]
  for i=0L,ncols-1L do begin		;{go through each column
    if icols[i] gt 0 then begin		;(read this column?
      val=data[i]
      
      cc=strlowcase(strtrim(val,2))
      case cc of
        'nan': if fmt eq 'R8' then val=!values.D_NAN else val=!values.F_NAN
        'inf': if fmt eq 'R8' then val=!values.D_INFINITY else val=!values.F_INFINITY
        'infinity': if fmt eq 'R8' then val=!values.D_INFINITY else val=!values.F_INFINITY
        else: begin
          if fmt eq 'I2' then val=(str_2_arr(strtrim(data[i],2)))[0]
          if fmt eq 'I4' then val=(str_2_arr(strtrim(data[i],2),/i4))[0]
          if fmt eq 'R4' then val=(str_2_arr(strtrim(data[i],2),/r4))[0]
          if fmt eq 'R8' then val=(str_2_arr(strtrim(data[i],2),/r8))[0]
        end
      endcase
      ;if fmtcode[i] eq 'I2' then val=(str_2_arr(strtrim(data[i],2)))[0]
      ;if fmtcode[i] eq 'I4' then val=(str_2_arr(strtrim(data[i],2),/i4))[0]
      ;if fmtcode[i] eq 'R4' then val=(str_2_arr(strtrim(data[i],2),/r4))[0]
      ;if fmtcode[i] eq 'R8' then val=(str_2_arr(strtrim(data[i],2),/r8))[0]

      if rdbcol gt mcols then message,'BUG: RDBCOL >? MCOLS!!'
      rdbtab.(rdbcol-1L)(kline-1L)=val & rdbcol=rdbcol+1L
	;remember, MCOLS is the number of columns actually being read

    endif				;ICOLS[I]>0)
  endfor				;I=0,NCOLS-1}

  rdnextln:				;the GOTOs come here)})})}
endwhile			;EOF(URDB)}
close,urdb & free_lun,urdb	;close RDB file)

;	make the output table
if n_tags(rdbtab) eq 0 then rdbtab=create_struct('HEAD',header) else $
	rdbtab=create_struct('HEAD',header,rdbtab)

if v ge 4 then help,rdbtab,/str

return,rdbtab
end
