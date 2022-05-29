function rd_gaiatsv,tsvfil,sep=sep,cols=cols,fmts=fmts,verbose=verbose,_extra=e
;+
;function	rd_gaiatsv
;	read in the Gaia DR2 data extracted into a local file with find_gaia_dr2.py
;	into a structure and return the structure
;	NOTE: find_gaia_dr2.py is part of python-cdsclient, get it from
;	http://cds.u-strasbg.fr/resources/doku.php?id=cdsclient
;
;syntax
;	gaiastr=rd_gaiatsv(tsvfil,sep=sep,cols=cols,fmts=fmts,verbose=verbose)
;
;parameters
;	tsvfil	[INPUT; required] name of input tab-separated-values file
;
;keywords
;	sep	[INPUT] 1-character separator (default is string(byte(9)), signifying a tab)
;	cols	[OUTPUT] column namesathe following column name changes are performed:
;		_r -> off_r
;		E(BR/RP) -> excess_bprp
;		BP-RP -> color_bprp
;		BP-G -> color_bpG
;		G-RP -> color_Grp
;		[Fe/H]temp -> FeH_template
;		E(BP-RP) -> reddening_bprp
;		b_E(BP-RP) -> lo_reddening_bprp
;		B_E(BP-RP) -> hi_reddening_bprp
;		b_Teff -> lo_Teff
;		B_Teff -> hi_Teff
;		b_AG -> lo_AG
;		B_AG -> hi_AG
;		b_Rad -> lo_Rad
;		B_Rad -> hi_Rad
;		b_Lum -> lo_Lum
;		B_Lum -> hi_Lum
;	fmts	[OUTPUTS] nominal formats for columns
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	find_gaia_dr2.py must have been run and its output redirected into a local file
;	no checks are performed to verify that the file being read in is proper
;	the following column name changes are performed:
;	requires subroutines
;		KILROY
;
;example
;	find_gaia_dr2.py --format tsv -m 1692919135 -r 30 "22 38 45.5747 -20 37 16.081" > gaia_FKAqr.tsv
;	gaiastr=rd_gaiatsv('gaia_FKAqr.tsv')
;
;history
;	Vinay Kashyap (Apr2018)
;	corrected for when first column had different number of dashes at beginning of
;	  data block (VK; Aug2018)
;-

;	check input
ok='ok' & np=n_params() & ntf=n_elements(tsvfil) & sztf=size(tsvfil,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if ntf eq 0 then ok='TSVFIL is undefined' else $
  if ntf gt 1 then ok='can only handle one file at a time' else $
   if sztf ne 7 then ok='TSVFIL must be a character string'
if ok ne 'ok' then begin
  print,'Usage: gaiastr=rd_gaiatsv(tsvfil,sep=sep,verbose=verbose)'
  print,'  read in the output of find_gaia_dr2.py and return as a structure'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
ss=string(byte(9b)) & if n_elements(sep) gt 0 then ss=sep[0] & nss=strlen(ss)
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

;	check that file exists
fil=file_search(tsvfil[0],count=nfil)
if nfil eq 0 then begin
  message,tsvfil[0]+': file not found',/informational & return,-1L
endif
if nfil gt 1 then begin
  print,fil
  message,tsvfil[0]+': too many files found',/informational & return,-1L
endif
gaiatsv=fil[0]

;	how many entries in the file?
spawn,'grep -v ^# '+gaiatsv+' | wc -l',tmp & nrow=long(tmp[0])-3L	;column names row, units row, and dashes row
;	less any blank lines
spawn,'grep ^$ '+gaiatsv+' | wc -l',tmp & nrow=nrow-long(tmp[0])

;	column names
spawn,"grep ^#Column "+gaiatsv+" | awk '{print $2}'",tmp & cols=strtrim(tmp,2)
ncol=n_elements(cols)
;	special cases
ok=where(cols eq '_r',mok) & if mok gt 0 then cols[ok[0]]='off_r'
ok=where(cols eq 'E(BR/RP)',mok) & if mok gt 0 then cols[ok[0]]='excess_bprp'
ok=where(cols eq 'BP-RP',mok) & if mok gt 0 then cols[ok[0]]='color_bprp'
ok=where(cols eq 'BP-G',mok) & if mok gt 0 then cols[ok[0]]='color_bpG'
ok=where(cols eq 'G-RP',mok) & if mok gt 0 then cols[ok[0]]='color_Grp'
ok=where(cols eq '[Fe/H]temp',mok) & if mok gt 0 then cols[ok[0]]='FeH_template'
ok=where(cols eq 'E(BP-RP)',mok) & if mok gt 0 then cols[ok[0]]='reddening_bprp'
ok=where(cols eq 'b_E(BP-RP)',mok) & if mok gt 0 then cols[ok[0]]='lo_reddening_bprp'
ok=where(cols eq 'B_E(BP-RP)',mok) & if mok gt 0 then cols[ok[0]]='hi_reddening_bprp'
ok=where(cols eq 'b_Teff',mok) & if mok gt 0 then cols[ok[0]]='lo_Teff'
ok=where(cols eq 'B_Teff',mok) & if mok gt 0 then cols[ok[0]]='hi_Teff'
ok=where(cols eq 'b_AG',mok) & if mok gt 0 then cols[ok[0]]='lo_AG'
ok=where(cols eq 'B_AG',mok) & if mok gt 0 then cols[ok[0]]='hi_AG'
ok=where(cols eq 'b_Rad',mok) & if mok gt 0 then cols[ok[0]]='lo_Rad'
ok=where(cols eq 'B_Rad',mok) & if mok gt 0 then cols[ok[0]]='hi_Rad'
ok=where(cols eq 'b_Lum',mok) & if mok gt 0 then cols[ok[0]]='lo_Lum'
ok=where(cols eq 'B_Lum',mok) & if mok gt 0 then cols[ok[0]]='hi_Lum'
;	extra special cases
lcolmax=long(alog10(ncol)+1) & iform="(i"+strtrim(lcolmax,2)+'.'+strtrim(lcolmax,2)+")"
for i=0L,ncol-1L do begin
  isvar=execute(cols[i]+'=0')
  if not isvar then cols[i]='col'+string(i,iform)
endfor

;	column formats
spawn,"grep ^#Column "+gaiatsv+" | awk '{print $3}' | sed 's,(,,' | sed 's,),,'",tmp & fmts=strtrim(tmp,2)
ok=where(fmts eq 'D12.5') & if mok gt 0 then fmts[ok[0]]='F12.5'
;
;	initialize outputs
cc='' & ccfmt=strupcase(strmid(fmts,0,1))
for i=0L,ncol-1L do begin
  cc=cols[i]
  if ccfmt[i] eq 'F' or ccfmt[i] eq 'E' or ccfmt[i] eq 'D' then cc=cc+'=dblarr('+strtrim(nrow,2)+')'
  if ccfmt[i] eq 'A' then cc=cc+'=strarr('+strtrim(nrow,2)+')'
  if ccfmt[i] eq 'I' then cc=cc+'=lonarr('+strtrim(nrow,2)+')'
  jnk=execute(cc)
endfor
cc='' & cfmt=strjoin('var'+string(lindgen(ncol),'(i3.3)'),',')+",form='("+strjoin(fmts,',a1')+")'"

header=[gaiatsv]

;	read in
datablock=0 & krow=0L
openr,ug,gaiatsv,/get_lun
while not eof(ug) do begin
  if keyword_set(datablock) then begin
    readf,ug,cc
    if strtrim(cc) ne '' then begin
      cc2=strtrim(strsplit(cc,ss,/extract,/preserve_null),2)
      for j=0L,ncol-1L do begin
        if cc2[j] ne '' and cc2[j] ne 'NOT_AVAILABLE' then begin
          if ccfmt[j] eq 'F' or ccfmt[j] eq 'E' or ccfmt[j] eq 'D' then jnk=execute(cols[j]+"["+strtrim(krow,2)+"]=double("+cc2[j]+")")
          if ccfmt[j] eq 'I' then jnk=execute(cols[j]+"["+strtrim(krow,2)+"]=long("+cc2[j]+")")
          if ccfmt[j] eq 'A' then jnk=execute(cols[j]+"["+strtrim(krow,2)+"]='"+cc2[j]+"'")
        endif
      endfor
      krow=krow+1L
      if krow eq long(krow/100)*100 then kilroy,dot=strtrim(nrow-krow,2)+' '
      ;jnk=execute("readf,ug,"+cfmt)
      ;for j=0L,ncol-1L do jnk=execute(cols[j]+"[i]=var"+string(j,'(i3.3)'))
    endif
  endif else begin
    readf,ug,cc
    print,cc
    header=[header,cc]
    if strpos(cc,'--------',0) ge 0 then datablock=1
    ;if datablock ne 0 then stop,datablock
  endelse
endwhile
close,ug & free_lun,ug

ccstr="gaiastr=create_struct('header',header,"+strjoin("'"+cols+"',"+cols,',')+')'
header=[header,ccstr]
jnk=execute(ccstr)

if vv gt 1000 then stop,'halting; type .CON to continue'

return,gaiastr
end
