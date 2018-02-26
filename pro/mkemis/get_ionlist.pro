function get_ionlist,z,dielec=dielec,chidir=chidir, _extra=e
;+
;function	get_ionlist
;	return a list of all the ions for a particular element that
;	CHIANTI has data for
;
;syntax
;	all_ions=get_ionlist(Z,dielec=dielec,chidir=chianti_topdir)
;
;parameters
;	z	[INPUT; required] atomic number
;		* must be a scalar integer!
;
;keywords
;	chidir	[INPUT] path to CHIANTI top directory
;		[default: /data/fubar/SCAR/CHIANTI/dbase]
;	dielec	[OUTPUT] integer array of same size as output, specifying
;		whether given ion refers to dielectronic recombination
;		lines (1) or not (0)
;	_extra	[JUNK] here only to prevent crashing program
;
;requires subroutines
;	SYMB2ZION
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Nov96; based on CHIANTI's SYNTHETIC.PRO)
;	added keyword DIELEC, modified to include CHIANTIv3 dielectronic
;	  recombination directories (VK; OctMM)
;	handle trailing whitespace in masterlist.ions (VK; Jun02)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;-

;	usage
ok='ok' & np=n_params() & nZ=n_elements(Z) & szZ=size(Z) & nszZ=n_elements(szZ)
if np eq 0 then ok='Insufficient parameters' else $
 if nZ gt 1 then ok='Cannot handle array of Z' else $
  if szZ[nszZ-2] eq 7 then ok='Z cannot be a string'
if ok ne 'ok' then begin
  print,'Usage: all_ions=get_ionlist(Z,dielec=dielec,chidir=chianti_topdir)'
  print,'  returns list of all ions with data available in CHIANTI'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;where is CHIANTI?
if not keyword_set(chidir) then chidir='/data/fubar/SCAR/CHIANTI/dbase'

;master list of ions
masterlist = chidir + '/masterlist/masterlist.ions'

;	initialize
dielec=0

openr,umi,masterlist,/get_lun				;open file
;							step through file
while not eof(umi) do begin
  line='' & readf,umi,line 		;read line
  ;	find first whitspace and/or first semi-colon
  cc=strtrim(line,2) & isp=strpos(cc,' ',0) & isc=strpos(cc,';',0)
  ic=isp > isc
  if isp ge 0 then ic=ic < isp
  if isc ge 0 then ic=ic < isc
  if ic ge 0 then line=strmid(cc,0,ic) else line=cc
  ;	now check whether it has a "d" appended
  cc=strmid(line,strlen(line)-1,1) & ide=0
  if strlowcase(cc) eq 'd' then begin
    name=strmid(line,0,strlen(line)-1) & ide=1
  endif else name=line
  symb2zion,name,iz,istate		;decode
  if iz eq z[0] then begin			;store
    if not is_keyword_set(ion) then ion=[istate] else ion=[ion,istate]
    if not is_keyword_set(dielec) then dielec=[ide] else dielec=[dielec,ide]
  endif
endwhile
close,umi,/all & free_lun,umi		;close file

nion=n_elements(ion)
if nion eq 0 then ion=-1L

return,ion
end
