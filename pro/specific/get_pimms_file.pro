function get_pimms_file,mission,instr,config,special=special,pdir=pdir
;+
;function	get_pimms_file
;	returns a string containing the name of the file containing the
;	effective areas for the given instrument.
;
;usage
;	filnam=get_pimms_file(mission,instr,config,special=special,pdir=pdir)
;
;parameters
;	mission	[INPUT; required] e.g., 'ROSAT', 'AXAF', etc.
;	instr	[INPUT] instrument (e.g., 'PSPC', 'ACIS', etc.)
;	config	[INPUT] instrument configuration (e.g., 'OPEN', 'BI', etc.)
;
;keywords
;	special	[INPUT] any special requests
;		* e.g., if MISSION='rosat', INSTR='pspc', and CONFIG='R4'
;		and SPECIAL='R7', the output will be 'rosat_pspc_r4tor7.area'
;	pdir	[INPUT; default='/data/fubar/SCAR/ardb/pimms/data/'] directory
;		containing the PIMMS data files
;
;history
;	vinay kashyap (Mar98)
;	changed PDIR default from /soft/prop_cli/config/pimms/data to
;	  /soft/pimms/data (VK; Feb03)
;	changed PDIR default from /soft/pimms/data to
;	  /soft/ciao/config/pimms/data/ (VK; Mar06)
;	changed PDIR default from /soft/ciao/config/pimms/data to
;	  !ARDB/pimms/data/ (VK; Mar13)
;	added explicit call to GETPOADEF() to get !ARDB (VK; May17)
;-

;	usage
nm=n_elements(mission) & ni=n_elements(instr)
if nm eq 0 then begin
  print,'Usage: filname=get_pimms_file(mission,instr,config,special=special,pdir=pdir)'
  print,'  returns PIMMS data file for specified instrument'
  return,'file not found'
endif

;	read in parameters...
tele=strlowcase(mission(0)) & if ni ne 0 then inst=strlowcase(instr(0))
if n_elements(config) ne 0 then cfg=strlowcase(config(0)) else cfg=''
;
;	...and keywords
if keyword_set(special) then spl=strlowcase(special) else spl=''
if keyword_set(pdir) then dirp=pdir else begin
  ardb=getpoadef('ARDB')
  ;dirp='/fubar/SCAR/ardb/pimms/data/'
  ;	check whether !ARDB is defined and look in there
  defsysv,'!ARDB',exists=iardb
  if iardb eq 1 then dirp=filepath('',root_dir=!ARDB,subdir=['pimms','data'])
endelse

;	help
if ni eq 0 then begin
  fils=findfile(dirp+'/'+tele+'*',count=nfils)
  if nfils eq 0 then begin
    message,'MISSION '+tele+' not understood',/info
  endif else begin
    print,'available files for this mission are:' & print,''
    for i=0,nfils-1 do print,fils(i)
    print,'' & print,'please specify the appropriate instrument!'
  endelse
  return,'file not found'
endif

;	first, all the files
fils=findfile(dirp+'/*',count=nfils)
if nfils eq 0 then begin
  message,'no files found in this directory!',/info
  return,'file not found'
endif

;	which of the files match the criteria?
ifils=lonarr(nfils)
for i=0,nfils-1 do begin
  file=strlowcase(fils(i)) & k=0
  j1=strpos(file,tele,0) & if j1 ge 0 then k=1		;matches MISSION
  j2=strpos(file,inst,(j1>0)) & if j2 lt 0 then k=0	;& matches INSTR
  if cfg ne '' then begin
    j3=strpos(file,cfg,(j2>0)) & if j3 lt 0 then k=0	;& matches CONFIG
  endif
  if spl ne '' then begin
    j4=strpos(file,spl,(j2>0)) & if j4 lt 0 then k=0	;& matches SPECIAL
  endif
  ifils(i)=k
endfor

;	identify the matched file
oo=where(ifils eq 1,moo)
case moo of
  0: arfile='file not found'
  1: arfile=fils(oo(0))
  else: begin					;{find the "nearest"
    del=fltarr(moo) & tmplt=tele+inst+cfg+spl
    ktmplt=total(float(byte(tmplt)))
    for i=0,moo-1 do begin			;{check each candidate file
      fil=fils(oo(i))
      kfil=total(float(byte(fil)))
      del(i)=abs(ktmplt-kfil)
    endfor					;I=0,MOO-1}
    tmp=min(del,imn) & arfile=fils(oo(imn))
  end						;found the nearest}
endcase

return,arfile
end
