;+
;MK_SCAR_HELP.PRO
;	create a HTML help file for all the SCAR routines
;	and make a backup copy and dump it in the home directory
;
;vinay kashyap (Jun98)
;  -- removed "obsolete" from list of directories (VK; 99May)
;  -- makes backup of all routines (VK; 99Nov)
;  -- check whether SCARDIR exists before blithely moving on (VK; MMJan)
;  -- write scar.html to $HOME/www.dmz/scar.html (VK; Jul02)
;-

;	initialize
if not keyword_set(scardir) then scardir='/data/fubar/SCAR/pro/'
tmp=findfile(scardir+'/*',count=count)
if count eq 0 then scardir='/data'+scardir
tmp=findfile(scardir+'/*',count=count)
if count eq 0 then message,'Directory '+scardir+' is missing?'
indir=scardir+$
  ['','external','misc','mkemis','scrypt','solar','specific','esempio','stat','timing']

;outfil='/home/kashyap/www.dmz/scar.html'
outfil=filepath('scar.html',root=getenv('HOME'),subdirectory='www.dmz')
title='<a href="http://hea-www.harvard.edu/PINTofALE/">PINTofALE: Stellar Coronal Analysis Routines</a>'
strict=1
verbose=1
qut=!quiet & !quiet=1

print,'INDIR: ',indir
print,'mk_html_help,indir,"'+outfil+'",title="'+title+$
	'",strict='+strtrim(strict,2)+',verbose='+strtrim(verbose,2)
mk_html_help,indir,outfil,title=title,strict=strict,verbose=verbose
spawn,'chmod a+r '+outfil,tmp
print,"=========================================="
print,outfil

!quiet=qut

;	make a backup copy while at it
;	(designed to only work for me.  -VK)
print,"=========================================="
print,"spawn,'cat /data/fubar/SCAR/notes/tar.com'"
spawn,'cat /data/fubar/SCAR/notes/tar.com'
print,"=========================================="
print,"spawn,'/data/fubar/SCAR/notes/tar.com',tmp"
spawn,'/data/fubar/SCAR/notes/tar.com',tmp

end
