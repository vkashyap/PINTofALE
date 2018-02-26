pro poarmission,poadir=poadir
;+
;procedure	poarmission
;	set permissions appropriately for SCAR and PINTofALE directories
;
;parameters	NONE
;
;keywords
;	poadir	[INPUT] the full path to the PoA directory
;		* if not set, assumes
;		  '/data/fubar/'+['SCAR','PINTofALE']
;		
;
;restrictions
;	highly site specific
;
;history
;	vinay kashyap (Apr2006)
;	changed dir pathnames (Aug2015)
;-

;	print usage in all cases
print,''
print,'************************************************************'
print,'Usage: poarmission,poadir=poadir'
print,'  set correct permissions for PoA directories'
print,'************************************************************'
print,''

;	figure out who I am
spawn,'whoami',username

;dirs='/data/fubar/'+['SCAR','PINTofALE']
dirs=['/fubar/SCAR','/data/astrostat/soft/PINTofALE','/proj/webhead/PINTofALE']
if keyword_set(poadir) then begin
  if size(poadir,/type) eq 7 then dirs=poadir else begin 
    print,poadir
    message,'PoAdir : keyword not understood, assuming default',$
	/informational
  endelse
endif
ndirs=n_elements(dirs)

;	make sure all directories have a+rx,g+w permissions
;	and all files have a+r,g+w permissions
for i=0,ndirs-1 do begin
  cmd='find '+dirs[i]+' -type d ! -name vaporware -user '+username+$
	' -exec chmod a+rx,g+w {} \;'
  cd,dirs[i]
  spawn,'ls '+dirs[i]
  print,cmd & spawn,cmd
  cmd='find '+dirs[i]+' -type f -user '+username+$
	' -exec chmod a+r,g+w {} \;'
  print,cmd & spawn,cmd
  cmd='find '+dirs[i]+' -type f -user '+username+$
	' -exec chgrp poawww {} \;'
  print,cmd & spawn,cmd
endfor

return
end
