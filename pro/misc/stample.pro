pro stample,package,nuthin=nuthin,noday=noday,notime=notime,nouser=nouser,$
	nopack=nopack,stacol=stacol,stasiz=stasiz,stathk=stathk,help=help,$
	_extra=e
;+
;procedure	stample
;	brand existing plot.  The standard format is to put
;	DATE (YYYY MON DD) TIME (HH:MM:SS) USER (name@host) PACKAGE
;	at the lower-left of a plot.
;
;	stample is a portmanteau (see Jabberwocky) of stamp and staple.
;
;syntax
;	stample,package,/nuthin,/noday,/notime,/nouser,/nopack,$
;	stacol=stacol,stasiz=stasiz,stathk=stathk,/help
;
;parameters
;	package	[INPUT] name of package
;		* default is 'PINTofALE <version>'
;
;keywords
;	nuthin	[INPUT] if set, the same as setting all of
;		NODAY,NOTIME,NOUSER,NOPACK
;	noday	[INPUT] if set, does not include the day
;	notime	[INPUT] if set, does not include the time
;	nouser	[INPUT] if set, does not include the user name
;	nopack	[INPUT] if set, does not include the package name
;	stacol	[INPUT] index of color in color-table to use
;		* default is !D.N_COLORS-10 or 1, whichever is greater
;	stasiz	[INPUT] character size, defaults to 1.
;	stathk	[INPUT] character thickness, defaults to 1.
;	help	[INPUT] if set, prints usage and quits
;	_extra	[JUNK] here only to avoid crashing the program
;
;common
;	cb_stample	{vname, vhost, vPoA}
;
;restrictions
;	requires subroutine SETSYSVAL
;
;history
;	vinay kashyap (DecMM)
;	added common block STAMPLE (VK; JanMMI)
;	added call to SETSYSVAL (VK; FebMMI)
;	back-compatibility fix for SYSTIME (VK FebMMI/A.Maggio)
;	bug corrected -- repeat use of "m" in call to CALDAT
;	  (VK Oct01/A.Maggio)
;	changed name of common block from STAMPLE to CB_STAMPLE
;	  (VK Oct02)
;-

;	usage
if (!D.WINDOW lt 0 and !D.NAME eq 'X') or keyword_set(help) then begin
  print,'Usage: stample,package,/nuthin,/noday,/notime,/nouser,/nopack,$'
  print,'       stacol=stacol,/help'
  print,'  brand the current plot'
  return
endif

;	common block to avoid multiple calls to define NAME@HOST
common cb_stample,vname,vhost,vPoA

;	initialize
PoA_version=2.99202301'		;until we change to FITS inputs and various other updadtes, the version number will stay as 2.99YYYYMM
mon=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
if float(strmid(!version.RELEASE,0,3)) lt 5.3 then begin
  ;	this suggested by Antonio Maggio for back compatibility with
  ;	older versions of IDL.  Not using the same for v5.3, because
  ;	the fewer the extra subroutines needed the better.
  get_juldate,jd		;this is an IDLASTRO procedure
  caldat,jd,month,d,y,h,minute,s
endif else caldat,systime(/julian),month,d,y,h,minute,s
cc=''
if n_elements(package) gt 0 then pack=strtrim(package[0],2) else begin
  if not keyword_set(vPoA) then setsysval,'PoA',vPoA,/getval
  if not keyword_set(vPoA) then vPoA=PoA_version
  pack='PINTofALE '+strtrim(string(vPoA,'(f5.2)'),2)
endelse

kol=!D.N_COLORS-10 > 1
csiz=0.6 & cthk=1.0
if keyword_set(stacol) then kol=long(stacol[0])>1
if keyword_set(stasiz) then csiz=float(stasiz[0])>0
if keyword_set(stathk) then cthk=float(stathk[0])>0

;	no user can be defined if not a unix system
userno=0 & if keyword_set(nouser) then userno=1
if !version.OS_FAMILY ne 'unix' then userno=1

if not keyword_set(noday) then cc=cc+string(y,'(i5)')+' '+$
	string(mon[month-1L])+' '+string(d,'(i2)')+' '
if not keyword_set(notime) then cc=cc+string(h,'(i2)')+':'+$
	string(minute,'(i2)')+':'+string(s,'(i2)')+' '
if not keyword_set(userno) then begin
  if not keyword_set(vname) then begin
    spawn,'whoami',c1
    if c1[0] ne '' then vname=c1[0] else vname=' '
  endif
  if not keyword_set(vhost) then begin
    c2=getenv('HOST')
    if c2[0] ne '' then vhost='@'+c2[0] else vhost=' '
  endif
  ;if n_elements(c1) eq 1 then cc=cc+c1[0]
  ;if c2[0] ne '' then cc=cc+'@'+c2[0]
  cc=cc+vname+vhost+' '
  cc=cc+' '
endif
if not keyword_set(nopack) then cc=cc+pack

xyouts,(!D.ORIGIN)[0],(!D.ORIGIN)[1],cc,/device,color=kol,$
	charsize=csiz,charthick=cthk

return
end
