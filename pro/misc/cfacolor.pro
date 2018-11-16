pro cfacolor,type,cmin=cmin,cmax=cmax,oldr=oldr,oldg=oldg,oldb=oldb, _extra=e
;+
;procedure	cfacolor
;	sets the color scale to one of the new (as of 14 Nov 2018) CfA branded schemes
;
;syntax
;	cfacolor,type,cmin=cmin,cmax=cmax,oldr=oldr,oldg=oldg,oldb=oldb
;
;parameters
;	type	[INPUT; required] specify how the colors are arranged
;		1 = 'CfA' : goes from CfA Red (1) to CfA Violet (255)
;		2 = 'invert' : goes from CfA Violet (1) to CfA Violet (255)
;		3 = 'red' : goes from CfA Red (CMIN) to white (CMAX)
;		4 = 'violet' : goes from CfA Violet (CMIN) to white (CMAX)
;		5 = 'blue' : goes from CfA Dark Blue (CMIN) to white (CMAX)
;
;keywords
;	cmin	[INPUT] minimum of color index to change (default=1; range=0:255)
;	cmax	[INPUT] maximum of color index to change (default=255; range=0:255)
;	oldR	[OUTPUT] old R array
;	oldG	[OUTPUT] old G array
;	oldB	[OUTPUT] old B array
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2018nov)
;-

;	initialize
dncol=!D.N_COLORS
if dncol gt 256 then begin
  message,'You have 24-color depth; setting DECOMPOSED=0',/informational
  device,decomposed=0,retain=2
endif
tvlct,rr,gg,bb,/get & oldr=rr & oldg=gg & oldb=bb
ctyp=['','CfA','invert','red','violet','blue'] & nctyp=n_elements(ctyp)
help=["	1 = 'CfA' : goes from CfA Red (1) to CfA Violet (255)",$
      "	2 = 'invert' : goes from CfA Violet (1) to CfA Violet (255)",$
      "	3 = 'red' : goes from CfA Red (CMIN) to white (CMAX)",$
      "	4 = 'violet' : goes from CfA Violet (CMIN) to white (CMAX)",$
      "	5 = 'blue' : goes from CfA Dark Blue (CMIN) to white (CMAX)"]

;	usage
ok='ok' & np=n_params() & ntyp=n_elements(type) & szt=size(type,/type)
if ntyp eq 0 then ok='Insufficient parameters' else $
 if ntyp gt 1 then ok='can only handle one TYPE at a time' else $
  if ntyp gt 4 and ntyp lt 7 then ok='TYPE must be char string or number'
;
if ok ne 'ok' then begin
  print,'Usage: cfacolor,type,cmin=cmin,cmax=cmax,oldr=oldr,oldg=oldg,oldb=oldb'
  print,'  load in a CfA branded color table'
  print,'  TYPE = '+strjoin(ctyp[1:*],',')
  if np ne 0 then message,ok,/informational
  for i=0L,n_elements(help)-1L do print,help[i]
  return
endif
ityp=long(type[0]) & if ityp ge nctyp then ityp=0

;	keywords
kmin=1b & if keyword_set(cmin) then kmin=(byte(cmin[0])>0b)<255b
kmax=255b & if keyword_set(cmax) then kmax=(byte(cmax[0])>0b)<255b
kk=kmax-kmin+1

;	definitions
redR=141 & redG=0  & redB=52
vioR=43  & vioG=53 & vioB=117
bluR=19  & bluG=26 & bluB=60

;	figure out the type
if szt eq 7 then begin
  cc=strupcase(type[0])
  if strpos(cc,'CFA') ge 0 then ityp=1 else $
   if strpos(cc,'INV') ge 0 then ityp=2 else $
    if strpos(cc,'RED') ge 0 then ityp=3 else $
     if strpos(cc,'VIO') ge 0 then ityp=4 else $
      if strpos(cc,'BLU') ge 0 then ityp=5 else $
       ityp=0
endif

;	make the new color tables
rc=rr & gc=rr & bc=rr
case ityp of
  1: begin
    rc[1:255]=byte(((255.-findgen(255))*redR+findgen(255)*vioR)/255.)
    gc[1:255]=byte(((255.-findgen(255))*redG+findgen(255)*vioG)/255.)
    bc[1:255]=byte(((255.-findgen(255))*redB+findgen(255)*vioB)/255.)
  end
  2: begin
    rc[1:255]=byte(((255.-findgen(255))*vioR+findgen(255)*redR)/255.)
    gc[1:255]=byte(((255.-findgen(255))*vioG+findgen(255)*redG)/255.)
    bc[1:255]=byte(((255.-findgen(255))*vioB+findgen(255)*redB)/255.)
  end
  3: begin
    rc[kmin:kmax]=byte(((kk-findgen(kk))*redR+findgen(kk)*256.)/float(kk))
    gc[kmin:kmax]=byte(((kk-findgen(kk))*redG+findgen(kk)*256.)/float(kk))
    bc[kmin:kmax]=byte(((kk-findgen(kk))*redB+findgen(kk)*256.)/float(kk))
  end
  4: begin
    rc[kmin:kmax]=byte(((kk-findgen(kk))*vioR+findgen(kk)*256)/float(kk))
    gc[kmin:kmax]=byte(((kk-findgen(kk))*vioG+findgen(kk)*256)/float(kk))
    bc[kmin:kmax]=byte(((kk-findgen(kk))*vioB+findgen(kk)*256)/float(kk))
  end
  5: begin
    rc[kmin:kmax]=byte(((kk-findgen(kk))*bluR+findgen(kk)*256)/float(kk))
    gc[kmin:kmax]=byte(((kk-findgen(kk))*bluG+findgen(kk)*256)/float(kk))
    bc[kmin:kmax]=byte(((kk-findgen(kk))*bluB+findgen(kk)*256)/float(kk))
  end
  else: begin
    message,'TYPE = ['+strtrim(ityp,2)+'] not understood',/informational
  end
endcase

;	load in the updated R,G,B
tvlct,rc,gc,bc

return
end
