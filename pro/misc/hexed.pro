function hexed,rr,gg,bb,rescale=rescale, _extra=e
;+
;function	hexed
;	converts numbers in the decimal system to hexadecimal strings
;
;syntax
;	hxstr=hexed(r,g,b,rescale=rescale)
;
;parameters
;	rr	[INPUT; required] (rgb) value to convert to hex
;		* can be a scalar or a triple or 2D array of the form [1,3]
;		* if input is byte, assumed to go from 0-255
;		* if input is integer or float, rescaled according to RESCALE
;		* if triple, assumed to be in sequence [r,g,b]
;	gg	[INPUT] (green) value to convert to hex
;		* used only if rr is a scalar
;	bb	[INPUT] (blue) value to convert to hex
;		* used only if rr is a scalar
;
;keywords
;	rescale	[INPUT] if set, and input is not byte, first divides the inputs by
;		this number and then multiplies by 255 and then converts to byte
;		* if input is float or double, default is 1.0
;		* if input is integer, default is 255
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	help,hexed(1,0,0.5,rescale=1)
;	#FF007F
;	print,hexed([1,0,0.5])
;	#FF007F
;	print,hexed(transpose([[255,0,0],[0,255,0]]))
;	#FF0000 #00FF00
;
;history
;	Vinay Kashyap (2015oct)
;	cleaned up; allow multiple values input as (N,3) array (VK; 2015nov)
;-

;	usage
ok='ok' & np=n_params()
nr=n_elements(rr) & ng=n_elements(gg) & nb=n_elements(bb)
szr=size(rr) & narr=nr & typecode=szr[n_elements(szr)-2L]
szg=size(gg) & szb=size(bb)
if np eq 0 then ok='Insufficient parameters' else $
 if nr eq 0 then ok='R is not given' else $
  if nr gt 1 then begin
    if szr[0] eq 1 then begin
      if nr ne 3 then ok='input not a triple -- use either array of size [3] or array of size ['+strtrim(nr,2)+',3]' else narr=1
    endif
    if szr[0] eq 2 then begin
      if szr[2] ne 3 then ok='input not an array of [N,3]' else narr=szr[1]
    endif
    if szr[0] gt 2 then ok='cannot handle 3D or higher arrays'
    if typecode ge 6 and typecode le 11 then ok='input must be byte, integer, or float only'
  endif
if ok ne 'ok' then begin
  print,'Usage: hexstr=hexed(r,g,b,rescale=rescale)
  print,'  converts decimal tuples to hexadecimal strings'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	rescaling factor
if typecode ne 1 then begin
  norm=1.0
  if typecode lt 4 or typecode gt 11 then norm=255.
endif
if n_elements(rescale) ne 0 then norm=rescale[0]+0.0

numarr=1L & hxstr=''
if szr[0] gt 0 then begin & if szr[0] eq 2 then numarr=szr[1] & hxstr=strarr(numarr) & endif

;	step through each row of the input
for i=0L,numarr-1L do begin
  if szr[0] eq 0 then begin
    zr=rr[i] & zg=0. & zb=0.
    if ng gt 0 then zg=gg[i]
    if nb gt 0 then zb=bb[i]
  endif else begin
    if szr[0] eq 1 and szr[1] eq 3 then begin
      zr=rr[0] & zg=rr[1] & zb=rr[2]
    endif
    if szr[0] eq 2 then begin
      zr=rr[i,0] & zg=rr[i,1] & zb=rr[i,2]
    endif
  endelse
  if typecode gt 1 then begin
    zr=byte(255.*(float(zr)/norm))
    zg=byte(255.*(float(zg)/norm))
    zb=byte(255.*(float(zb)/norm))
  endif
  cr=string(zr,form='(Z0)') & if strlen(cr) eq 1 then cr='0'+cr
  cg=string(zg,form='(Z0)') & if strlen(cg) eq 1 then cg='0'+cg
  cb=string(zb,form='(Z0)') & if strlen(cb) eq 1 then cb='0'+cb
  hxstr[i]='#'+cr+cg+cb
endfor

return,hxstr
end
