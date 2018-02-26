function roofn,num,base,floor=floor,pwr=pwr, _extra=e
;+
;function	roofn
;	returns the smallest number greater than a specified number
;	that is a power of the specified base.
;
;	output will be of same type as BASE.
;
;syntax
;	nroof=roofn(num,base,floor=floor,pwr=pwr)
;
;parameters
;	num	[INPUT; required] number to "round up"
;	base	[INPUT; default=2] base
;
;keywords
;	floor	[INPUT] if set, returns the largest integer LESS THAN M
;		that is a power of N.
;	pwr	[OUTPUT] the appropriate index [N^(PWR)]
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (April 1995)
;	added keyword _EXTRA (VK; Mar99)
;	cleaned up a bit (VK; Apr99)
;	cleaned up some more (VK; Dec07)
;	completely rewritten to be more flexible and faster (VK; Jan13)
;-

compile_opt idl2

;	usage
ok='ok' & np=n_params()
nm=n_elements(num) & sm=size(num,/type)
nb=n_elements(base) & sb=size(base,/type)
if np lt 2 then ok='insufficient parameters' else $
 if nm eq 0 then ok='NUM is missing' else $
  if (sm ge 6 and sm le 11) then ok='NUM is not a number' else $
   if nb eq 0 then ok='BASE is missing' else $
    if nb gt 1 then ok='BASE cannot be an array' else $
     if (sb ge 6 and sb le 11) then ok='BASE is not a number' else $
      if base[0] le 0 then ok='BASE cannot be .le. 0'
if ok ne 'ok' then begin
  print, 'Usage: nroof=roofn(num,base,pwr=pwr)'
  print,'        nfloor=roofn(num,base,/floor,pwr=pwr)'
  print, '  returns the nearest integer (higher or lower) to a given number'
  print, '  NUM that is a power of a specified base N'
  if np ne 0 then message,ok,/info
  return,0L
endif

;	initialize output
PWR=dblarr(nm) & out=1L

;	first check to see that NUM is not 0 or negative
if nm eq 1 then begin	;(if scalar, break this out to preserve output type
  if num[0] le 0 then $
   if keyword_set(floor) then return,0 else return,1

  pp=alog(num[0])/alog(base[0])

endif else begin	;NM=1)(NM>1
  o0=where(NUM le 0,mo0,complement=o1,ncomplement=mo1)
  if mo0 gt 0 then $
   if keyword_set(floor) then PWR[o0]=0 else PWR[o0]=1

  pp=pwr & if mo1 gt 0 then pp[o1]=alog(num[o1])/alog(base[0])

endelse			;NUM is vector)

if keyword_set(floor) then pwr=long(pp) else pwr=ceil(pp)

return,base[0]^pwr
end

;	this was the original code, left here just for laffs
;mm=abs(m)
;
;;	check parameters
;if n_n eq 0 then begin
;  message,'assuming base 2',/info & n=2L
;endif
;
;mm = long(abs(m))
;nn = long(abs(n))
;
;m0 = 1L & k = 0L & diff = mm-m0
;
;while diff gt 0 do begin
;  k=k+1L & m0 = m0*nn & diff = mm-m0
;endwhile
;pwr = k
;
;if keyword_set(floor) then m0 = m0/nn
;if keyword_set(floor) then pwr = pwr-1
;
;return,m0
;end
