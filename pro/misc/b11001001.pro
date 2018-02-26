function b11001001,b,str=str,otto=otto
;+
;function	b11001001
;	returns the 0s and 1s of a byte variable
;
;syntax
;	d=b11001001(b,/str,/otto)
;
;parameters
;	b	[INPUT; required] scalar or array of bytes
;
;keywords
;	str	[INPUT] if set, returns the 0s and 1s in a string
;		* overrides OTTO
;	otto	[INPUT] if set, returns an 8xN(B) byte array as output
;
;history
;	vinay kashyap (Jul97)
;-

;	usage
nb=n_elements(b)
if nb eq 0 then begin
  print,'Usage: d=b11001001(b,/str,/otto)'
  print,'  returns 0s and 1s for byte variable'
  return,0
endif

;	initialize
bb=[b]
if keyword_set(str) then d=strarr(nb) else $
  if keyword_set(otto) then d=bytarr(8,nb) else d=lonarr(nb)

for i=8,1,-1 do begin
  bn=bb-2^(i-1)
  op=where(bn ge 0,mop) & om=where(bn lt 0,mom)
  if mop gt 0 then if keyword_set(str) then d(op)=d(op)+'1' else $
    if keyword_set(otto) then d(i-1,op)=byte(1) else d(op)=d(op)+10L^(i-1)
  if mom gt 0 and keyword_set(str) then d(om)=d(om)+'0'
  if mop gt 0 then bb(op)=bn(op)
endfor

return,d
end
