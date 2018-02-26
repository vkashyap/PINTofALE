function lat2arab,latin
;+
;function 	lat2arab
;		return arabic (decimal) equivalent of a latin numeral
;
;syntax
;	arabic=lat2arab(latin)
;
;parameters
;	latin	[INPUT; required] the latin numeral
;
;keywords	NONE
;
;restrictions
;	* spaces, "."s, and anything that is not part of the
;	  roman numeral system (I V X L C D M) is ignored.
;	* there is no reverse function, to go from arabic
;	  to latin, because that is not a unique transformation.
;	  e.g., 1995 could be either MDCCCCLXXXXV (the BBC way),
;	  or MVM, or all steps in between.
;	* guaranteed to work only in simple cases.
;	  DON'T push the envelope!
;	* CAVEATVSVATOR
;
;history
;	vinay kashyap (Nov95)
;-

np=n_params(0)
if np eq 0 then begin
  print,'Usage: arabic=lat2arab(latin)'
  print,'  returns decimal equivalent of latin numeral'
  return,0 ;->
endif

num=strtrim(strupcase(latin),2) & lnum=strlen(num)
larr=strarr(lnum) & narr=lonarr(lnum)

for i=0,lnum-1 do begin
  c1=strmid(num,i,1) & larr(i)=c1
  case c1 of
    'I': narr(i)=1
    'V': narr(i)=5
    'X': narr(i)=10
    'L': narr(i)=50
    'C': narr(i)=100
    'D': narr(i)=500
    'M': narr(i)=1000
    else: narr(i)=0
  endcase
endfor

;initialize
arab=0 & tmp=narr(lnum-1)
for i=lnum-2,0,-1 do begin
  if narr(i) ge narr(i+1) then begin
    arab=arab+narr(i)+tmp & tmp=0		;add up
  endif
  if narr(i) lt narr(i+1) then begin
    tmp=tmp-narr(i)				;accumulate subtractions
  endif
endfor
arab=arab+tmp					;final

return,arab
end
