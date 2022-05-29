function modalpoint,array,eps=eps,verbose=verbose, _extra=e
;+
;function modalpoint
;	returns the mode of distribution defined by set of unbinned array values
;
;syntax
;	arraymode=modalpoint(array,eps=eps,verbose=verbose)
;
;parameters
;	array	[INPUT; required] array of values for which mode must be found
;
;keywords
;	eps	[INPUT; default=1e-6] a small number
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	sort the array, divide into two sets around midpoint of range,
;	choose the set that contains more elements, and repeat until
;	the range is sufficiently small, and declare the midpoint of
;	the range to be the mode
;
;warning
;	may fail for multimodal distributions
;
;example
;	for i=0,20 do print,i,modalpoint(randomn(seed,10000L)+i)
;
;history
;	translated to IDL by Vinay Kashyap from C code written for BEHR
;	  by Taeyoung Park c.2003 (MarMMVI)
;	now correctly handles case when input is quantized (VK; SepMMVI)
;	added edge case where if split is even, stops right there (VK; MayMMXXII)
;-

;	usage
ok='ok' & np=n_params() & na=n_elements(array)
if np eq 0 then ok='Insufficient parameters' else $
 if na eq 0 then ok='Input array is undefined' else $
  if na lt 2 then ok='Array must have at least 2 elements'
if ok ne 'ok' then begin
  print,'Usage: arraymode=modalpoint(array,eps=eps,verbose=verbose)'
  print,'  return mode of array'
  if np ne 0 then message,ok,/informational
  return,!values.F_NAN
endif

;	inputs and some special cases
if na lt 3 then return,mean(array)
ok=where(finite(array) ne 0,mok)
if mok eq 0 then return,!values.F_NAN
arr=array[ok] & os=sort(arr) & arr=arr[os]
;
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if keyword_set(eps) then epsilon=double(eps[0]) else epsilon=1d-6

;	step through the array and find the mode
go_on=1
narr=n_elements(arr) & amax=max(arr,min=amin,/nan)
while go_on do begin
  if vv gt 10 then print,strtrim(narr,2)+'.. ',format='($,a)'
  o1=where(arr gt 0.5*(amin+amax),mo1)
  if mo1 eq 0 or mo1 eq narr then begin
    ;message,'BUG?'
    ;not a bug, this means that they are all identical, so quit right here
    return,median(arr)
  endif
  if o1[0] gt narr/2 then tmparr=arr[0:o1[0]-1L] else tmparr=arr[o1]
  if o1[0] eq narr/2 then begin
    if vv gt 0 then message,'evenly split, might be multimodal; cannot deal',/informational
    tmparr=arr
    go_on=0	;quit at this stage
  endif
  if vv gt 100 then print,narr/2,mo1,o1[0],min(tmparr),max(tmparr)
  arr=tmparr
  narr=n_elements(arr) & amax=max(arr,min=amin,/nan)
  if narr eq 1 then go_on=0	;stop when there is only one element
  if amax-amin lt epsilon then go_on=0	;stop when range gets too small
endwhile

return,0.5*(amin+amax)
end
