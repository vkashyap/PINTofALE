pro ogipzip,arr,carr,ngrp,fchan,nchan,eps=eps,chan0=chan0,cchan=cchan,$
	_extra=e
;+
;procedure	ogipzip
;	compress a numerical array to return only those elements
;	which are above some threshold, and the means to reconstruct
;	the full fledged array.  basically this is how OGIP RMF matrix
;	values are compressed and stored.
;
;syntax
;	ogipzip,arr,carr,ngrp,fchan,nchan,eps=eps,chan0=chan0,cchan=cchan
;
;parameters
;	arr	[INPUT; required] input 1D array to be compressed
;	carr	[OUTPUT; required] compressed array
;	ngrp	[OUTPUT; required] number of groups into which ARR is broken
;	fchan	[OUTPUT; required] indices of first elements of each group
;	nchan	[OUTPUT; required] number of elements in each group
;
;keywords
;	eps	[INPUT] a small number, below which to ignore all values
;		* default is 1e-6
;	chan0	[INPUT] index number of first channel, to determine
;		whether FCHAN should be 0-based or 1-based
;		* default is 1-based, because it appears that that is
;		  how Chandra RMFs are being made.  set CHAN0 explicitly
;		  to 0 to make it 0-based.
;	cchan	[OUTPUT] channels corresponding to CARR
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Nov2001)
;-

;	usage
ok='ok' & np=n_params() & na=n_elements(arr)
if np lt 5 then ok='Insufficient parameters' else $
 if na eq 0 then ok='ARR is undefined'
if ok ne 'ok' then begin
  print,'Usage: ogipzip,arr,carr,ngrp,fchan,nchan,eps=eps,chan0=chan0'
  print,'  OGIP RMF style array compression'
  if np ne 0 then message,ok,/info
  return
endif

;	recast input
a=[arr[*]]

;	check keywords
thr=1e-6 & if keyword_set(eps) then thr=eps[0]
firstchan=1 & if n_elements(chan0) ne 0 then firstchan=long(chan0[0])

;	compress away..
oo=where(a gt thr,moo)
cchan=oo+firstchan

;	special cases
if moo eq 0 then begin
  carr=-1L & ngrp=0 & fchan=0+firstchan & nchan=0 & return
endif
if moo eq na then begin
  carr=arr & ngrp=1 & fchan=0+firstchan & nchan=na & return
endif

;	non-special case
carr=arr[oo] & ngrp=1 & fchan=[oo[0]+firstchan] & nchan=[1L]
for i=1L,moo-1L do begin
  if oo[i]-oo[i-1] eq 1 then nchan[ngrp-1L]=nchan[ngrp-1L]+1 else begin
    ngrp=ngrp+1 & fchan=[fchan,oo[i]+firstchan] & nchan=[nchan,1]
  endelse
endfor

return
end
