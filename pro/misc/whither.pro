function whither,idx1,idx2,n12,y1n2=y1n2,n1y2=n1y2,verbose=verbose, _extra=e
;+
;procedure	whither
;	returns the indices that are common to two sets of array indices,
;	assuming that both are subsets of a larger set.
;
;syntax
;	y1y2=whither(idx1,idx2,n12,y1n2=y1n2,n1y2=n1y2,verbose=verbose)
;
;parameters
;	idx1	[INPUT; required] indices of array 1
;	idx2	[INPUT; required] indices of array 2
;	n12	[OUTPUT] number of elements in output
;
;keywords
;	y1n2	[OUTPUT] indices in IDX1 but not in IDX2
;	n1y2	[OUTPUT] indices not in IDX1 but in IDX2
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	a=randomn(seed,1000) & idx1=where(a lt 0) & idx2=where(abs(a) lt 0.3)
;	plot,a,psym=1 & setkolor,'red',1 & setkolor,'green',2
;	oplot,idx1,a[idx1],psym=1,col=1 & oplot,idx2,a[idx2],psym=4,col=2
;	y1y2=whither(idx1,idx2,y1n2=y1n2,n1y2=n1y2) & oplot,y1y2,a[y1y2],psym=1,col=3
;	setkolor,'blue',3 & setkolor,'yellow',4 & setkolor,'orange',5
;	oplot,y1n2,a[y1n2],psym=2,col=4 & oplot,n1y2,a[n1y2],psym=2,col=5
;
;history
;	vinay kashyap (Apr02)
;-

;	default outputs (assume nothing matches and nothing exists)
y1y2=-1L & n12=0L & y1n2=-1L & n1y2=-1L

;	usage
ok='ok' & np=n_params() & n1=n_elements(idx1) & n2=n_elements(idx2)
if np lt 2 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='IDX1 is undefined' else $
  if n2 eq 0 then ok='IDX2 is undefined'
if ok ne 'ok' then begin
  print,'Usage: y1y2=whither(idx1,idx2,n12,y1n2=y1n2,n1y2=n1y2,verbose=verbose)'
  print,'  returns common and non-common elements of 2 sub arrays'
  if np ne 0 then message,ok,/info
  return,y1y2
endif

;	keywords
vv = 0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1

;	default outputs (assume nothing matches)
y1n2=idx1 & n1y2=idx2

;	some rudimentary checks
if idx1[0] lt 0 then begin
  if vv gt 0 then message,'IDX1 is <empty>',/info
  return,y1y2
endif
if idx2[0] lt 0 then begin
  if vv gt 0 then message,'IDX2 is <empty>',/info
  return,y1y2
endif

;	what must the superset be?
mx1=max(idx1)+1L & mx2=max(idx2)+1L
n=mx1 > mx2
ii=bytarr(n)

;	which elements belong where?
ii[idx1]=1 & ii[idx2]=ii[idx2]+2

;	present in IDX1 and IDX2
y1y2=where(ii eq 3,n12)

;	present in IDX1 but not in IDX2
if arg_present(y1n2) then y1n2=where(ii eq 1)

;	present in IDX2 but not in IDX1
if arg_present(n1y2) then n1y2=where(ii eq 2)

return,y1y2
end
