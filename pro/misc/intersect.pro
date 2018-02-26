pro intersect,lis1,lis2,y1y2,y1n2,n1y2,val=val,$
	ny1y2=ny1y2,ny1n2=ny1n2,nn1y2=nn1y2,verbose=verbose,$
	_extra=e
;+
;procedure	intersect
;	procedure to obtain common elements of two numeric arrays
;	sorted in ascending order and with no duplicates
;
;syntax
;	intersect,lis1,lis2,y1y2,y1n2,n1y2,/val,$
;	ny1y2=ny1y2,ny1n2=ny1n2,nn1y2=nn1y2,verbose=verbose
;
;parameters
;	lis1	[INPUT; required] input list #1, array
;	lis2	[INPUT; required] input list #2, array
;	y1y2	[OUTPUT; required] elements in LIS1 also in LIS2
;	y1n2	[OUTPUT] elements in LIS1 that are not in LIS2
;	n1y2	[OUTPUT] elements not in LIS1 that are in LIS2
;
;keywords
;	val	[INPUT] if set, returns the actual values in x1x2 rather
;		than the array position numbers
;	ny1y2	[OUTPUT] number of elements in Y1Y2
;	ny1n2	[OUTPUT] number of elements in Y1N2
;	nn1y2	[OUTPUT] number of elements in N1Y2
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	Y1Y2 contains indices of elements of LIS1 that are common to both LIS1 and LIS2.
;	Y1N2 contains indices of elements of LIS1 that are not in LIS2.
;	N1Y2 contains indices of elements of LIS2 that are not in LIS1.
;	if Y1Y2/Y1N2/N1Y2 is empty, it is set equal to -1L
;	if keyword VAL is set, actual array elements are returned instead of indices,
;	and if Y1Y2/Y1N2/N1Y2 is empty, it is set equal to !values.F_NAN
;
;history
;	algorithm by mala viswanath
;	idl version (but not optimized for it) by vinay kashyap (8/27/1992)
;	added keyword VAL (3/26/93)
;	added keyword VERBOSE (4/14/02)
;	added keywords NY1Y2,NY1N2,NN1Y2, cleaned up edge case outputs, check for sorted arrays (VK; 7/18/15)
;-

;	usage
ok='ok' & np=n_params() & n1=n_elements(lis1) & n2=n_elements(lis2)
if np lt 2 then ok='Insufficient input parameters' else $
 if n1 eq 0 then ok='LIS1 is undefined' else $
  if n2 eq 0 then ok='LIS2 is undefined' else $
   if n1 eq 1 then ok='LIS1 has only one element, use something else' else $
    if n2 eq 1 then ok='LIS2 has only one element, use something else' else $
     if not arg_present(y1y2) then ok='Y1Y2 absent, why bother running this?'
if ok ne 'ok' then begin
  print,'Usage: intersect,lis1,lis2,y1y2,y1n2,n1y2,/val,$'
  print,'       ny1y2=ny1y2,ny1n2=ny1n2,nn1y2=nn1y2,verbose=verbose'
  print,'  returns common and non-common elements of 2 numeric lists sorted in ascending order'
  if np ne 0 then message,ok,/info
  return
endif

;	keywords
getval = 0 & if keyword_set(val) then getval = 1
vv = 0 & if keyword_set(verbose) then vv=long(verbose(0)) > 1

y1y2 = -1L & y1n2 = y1y2 & n1y2 = y1y2
if getval then y1y2 = y1y2*!values.F_NAN
if getval then y1n2 = y1n2*!values.F_NAN
if getval then n1y2 = n1y2*!values.F_NAN

i = 0 & j = 0 & k = 0 & i1 = 0 & i2 = 0

if vv gt 0 then begin
  print,'# elements in LIS1: '+strtrim(n1,2)
  print,'# elements in LIS2: '+strtrim(n2,2)
endif

z1=lis1 & z2=lis2
;	check that the arrays are sorted, and if not, quit with error
dtmp1=z1[1:*]-z1 & dtmp2=z2[1:*]-z2
oy1=where(dtmp1 lt 0,moy1) & oy2=where(dtmp2 lt 0,moy2)
if moy1 gt 0 then begin & message,'LIS1 must be sorted in ascending order',/informational & return & endif
if moy2 gt 0 then begin & message,'LIS2 must be sorted in ascending order',/informational & return & endif
;	(some day, implement local sorting to avoid above problem)
;ii1=lindgen(n1) & ii2=lindgen(n2)
;os1=sort(lis1) & z1=lis1[os1] & jj1=ii1[os1] & ks1=sort(jj1) & k1=(ii1[jj1])[ks1]
;os2=sort(lis2) & z2=lis2[os2] & jj2=ii2[os2] & ks2=sort(jj2) & k2=(ii2[jj2])[ks2]

;first, get y1y2
i = 0 & i1 = 0 & i2 = 0
while i1 lt n1 and i2 lt n2 do begin
  if vv gt 5 and i eq 1000L*long(i/1000) then print,form="($,a)",'.'
  x = z1(i1) - z2(i2)
  case 1 of
    x lt 0: i1=i1+1L
    x gt 0: i2=i2+1L
    else: begin
      if i eq 0 and getval eq 0 then y1y2 = [ i1 ]
      if i eq 0 and getval eq 1 then y1y2 = [ z1(i1) ]
      if i gt 0 and getval eq 0 then y1y2 = [y1y2,i1]
      if i gt 0 and getval eq 1 then y1y2 = [y1y2,z1(i1)]
      i1=i1+1L & i2=i2+1L & i=i+1L
    end
  endcase
endwhile
ny1y2=n_elements(y1y2)
if ny1y2 eq 1 then begin
  if getval and finite(y1y2[0]) eq 0 then ny1y2=0L
  if not getval and y1y2[0] lt 0 then ny1y2=0L
endif
if vv gt 2 then message,'Y1Y2: done',/info

if np eq 3 then return
;now, get y1n2
i=0 & i1 = 0 & i2 = 0 & flag = 1
while flag do begin
  if vv gt 5 and i eq 100L*long(i/100) then print,form="($,a)",strtrim(i,2)
  x = z1(i1) - z2(i2)
  case 1 of
    x eq 0: begin
      i1=i1+1L & i2=i2+1L
    end
    x lt 0: begin
      if i eq 0 and getval eq 0 then y1n2 = [ i1 ]
      if i eq 0 and getval eq 1 then y1n2 = [ z1(i1) ]
      if i gt 0 and getval eq 0 then y1n2 = [y1n2,i1]
      if i gt 0 and getval eq 1 then y1n2 = [y1n2,z1(i1)]
      i1=i1+1L & i=i+1L
    end
    else: i2=i2+1L
  endcase
  if i2 eq n2 and i1 lt n1 then begin
    flag = 0
    if i eq 0 and getval eq 0 then y1n2 = indgen(n1)
    if i eq 0 and getval eq 1 then y1n2 = z1
    if i gt 0 and getval eq 0 then y1n2 = [y1n2,indgen(n1-i1)+i1]
    if i gt 0 and getval eq 1 then y1n2 = [y1n2,z1(i1:*)]
  endif
  if i1 eq n1 then flag = 0
endwhile
ny1n2=n_elements(y1n2)
if ny1n2 eq 1 then begin
  if getval and finite(y1n2[0]) eq 0 then ny1n2=0L
  if not getval and y1n2[0] lt 0 then ny1n2=0L
endif
if vv gt 2 then message,'Y1N2: done',/info

if np eq 4 then return
;finally, get n1y2
i=0 & i1 = 0 & i2 = 0 & flag = 1
while flag do begin
  if vv gt 5 and i eq 100L*long(i/100) then print,form="($,a)",strtrim(i,2)
  x = z2(i2) - z1(i1)
  case 1 of
    x eq 0: begin
      i2=i2+1L & i1=i1+1L
    end
    x lt 0: begin
      if i eq 0 and getval eq 0 then n1y2 = [ i2 ]
      if i eq 0 and getval eq 1 then n1y2 = [ z2(i2) ]
      if i gt 0 and getval eq 0 then n1y2 = [n1y2,i2]
      if i gt 0 and getval eq 1 then n1y2 = [n1y2,z2(i2)]
      i2=i2+1L & i=i+1L
    end
    else: i1=i1+1L
  endcase
  if i1 eq n1 and i2 lt n2 then begin
    flag = 0
    if i eq 0 and getval eq 0 then n1y2 = indgen(n1)
    if i eq 0 and getval eq 1 then n1y2 = z2
    if i gt 0 and getval eq 0 then n1y2 = [n1y2,indgen(n2-i2)+i2]
    if i gt 0 and getval eq 1 then n1y2 = [n1y2,z2(i2:*)]
  endif
  if i2 eq n2 then flag = 0
endwhile
nn1y2=n_elements(n1y2)
if nn1y2 eq 1 then begin
  if getval and finite(n1y2[0]) eq 0 then nn1y2=0L
  if not getval and n1y2[0] lt 0 then nn1y2=0L
endif
if vv gt 2 then message,'N1Y2: done',/info

return
end
