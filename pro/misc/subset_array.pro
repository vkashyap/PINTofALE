function subset_array,arr,kd,idx,verbose=verbose,_extra=e
;+
;function	subset_array
;	extract a subarray from a multi-dimensional array that includes
;	only the elements at specified indices at the specified dimension
;
;syntax
;	subarray=subset_array(arr,kd,idx,verbose=verbose)
;
;parameters
;	arr	[INPUT; required] a multi-dimensional array from which to
;		extract a subset
;	kd	[INPUT; required] the dimension/axis on which to filter
;	idx	[INPUT; required] indices along dimension KD to pick out
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;warning
;	because IDL throws away trailing dimensions of size 1, permanently
;	reducing the dimensionality of the array, the input cannot have any
;	columns with size 1.
;
;example
;	.run subset_array
;
;history
;	Vinay Kashyap (2019feb)
;	  Michael Stevens and Kathy Reeves posed this question on 2018-jul-23,
;	  to extract a subset of an array sliced at a given dimension without
;	  using EXECUTE().  It always seemed like it has to be possible because
;	  array indices are well behaved.  The trick is to note that it is easy
;	  to extract the indices for the first dimension.  So we can transpose
;	  the N-dimensional matrix such that the dimension of interest is now
;	  the first one in the list, do the extraction, then transpose back.
;	  Et voila.
;	added warning about columns of size 1 (VK; 2019mar)
;-

;	usage
ok='ok' & np=n_params() & na=n_elements(arr) & sza=size(arr) & nd=n_elements(kd) & ni=n_elements(idx)
if np lt 3 then ok='Insufficient parameters' else $
 if na eq 0 then ok='ARR is not defined' else $
  if nd eq 0 then ok='KD: Dimension slice is not defined' else $
   if ni eq 0 then ok='IDX: subset indices are not defined' else $
    if na eq 1 then ok='ARR is not an array' else $
     if nd gt 1 then ok='can only handle one slice at a time' else $
      if kd[0] gt sza[0]-1 then ok='KD>size(ARR): asking for non-existant dimensional slices' else $
       if min(idx) lt 0 then ok='IDX has -ves: does not appear to be indices' else $
        if max(idx) gt sza[kd+1L]-1L then ok='IDX overflows size of input ARR'
if ok ne 'ok' then begin
  print,'Usage: subarray=subset_array(array,kd,idx,verbose=verbose)'
  print,'  extract subarray that includes specified indices from specified dimensional slice'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

zz=long(0*arr)-1L	;first make a holding index array

dims=sza[1:sza[0]]	;what are the dimensions we have to deal with?
o1=where(dims eq 1,mo1) & if mo1 ne 0 then message,$
	'CAVEAT VSVATOR: some columns have size 1.  They may inadvertantly get compressed.',/informational
newdims=dims & newdims[kd]=ni	;what are the output dimensions?
idims=lindgen(sza[0])	;just the original sequence of dimensions
kdims=shift(idims,-kd)	;shift like this such that the KDth dimension becomes the first one
arr2=transpose(arr,kdims)	;do the shift
sz2=size(arr2) & nn=1L & for i=1L,sza[0]-1L do nn=nn*sz2[i+1L]	;skipping size
for i=0L,ni-1L do zz[idx[i]+lindgen(nn)*sz2[1]]=idx[i]+lindgen(nn)*sz2[1]	;fill in the good indices in ZZ
ok=where(zz gt -1L,mok)	;keep these parts of ZZ
if mok gt 0 then begin
  arr2=arr2[ok]		;filter out the unneeded elements
endif else begin
  message,'BUG! no elements selected',/informational
  return,-1L
endelse
newarr=make_array(shift(newdims,-kd),value=arr[0])	;make a new placeholder array of the right (shifted) dimensions
newarr[0]=arr2						;and push everything into it

subarr=transpose(newarr,sort(kdims))			;shift back to original order, this is the output

if vv gt 1000 then stop,'halting; type .CON to continue'

return,subarr
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	example usage
subarr=subset_array()

;	try this case
if not keyword_set(verbose) then verbose=10

print,'from [6,4,5,3] extract [6,4,[0,3],3]'
arr=360-indgen(6,4,5,3) & kd=2 & idx=[0,3]
arr=reform(randomn(seed,360),[6,4,5,3]) & kd=2 & idx=[0,3]

;	what should we extract?
print,'arr[*,*,0,*] ==>'
print,arr[*,*,0,*]
print,'arr[*,*,3,*] ==>'
print,arr[*,*,3,*]
;print,'ARR=indgen(6,4,5,3) & KD=2 & IDX=[0,3]'
print,'ARR=reform(randomn(seed,360),[6,4,5,3]) & KD=2 & IDX=[0,3]'

;	extrac subarray
subarr=subset_array(arr,kd,idx,verbose=verbose)

;	print output, compare above
print,'subarr[*,*,0,*] (compare arr[*,*,0,*]) ==>'
print,subarr[*,*,0,*]
print,'subarr[*,*,1,*] (compare arr[*,*,3,*]) ==>'
print,subarr[*,*,1,*]

end
