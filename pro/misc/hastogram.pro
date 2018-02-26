function hastogram,list,x,wts=wts,verbose=verbose, _extra=e
;+
;function	hastogram
;	returns a frequency histogram over a specified grid,
;	calculated in a fast and clever way.  really.
;	works in haste, without waste.
;
;	works best when the grid is sparsely populated.
;	if there are more items than bins, use IDL's
;	built-in HISTOGRAM() function instead.
;
;syntax
;	f=hastogram(list,x,wts=wts,verbose=verbose)
;
;parameters
;	list	[INPUT; required] list of numbers to bin
;	x	[I/O] the required binning scheme
;		* if not set, assumes a linear grid of size 100 from
;		  min(LIST) to max(LIST)
;		* if scalar or 1-element vector, assumes this to be the
;		  number of bins
;		  * if -ve, log grid, else linear
;		* if vector with more than 1 element, assumes it to be
;		  the bin boundaries.
;		* WARNING: if input as scalar, will overwrite with the
;		  adopted grip on output
;
;keywords
;	wts	[INPUT] if given, appropriately weights each element of LIST
;		* default is unity
;		* if size does not match size(LIST), will be appropriately
;		  interpolated
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	the problem is the following: if the number of bins is large,
;	accumulating a frequency histogram by linear search takes too
;	long, esp. in IDL if one uses for-loops.
;	so, first create a monster array containing both the list elements
;	and the grid values.  then sort this array.  this results in an
;	array where the list elements are all nicely placed within the
;	correct bins.  now, if we've been keeping track of the positions
;	of the list elements, it's an easy job to count up the number in
;	each bin of the grid.  to do the latter, we first create a new
;	array made up of -1s in positions of list elements and position
;	indices for grid values and reorder according to the above sort.
;	then, replace each -1 by the nearest non-(-1) from the left.  now
;	each list element is assigned the correct bin number.  voila!
;
;restrictions
;	X must be sorted in increasing order.
;	(no longer, as of Jun02)
;
;subroutines
;	KILROY
;
;history
;	vinay kashyap (Jan98)
;	added kludge to speed up in case max(X) < max(LIST) (VK; Sep98)
;	added quit in case X is not sorted in ascending order (VK; JanMMI)
;	added keyword VERBOSE, and removed restriction on order of X
;	  (VK; Jun02)
;	changed ints to longs (VK; Aug08)
;	converted bracket notation to IDL5+ and fixed bug with inputs being
;	  array of size [1,N] (VK;Mar10)
;-

;	usage
nn=n_elements(list)
if nn eq 0 then begin
  print,'Usage: f=hastogram(list,x,wts=wts,verbose=verbose)'
  print,'  returns frequency histogram accumulated over irregular grid'
  return,-1L
endif

;	verbosity
vv=0L & if keyword_set(verbose) then vv=long(verbose[0]) > 1

;	decode X
nx=n_elements(x)
case nx of					;{decode X
  0: begin				;(use hardcoded default
     xmin=min(list,max=xmax) & nbin=100L & dx=float(xmax-xmin)/float(nbin)
     xx=[xmin+findgen(nbin)*dx,xmax]
  end					;NX=0)
  1: begin				;(number of bins given
     nbin=x(0)
     if nbin eq 0 then nbin=100L
     if nbin lt 0 then begin		;(log grid
       xmin=min(alog10(abs(list)),max=xmax)
       dx=float(xmax-xmin)/float(nbin)
       xx=[xmin+findgen(nbin)*dx,xmax] & xx=10.^(xx)
     endif else begin			;)(linear grid
       xmin=min(list,max=xmax)
       dx=float(xmax-xmin)/float(nbin)
       xx=[xmin+findgen(nbin)*dx,xmax]
     endelse				;)
  end					;NX=1)
  else: begin
    xx=x[*]				;all is well
    ;	sort if required
    if total((xx[1:*]-xx)<0) ne 0 then begin
      message,'resorting X to be in ascending order',/info
      xx=x(sort(x))
    endif
  end
endcase						;X}
nx=n_elements(xx)

if nn gt 2L*nx and vv gt 0 then begin
  message,"There are more list elements than bins; consider",/informational
  message,"using IDL's built-in HISTOGRAM() function instead",/informational
endif

;;	check to see that XX are sorted in increasing order
;if xx[1] lt xx[0] then begin
;  message,'binning grid must be sorted in increasing order',/info
;  return,-1L
;endif

;	poundage
nw=n_elements(wts) & ww=lonarr(nn)+1		;default=1
if nw eq 1 then ww[*]=wts[0]			;all set to same
if nw gt 1 and nw ne nn then $			;interpolate over range
	ww=interpol(wts,findgen(nw),findgen(nn)*(nw-1.)/(nn-1.))
if nw eq nn then ww=wts[*]			;specified

;	declare the output
h=lonarr(nx-1) +$			;regular frequency histogram
	0*ww[0]				;changed to appropriate type

;	find the indices
allX=[xx[0],list,xx[1:*]]	;the large array containing LIST and grid
	;NOTE: XX is split so that min(LIST)=xx[0] is not lost in the sort
oX=sort(allX)		;sort above to move list elements into correct bins
allI=[0L,lonarr(nn)-1L,lindgen(nx-1L)+1L]	;the position remembering array
allW=[0L,ww,lonarr(nx-1L)]		;and the weight remembering array
allI=allI[oX]			;positions have been reordered per above sort
allW=allW[oX]				;same with the weights
    ;kludge in case max(X) < max(LIST)
    ugh=where(allI ge 0,mugh)
    if ugh[mugh-1L] lt nn+nx-1 then allI[ugh[mugh-1L]+1:*]=allI[ugh[mugh-1L]]
oI=where(allI lt 0,moI)			;remember where the LIST elements are
oJ=oI					;copy of same, will be modified
while moI gt 0 do begin		;{while there are list elements w/o bin indices
  if vv gt 1 then kilroy; was here
  allI[oJ]=allI[oJ-1L]			;assign bin value
  moJ=moI & oJ=where(allI lt 0,moI)	;any left?
  if moI eq moJ then moI=0		;in case XMIN > min(LIST)
endwhile			;MOI>0}

;	extract the indices
oo=allI[oI]			;these are the bin indices of LIST elements
f=allW[oI]			;extracting the weights in the right order
;reverse_indices=	;mapping back each LIST element to a bin

;	now make the frequency histogram
;	we expect to use this routine when the number of bins are much
;	larger than the number of LIST elements, i.e., when very few of
;	the bins are actually populated.  so minimal qualms about stepping
;	through only the populated bins using a (ugh) for-loop
oq=oo[uniq(oo,sort(oo))] & noq=n_elements(oq)	;all the populated bins
if vv gt 1 then begin
  if 2*noq lt n_elements(list) then begin
    message,'Number of bins < number of list elements',/informational
    message,'expect slow down; consider using HISTOGRAM()',/informational
  endif
endif
for i=0L,noq-1L do begin		;{for each populated bin
  k=oq[i] & ok=where(oo eq k)		;find the LIST elements
  if k ge 0 and k lt nx-1 then h[k]=total(f[ok])	;and approp. weights
endfor					;I=0,NOQ-1}

;	outputs
x=xx
return,h

end
