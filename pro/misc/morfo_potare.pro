function morfo_potare,skel,minpix,verbose=verbose, _extra=e
;+
;function	morfo_potare
;	prunes the little tributaries from a morphological skeleton
;	of an image and returns a cleaner image of the skeleton.
;
;	this pruning works as follows: first, branch points in the
;	skeleton are identified, and they and their neighboring points
;	are extracted and set aside; next, the remaining points are
;	segmented into distinct and separate lines and then all segments
;	which are smaller than some threshold are discarded; and finally
;	the branch points and their neighborhoods are added back on.
;
;	this is not the standard algorithm for pruning, as far as I know.
;
;syntax
;	pruned_skeleton=morfo_potare(skel,minpix,verbose=verbose)
;
;parameters
;	skel	[INPUT; required] the input 2D bitmap image that
;		contains a morphological skeleton
;		* NOTE: input _must_ be a skeleton, otherwise the
;		  program will fail miserably
;	minpix	[INPUT] minimum number of pixels that must exist
;		in a given branch for it to be kept on
;		* default is 4
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	MORFO_SEGMENTO
;
;history
;	vinay kashyap (Jul2006)
;	added _EXTRA keyword; changed call from MORFO_FILTRI to
;	  MORFO_SEGMENTO (VK; May2007)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(skel) & szs=size(skel)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SKEL is undefined' else $
  if szs[1] lt 4 then ok='SKEL x-axis is collapsed' else $
   if szs[2] lt 4 then ok='SKEL y-axis is collapsed'
if ok ne 'ok' then begin
  print,'Usage: pruned_skeleton=morfo_potare(skel,minpix,verbose=verbose)'
  print,'  prunes tributaries from morphological skeleton'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
mpix=4 & if keyword_set(minpix) then mpix=long(minpix[0])>1

;	first, make sure it is a bitmap
AA = skel gt 0

;	second, find points which have more than 2 neighbors
bb=fltarr(3,3)+1 & bAA=AA*convol(AA,bb,/edge_truncate)
ox=where(bAA gt 3,mox)
if mox eq 0 then begin
  message,'No intersecting points found',/information
  return,skel
endif

;	remove the intersecting points and all its neighbors
xAA=0*AA & xAA[ox]=1 & yAA = AA and 1-xAA

;	find the interconnected regions and the number of pixels in them
zAA=morfo_segmento(yAA,areas,verbose=vv)

;	remove all regions which have less than the required number of pixels
oz=where(areas le mpix,moz)
if moz eq 0 then begin
  message,'Nothing to prune',/informational
  return,skel
endif
for i=0L,moz-1L do begin
  ii=oz[i]
  oo=where(zAA eq ii,moo) & if moo eq 0 then message,'BUG!'
  zAA[oo]=0
endfor

;	convert back to bitmap
yAA = zAA gt 0

;	put the intersecting points back in
yAA[ox]=1

;	close the gaps by dilating along the original skeleton
npix=total(yAA) & nnpix=0
while npix ne nnpix do begin
  if vv gt 0 then kilroy
  nnpix=npix
  img=dilate(yAA,bb) and AA
  npix=total(img)
endwhile

return,img
end
