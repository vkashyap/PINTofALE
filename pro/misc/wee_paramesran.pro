function wee_paramesran,img,firstk=firstk,verbose=verbose,$
	_extra=e
;+
;function	wee_paramesran,img
;	computes and returns the sharpness of an image
;	(see Wee, C.-Y., and Paramesran, R., 2008, ICSP 2008)
;	http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&tp=&arnumber=4697259
;
;parameters
;	img	[INPUT; required] 2D image
;
;keywords
;	firstk	[INPUT] how many of the eigenvalues to include in
;		measure of sharpness
;		* default is to include all of them
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here to prevent program from crashing
;
;history
;	vinay kashyap (2014jul)
;-

;usage
ok='ok' & np=n_params() & ni=n_elements(img) & szi=size(img)
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMG not defined' else $
  if szi[0] ne 2 then ok='IMG must be 2D array' else $
   if szi[1] le 1 then ok='IMG must have more than 1 row' else $
    if szi[2] le 1 then ok='IMG must have more than 1 column'
if ok ne 'ok' then begin
  print,'Usage: sharpness=wee_paramesran(img,verbose=verbose)'
  print,'  computes and returns the sharpness of an image'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
kfirst=szi[1]>szi[2]
if keyword_set(firstk) then kfirst=firstk

;(1) normalize input img by its energy
g=img/sqrt(total(double(img)^2))

;(2) compute mean operator
mu=mean(g,/nan)

;(3) compute mean deviation covariance matrix
gg = g-mu
nn = (szi[1]+szi[2])/2.
sg = (gg # transpose(gg))/(nn-1L)

;() compute singular values matrix
svdc,sg,ww,aa,vv
nw=n_elements(ww)
if kfirst gt nw then kfirst=nw

;(4) compute eigenvalues lambda
lambda = ww^2

;(5) compute trace of eigenvalues
sharpness=total(lambda[0:kfirst-1L])

return,sharpness
end
