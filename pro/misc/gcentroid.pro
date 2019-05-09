function gcentroid,img,xsig=xsig,ysig=ysig,itermax=itermax,eps=eps,$
	verbose=verbose, _extra=e
;+
;function	gcentroid
;	Interpolates in pixels to find the location of a peak, assuming
;	that the peak is similar to a Gaussian.  Based on SDSS memo by Lupton &
;	Gunn (1998; https://www.astro.princeton.edu/~rhl/photomisc/centroiding.ps),
;	returns the (i,j) location indices of the peak
;
;syntax
;	ipos=gcentroid(img,xsig=xsig,ysig=ysig,itermax=itermax,eps=eps,verbose=verbose)
;
;parameters
;	img	[INPUT; required] 2D array in which to find the location of the peak
;
;keywords
;	xsig	[OUTPUT] computed width along X axis (not implemented yet)
;	ysig	[OUTPUT] computed width along Y axis (not implemented yet)
;	itermax	[INPUT; default=100] max number of iterations to compute
;		intersection point of quadratic curves
;	eps	[INPUT; default=1e-6] a small number
;		* used here at 10*EPS to decide when iterations have converged
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;notes
;	- to see an example usage case, do
;	  .run gcentroid
;	- slight difference from Lupton & Gunn is that the quadratic solution iterates
;	  until an asymptotic limit is reached, not just once.
;
;history
;	Vinay Kashyap (MMXIX.V)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(img) & szi=size(img)
nx=0L & ny=0L & if szi[0] gt 1 then begin & nx=szi[1] & ny=szi[2] & endif
if np lt 1 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IMG is undefined' else $
  if szi[0] ne 2 then ok='can only handle 2D images' else $
   if nx lt 3 or ny lt 3 then ok='IMG too small'
if ok ne 'ok' then begin
  print,'Usage: pos=gcentroid(img,xsig=xsig,ysig=ysig,itermax=itermax,eps=eps,verbose=verbose)'
  print,'  finds peak location by assuming Gaussian shape and returns it'
  if np ne 0 then message,ok,/informational
  return, -1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
maxiter=100L & if keyword_set(itermax) then maxiter=long(itermax[0])>1L
;
epsilon=1e-6 & if keyword_set(eps) then epsilon=eps[0]

;	first smooth the input image and find the peak
kern=fltarr(3,3)+1./9.
simg=convol(float(img),kern,/edge_wrap)
imgmax=max(simg,imx) & ii=array_indices(img,imx) & i0=ii[0] & j0=ii[1]
pos=[i0,j0]	;base solution
if vv gt 1 then print,'max @ ('+strjoin(strtrim(pos,2),',')+')'

;	make a 3x3 postage stamp around the maximum
ki=1 & kj=1
iA=i0-1L & iB=i0+1L & jA=j0-1L & jB=j0+1L
;	{this block just in case the maximum is flush along an edge
if iA lt 0 then begin & iA=0L & iB=iB+1L & ki=0 & endif
if jA lt 0 then begin & jA=0L & jB=jB+1L & kj=0 & endif
if iB eq nx then begin & iB=iB-1L & iB=nx-1L & ki=2 & endif
if jB eq ny then begin & jB=jB-1L & jB=ny-1L & kj=2 & endif
;	}
i3x3=simg[iA:iB,jA:jB]

;	compute difference measures
sxm=i3x3[0,0]-i3x3[2,0] & d2xm=2.*i3x3[1,0]-i3x3[0,0]-i3x3[2,0]
sx0=i3x3[0,1]-i3x3[2,1] & d2x0=2.*i3x3[1,1]-i3x3[0,1]-i3x3[2,1]
sxp=i3x3[0,2]-i3x3[2,2] & d2xp=2.*i3x3[1,2]-i3x3[0,2]-i3x3[2,2]
sym=i3x3[0,0]-i3x3[0,2] & d2ym=2.*i3x3[0,1]-i3x3[0,0]-i3x3[0,2]
sy0=i3x3[1,0]-i3x3[1,2] & d2y0=2.*i3x3[1,1]-i3x3[1,0]-i3x3[1,2]
syp=i3x3[2,0]-i3x3[2,2] & d2yp=2.*i3x3[2,1]-i3x3[2,0]-i3x3[2,2]
x0m=sxm/2./d2xm
x00=sx0/2./d2x0
x0p=sxp/2./d2xp
y0m=sym/2./d2ym
y00=sy0/2./d2y0
y0p=syp/2./d2yp

;	first order solution
x0=i0 & y0=j0
case ki of
  0: if d2xm gt 0 then x0=x0m+i0 
  1: if d2x0 gt 0 then x0=x00+i0
  2: if d2xp gt 0 then x0=x0p+i0
endcase
case kj of
  0: if d2ym gt 0 then y0=y0m+j0 
  1: if d2y0 gt 0 then y0=y00+j0
  2: if d2yp gt 0 then y0=y0p+j0
endcase
pos=[x0,y0]
if vv gt 1 then print,'first order pos @ ('+strjoin(strtrim(pos,2),',')+')'

;	quadratic solution
if d2xm gt 0 and d2x0 gt 0 and d2xp gt 0 and d2ym gt 0 and d2y0 gt 0 and d2yp gt 0 then begin
  xx=findgen(301)*0.01-1.5 & yy=xx
  xcur0=x0-i0 & ycur0=y0-j0
  ycur=y0-j0 & xcur=x00+ycur*(x0p+x0m)/2.-ycur^2*(2.*x00-x0p-x0m)
  ycur=y00+xcur*(y0p+y0m)/2.-xcur^2*(2.*y00-y0p-y0m)
  dcur0=sqrt((xcur-xcur0)^2+(ycur-ycur0)^2)
  go_on=1 & iter=0L
  while go_on do begin
    xcur0=xcur & ycur0=ycur & iter=iter+1L
    xcur=x00+ycur*(x0p+x0m)/2.-ycur^2*(2.*x00-x0p-x0m)
    ycur=y00+xcur*(y0p+y0m)/2.-xcur^2*(2.*y00-y0p-y0m)
    dcur=sqrt((xcur-xcur0)^2+(ycur-ycur0)^2)
    if dcur0-dcur lt epsilon then go_on=0
    if iter gt maxiter then go_on=0
    dcur0=dcur
    if vv gt 5 then kilroy
    xcur = (xcur>(-3))<3
    ycur = (ycur>(-3))<3
  endwhile
  if vv gt 10 then print,xcur,ycur
  pos=[xcur+i0,ycur+j0]
  if vv gt 1 then print,'quadratic pos @ ('+strjoin(strtrim(pos,2),',')+')'
  if vv gt 5 then print,'iterations = '+strtrim(iter,2)
endif else if vv gt 1 then message,'no quadratic solution found; returning first-order',/informational

if vv gt 1000 then stop,'halting; type .CON to continue'

return,pos
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example usage case
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	syntax
jnk=gcentroid()

;	initialize
if not keyword_set(nsim) then nsim=100000L
if not keyword_set(verbose) then verbose=100

;	make an image of a Gaussian 
rx=randomn(seed,nsim) & ry=randomn(seed,nsim)
img=hist_2d(rx,ry,min1=-3,min2=-3,max1=3,max2=3,bin1=0.1,bin2=0.1)
szi=size(img) & icen=szi[1]/2 & jcen=szi[2]/2
help,img

;	find the peak of the image
pos=gcentroid(img,verbose=verbose)
print,'theoretical location = ('+strtrim(icen,2)+','+strtrim(jcen,2)+')'
print,'Gaussian centroid location = ('+strtrim(pos[0],2)+','+strtrim(pos[1],2)+')'
print,'NOTE: actual location of peak can vary because of randomness in construction'

;	show
contour,smooth(img,3),/xs,/ys,/downhill,levels=max(smooth(img,3))*[0.1,0.5,0.8,0.95,0.99]
plots,pos,psym=4,symsize=4,thick=2
oplot,icen*[1,1],!y.crange,line=1
oplot,!x.crange,jcen*[1,1],line=1

end
