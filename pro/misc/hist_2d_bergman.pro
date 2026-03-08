function hist_2d_bergman,vx,vy,xgrid,ygrid, verbose=verbose, _extra=e
;+
;function	hist_2d_bergman
;	wrapper to hist_2d() to construct an image in grid from a sample of coordinates
;	(because HIST_2D is a bit too finicky)
;
;syntax
;	img=hist_2d_bergman(vx,vy,xgrid,ygrid, verbose=verbose)
;
;parameters
;	VX	[INPUT; required] samples of X-coordinate values
;	VY	[INPUT; required] samples of Y-coordinate values
;	XGRID	[INPUT; required] X-axis grid
;	YGRID	[INPUT; required] Y-axis grid
;		* values of VX,VY that fall outside the range defined by XGRID,YGRID are discarded
;
;keywords
;	VERBOSE	[INPUT; default=0] controls chatter
;
;history
;	Vinay Kashyap (2026-Mar-8)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(VX) & ny=n_elements(VY) & nxg=n_elements(XGRID) & nyg=n_elements(YGRID)
if np lt 4 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='VX is not defined' else $
  if ny eq 0 then ok='VY is not defined' else $
   if nx ne ny then ok='VX and VY are not compatible' else $
    if nxg eq 0 then ok='XGRID is not defined' else $
     if nyg eq 0 then ok='YGRID is not defined' else $
      if nxg eq 1 then ok='XGRID must have at least 2 bounds' else $
       if nyg eq 1 then ok='YGRID must have at least 2 bounds'
if ok ne 'ok' then begin
  print,'Usage: img=hist_2d_bergman(vx,vy,xgrid,ygrid, verbose=verbose)'
  print,'  make an image in grid over specified axis grids'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

;	the range
xmin=min(xgrid,max=xmax)
ymin=min(ygrid,max=ymax)
if xmin eq xmax then ok='X range is 0'
if ymin eq ymax then ok='Y range is 0'
if ok ne 'ok' then begin & message,ok,/informational & return,-1L & endif

;	check the grid order
gxord=1 & if xgrid[0] gt xgrid[nxg-1L] then gxord=-1
gyord=1 & if ygrid[0] gt ygrid[nyg-1L] then gyord=-1
;	pick out the bin lower bounds
xlo=xgrid[0L:nxg-2L] & if gxord lt 0 then xlo=xgrid[1:*]
ylo=ygrid[0L:nyg-2L] & if gyord lt 0 then ylo=ygrid[1:*]

if vv gt 10 then begin
  print,'x,y grid orders',gxord,gyord
  print,'x range',xmin,xmax
  print,'y range',ymin,ymax
endif

ivx=lonarr(nx) & ivy=lonarr(ny)
if gxord eq 1 then begin & kx0=0L & kx1=nxg-1L & dkx=1L & endif else begin & kx0=nxg-1L & kx1=0L & dkx=-1L & endelse
if gyord eq 1 then begin & ky0=0L & ky1=nyg-1L & dky=1L & endif else begin & ky0=nyg-1L & ky1=0L & dky=-1L & endelse

for k=kx0,kx1,dkx do begin & oo=where(vx ge xgrid[k],moo) & if moo gt 0 then ivx[oo]=k & endfor
for k=ky0,ky1,dky do begin & oo=where(vy ge ygrid[k],moo) & if moo gt 0 then ivy[oo]=k & endfor
oo=where(ivx eq kx1,moo) & if moo gt 0 then ivx[oo]=!values.F_NAN
oo=where(ivy eq ky1,moo) & if moo gt 0 then ivy[oo]=!values.F_NAN

if gxord eq 1 then begin & kx0=0L & kx1=nxg-2L & endif else begin & kx0=1 & kx1=nxg-1L & endelse
if gyord eq 1 then begin & ky0=0L & ky1=nyg-2L & endif else begin & ky0=1 & ky1=nyg-1L & endelse
img=hist_2d(ivx,ivy,min1=kx0,min2=ky0,max1=kx1,max2=ky1,bin1=1,bin2=1)

if vv gt 1000 then stop,'halting; type .CON to continue'

return,img
end
