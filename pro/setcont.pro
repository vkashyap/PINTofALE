function setcont,x,y,ysig,ycont=ycont,xcont=xcont,const=const,$
	xrange=xrange,yrange=yrange,lock=lock,ccol=ccol, _extra=e
;+
;function	setcont
;	interactively determine a continuum
;
;syntax
;	cont=setcont(x,y,ysig,ycont=ycont,xcont=xcont,/const,$
;	xrange=xrange,yrange=yrange,ccol=ccol,/lock)
;
;parameters
;	x	[INPUT; required] abscissa values
;		* must be an array with at least 2 elements
;	y	[INPUT; required] Y(X)
;		* length must match X
;	ysig	[OUTPUT] errors on returned continuum (see description below)
;
;keywords
;	ycont	[I/O] if set on input, begins work from here
;		* overwritten by new values on output
;	xcont	[I/O] if N(YCONT) .NE. N(X), then look here for the
;		actual abscissae values.. YCONT=YCONT(XCONT)
;		* if N(XCONT) .NE. N(YCONT), both are ignored
;	const	[I/O] if set, returns piecewise >constant< continum,
;		as opposed to piecewise >linear< curve
;		* can be toggled by clicking outside the plot area
;	xrange	[INPUT] xrange of plot -- passed with no comment to PLOT
;	yrange	[INPUT] yrange of plot -- passed with no comment to PLOT
;	ccol	[INPUT] color of continuum line
;		* default is 150*!D.N_COLORS/256
;	lock	[INPUT] if set, locks up the relative positions of YCONT
;		* has an effect ONLY if YCONT is set
;		* side-effects:
;		  -- does not increase number of N({XY}CONT)
;		  -- left or middle click moves the entire curve up and down
;	_extra	[JUNK] here only to prevent crashing the program
;
;error calculation
;	since the continuum is being set interactively, there is no meaningful
;	way to define a statistical error on it.
;
;restrictions
;	works in X.  only?
;
;history
;	vinay kashyap (Oct98)
;	increased default plotting YRANGE (VK; Nov98)
;	added keyword LOCK (VK; Dec98)
;	cursor changes shape when mouse button is pressed (VK; FebMM)
;	allowed return of zero errors, cleaned up error handling, added CCOL
;	  cleaned up LOCK (VK; MarMM)
;	undid the cursor shape changes (VK; AprMM)
;	bug: wasn't holding the poisson v/s zero error toggle (VK; JulMM)
;	improved color-scale setting for 24-bit consoles (VK; FebMMI)
;	changed default to produce errors=0 rather than Poisson errors
;	  (VK; Aug03)
;-

;	usage
nx=n_elements(x) & ny=n_elements(y)
if nx lt 2 or nx ne ny then begin
  print,'Usage: cont=setcont(x,y,ysig,ycont=ycont,xcont=xcont,/const,$'
  print,'       xrange=xrange,yrange=yrange,ccol=ccol,/lock)'
  print,'  return interactively determined continuum'
  if nx eq 0 then return,[-1L] else return,0*x
endif

;	sundry initializations
xmin=min(x,max=xmax) & ymin=min(y,max=ymax)
draglim=1e-4*sqrt((xmax-xmin)^2+(ymax-ymin)^2)
dncolors=256. > !D.N_COLORS	;24-bit color screen temporary fix
if not keyword_set(ccol) then ccol=fix(150.*float(!d.n_colors)/dncolors+1)
nc=0L

;	is there an initial guess?
nyc=n_elements(ycont) & nxc=n_elements(xcont)
if nyc eq nx then begin & xcont=x & nxc=nx & endif
if nyc ne 0 and nyc eq nxc then begin
  nc=nyc & xc=xcont & yc=ycont & yce=0*(sqrt(abs(yc)+0.75)+1.)
endif

;	keyword check
if n_elements(xrange) lt 2 then xrange=[xmin,xmax]
if n_elements(yrange) lt 2 then begin
  dy=ymax-ymin
  yrange=[ymin-0.1*dy,ymax+0.1*dy]
endif
if not keyword_set(const) then const=0

;	are the continuum points in lockstep?
lockstep=0
if keyword_set(lock) then begin
  if nc gt 0 then lockstep=1
endif

;	until otherwise stated, go and mark points for the continuum
go_on=1
if not keyword_set(lockstep) then begin
  print,'LEFT CLICK to mark, LEFT DRAG to move'
  print,'MIDDLE CLICK/DRAG to delete, RIGHT CLICK/DRAG to exit'
  print,'CLICK ABOVE plot area to toggle piecewise constant curve'
  print,'CLICK BELOW plot area to toggle Poisson errors v/s zero errors'
endif else begin
  print,'LEFT to move, MIDDLE to undo, RIGHT to exit'
endelse

while go_on eq 1 do begin			;{endless loop

  ;	plot
  plot,x,y,/xs,/ys,psym=10,xrange=xrange,yrange=yrange
  ;	reset xmin,xmax,ymin,ymax from the plot
  xmin=min(!X.CRANGE) & xmax=max(!X.CRANGE)
  ymin=min(!Y.CRANGE) & ymax=max(!Y.CRANGE)

  ;	is there a continuum to overplot?
  if lockstep then oplot,xc,yc,col=ccol else begin
    nxc=n_elements(xc) & nyc=n_elements(yc) & nyce=n_elements(yce)
    if nyce eq 0 and nyc ne 0 then yce=0*yc
    if nxc eq 0 then begin & ycc=0*y & ycce=0*y & endif
    if nxc eq 1 then begin & ycc=0*y+yc(0) & ycce=0*y+yce(0) & endif
    if nxc gt 1 then begin
      if keyword_set(const) then begin
        xcc=fltarr(2L*nc) & ycc=xcc & ycce=xcc
	dxc=xc(1:*)-xc
        for i=0,nc-1 do begin
	  if i eq 0 then xcc(2*i)=xmin else xcc(2*i)=xc(i)-0.5*dxc(i-1)
	  if i lt nc-1 then xcc(2*i+1)=xc(i)+0.5*dxc(i) else xcc(2*i+1)=xmax
	  ycc([2*i,2*i+1])=yc(i)
	  ycce([2*i,2*i+1])=yce(i)
        endfor
        ycc=interpol(ycc,xcc,x)
        ycce=interpol(ycce,xcc,x)
      endif else begin
	ycc=interpol(yc,xc,x) & ycce=interpol(yce,xc,x)
      endelse
    endif
    oplot,x,ycc,col=ccol
    oplot,x,ycc+ycce,col=ccol,line=1 & oplot,x,ycc-ycce,col=ccol,line=1
    if nc gt 0 then oplot,xc,yc,psym=1
  endelse

  ;	mark
  cursor,x0,y0,/down,/data
  mbutton=!mouse.button
  if !D.NAME eq 'X' then begin
    ;if mbutton eq 1 then device,cursor_standard=74
    ;if mbutton eq 2 then device,cursor_standard=82
    ;if mbutton eq 4 then device,cursor_standard=100
  endif
  cursor,x1,y1,/up,/data
  ;if !D.NAME eq 'X' then device,/cursor_original
  drag=sqrt((x1-x0)^2+(y1-y0)^2)
  if drag lt draglim then drag=0.
  y0e=0*(sqrt(abs(y0)+0.75)+1.)
  y1e=0*(sqrt(abs(y1)+0.75)+1.)

  outclik=0		;click was inside plot area?
  if x0 gt xmax then outclik=outclik+2^0		;right of plot area
  if y0 gt ymax then outclik=outclik+2^1		;above plot area
  if x0 lt xmin then outclik=outclik+2^2		;left of plot area
  if y0 lt ymin then outclik=outclik+2^3		;below plot area
  ;
  if keyword_set(lockstep) then outclik=0
  if outclik ne 0 then begin
    mbutton=0	;mark/delete no points
    ;	toggle piecewise constant if click is above plot area
    if outclik ge 2 and outclik le 6 then begin
      const=1-const
      c1='changing from piecewise '
      if const eq 0 then c1=c1+'constant to piecewise linear' else $
	c1=c1+'linear to piecewise constant' & message,c1,/info
    endif
    ;	toggle poisson errors v/s zero errors if click is below plot area
    if outclik ge 8 then begin
      if n_elements(yc) gt 0 then begin
        c1='changing errors from '
        if n_elements(yce) eq 0 then yce=0*yc
        if yce(0) eq 0 then begin
	  message,c1+'zero to Poisson',/info
	  yce=sqrt(abs(yc)+0.75)+1.
        endif else begin
	  message,c1+'Poisson to zero',/info
	  yce(*)=0.
        endelse
      endif else message,'Nothing to make an error out of yet!',/info
    endif
  endif

  if mbutton eq 1 then begin					;(LEFT
    if drag eq 0 then begin				;(click
      ;	mark this point
      if nc eq 0 then begin
	xc=[x0] & yc=[y0] & yce=[y0e]
      endif else begin
	if lockstep then begin		;(just move the whole curve
	  ;find nearest point to X0 and peg it to Y0
	  tmp=min(abs(xc-x0),ic) & dely=y0-yc(ic)
	  yc=yc+dely
	endif else begin		;)(else add a new point
	  xc=[xc,x0] & yc=[yc,y0]
	  if yce(0) ne 0 then yce=[yce,y0e] else yce=[yce,0]
	  ox=sort(xc) & xc=xc(ox) & yc=yc(ox) & yce=yce(ox)
	endelse				;lockstep?)
      endelse
    endif else begin					;)(click+drag
      ;	find nearest point and move it to new position
      if nc gt 0 then begin
	if lockstep then begin		;(move curve as a whole
	  tmp=min(abs(xc-x0),ic) & dely=y1-yc(ic)
	  yc=yc+dely
	endif else begin			;)(move just this point
          dd=sqrt((xc-x0)^2+(yc-y0)^2) & tmp=min(dd,ic)
	  xc(ic)=x1 & yc(ic)=y1 & if yce(0) ne 0 then yce(ic)=y1e
	  ox=sort(xc) & xc=xc(ox) & yc=yc(ox) & yce=yce(ox)
	endelse				;lockstep?)
      endif
    endelse						;drag)
    nc=n_elements(xc)
  endif								;LEFT)

  if mbutton eq 2 then begin					;(MIDDLE
    if lockstep then begin		;(undo previous shift
      if not keyword_set(dely) then dely=0.
      yc=yc-dely
      dely=-dely
    endif else begin			;)(delete nearest point
      ;	delete nearest point
      if nc gt 0 then begin
        dd=sqrt((xc-x0)^2+(yc-y0)^2) & tmp=min(dd,ic)
        ii=lindgen(nc) & oi=where(ii ne ic,moi)
        if moi gt 0 then begin
	  xc=xc(oi) & yc=yc(oi) & yce=yce(oi) & nc=moi
        endif else nc=0L
      endif
    endelse				;lockstep?)
  endif								;MIDDLE)

  if mbutton eq 4 then begin					;(RIGHT
    go_on=0					;exit
  endif								;RIGHT)

endwhile					;end of endless loop}

;	outputs
nxc=n_elements(xc) & nyc=n_elements(yc)
if nxc gt 0 and nxc eq nyc then begin
  xcont=xc & ycont=yc & ysig=yce
endif
;
if nxc eq 0 then begin & cont=0*y & ysig=0*y & endif
if nxc eq 1 then begin & cont=0*y+yc(0) & ysig=0*y+yce(0) & endif
if nxc gt 1 then begin				;(interpolate YC onto X
  if keyword_set(const) then begin
    xcc=fltarr(2L*nc) & cont=xcc & ysig=xcc
    dxc=xc(1:*)-xc
    for i=0,nc-1 do begin
      if i eq 0 then xcc(2*i)=xmin else xcc(2*i)=xc(i)-0.5*dxc(i-1)
      if i lt nc-1 then xcc(2*i+1)=xc(i)+0.5*dxc(i) else xcc(2*i+1)=xmax
      cont([2*i,2*i+1])=yc(i)
      ysig([2*i,2*i+1])=yce(i)
    endfor
    cont=interpol(cont,xcc,x)
    ysig=interpol(ysig,xcc,x)
  endif else begin
    cont=interpol(yc,xc,x) & ysig=interpol(yce,xc,x)
  endelse
endif						;YC(XC))

;;	error-bars
;if ysig(0) ne 0 then ysig=sqrt(abs(cont)+0.75)+1.

;	plot
plot,x,y,/xs,/ys,psym=10,xrange=xrange,yrange=yrange
oplot,x,cont,col=ccol
oplot,x,cont+ysig,col=ccol,line=1 & oplot,x,cont-ysig,col=ccol,line=1
if nc gt 0 then oplot,xc,yc,psym=1

return,cont
end
