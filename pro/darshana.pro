pro darshana,x,y,ststr,wid=wid,pid=pid,xsize=xsize,ysize=ysize,colmax=colmax,$
	stretch=stretch, _extra=e
;+
;procedure	darshana
;	display 1D curve with annotation
;
;	"darshana" is the verb form of "vision", the Joan kind.
;	pronounced "the-r-shun-ah"
;
;syntax
;	darshana,x,y,ststr,wid=wid,pid=pid,xsize=xsize,ysize=ysize,$
;	colmax=colmax,stretch=stretch, PLOT_KEYWORDS
;
;parameters
;	x	[INPUT; required] points where the curve is defined
;	y	[INPUT; required] the curve Y[X] (usually a spectrum)
;		* size of Y must match that of X
;	ststr	[INPUT; required] the state structure containing all
;		the necessary information on what labels to put where.
;		* see description in KALPANA.PRO
;
;keywords
;	wid	[INPUT] window ID for plots (ignored if device is not X)
;	pid	[INPUT] if set, plots only the PIDth plot
;	xsize	[INPUT] window size (goes with WID)
;	ysize	[INPUT] window size (goes with WID)
;	colmax	[INPUT] maximum available color (!D.N_COLORS, we need
;		this hack to handle postscript plots)
;	stretch	[INPUT] set to float number to stretch all y-positions
;		by constant amount
;		* default is 1
;		* 0, etc. are ignored
;	_extra	[INPUT ONLY] pass defined keywords to PLOT
;		(anything except the stuff provided for in STSTR.WINDOW,
;		TITLE,[XY]TITLE,[XY]STYLE,[XY]LOG,CHARSIZE,XRANGE,YRANGE)
;
;side-effects
;	makes plots
;
;history
;	vinay kashyap (SepMIM)
;	added _EXTRA to PLOT (VK; AugMM)
;-

;	usage
ok='ok'
nx=n_elements(x) & ny=n_elements(y) & np=n_params()
if nx eq 0 then ok='Missing X' else $
 if ny eq 0 then ok='Missing Y(X)' else $
  if nx ne ny then ok='Y(X) and X are incompatible' else $
   if np lt 3 then ok='insufficiant parameters'
if ok ne 'ok' then begin
  print,'Usage: darshana,x,y,ststr,wid=wid,pid=pid,xsize=xsize,ysize=ysize,$'
  print,'       colmax=colmax,stretch=stretch, PLOT_KEYWORDS'
  print,'  display curve with annotation'
  message,ok,/info & return
endif

;	check keywords
if not keyword_set(xsize) then xsize=!d.x_size
if not keyword_set(ysize) then ysize=!d.y_size
if n_elements(wid) ne 0 then begin
  if !D.NAME eq 'X' then window,wid,xsize=xsize,ysize=ysize
endif
if not keyword_set(colmax) then colmax=!D.N_COLORS
yscale=1.0 & if keyword_set(stretch) then yscale=float(stretch)
if yscale eq 0.0 then yscale=1.0
;
idplot=0 & multiplot=1	;NOTE: just coz this is set does not mean multi plots
if n_elements(pid) gt 0 then begin
  multiplot=0 & idplot=long(pid(0)) > 0
endif else multiplot=1	

;	is state structure defined?
ns=n_elements(ststr) & ms=n_tags(ststr) & ok='ok'
if ns eq 0 then ok='State undefined.' else $
 if ms eq 0 then ok='State variable is junk.' else begin
   stnam=tag_names(ststr)
   if stnam(0) ne 'WINDOW' then ok='incorrect field.' else $
    if stnam(1) ne 'LOC' then ok='incorrect field.' else begin
      nl=n_elements(ststr.(1).(0))
      if ms ne nl+2L then ok='insufficient fields.'
    endelse
 endelse
if ok ne 'ok' then begin
  message,'Input state error: '+ok,/info
  return
endif

;	break into components
wstr=ststr.WINDOW & xloc=ststr.LOC.X & yloc=ststr.LOC.Y & iloc=ststr.LOC.GROUP
;
nrow=(ststr.WINDOW.multi)(1) > 1
ncol=(ststr.WINDOW.multi)(2) > 1
nplot=n_elements(ststr.WINDOW.XMIN)
;nplot=nrow*ncol < n_elements(ststr.WINDOW.XMIN)
pmulti=!p.multi
if multiplot ne 0 then !p.multi=wstr.MULTI else !p.multi=[0,1,1,0,0]
;
title=wstr.TITLE
stitle=strarr(nplot) & titlex=stitle & titley=stitle
stylex=intarr(nplot) & styley=stylex & logx=stylex & logy=stylex
chars=fltarr(nplot) & xmin=chars & xmax=chars & ymin=chars & ymax=chars
;
subt=wstr.SUBTITLE & n=n_elements(subt) & m=n<nplot
	stitle(*)=subt(0) & stitle(0:m-1L)=subt(0:m-1L)
xtit=wstr.XTITLE & n=n_elements(xtit) & m=n<nplot
	titlex(*)=xtit(0) & titlex(0:m-1L)=xtit(0:m-1L)
ytit=wstr.YTITLE & n=n_elements(ytit) & m=n<nplot
	titley(*)=ytit(0) & titley(0:m-1L)=ytit(0:m-1L)
xsty=wstr.XSTYLE & n=n_elements(xsty) & m=n<nplot
	stylex(*)=xsty(0) & stylex(0:m-1L)=xsty(0:m-1L)
ysty=wstr.YSTYLE & n=n_elements(xsty) & m=n<nplot
	styley(*)=xsty(0) & styley(0:m-1L)=ysty(0:m-1L)
xlg=wstr.XLOG & n=n_elements(xlg) & m=n<nplot
	logx(*)=xlg(0) & logx(0:m-1L)=xlg(0:m-1L)
ylg=wstr.YLOG & n=n_elements(ylg) & m=n<nplot
	logy(*)=ylg(0) & logy(0:m-1L)=ylg(0:m-1L)
char=wstr.CHARS & n=n_elements(char) & m=n<nplot
	chars(*)=char(0) & chars(0:m-1L)=char(0:m-1L)
xmn=wstr.XMIN & n=n_elements(xmn) & m=n<nplot
	xmin(*)=xmn(0) & xmin(0:m-1L)=xmn(0:m-1L)
xmx=wstr.XMAX & n=n_elements(xmx) & m=n<nplot
	xmax(*)=xmx(0) & xmax(0:m-1L)=xmx(0:m-1L)
ymn=wstr.YMIN & n=n_elements(ymn) & m=n<nplot
	ymin(*)=ymn(0) & ymin(0:m-1L)=ymn(0:m-1L)
ymx=wstr.YMAX & n=n_elements(ymx) & m=n<nplot
	ymax(*)=ymx(0) & ymax(0:m-1L)=ymx(0:m-1L)

;	display
for i=0L,nplot-1L do begin		;{now plot
  if multiplot eq 0 and i ne idplot then goto,endplot	;{yeah, a goto.

  x0=xmin(i) & x1=xmax(i) & y0=ymin(i) & y1=ymax(i) & tt=stitle(i) & subt=''
  minx=x0 < x1
  maxx=x1 > x0
  ;	if YMAX = YMIN, float YMAX
  oo=where(x ge minx and x le maxx,moo)
  if y0 eq y1 and moo gt 0 then y1=max(y(oo))*1.5
  miny=y0 < y1
  maxy=y1 > y0

  if multiplot eq 0 or nplot eq 1 then begin
    tt=title & subt=stitle(0)
  endif
  plot,x,y,psym=10,xrange=[x0,x1],yrange=[y0,y1],$
	xstyle=stylex(i),ystyle=styley(i),xlog=logx(i),ylog=logy(i),$
	xtitle=titlex(i),ytitle=titley(i),title=tt,subtitle=subt,$
	charsize=chars(i),/nodata,color=colmax-1L, _extra=e
  oplot,x,y,psym=10,color=colmax-2L

  oi=where(xloc ge minx and xloc le maxx,moi)
  for j=0L,moi-1L do begin			;{plot labels
    k=oi(j) & kg=iloc(k)
    ll=ststr.(k+2) & llg=ststr.(kg+2)
    pos=ll.POS & label=ll.LABEL & xpath=ll.XPATH & ypath=ll.YPATH
    align=llg.ALIGN & orient=llg.ORIENT & size=llg.SIZE & thick=llg.THICK
    labc=llg.LABCOLOR & linc=llg.LINCOLOR & arr=llg.ARRANGE
    uline=llg.UNDERLINE & sline=llg.SIDELINE

    if ll.SIZE le 0 or size le 0 then goto,nextlabl	;{yeah, another goto
    	;this is how a label can be hidden or displayed as needed

    ;	if label position is at or below the minimum yrange, float it up
    if pos(1) le miny then begin
      tmp=min(abs(pos(0)-x),imn) & pos(1)=y(imn)+0.1*(maxy-y(imn))
    endif

    ;	some initializations
    xchar=!d.x_ch_size*size & ychar=!d.y_ch_size*size

    ;		define box in which to put label
    ;	break label down into fields
    c1=label & cc=str_sep(c1,'|') & ncc=n_elements(cc)
    if strupcase(arr) eq 'ROW' then begin
      ;	get rid of the field separator
      c1=' ' & for ic=0L,ncc-1L do c1=c1+cc(ic)+' '
    endif else begin
      c1=' ' & for ic=0L,ncc-2L do c1=c1+cc(ic)+'!C ' & c1=c1+cc(ncc-1L)
    endelse
    label=c1
    ;	find maximum width
    bwidth=1
    for ic=0L,ncc-1L do begin
      c1=str_sep(cc(ic),'!') & nc1=(n_elements(c1)-1L)
      bw=strlen(cc(ic))-2*nc1+1
      if strupcase(arr) eq 'ROW' then bwidth=bwidth+bw else $
      	bwidth = bwidth > bw
    endfor
    if strupcase(arr) eq 'ROW' then ncc=1L
    ;	pad and convert to pixel sizes
    bwidth=(bwidth+1.)*xchar & bheight=(ncc+1.)*ychar
    ;	convert POS into device coords
    pd=convert_coord(pos,/data,/to_dev)
    ;	set the bounds
    dw=0.5*bwidth & dh=0.5*bheight
    xll=-dw & yll=-dh & xlr=dw & ylr=-dh & xul=-dw & yul=dh & xur=dw &
    yur=dh
    ;	rotate the box
    sinO=sin(-orient*!pi/180.) & cosO=cos(-orient*!pi/180.)
    xllr=xll*cosO+yll*sinO & yllr=-xll*sinO+yll*cosO
    xlrr=xlr*cosO+ylr*sinO & ylrr=-xlr*sinO+ylr*cosO
    xurr=xur*cosO+yur*sinO & yurr=-xur*sinO+yur*cosO
    xulr=xul*cosO+yul*sinO & yulr=-xul*sinO+yul*cosO
    ;
    dropx=(-0.5*xchar)*cosO+(1.5*ychar)*sinO
    dropy=-(-0.5*xchar)*sinO+(1.5*ychar)*cosO
    ;	find the base position on the box
    case strupcase(align) of
      'LEFT': begin
	xll=xulr & xlr=xllr & xur=xlrr & xul=xurr
	yll=yulr & ylr=yllr & yur=ylrr & yul=yurr
	xlab=xll-dropx & ylab=yll-dropy
      end
      'RIGHT': begin
	xll=xlrr & xlr=xurr & xur=xulr & xul=xllr
	yll=ylrr & ylr=yurr & yur=yulr & yul=yllr
	xlab=xur-dropx & ylab=yur-dropy
      end
      'UP': begin
	xll=xurr & xlr=xulr & xur=xllr & xul=xlrr
	yll=yurr & ylr=yulr & yur=yllr & yul=ylrr
	xlab=xlr-dropx & ylab=ylr-dropy
      end
      else: begin
	xll=xllr & xlr=xlrr & xur=xurr & xul=xulr
	yll=yllr & ylr=ylrr & yur=yurr & yul=yulr
	xlab=xul-dropx & ylab=yul-dropy
      end
    endcase
    xbase=0.5*(xll+xlr) & ybase=0.5*(yll+ylr)
    ;	move the box
    xmid=pd(0)-xbase & ymid=pd(1)-ybase
    xll=xll+xmid & xlr=xlr+xmid & xur=xur+xmid & xul=xul+xmid
    yll=yll+ymid & ylr=ylr+ymid & yur=yur+ymid & yul=yul+ymid
    xlab=xlab+xmid & ylab=ylab+ymid
    ;	place the label
    xyouts,xlab,ylab,label,color=labc,charsize=size,$
	orient=orient,/clip,/dev

    ;	define the underline and corresponding pincers
    uln=abs(uline) & if uln le 1 then fac=uln else fac=1.-(uln-fix(uln))
    xwdt=fac*0.5*(xlr-xll) & ywdt=fac*0.5*(ylr-yll)
    midx=0.5*(xll+xlr) & midy=0.5*(yll+ylr)
    if uline lt 0 then begin
      midx=0.5*(xul+xur) & midy=0.5*(yul+yur)
    endif
    pincyl=(uline/100.)*(yul-yll) & pincyr=(uline/100.)*(yur-ylr)
    pincxl=(uline/100.)*(xul-xll) & pincxr=(uline/100.)*(xur-xlr)
    ux=[midx-xwdt+pincxl,midx-xwdt,midx+xwdt,midx+xwdt+pincxr]
    uy=[midy-ywdt+pincyl,midy-ywdt,midy+ywdt,midy+ywdt+pincyr]
    if uline ne 0 then plots,ux,uy,/dev,thick=thick,color=linc

    ;	define the sideline and corresponding pincers
    sln=abs(sline) & if sln le 1 then fac=sln else fac=1.-(sln-fix(sln))
    xwdt=fac*0.5*(xul-xll) & ywdt=fac*0.5*(yul-yll)
    midx=0.5*(xll+xul) & midy=0.5*(yll+yul)
    if sline lt 0 then begin
      midx=0.5*(xlr+xur) & midy=0.5*(ylr+yur)
    endif
    pincxu=(sline/100.)*(xur-xul) & pincxd=(sline/100.)*(xlr-xll)
    pincyu=(sline/100.)*(yur-yul) & pincyd=(sline/100.)*(ylr-yll)
    sx=[midx-xwdt+pincxu,midx-xwdt,midx+xwdt,midx+xwdt+pincxd]
    sy=[midy-ywdt+pincyu,midy-ywdt,midy+ywdt,midy+ywdt+pincyd]
    if sline ne 0 then plots,sx,sy,/dev,thick=thick,color=linc

    ;	define the path FROM label TO feature
    if yloc(k) eq 0 then begin
      tmp=min(abs(xloc(k)-x),iy) & yloc(k)=y(iy)
    endif
    xp=[pos(0),xpath,xloc(k)] & yp=[pos(1),ypath,yloc(k)]
    oplot,xp,yp,thick=thick,color=linc

    nextlabl:	;the goto comes here}
  endfor					;J=0,MOI-1}

  endplot:	;the goto comes here}
endfor					;I=0,NPLOT-1}

;	overplot main title
if multiplot ne 0 and nplot gt 1 then begin
  !p.multi=[1,0,0]
  plot,[0,1],[0,1],/nodata,xstyle=4,ystyle=4
  xyouts,0.,1.01,title,charsize=max(ststr.WINDOW.CHARS)
endif

;	reset graphics
!p.multi=pmulti

return
end
