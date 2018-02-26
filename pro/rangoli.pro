pro rangoli,ststr,pstyle,x,y,pnum=pnum,lnum=lnum,group=group,$
	ytop=ytop,ybot=ybot,oststr=oststr, _extra=e
;+
;procedure	rangoli
;	update the location of the labels according to specified style
;
;	"rangoli" is the intricate pattern people draw as decorations
;	during festivals and such.
;	pronounced "run-go-li"
;
;syntax
;	rangoli,ststr,pstyle,pnum=pnum,lnum=lnum,group=group
;
;parameters
;	ststr	[I/O; required] the structure describing the state
;		of the plot.  see MUDRA for a description.
;	pstyle	[INPUT] the plotting style
;		0	no change
;		1.F	align vertically to hang from ceiling
;			"F" is the length of YPATH as fraction of YRANGE
;		2.F	as (1), but stagger along X to avoid overlaps
;		3.F	align vertically starting from YTOP
;		4.F	as (2), but stagger along X
;		5.F	align horizontally along YTOP
;		6.F	as (2), but stagger along Y
;		-1	connect (LOC.X,YBOT) to L#.POS with straight line
;		-2	connect with half-T
;		-3	connect with inverted half-T
;		-4.F	connect with half-Y
;			"F" is length of straight leg, as fraction of YRANGE
;		-5.F	connect with inverted half-Y
;		-6.F	as (5), with another straight leg at the end
;	X	[INPUT] points where the curve is defined
;	Y	[INPUT] the curve Y[X] (usually a spectrum)
;		* if N(Y).NE.N(X), both are ignored
;		* in fact, they're ignored anyway unless LOC.Y or the
;		  window ranges are placeholders
;
;keywords
;	pnum	[INPUT; default=0] target plotting window
;	lnum	[INPUT; default=all] array of label indices to reposition
;	group	[INPUT] if given, makes the grouped arrangement "permanent"
;	ytop	[INPUT] 2-element array determines the location of the
;		label, YPOS=YMAX-(YTOP(0)+YTOP(1)*(YMAX-YMIN))
;		* YTOP(1)=0 if not specified
;		* special case: if single element && <0.5, then
;		  YTOP(1) <-- YTOP(0) & YTOP(0)=0.
;	ybot	[INPUT] 2-element array determines the location where the
;		"arrows" end, YLOC=YBOT(0)+YBOT(1)*Y(XLOC)
;		* if scalar or 1-element array, then YBOT(1)=1.
;	oststr	[OUTPUT] old STSTR, the unchanged version
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (MIMSep)
;-

;	usage
ns=n_tags(ststr)
if ns eq 0 then begin
  print,'Usage: rangoli,ststr,pstyle,x,y,pnum=pnum,lnum=lnum,group=group,$'
  print,'       ytop=ytop,ybot=ybot,oststr=oststr'
  print,'  update the label locations according to specified style'
  return
endif

;	check input
nloc=n_tags(ststr.(1))
if nloc eq 0 then begin
  message,'input not understandable',/info & return
endif
nlabel=n_elements(ststr.(1).(0))
if nlabel ne ns-2L then begin
  message,'input messed up',/info & return
endif
oststr=ststr
winstr=ststr.(0)
locstr=ststr.(1) & xloc=locstr.X & yloc=locstr.Y & yxloc=yloc
;
np=n_elements(pstyle) & pp=0. & if np ne 0 then pp=float(pstyle(0))
;
nx=n_elements(x) & ny=n_elements(y)
if nx gt 0 and nx eq ny then begin
  minx=min(xloc,max=maxx)
  for i=0L,nlabel-1L do begin
    tmp=min(abs(xloc(i)-x),iy) & yxloc(i)=y(iy)
  endfor
  miny=min(y)/2. & maxy=max(yxloc)*2.
endif else begin
  minx=0.0 & maxx=1.0 & miny=0.0 & maxy=1.0
endelse

;	check keywords
ip=0L & if keyword_set(pnum) then ip=long(pnum)
;
ll=lindgen(nlabel) & if n_elements(lnum) gt 0 then ll=[long(lnum(*))]
;
gg=locstr.GROUP & ng=n_elements(group)
for i=0L,ng-1L do if i lt nlabel then gg(i)=group(i)
;
topy=[0.,0.2] & nt=n_elements(ytop)
case nt of
  0:	;use default
  1: begin
    if ytop(0) lt 0.5 then topy(1)=ytop(0) else topy=[ytop(0),0.]
  end
  else: topy=[ytop(0),ytop(1)]
endcase
;
boty=[0.,1.1] & nb=n_elements(ybot)
case nb of
  0:	;use default
  1: boty=[ybot(0),1.1]
  else: boty=[ybot(0),ybot(1)]
endcase

;	double check keywords
nrow=(winstr.MULTI)(1) & ncol=(winstr.MULTI)(1) & nplot=nrow*ncol > 1
if ip lt 0 or ip ge nplot then ip=0L
ol=where(ll ge 0 and ll lt nlabel,mol)
if mol eq 0 then begin
  message,'say which labels again?',/info & return
endif
ll=ll(ol)

;	get ranges for window of choice
xmin=winstr.XMIN(ip) & xmax=winstr.XMAX(ip)
ymin=winstr.YMIN(ip) & ymax=winstr.YMAX(ip)
;
;	if the ranges are meaningless, get em from the spectrum
if xmin ge xmax then begin
  ststr.WINDOW.XMIN(ip)=minx & ststr.WINDOW.XMAX(ip)=maxx
endif
if ymin ge ymax then begin
  ststr.WINDOW.YMIN(ip)=miny & ststr.WINDOW.YMAX(ip)=maxy
endif

;	handle the styles
for i=0L,mol-1L do begin		;{step through the labels
  il=ll(i)			;current label index
  ff=abs(pp-fix(pp)) & pp=fix(pp)	;basic style
  case pp of				;{step through the styles
    1: begin
      c=ststr.(il+2).LABEL & ncom=2*(n_elements(str_sep(c,'!'))-1L)
      ;	find the length of string
      clen=strlen(c)-ncom & csiz=ststr.(il+2).SIZE & chars=!d.x_ch_size
      x1=xloc(il) & y1=ymax & tmp=convert_coord(x1,y1,/data,/to_dev)
      x1d=tmp(0) & y1d=tmp(1) & x0d=x1d & y0d=y1d-clen*csiz*chars
      ;	and hang it from the top
      tmp=convert_coord(x0d,y0d,/dev,/to_data) & x0=tmp(0) & y0=tmp(1) > ymin
      if ff eq 0 then ff=0.1
      dropy=ff*(ymax-ymin) & yloc(il)=y0-dropy
      ;	and connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      xpath(*)=xloc(il) & ypath=y0-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=yloc(il)
      ststr.(il+2).ALIGN='LEFT' & ststr.(il+2).ARRANGE='ROW'
      ststr.(il+2).ORIENT=90.0
      ststr.(il+2).POS=[xloc(il),y0]
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
      ;ststr.(il+2).UNDERLINE=0.0 & ststr.(il+2).SIDELINE=0.0
    end
    2: begin
      message,'sorry, not implemented yet',/info
    end
    3: begin
      ;	where to place label along Y?
      y0=ymax-(topy(0)+topy(1)*(ymax-ymin))
      if ff eq 0 then ff=0.1
      dropy=ff*(ymax-ymin) & yloc(il)=y0-dropy
      ;	and connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      xpath(*)=xloc(il) & ypath=y0-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=yloc(il)
      ststr.(il+2).ALIGN='LEFT' & ststr.(il+2).ARRANGE='ROW'
      ststr.(il+2).ORIENT=90.0
      ststr.(il+2).POS=[xloc(il),y0]
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
      ststr.(il+2).UNDERLINE=0.0 & ststr.(il+2).SIDELINE=0.0
    end
    4: begin
      message,'sorry, not implemented yet',/info
    end
    5: begin
      ;	where to place it along Y?
      y0=ymax-(topy(0)+topy(1)*(ymax-ymin))
      if ff eq 0 then ff=0.5
      dropy=ff*(ymax-ymin) & yloc(il)=y0-dropy
      ;	and connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      xpath(*)=xloc(il) & ypath=y0-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=yloc(il)
      ststr.(il+2).ALIGN='DOWN' & ststr.(il+2).ARRANGE='ROW'
      ststr.(il+2).ORIENT=0.0
      ststr.(il+2).POS=[xloc(il),y0]
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
    end
    6: begin
      message,'sorry, not implemented yet',/info
    end
    -1: begin					;(straight line
      ;	figure out end point
      y0=boty(0)+boty(1)*yxloc(il)
      ;	and the label position is, just for reference --
      pos=ststr.(il+2).POS
      ;	connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      dropx=pos(0)-xloc(il) & dropy=pos(1)-y0
      xpath=xloc(il)-dropx*fpath & ypath=y0-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=y0
      if pos(0) gt xloc(il) then ststr.(il+2).ALIGN='LEFT'
      if pos(0) lt xloc(il) then ststr.(il+2).ALIGN='RIGHT'
      ststr.(il+2).ORIENT=0.0
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
    end						;straight line)
    -2: begin					;(half-T
      ;	figure out end point
      y0=boty(0)+boty(1)*yxloc(il)
      ;	and the label position is, just for reference --
      pos=ststr.(il+2).POS
      ;	connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      dropy=pos(1)-y0
      xpath(*)=xloc(il) & ypath=pos(1)-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=y0
      if pos(0) gt xloc(il) then ststr.(il+2).ALIGN='LEFT' else $
      	ststr.(il+2).ALIGN='RIGHT'
      ststr.(il+2).ORIENT=0.0
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
    end						;half-T)
    -3: begin					;(inverted half-T
      ;	figure out end point
      y0=boty(0)+boty(1)*yxloc(il)
      ;	and the label position is, just for reference --
      pos=ststr.(il+2).POS
      ;	connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      dropx=pos(0)-xloc(il)
      xpath=pos(0)-dropx*fpath & ypath(*)=yloc(il)
      ;	update
      ststr.LOC.Y(il)=y0
      if pos(0) gt xloc(il) then ststr.(il+2).ALIGN='DOWN'
      if pos(0) lt xloc(il) then ststr.(il+2).ALIGN='UP'
      ststr.(il+2).ORIENT=0.0
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
    end						;inverted half-T)
    -4: begin					;(half-Y
      ;	figure out end point
      y0=boty(0)+boty(1)*yxloc(il)
      ;	and the label position is, just for reference --
      pos=ststr.(il+2).POS
      ;	and the location of the vertex is
      if ff eq 0 then ff=0.1
      xv=xloc(il) & yv=y0+ff*(ymax-ymin)
      ;	connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      dropy=yv-y0
      xpath(*)=xv & ypath=yv-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=y0
      ststr.(il+2).ALIGN='DOWN'
      ststr.(il+2).ORIENT=0.0
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
    end						;half-Y)
    -5: begin					;(inverted half-Y
      ;	figure out end point
      y0=boty(0)+boty(1)*yxloc(il)
      ;	and the label position is, just for reference --
      pos=ststr.(il+2).POS
      ;	and the location of the vertex is
      if ff eq 0 then ff=0.1
      xv=pos(0) & yv=y0+ff*(ymax-ymin)
      ;	connect the dots
      xpath=ststr.(il+2).XPATH & ypath=ststr.(il+2).YPATH
      npath=n_elements(xpath) & fpath=findgen(npath)/npath
      dropx=xv-xloc(il) & dropy=yv-y0
      xpath=xv-dropx*fpath & ypath=yv-dropy*fpath
      ;	update
      ststr.LOC.Y(il)=y0
      ststr.(il+2).ALIGN='DOWN'
      ststr.(il+2).ORIENT=0.0
      ststr.(il+2).XPATH=xpath & ststr.(il+2).YPATH=ypath
    end						;inverted half-Y)
    -6: begin					;(sideways Z
      message,'sorry, not implemented yet',/info
    end						;sideways Z)
    else: message,'plotting style '+strtrim(pp,2)+' is unknown',/info ;nothing
  endcase				;PP}
endfor					;I=0,MOL-1}

return
end
