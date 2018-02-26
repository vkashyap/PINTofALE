pro new_plotid,line,lamda,spec,way=way,finger=finger,_extra=e
;+
;procedure	new_plotid
;	easily said: uses the output of LINEID to generate publication
;	quality graph with all the bells and whistles.
;
;syntax
;	new_plotid,line_id,lamda,spec,way=way,/finger, PLOT_KEYWORDS
;
;parameters
;	line	[INPUT; required] Line ID information
;		* type
;		  IDL> help,lineid(),/str
;		  for description and format
;	lamda	[INPUT; required] wavelengths at which spectrum is defined
;	spec	[INPUT; required] the spectrum
;
;keywords
;	way	[I/O] a structure that contains all the pertinent info
;		for the labels
;		* the actual labels (LINE will be generally left unused),
;		  label positions, label orientation, line types, line
;		  thickness, and character sizes
;		* WAY may be overwritten by setting FINGER
;	finger	[INPUT] if set, allows interactive placement of line labels
;	_extra	[INPUT] allows this routine to act as wrapper to PLOT
;
;restrictions
;	* requires heavy graphical capability
;	* LINEID must have been run previously
;	* uses subroutines LABLPLOT, WMENU
;
;history
;	vinay kashyap (Jan96)
;	button press status now stored in !MOUSE, not !ERR (VK; Apr09)
;-

;	usage
if n_params() lt 3 then begin
  print,'Usage: new_plotid,line_id,lamda,spec,way=way,/finger'
  print,'  plots spectrum with overlaid line IDs'
  return
endif

;	check inputs for consistency
nope=0 & szl=size(line) & nszl=n_elements(szl)
nx=n_elements(lamda) & ny=n_elements(spec) & comm=(tag_names(line))(1)
if szl(nszl-2) ne 8 then nope=1			;LINE not a structure
if comm eq 'WVL_COMMENT' then nope=2		;LINEID_HELP
if nx ne ny then nope=3				;incompatible array sizes
if nx eq 0 or ny eq 0 then nope=4		;missing array
if nx eq 1 then nope=5				;array has 1 element
case nope of
  1: message,'Line ID info not in a structure',/info
  2: message,'Line ID info not available',/info
  3: message,'spectral array size mismatch',/info
  4: message,'missing spectral array',/info
  5: message,'need at least 2 points in spectral array',/info
  else:	;no problemo
endcase
if nope gt 0 then return

;	initialize
x=lamda & y=spec & wvl=line.wvl & nid=n_elements(wvl) & label=strarr(nid)

;	set up defaults
xp=wvl					;x-position of labels
yp=fltarr(nid)+1.3*(max(y)-min(y))	;y-position of labels
lbi=intarr(nid)				;label type
lbl=strarr(nid)				;labels
ori=intarr(nid)				;label orientation [H,V,I,>]
fnt=intarr(nid)+3			;label font number
lsz=fltarr(nid)+1.			;label size
uln=intarr(nid)				;under/side line (1) or not (0)?
clt=intarr(nid)				;connecting line type [|,/,/|,L,T,N]
lln=fltarr(nid)+1			;2nd leg lengths, fraction of default
dot=intarr(nid)				;linestyle, for plot
lwd=fltarr(nid)+1.			;line sizes

;	generate appropriate label
for i=0,nid-1 do lbl(i)=genlabel(line.(i+1),lbi(i),wvl(i))

;	check keywords
if keyword_set(way) then begin		;decipher label style
  szw=size(way) & nszw=n_elements(szw)
  if n_tags(way) eq 12 then begin	;style is in right structure
    xp=way.(0) & yp=way.(1) & lbi=way.(2) & lbl=way.(3)
    ori=way.(4) & fnt=way.(5) & lsz=way.(6) & uln=way.(7)
    clt=way.(8) & lln=way.(9) & dot=way.(10) & lwd=way.(11)
  endif
endif

;	plot
labl_plot,x,y,wvl,lbl,ori,fnt,lsz,uln,clt,lln,dot,lwd,$
	psym=10,/xstyle,/ystyle,xtitle='!3['+string(byte(197))+']',_extra=e

;	are we done?  then quit!
nope=0
if not keyword_set(finger) then nope=1
if !d.name ne 'X' then nope=2
if nope ne 0 then return

;	interactive stuff
go_on=1 & storxp=xp & storyp=yp & storsty=sty

shelp=[	'	Use mouse buttons to select label to reformat',$
	'LEFT: reposition label',$
	'MIDDLE: allow label restyle via keyboard input',$
	'RIGHT: allow connecting line restyle via keyboard input',$
	'	keyboard options',$
	'q,Q,x,X: quit',$
	'c: cycle possible labels	C: new label',$
	'd: toggle dotted/solid lines',$
	'f: cycle font			F: specify vector font number',$
	'u: undo previous positioning	U: undo previous styling',$
	'z,Z: zero out x-offset',$
	'"_": toggle underline']

print,'		Use mouse buttons to select label'
print,'LEFT/MIDDLE allows repositioning/restyling; RIGHT allows reset/quit'
print,'	L/M: left click or drag to position label'
print,'	     middle click to cycle style, drag to select style'
print,'	     right click or drag to go to next label'
print,'	R: left click to reset initial position, left drag to zero x-offset'
print,'	   middle click or drag to reset initial label style'
print,'	   right click or drag to QUIT'

idlvers=float(!VERSION.release) & if idlvers le 0 then idlvers=7
while go_on eq 1 do begin			;{reposition/restyle labels
  mbutton=0 & ohy=0L

  repeat begin				;(select label by clicking
    cursor,xc,yc,/change,/data
    if idlvers lt 5 then mbutton=!err else mbutton=!MOUSE.button
    oo=where(abs(wvl-xc) eq min(abs(wvl-xc)))
    if ohy ne oo(0) then begin
      xyouts,wvl(ohy),max(y),'!95!3',charsize=3,color=0,align=0.5	;erase
      xyouts,wvl(oo(0)),max(y),'!95!3',charsize=3,align=0.5		;mark
      ohy=oo(0)
    endif
  endrep until (mbutton ne 0)		;label selected)
  cursor,xc,yc,/up,/data		;clear left-over button transitions

  lablplot,x,y,wvl,label,xp,yp,sty,psym=10,hylyt=ohy+1,$
	xtitle='!3['+string(byte(197))+']',_extra=e 	;replot

  if mbutton eq 4 then begin			;{right button
    cursor,x0,y0,/down,/data
    if idlvers lt 5 then mbt=!err else mbt=!MOUSE.button
    cursor,x1,y1,/up,/data
    drag=abs(x1-x0)+abs(y1-y0)
    case mbt of
      1: begin				;left button
	if drag eq 0 then begin		;click
	  xp(ohy)=storxp(ohy)
	  yp(ohy)=storyp(ohy)
	endif else begin		;drag
	  xp(ohy)=wvl(ohy)
	endelse
      end
      2: begin				;middle button
	sty(ohy)=storsty(ohy)
      end
      4: go_on=0			;QUIT
      else: message,'unknown button transition',/info
    endcase
  endif						;right button}

  mbt=0 & if mbutton eq 4 then mbt=4
  while mbt ne 4 do begin			;{left/middle button
    mbt=0 & cursor,x0,y0,/nowait	;to figure out if there was a drag
    repeat begin			;to print out coordinates
      cursor,x0,y0,/change,/data
      if idlvers lt 5 then mbt=!err else mbt=!MOUSE.button
      print,form='($,g10.4,g10.4,a)',x0,y0,string("15b)
    endrep until (mbt ne 0)
    cursor,x1,y1,/up,/data	;at the button-up transition
    drag=abs(x1-x0)+abs(y1-y0)
    case mbt of
      1: begin				;left button
	xp(ohy)=x1 & yp(ohy)=y1	;reposition
	lablplot,x,y,wvl,label,xp,yp,sty,psym=10,hylyt=ohy+1,$
	  xtitle='!3['+string(byte(197))+']',_extra=e	;replot
      end
      2: begin				;middle button
	if drag eq 0 then begin		;click
	  sty(ohy)=sty(ohy)+1 & if sty(ohy) gt maxsty then sty(ohy)=0
				;cycle through styles
	  lablplot,x,y,wvl,label,xp,yp,sty,psym=10,hylyt=ohy+1,$
	    xtitle='!3['+string(byte(197))+']',_extra=e	;replot
	endif else begin		;drag
	  isty=wmenu(stylist,init=sty(ohy))
	  sty(ohy)=isty		;select style
	  lablplot,x,y,wvl,label,xp,yp,sty,psym=10,hylyt=ohy+1,$
	    xtitle='!3['+string(byte(197))+']',_extra=e	;replot
	endelse
      end
      else: print,form='($,21a)',replicate(' ',20),string("15b)
    endcase
  endwhile					;left/middle button}

  lablplot,x,y,wvl,label,xp,yp,sty,psym=10,xtitle='!3['+string(byte(197))+']',$
	_extra=e	;replot
endwhile					;go_on=1}

xpos=xp & ypos=yp & style=sty

return
end
