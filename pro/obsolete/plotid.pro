pro plotid,line,lamda,spec,xpos=xpos,ypos=ypos,style=style,finger=finger,$
	nolam=nolam,desig=desig,_extra=e
;+
;procedure	plotid
;	easily said: uses the output of LINEID to generate publication
;	quality graph with all the bells and whistles.
;
;obsolescence
;	(to be) supplanted by KALPANA
;
;syntax
;	plotid,line,lamda,spec,xpos=xpos,ypos=ypos,style=style,$
;	/finger,/nolam,/desig, PLOT_KEYWORDS
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
;	xpos	[I/O] x-axis positions of line labels
;	ypos	[I/O] y-axis positions of line labels
;		* if supplied (and consistent with LAMDA/SPEC),
;		  overrides positions calculated in situ
;		* use this to recreate graphs generated interactively!
;	style	[I/O] style of labeling, integer array of same size as [XY]POS
;		* type "lablplot,/help" for a description
;	finger	[INPUT] if set, allows interactive placement of line labels
;	nolam	[INPUT] if set, labels do NOT include wavelength matches
;	desig	[INPUT] if set, labels WILL include transition levels
;	_extra	[INPUT] allows this routine to act as wrapper to PLOT
;
;restrictions
;	* requires heavy graphical capability
;	* LINEID must have been run previously
;	* uses subroutines LABLPLOT, WMENU
;	* to be superseded by KALPANA
;
;history
;	vinay kashyap (Jan97)
;	added IDL5 features (VK; Jan98)
;	hid IDL4 calling sequence 'feature' (VK; 9/9/99)
;	more plotting styles added to LABLPLOT (VK; Jan01)
;-

message,'OBSOLETE -- use KALPANA if available',/info

;	usage
if n_params() lt 3 then begin
  print,'Usage: plotid,line_id,lamda,spec,xpos=xp,ypos=yp,style=sty,/finger'
  print,'  plots spectrum with overlaid line IDs'
  return
endif

;	IDL5 compatibility line
forward_function arg_present

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
rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[rom,'X'+rom,'XX'+rom,'XXX'+rom,'']	;roman numerals from 1-40, plus unknown
atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe',$
	'Co','Ni','Cu','Zn','?']	;elements from 1-30, plus unknown
stylist=['0: sideways+little vertical line',$
	 '1: sideways+slanted line',$
	 '2: sideways+inverted half-Y',$
	 '3: horizontal+slanted line',$
	 '4: horizontal+half-T',$
	 '5: horizontal+inverted half-T',$
	 '6: vertical+slanted line',$
	 '7: vertical+half-T',$
	 '8: vertical+inverted half-T',$
	 '9: sideways+inverted half-Y',$
	 '10: horizontal+half-Y',$
	 '11: vertical+half-Y',$
	 '12: sideways+no line',$
	 '13: horizontal+no line',$
	 '14: vertical+no line'] & maxsty=n_elements(stylist)-1

;	check keywords
if keyword_set(xpos) then xp=xpos
if keyword_set(ypos) then yp=ypos
if n_elements(xp) ne nid then xp=wvl
if n_elements(yp) ne nid then yp=fltarr(nid)+1.3*(max(y)-min(y))
;
if keyword_set(style) then sty=fix(style)
if n_elements(sty) ne nid then sty=intarr(nid)

;	generate labels
for i=1,nid do begin
  ;	unpack LINE
  z=line.(i).z & ion=line.(i).ion & labl=line.(i).labl
  if labl(0) eq 'Unknown' then begin
    labl=['Un','known'] & z=n_elements(atom) & ion=n_elements(rom)
  endif
  ll=strtrim(string(abs(line.(i).wvl),'(f8.2)'),2)
  oo=where(line.(i).wvl lt 0) & if oo(0) ne -1 then ll(oo)=ll(oo)+'*'
	;the '*' denotes theoretical wavelength (cf. CHIANTI)
  lid=n_elements(z) & cc=''
  for j=0,lid-1 do begin
    cc=cc+atom(z(j)-1)+rom(ion(j)-1)	;ELEMENT IONIC_STATE
    if keyword_set(desig) then begin	;(level designation
      case desig of
      1: begin				;only level designations
		;for no reason whatsoever, I'm putting a "(" here
	      c1=strtrim((str_sep(labl(1,j),')'))(1),2) & lc1=strlen(c1)
	      cc=cc+' [!U'+strmid(c1,0,1)+'!N'+strmid(c1,1,1)+$
		 '!D'+strmid(c1,3,lc1-3)+'!N -> !U'
		;for no reason whatsoever, I'm putting a "(" here
	      c1=strtrim((str_sep(labl(0,j),')'))(1),2) & lc1=strlen(c1)
	      cc=cc+strmid(c1,0,1)+'!N'+strmid(c1,1,1)+'!D'+$
		 strmid(c1,3,lc1-3)+'!N]'
	;c1=strtrim(((str_sep(labl(1,j),')'))(1),2) & lc1=strlen(c1)
	;cc=cc+' [!U'+strmid(c1,0,1)+'!N'+strmid(c1,1,1)+$
	;	'!D'+strmid(c1,3,lc1-3)+'!N -> !U'
	;c1=strtrim(((str_sep(labl(0,j),')'))(1),2) & lc1=strlen(c1)
	;cc=cc+strmid(c1,0,1)+'!N'+strmid(c1,1,1)+'!D'+$
	;	strmid(c1,3,lc1-3)+'!N]'
      end
      else: cc=cc+' ['+labl(1,j)+' -> '+labl(0,j)+']'	;everything
      endcase
    endif				;DESIG)
    if not keyword_set(nolam) then cc=cc+' !4k!3'+ll(j)	;matched wavelength
    cc=cc+'|'				;field seperator
  endfor
  label(i-1)=strmid(cc,0,strlen(cc)-1)	;remove the last "|"
endfor

;	plot
lablplot,x,y,wvl,label,xp,yp,sty,psym=10,/xstyle,/ystyle,$
	xtitle='!3['+string(byte(197))+']',_extra=e

;	are we done?  then quit!
nope=0
if not keyword_set(finger) then nope=1
if !d.name ne 'X' then nope=2
if nope ne 0 then begin
  xpos=xp & ypos=yp & style=sty & return
endif

;	interactive stuff
go_on=1 & storxp=xp & storyp=yp & storsty=sty

;	but no point in continuing if output will not be saved, eh?
if float(!version.release) lt 5 then begin
  if not keyword_set(xpos) or not keyword_set(ypos) or not keyword_set(style) then begin
    help,xpos,ypos,style
    message,'Check calling sequence! Type .CON to continue',/info
    stop
  endif
endif
if float(!version.release) ge 5. then begin
  if not arg_present(xpos) then begin
    message,'XPOS will not be returned',/info & nope=3
  endif
  if not arg_present(ypos) then begin
    message,'YPOS will not be returned',/info & nope=4
  endif
  if not arg_present(style) then begin
    message,'STYLE will not be returned',/info & nope=5
  endif
  if nope ge 3 then message,'type .SKIP and .CON to continue'
endif

print,'		Use mouse buttons to select label'
print,'LEFT/MIDDLE allows repositioning/restyling; RIGHT allows reset/quit'
print,'	L/M: left click or drag to position label'
print,'	     middle click to cycle style, drag to select style'
print,'	     right click or drag to go to next label'
print,'	R: left click to reset initial position, left drag to zero x-offset'
print,'	   middle click or drag to reset initial label style'
print,'	   right click or drag to QUIT'

while go_on eq 1 do begin			;{reposition/restyle labels
  mbutton=0 & ohy=0L

  repeat begin				;(select label by clicking
    cursor,xc,yc,/change,/data & mbutton=!err
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
    cursor,x0,y0,/down,/data & mbt=!err & cursor,x1,y1,/up,/data
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
      cursor,x0,y0,/change,/data & mbt=!err
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
