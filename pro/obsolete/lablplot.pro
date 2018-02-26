pro lablplot,x,y,z,label,xpos,ypos,style,hylyt=hylyt,help=help,$
	room=room,drop=drop,bars=bars,gaps=gaps, _extra=e
;+
;procedure	lablplot
;	plots labels of the appropriate style at the specified positions
;	on a spectrum
;
;planned obsolescence
;	(to be) superseded by DARSHANA
;
;syntax
;	lablplot,x,y,z,label,xpos,ypos,style,hylyt=hylyt,help=help,$
;	room=room,drop=drop,bars=bars,gaps=gaps, PLOT_KEYWORDS
;
;parameters
;	x	[INPUT; required] x-axis values
;	y	[INPUT; required] y-axis values
;	z	[INPUT; required] x-axis positions labels are identified with
;	label	[INPUT; required] said labels
;		* if there are multiple fields (see STYLE), each field must
;		  be separated by "|"
;	xpos	[INPUT; required] x-axis positions of line labels
;	ypos	[INPUT; required] y-axis positions of line labels
;	style	[INPUT; required] style of labeling (see array HLP below)
;
;keywords
;	hylyt	[INPUT] if set, "highlights" LABEL(HYLYT-1)
;	help	[INPUT] if set, prints out usage and STYLE help
;
;	room	[INPUT; default=1.] how much room to provide for the labels
;		as a fraction of the dynamic range in Y
;	drop	[INPUT; default=1.5] spacing in the case of multiple fields
;	bars	[INPUT; default=0.25] length of little vertical bars for
;		some special styles
;	gaps	[INPUT; default=0.1] gap between Y(Z) and line to LABEL(Z)
;
;	_extra	[INPUT] allows this routine to act as wrapper to PLOT
;
;restrictions
;	* requires graphical capability
;	* best called from PLOTID
;	* ALL parameters are needed
;	* X and Y must be of the same size, as must Z,LABEL,XPOS,YPOS, & STYLE
;
;history
;	vinay kashyap (Jan97)
;	added styles 9-14 (VK; Jan01)
;-

message,'OBSOLETE!',/informational

;	help
hlp=['Usage: lablplot,x,y,z,label,xpos,ypos,style,hylyt=hylyt,/help$',$
  '        room=room,drop=drop,bars=bars,gaps=gaps',$
  '  plots y(x) with appropriate labels label(z) at (xpos,ypos)',$
  '',$
  'label fields may be arranged',$
  'horizontally (H), vertically (V), or horizontal, but turned on side (I)',$
  '',$
  'lines pointing the label to the appropriate position may be drawn as',$
  '-- little vertical lines (|)',$
  '-- no line at all (.)',$
  '-- straight lines (/)',$
  '-- little vertical lines plus straight lines (/|)',$
  '-- inversion of above, half-Ys (Y)',$
  '-- straight up, then left or right (T)',$
  '-- sideways, then straight up (_L)',$
  '',$
  'only the following are encoded:',$
  '	STYLE=0: I|	STYLE=1: I/	STYLE=2: I/|',$
  '	STYLE=3: H/	STYLE=4: HT	STYLE=5: H_L',$
  '	STYLE=6: V/	STYLE=7: VT	STYLE=8: V_L',$
  '	STYLE=9: I/| (with short /)',$
  '	STYLE=10: HY	STYLE=11: VY',$
  '	STYLE=12: I.	STYLE=13: H.	STYLE=14+: V.']
if keyword_set(help) then for i=0,n_elements(hlp)-1 do print,hlp(i)

;	usage
if n_params() lt 7 then begin
  if not keyword_set(help) then begin
    print,hlp(0) & print,hlp(1) & print,hlp(2)
  endif
  return
endif

;	check inputs for consistency
nope=0
nx=n_elements(x) & ny=n_elements(y)
nz=n_elements(z) & nl=n_elements(label)
mx=n_elements(xpos) & my=n_elements(ypos) & ns=n_elements(style)
if (nx eq 0) or (ny eq 0) then nope=nope+1
if (nx ne ny) then nope=nope+2
if nz eq 0 or nl eq 0 or mx eq 0 or my eq 0 or ns eq 0 then nope=nope+4
if (nz gt nx) then nope=nope+8
if (nz ne nl) or (nz ne mx) or (nz ne my) or (nz ne ns) then nope=nope+16
if (nl ne mx) or (nl ne my) or (nl ne ns) then nope=nope+32
if (mx ne my) or (mx ne ns) then nope=nope+64
if (my ne ns) then nope=nope+128
if nope ne 0 then begin
  message,'	Error code: '+strtrim(nope,2),/info & return
endif

;	initialize
if n_elements(room) eq 0 then room=1.
if n_elements(drop) eq 0 then drop=1.5
if n_elements(bars) eq 0 then bars=0.25
if n_elements(gaps) eq 0 then gaps=0.1
xdrop=drop*!d.x_ch_size & ydrop=drop*!d.y_ch_size 
oldcol=!p.color & oldthk=!p.thick
xr=[min(x),max(x)] & yr=[min(y),max(y)] & dy=(yr(1)-yr(0))
yr=[yr(0),yr(1)+room*dy]

;	plot
plot,x,y,psym=10,xr=xr,yr=yr,/xs,/ys,_extra=e

;	step through the labels and put 'em on da plot
for i=0,nl-1 do begin
  c1=label(i) & cc=str_sep(c1,'|') & ncc=n_elements(cc)	;multiple fields?
  sty=style(i)						;label style

  ;	where to position the label?
  x0=z(i) & x1=xpos(i) & oo=where(abs(x-x0) eq min(abs(x-x0)))
  y1=ypos(i) & y0=y(oo(0))+gaps*dy

  ;	in device coords, these would be
  tmp=convert_coord([x0,x1],[y0,y1],/data,/to_dev)
  x0d=tmp(0,0) & x1d=tmp(0,1) & y0d=tmp(1,0) & y1d=tmp(1,1)

  ;	how many positioning commands/field?
  tmp=str_sep(c1,'!') & upad=n_elements(tmp)-1 & upad=upad*2/ncc

  ;	position of the underline
  ;in horizontal mode
  hxa=x1d-(strlen(c1)-upad*ncc)*!d.x_ch_size/2.
  hxb=x1d+(strlen(c1)-upad*ncc)*!d.x_ch_size/2.
  hya=y1d & hyb=hya+0.4*!d.y_ch_size
  ;in vertical mode
  vxa=x1d-max(strlen(cc)-upad)*!d.x_ch_size/2.
  vxb=x1d+max(strlen(cc)-upad)*!d.x_ch_size/2.
  vya=y1d-ncc*ydrop & vyb=vya+0.4*!d.y_ch_size
  ;in sideways mode
  sxa=x1d-(ncc+1)*xdrop/2. & sxb=x1d+(ncc+1)*xdrop/2.
  sya=y1d & syb=sya+0.4*!d.y_ch_size

  if keyword_set(hylyt) then begin
    if hylyt-1 eq i then begin
      !p.color=oldcol*0.8 & !p.thick=oldthk*1.5		;highlight
    endif
  endif

  case sty of
    0: begin							;I| (0)
      x1=x0 & x1d=x0d				;override x-pos offsets
      ya=y1-bars*dy>y0 & yb=y1-gaps*dy>ya
      oplot,[x0,x0],[ya,yb]			;vertical bar
      for j=0,ncc-1 do xyouts,x1d+j*xdrop,y1d,cc(j),$
	orientation=90,/dev			;label
    end
    1: begin							;I/
      ya=y0 & yb=y1
      oplot,[x0,x1],[ya,yb]			;straight line
      for j=0,ncc-1 do xyouts,sxa+(j+1.5)*xdrop,sya+0.2*abs(ydrop),cc(j),$
	orientation=90,/dev			;label
      plots,[sxa,sxa,sxb,sxb],[syb,sya,sya,syb],/dev	;underline
    end
    2: begin							;I/|
      ya=y1-bars*dy>y0 & yb=y1
      oplot,[x1,x1],[ya,yb]			;vertical bar
      oplot,[x0,x1],[y0,ya]			;straight line
      for j=0,ncc-1 do xyouts,sxa+(j+1.5)*xdrop,sya+0.2*abs(ydrop),cc(j),$
	orientation=90,/dev			;label
      plots,[sxa,sxa,sxb,sxb],[syb,sya,sya,syb],/dev	;underline
    end
    3: begin							;H/
      oplot,[x0,x1],[y0,y1]			;straight line
      c2='' & for j=0,ncc-2 do c2=c2+cc(j)+' ' & c2=c2+cc(ncc-1)
      xyouts,x1d,y1d+0.5*abs(ydrop),c2,align=0.5,/dev	;label
      plots,[hxa,hxa,hxb,hxb],[hyb,hya,hya,hyb],/dev	;underline
    end
    4: begin							;HT
      oplot,[x0,x0],[y0,y1]			;vertical line
      oplot,[x0,x1],[y1,y1]			;horizontal line
      c2='' & for j=0,ncc-2 do c2=c2+cc(j)+' ' & c2=c2+cc(ncc-1)
      xx=x1d+!d.x_ch_size & if x1 lt x0 then xx=x1d-(hxb-hxa)
      xyouts,xx,y1d-0.5*!d.y_ch_size,c2,/dev	;label
    end
    5: begin							;H_L
      oplot,[x0,x1],[y0,y0]			;horizontal line
      oplot,[x1,x1],[y0,y1]			;vertical line
      c2='' & for j=0,ncc-2 do c2=c2+cc(j)+' ' & c2=c2+cc(ncc-1)
      xyouts,x1d,y1d+0.5*abs(ydrop),c2,align=0.5,/dev	;label
      plots,[hxa,hxa,hxb,hxb],[hyb,hya,hya,hyb],/dev	;underline
    end
    6: begin							;V/
      plots,[x0d,0.5*(vxa+vxb)],[y0d,vya],/dev	;straight line
      for j=0,ncc-1 do xyouts,x1d,y1d-j*ydrop,cc(j),/dev,align=0.5	;label
      plots,[vxa,vxa,vxb,vxb],[vyb,vya,vya,vyb],/dev	;underline
    end
    7: begin							;VT
      oplot,[x0,x0],[y0,y1]				;vertical line
      oplot,[x0,x1],[y1,y1]				;horizontal line
      xx=x1d+!d.x_ch_size & if x1 lt x0 then xx=x1d-(vxb-vxa)
      for j=0,ncc-1 do xyouts,xx,y1d-j*ydrop,cc(j),/dev	;label
    end
    8: begin							;V_L
      oplot,[x0,x1],[y0,y0]				;horizontal line
      oplot,[x1,x1],[y0,y1]				;vertical line
      for j=0,ncc-1 do xyouts,x1d,y1d-j*ydrop,cc(j),/dev,align=-0.1	;label
    end
    9: begin							;I/|
      ya=y1-bars*dy>y0 & yb=y1
      oplot,[x1,x1],[ya,yb]				;vertical bar
      oplot,[x0,x1],[ya-bars*dy>y0,ya]			;straight line
      for j=0,ncc-1 do xyouts,sxa+(j+1.5)*xdrop,sya+0.2*abs(ydrop),cc(j),$
	orientation=90,/dev				;label
      plots,[sxa,sxa,sxb,sxb],[syb,sya,sya,syb],/dev	;underline
    end
    10: begin						;HY
      ya=y1-bars*dy>y0 & yb=y1
      oplot,[x0,x0],[y0,ya]				;vertical line
      oplot,[x0,x1],[ya,yb]				;slanted line
      c2='' & for j=0,ncc-2 do c2=c2+cc(j)+' ' & c2=c2+cc(ncc-1)
      xyouts,x1d,y1d+0.5*abs(ydrop),c2,align=0.5,/dev	;label
    end
    11: begin						;VY
      ya=y1-bars*dy>y0 & yb=y1
      oplot,[x0,x0],[y0,ya]				;vertical line
      oplot,[x0,x1],[ya,yb]				;slanted line
      xx=x1d+!d.x_ch_size & if x1 lt x0 then xx=x1d-(vxb-vxa)
      for j=0,ncc-1 do xyouts,xx,y1d-j*ydrop+ydrop,cc(j),$
	/dev,orientation=90	;label
    end
    12: begin						;I.
      for j=0,ncc-1 do xyouts,x0d+j*xdrop,y1d,cc(j),orientation=90,/dev	;label
    end
    13: begin						;H.
      c2='' & for j=0,ncc-2 do c2=c2+cc(j)+' ' & c2=c2+cc(ncc-1)
      xyouts,x0d,y1d+0.5*abs(ydrop),c2,align=0.5,/dev	;label
    end
    else: begin						;V.
      xx=x0d+!d.x_ch_size & if x1 lt x0 then xx=x0d-(vxb-vxa)
      for j=0,ncc-1 do xyouts,x0d-j*xdrop,y1d,cc(j),/dev,$
	orientation=90,align=0.5	;label
    end
  endcase

  if keyword_set(hylyt) then begin
    !p.color=oldcol & !p.thick=oldthk			;reset
  endif
endfor

return
end
