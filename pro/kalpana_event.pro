;+
;procedure	kalpana_event
;	widget event handler subroutine for KALPANA.  see that routine
;	for a description of variables, etc.  only rudimentary consistency
;	checks are carried out here.
;	uses subroutines PICKLABL and MOVELABL, which are included in
;	this file.
;syntax
;	kalpana_event,x,y,ststr,widg,oststr=oststr,wshow=wshow,wwork=wwork,$
;		outroot=outroot
;parameters
;	X,Y,STSTR,WIDG	[INPUT; required]
;keywords
;	OSTSTR,WSHOW,WWORK	[INPUT]
;	OUTROOT	[INPUT] filename root for output
;
;procedure	PICKLABL
;	interactively select label indices
;usage
;	picklabl,x,y,xloc,xlab,ylab,ilab
;parameters
;	X	[INPUT] as in KALPANA
;	Y	[INPUT] as in KALPANA
;	XLOC	[INPUT] STSTR.LOC.X
;	XLAB	[INPUT] STSTR.L#.POS(0)
;	YLAB	[INPUT] STSTR.L#.POS(1)
;	ILAB	[OUTPUT] the "#"-1 in STSTR.L#
;
;procedure	MOVELABL
;	interactively move label locations and connecting lines
;usage
;	movelabl,x,y,ststr,pnum,lnum,coarse=coarse,dcor=dcor,springy=springy
;parameters
;	X	[INPUT] as in KALPANA
;	Y	[INPUT] as in KALPANA
;	STSTR	[I/O] as in KALPANA
;	PNUM	[INPUT] plot window
;	LNUM	[I/O] selected label -- can be changed via keyboard control
;keywords
;	COARSE	[INPUT] if set, forces intermediate points to line up
;	DCOR	[INPUT; default=0.1] "fineness" of the coarseness
;	SPRINGY	[I/O] if set, keeps XPATH/YPATH relative positions as is
;		if end positions are changed
;		* if <0, then moves whole system bodily
;
;subroutines
;	PICKLABL
;	MOVELABL
;
;history
;	vinay kashyap (SepMIM)
;	handle color tables in 24-bit displays (VK; Jun02)
;	changed call from STR2ARR to STR_2_ARR (VK; Apr05)
;	button press status now stored in !MOUSE, not !ERR (VK; Apr09)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;{
pro picklabl,x,y,xloc,xlab,ylab,ilab
;(plus)
;procedure	picklabl
;		see above for description
;(minus)

;	move cursor to pick out a label, left click to select,
;	middle click to deselect, right click to exit

;	initialize
ok='ok' & ilab=-1L
ymx=!y.crange(1) & ymn=!y.crange(0) & ymark=ymx-0.05*(ymx-ymn)
idlvers=float(!VERSION.release) & if idlvers le 0 then idlvers=7

print,'	Left to Select, Middle to Deselect, Right to Exit'
while ok eq 'ok' do begin		;{endlessly loop
  mbutton=0 & ojlab=0L

  repeat begin			;{select label by clicking
    cursor,xc,yc,/change,/data
    if idlvers lt 5 then mbutton=!err else mbutton=!MOUSE.button
    tmp=min(sqrt((xc-xlab)^2+(yc-ylab)^2),jlab)
    if jlab ne ojlab then begin
      plots,xlab(ojlab),ylab(ojlab),color=0,psym=4,symsize=3,/data
      xyouts,xlab(ojlab),ymark,'!95!3',charsize=3,color=0,align=0.5
      ojlab=jlab
    endif
    plots,xlab(jlab),ylab(jlab),color=!d.n_colors-1,psym=4,symsize=3,/data
    xyouts,xlab(ojlab),ymark,'!95!3',charsize=3,color=!d.n_colors-1,align=0.5
  endrep until (mbutton ne 0)	;label selected}

  case mbutton of
    1: begin				;left
      if ilab(0) eq -1 then ilab=[jlab] else ilab=[ilab,jlab]
      plots,xlab(jlab),ylab(jlab),color=!d.n_colors-1,psym=2,symsize=3,/data
    end
    2: begin				;middle
      if ilab(0) ne -1 then begin
	oo=where(ilab ne jlab,moo)
	if moo ne 0 then ilab=ilab(oo) else ilab=-1L
        plots,xlab(jlab),ylab(jlab),color=0,psym=2,symsize=3,/data
      endif
    end
    4: ok='quit'			;right
  endcase

endwhile				;OK}

;	return only the unique elements
ilab0=ilab(0) & oi=where(ilab ne ilab0,moi)
if moi gt 0 then begin
  ilab=ilab(oi) & ilab=ilab(uniq(ilab,sort(ilab))) & ilab=[ilab0,ilab]
endif else ilab=[ilab0]

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;}{
pro movelabl,x,y,ststr,pnum,lnum,coarse=coarse,dcor=dcor,springy=springy,$
	ytop=ytop,ybot=ybot
;(plus)
;procedure	movelabl
;		see above for description
;(minus)

;	error check
if n_tags(ststr) eq 0 then return
ncol=(ststr.WINDOW.MULTI)(1) > 1
nrow=(ststr.WINDOW.MULTI)(2) > 1
nplot=n_elements(ststr.WINDOW.XMIN)
if pnum lt 0 or pnum ge nplot then return
nlabel=n_elements(ststr.LOC.X)
if lnum lt 0 or lnum ge nlabel then return

;	initialize
labcolor=ststr.(lnum+2).LABCOLOR
ymx=!y.crange(1) & ymn=!y.crange(0) & ymark=ymx-0.05*(ymx-ymn)
if labcolor ne !d.n_colors-1 then ststr.(lnum+2).LABCOLOR=!d.n_colors-1 else $
  ststr.(lnum+2).LABCOLOR=labcolor-1
wstr=ststr.WINDOW & lstr=ststr.LOC
xloc=(lstr.X)(lnum) & yloc=(lstr.Y)(lnum)
if yloc eq 0 then begin
  tmp=min(abs(xloc-x),iy) & yloc=y(iy)
  lstr.Y(lnum)=yloc
endif
istr=ststr.(lnum+2) & pos=istr.POS & xpath=istr.XPATH & ypath=istr.YPATH
posx=pos(0) & posy=pos(1)
ymn=wstr.YMIN(pnum) & ymx=wstr.YMAX(pnum)
if posy le ymn then posy=yloc+0.1*(ymx-yloc)
xx=[posx,xpath,xloc] & yy=[posy,ypath,yloc] & istr.POS=[posx,posy]
nfld=n_tags(istr) & namfld=tag_names(istr)
;
opos=pos & oxx=xx & oyy=yy & oxpath=xpath & oypath=ypath & oyloc=yloc
;
if not keyword_set(yfrac) then yfrac=0.1
if not keyword_set(dcor) then coarseness=0.2 else coarseness=float(dcor(0))
if not keyword_set(springy) then springy=0
;
;menu option description for CW_PDMENU
menu_opts=['oops','QUIT',$
	'Label Styles','hang/top','hang/top/slide',$
	  'vert/top','vert/top/slide','flat/top','flat/top/slide',$
	  'straight line','half-T','inverted half-T','half-Y',$
	  'inverted half-Y','drop+half-Y',$
	'Rigidity','Springy','Hingy','Sticky',$
	'Reorient','0','90','+90','-90',$
	'Realign','Left','Right','Up','Down',$
	'New Label','Next','Previous','Keyboard']
menu_type=[ 0,0,	1,0,0, 0,0,0,0, 0,0,0,0, 0,2,	1,0,0,2,$
	1,0,0,0,2,	1,0,0,0,2,	3,0,0,2]
menu_desc=strtrim(menu_type,2)+'\'+menu_opts

print,'	Left Click+Drag to reposition, Middle Click+Drag to undo'
print,'	Left Click to add point, Middle Click to delete point'
print,'	Right Click to exit, Right Click+Drag for more options'

ok='ok'
while ok eq 'ok' do begin		;{endlessly loop
  xyouts,xloc,ymark,'!95!3',charsize=3,color=!d.n_colors-1,align=0.5
  mbutton=0 & ojlab=0L
  cursor,x0,y0,/down,/data
  if idlvers lt 5 then mbutton=!err else mbutton=!MOUSE.button
  cursor,x1,y1,/up,/data
					;click and/or drag
  cdev=convert_coord([x0,x1],[y0,y1],/data,/to_dev)
  drag=long(sqrt((cdev(0,0)-cdev(0,1))^2+(cdev(1,0)-cdev(1,1))^2))
  					;if drag = 0, it's a click
  minx=min(!X.CRANGE) & maxx=max(!X.CRANGE)
  miny=min(!Y.CRANGE) & maxy=max(!Y.CRANGE)
  if x0 lt minx or x0 gt maxx then mbutton=0	;ignore if outside range
  if y0 lt miny or y0 gt maxy then mbutton=0	;ignore if outside range

  if mbutton eq 4 then begin			;(right
    if drag eq 0 then ok='quit' else begin
      base_menu=widget_base(title='MENU')
      menu_drop=cw_pdmenu(base_menu,menu_desc,/return_name)
      widget_control,base_menu,/realize
      event=widget_event(base_menu) & val=strlowcase(event.value)
      widget_control,base_menu,/destroy
      pstyle=0. & newlnum=lnum
      case val of
	'oops':	;nothing
	'quit': ok='quit'
	'hang/top': pstyle=1.+yfrac
	'hang/top/slide': pstyle=2.+yfrac
	'vert/top': pstyle=3.+yfrac
	'vert/top/slide': pstyle=4.+yfrac
	'flat/top': pstyle=5.+yfrac
	'flat/top/slide': pstyle=6.+yfrac
	'straight line': pstyle=-1
	'half-t': pstyle=-2
	'inverted half-t': pstyle=-3
	'half-y': pstyle=-4-yfrac
	'inverted half-y': pstyle=-5-yfrac
	'drop+half-y': pstyle=-6-yfrac
	'springy': springy=1
	'hingy': springy=0
	'sticky': springy=-1
	'next': newlnum=lnum+1L
	'previous': newlnum=lnum-1L
	'keyboard': begin
	  c1='' & print,prompt='type label index ['+strtrim(lnum,2)+']',c1
	  if strtrim(c1,2) ne '' then newlnum=long(c1)
	end
	'0': istr.ORIENT=0.
	'90': istr.ORIENT=90.
	'+90': istr.ORIENT=istr.ORIENT+90.
	'-90': istr.ORIENT=istr.ORIENT-90.
	'left': istr.ALIGN='LEFT'
	'right': istr.ALIGN='RIGHT'
	'up': istr.ALIGN='UP'
	'down': istr.ALIGN='DOWN'
	else: message,'Option not yet implemented: '+val,/info
      endcase
      if pstyle ne 0 then begin
	rangoli,ststr,pstyle,x,y,pnum=pnum,lnum=lnum,ytop=ytop,ybot=ybot
	reinit=1
      endif
      if newlnum ne lnum then begin
	;	reset
	ststr.(lnum+2).LABCOLOR=labcolor
	lnum=newlnum mod nlabel > 0
	labcolor=ststr.(lnum+2).LABCOLOR
	if labcolor ne !d.n_colors-1 then ststr.(lnum+2).LABCOLOR=$
		!d.n_colors-1 else ststr.(lnum+2).LABCOLOR=labcolor-1
	reinit=1
      endif
      if keyword_set(reinit) then begin
	reinit=0
	wstr=ststr.WINDOW & lstr=ststr.LOC
	xloc=(lstr.X)(lnum) & yloc=(lstr.Y)(lnum)
	if yloc eq 0 then begin
	  tmp=min(abs(xloc-x),iy) & yloc=y(iy)
	  lstr.Y(lnum)=yloc
	endif
	istr=ststr.(lnum+2)
	pos=istr.POS & xpath=istr.XPATH & ypath=istr.YPATH
	posx=pos(0) & posy=pos(1)
	xx=[posx,xpath,xloc] & yy=[posy,ypath,yloc] & istr.POS=[posx,posy]
	nfld=n_tags(istr) & namfld=tag_names(istr)
	opos=pos & oxx=xx & oyy=yy & oxpath=xpath & oypath=ypath & oyloc=yloc
      endif
    endelse
  endif						;right)

  if drag eq 0 then begin	;(click
    if mbutton eq 1 then begin			;(left
      ;	add this point to the PATHs
      if keyword_set(coarse) then begin		;(line it up
        xmin=min(xx,max=xmax) & ymin=min(yy,max=ymax)
	delx=(xmax-xmin) & dely=(ymax-ymin)
	offx=min(abs(x0-xx),ix) & offy=min(abs(y0-yy),iy)
	if offx/delx lt coarseness then x0=xx(ix)
	if offy/dely lt coarseness then y0=yy(iy)
      endif					;COARSE)
      xpath=[xpath,x0] & ypath=[ypath,y0]	;append
    endif					;left)
    if mbutton eq 2 then begin			;(middle
      tmp=min(sqrt((x0-xpath)^2+(y0-ypath)^2),ipt)
      npath=n_elements(xpath)
      if npath gt 1 then begin
	oi=where(lindgen(npath) ne ipt,moi)
	xpath=xpath(oi) & ypath=ypath(oi)	;delete
      endif else print,'There is nothing left to delete!'
    endif					;middle)
    for ifld=0L,nfld-1L do begin	;{now restructure ISTR
      if ifld eq 0 then jstr=create_struct(namfld(ifld),istr.(ifld)) else begin
	if namfld(ifld) eq 'XPATH' then jstr=$
	  create_struct(jstr,'XPATH',xpath) else $
	    if namfld(ifld) eq 'YPATH' then jstr=$
	      create_struct(jstr,'YPATH',ypath) else $
	        jstr=create_struct(jstr,namfld(ifld),istr.(ifld))
      endelse
    endfor				;IFLD=0,NFLD-1}
    istr=jstr
    xx=[pos(0),xpath,xloc] & yy=[pos(1),ypath,yloc]
  endif else begin		;click)(click+drag
    if mbutton eq 1 then begin		;(left
      tmp=min(sqrt((x0-xx)^2+(y0-yy)^2),ipt)
      nn=n_elements(xx)
      opos=pos & oxx=xx & oyy=yy & oxpath=xpath & oypath=ypath
      case ipt of		;{which point?
	0: begin		;label POSition
	  if keyword_set(coarse) then begin
            xmin=min(xx,max=xmax) & ymin=min(yy,max=ymax)
	    if abs((x1-lstr.X(lnum))/(xmax-xmin)) lt coarseness then $
		x1=lstr.X(lnum)
	  endif
	  if keyword_set(springy) then begin
	    if springy gt 0 then begin
	      xden=pos(0)-lstr.X(lnum) & yden=pos(1)-lstr.Y(lnum)
	      fxpath=(xpath-lstr.X(lnum)) & fypath=(ypath-lstr.Y(lnum))
	      if xden ne 0 then fxpath=fxpath/xden else fxpath=0*fxpath
	      if yden ne 0 then fypath=fypath/yden else fypath=0*fypath
	      xpath=fxpath*(x1-lstr.X(lnum))+lstr.X(lnum)
	      ypath=fypath*(y1-lstr.Y(lnum))+lstr.Y(lnum)
	    endif
	    if springy lt 0 then begin
	      delx=pos(0)-x1 & dely=pos(1)-y1
	      xpath=xpath-delx & ypath=ypath-dely
	      lstr.Y(lnum)=lstr.Y(lnum)-dely
	    endif
	    istr.XPATH=xpath & istr.YPATH=ypath
	  endif
	  pos(0)=x1 & pos(1)=y1
	  istr.POS=pos
	end
	nn-1: begin		;LOCation of end point
	  if keyword_set(springy) then begin
	    if springy gt 0 then begin
	      fypath=(ypath-pos(1))/(lstr.Y(lnum)-pos(1))
	      ypath=fypath*(y1-pos(1))+pos(1)
	    endif
	    if springy lt 0 then begin
	      dely=lstr.Y(lnum)-y1
	      ypath=ypath-dely & pos(1)=pos(1)-dely
	    endif
	    istr.YPATH=ypath & istr.POS=pos
	  endif
	  lstr.Y(lnum)=y1
	end
	else: begin
          if keyword_set(coarse) then begin		;(line it up
            xmin=min(xx,max=xmax) & ymin=min(yy,max=ymax)
	    delx=(xmax-xmin) & dely=(ymax-ymin)
	    offx=min(abs(x1-xx),ix) & offy=min(abs(y1-yy),iy)
	    if offx/delx lt coarseness then x1=xx(ix)
	    if offy/dely lt coarseness then y1=yy(iy)
          endif					;COARSE)
	  xpath(ipt-1L)=x1 & ypath(ipt-1L)=y1
	  istr.XPATH=xpath & istr.YPATH=ypath
	end
      endcase			;IPT}
    endif				;left)
    if mbutton eq 2 then begin		;(middle
      pos=opos & yloc=oyloc & xpath=oxpath & ypath=oypath
      lstr.Y=yloc
      istr.POS=pos
      istr.XPATH=xpath
      istr.YPATH=ypath
    endif				;middle)
  endelse			;click+drag)

  ;	now put the state structure back together and give darshana
  tstr=create_struct('WINDOW',wstr,'LOC',lstr)
  for i=0L,nlabel-1L do begin
    if i eq lnum then tstr=create_struct(tstr,'L'+strtrim(i,2),istr) else $
	tstr=create_struct(tstr,'L'+strtrim(i,2),ststr.(i+2))
  endfor
  darshana,x,y,tstr,pid=pnum

  ;	reinitialize
  ststr=tstr & wstr=tstr.WINDOW & lstr=tstr.LOC
  xloc=(lstr.X)(lnum) & yloc=(lstr.Y)(lnum)
  if yloc eq 0 then begin
    tmp=min(abs(xloc-x),iy) & yloc=y(iy)
    lstr.Y(lnum)=yloc
  endif
  istr=tstr.(lnum+2) & pos=istr.POS & xpath=istr.XPATH & ypath=istr.YPATH
  xx=[pos(0),xpath,xloc] & yy=[pos(1),ypath,yloc]
  nfld=n_tags(istr) & namfld=tag_names(istr)
  ;
  opos=pos & oxx=xx & oyy=yy & oxpath=xpath & oypath=ypath & oyloc=yloc

endwhile				;OK}

ststr.(lnum+2).LABCOLOR=labcolor

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;}

pro kalpana_event,x,y,ststr,widg,oststr=oststr,wshow=wshow,wwork=wwork,$
	outroot=outroot, _extra=e
;(plus)
;procedure	kalpana_event
;		see above for description
;(minus)

;	rudimentary checking
nx=n_elements(x) & ny=n_elements(y) & ns=n_tags(ststr) & nw=n_tags(widg)
ok='ok'
if nx eq 0 then ok='undefined input' else $
 if ny eq 0 then ok='undefined spectrum' else $
  if nx ne ny then ok='ill-defined spectrum' else $
   if ns eq 0 then ok='unknown state' else $
    if nw eq 0 then ok='cannot understand widget'
if ok ne 'ok' then begin
  print,'Usage: kalpana_event,x,y,ststr,widg'
  print,'  widget event handler for kalpana'
  message,ok,/info & return
endif

;	initialize stuff
;number of plots in window
oncol=(ststr.window.multi)(1) > 1
onrow=(ststr.window.multi)(2) > 1
onplot=n_elements(ststr.window.xmin)
ncol=oncol & nrow=onrow & nplot=onplot
zplot=lindgen(nplot)
;number of labels
nlabel=n_elements(ststr.LOC.X) & zlabel=lindgen(nlabel)
;these variables control which plot, and which label
pnum=0L & lnum=0L
;what type of output?
if not keyword_set(outroot) then outroot='root'
zlandscape=0 & zencap=0 & zcolor=0
;label arrangements
zalign=['LEFT','RIGHT','UP','DOWN'] & zarrange=['COLUMN','ROW']
zhide=['DISPLAY','HIDE']
;for UNDO's sake, save the state
savstr=ststr

;	disentangle WIDG
base=widg.base
;
base_col1=widg.col1			;	first column
k=0       & col1_menu=widg.menu(k)		;first row
  k=k+1   & menu_show=widg.menu(k)
  k=k+1   & menu_quit=widg.menu(k)
  k=k+1   & menu_save=widg.menu(k)
  k=k+1   & menu_bkup=widg.menu(k)
  k=k+1   & menu_help=widg.menu(k)
  k=k+1   & menu_prnt=widg.menu(k)
    k=k+1 & prnt_post=widg.menu(k)
    k=k+1 & prnt_gifs=widg.menu(k)
  k=k+1   & menu_file=widg.menu(k)
  k=k+1   & menu_dpar=widg.menu(k)
k=0       & col1_txts=widg.txts(k)		;second row
  k=k+1   & txts_stts=widg.txts(k)
k=0       & col1_wndo=widg.wind(k)		;third row
  k=k+1   & wndo_wind=widg.wind(k)
    k=k+1 & wind_ncol=widg.wind(k)
    k=k+1 & wind_nrow=widg.wind(k)
    k=k+1 & wind_titl=widg.wind(k)
  k=k+1   & wndo_cols=widg.wind(k)
    k=k+1 & cols_pltn=widg.wind(k)
    k=k+1 & cols_frmn=widg.wind(k)
    k=k+1 & cols_bkgn=widg.wind(k)
k=0       & col1_plot=widg.plot(k)		;fourth row
  k=k+1   & plot_row1=widg.plot(k)
  k=k+1   & row1_pnum=widg.plot(k)
  k=k+1   & row1_subt=widg.plot(k)
  k=k+1   & row1_csiz=widg.plot(k)
  k=k+1   & row1_dstr=widg.plot(k)
    k=k+1 & dstr_plot=widg.plot(k)
  k=k+1   & plot_rowx=widg.plot(k)
    k=k+1 & rowx_xtit=widg.plot(k)
    k=k+1 & rowx_xrng=widg.plot(k)
    k=k+1 & rowx_xsty=widg.plot(k)
    k=k+1 & rowx_xlog=widg.plot(k)
  k=k+1   & plot_rowy=widg.plot(k)
    k=k+1 & rowy_ytit=widg.plot(k)
    k=k+1 & rowy_yrng=widg.plot(k)
    k=k+1 & rowy_ysty=widg.plot(k)
    k=k+1 & rowy_ylog=widg.plot(k)
k=0       & col1_labl=widg.labl(k)		;fifth row
  k=k+1   & labl_row0=widg.labl(k)
    k=k+1 & row0_sele=widg.labl(k)
    k=k+1 & row0_next=widg.labl(k)
    k=k+1 & row0_prev=widg.labl(k)
    k=k+1 & row0_labl=widg.labl(k)
    k=k+1 & row0_hide=widg.labl(k)
    k=k+1 & row0_newl=widg.labl(k)
    k=k+1 & row0_dele=widg.labl(k)
      k=k+1 & dele_labl=widg.labl(k)
  k=k+1   & labl_row1=widg.labl(k)
    k=k+1 & row1_ltxt=widg.labl(k)
  k=k+1   & labl_row2=widg.labl(k)
    k=k+1 & row2_styl=widg.labl(k)
    k=k+1 & row2_algn=widg.labl(k)
    k=k+1 & row2_arng=widg.labl(k)
  k=k+1   & labl_row3=widg.labl(k)
    k=k+1 & row3_ornt=widg.labl(k)
    k=k+1 & row3_size=widg.labl(k)
    k=k+1 & row3_thck=widg.labl(k)
  k=k+1   & labl_row4=widg.labl(k)
    k=k+1 & row4_labc=widg.labl(k)
    k=k+1 & row4_labn=widg.labl(k)
    k=k+1 & row4_linc=widg.labl(k)
    k=k+1 & row4_linn=widg.labl(k)
  k=k+1   & labl_row5=widg.labl(k)
    k=k+1 & row5_ulin=widg.labl(k)
    k=k+1 & row5_slin=widg.labl(k)
;
base_col2=widg.col2			;	second column
k=0   & col2_move=widg.tweak(k)
k=k+1 & col2_undo=widg.tweak(k)
k=k+1 & col2_grps=widg.tweak(k)
  k=k+1 & grps_alle=widg.tweak(k)
  k=k+1 & grps_sele=widg.tweak(k)
  k=k+1 & grps_dele=widg.tweak(k)
  k=k+1 & grps_undo=widg.tweak(k)
k=k+1 & col2_alny=widg.tweak(k)
  k=k+1 & alny_alle=widg.tweak(k)
  k=k+1 & alny_grps=widg.tweak(k)
  k=k+1 & alny_sele=widg.tweak(k)
k=k+1 & col2_ornt=widg.tweak(k)
  k=k+1 & ornt_alle=widg.tweak(k)
  k=k+1 & ornt_grps=widg.tweak(k)
  k=k+1 & ornt_sele=widg.tweak(k)
k=k+1 & col2_lsiz=widg.tweak(k)
  k=k+1 & lsiz_alle=widg.tweak(k)
  k=k+1 & lsiz_grps=widg.tweak(k)
  k=k+1 & lsiz_sele=widg.tweak(k)
k=k+1 & col2_thck=widg.tweak(k)
  k=k+1 & thck_alle=widg.tweak(k)
  k=k+1 & thck_grps=widg.tweak(k)
  k=k+1 & thck_sele=widg.tweak(k)
k=k+1 & col2_arng=widg.tweak(k)
  k=k+1 & arng_alle=widg.tweak(k)
  k=k+1 & arng_grps=widg.tweak(k)
  k=k+1 & arng_sele=widg.tweak(k)
k=k+1 & col2_ulin=widg.tweak(k)
  k=k+1 & ulin_alle=widg.tweak(k)
  k=k+1 & ulin_grps=widg.tweak(k)
  k=k+1 & ulin_sele=widg.tweak(k)
k=k+1 & col2_slin=widg.tweak(k)
  k=k+1 & slin_alle=widg.tweak(k)
  k=k+1 & slin_grps=widg.tweak(k)
  k=k+1 & slin_sele=widg.tweak(k)
k=k+1 & col2_posn=widg.tweak(k)
k=k+1 & col2_ties=widg.tweak(k)
k=k+1 & col2_coar=widg.tweak(k)
k=k+1 & col2_ytop=widg.tweak(k)
k=k+1 & col2_ybot=widg.tweak(k)
k=k+1 & col2_yfrc=widg.tweak(k)
k_help=widg.help & nhelp=n_elements(k_help)	;help

ok=''
while ok ne 'quit' do begin			;{endless loop
  
  ;catch a falling event
  event=widget_event(base) & widget_control,event.id,get_uvalue=ev
  widget_control,txts_stts,set_value=ev

  replot=0			;flag to replot (1) or no (0)
  rewstr=0			;flag to reset STSTR.WINDOW (1) or no (0)

  ;	label location info available all the time
  xloc=ststr.LOC.X & nlabel=n_elements(xloc) & xlab=0*xloc & ylab=xlab
  for ilab=0L,nlabel-1L do begin
    xlab(ilab)=(ststr.(ilab+2).POS)(0)
    ylab(ilab)=(ststr.(ilab+2).POS)(1)
  endfor

  case ev of					;{handle the events

    ;{	first row, first column
    'menu_show': begin				;display
      widget_control,txts_stts,set_value='Displaying plot(s)'
      replot=1
    end
    'menu_quit': ok='quit'			;exit
    'menu_save': begin				;save
      widget_control,txts_stts,set_value='Overwriting original input state'
      oststr=ststr
    end
    'menu_bkup': begin				;restore
      widget_control,txts_stts,set_value='Restoring from original input state'
      ststr=oststr
      ;	and various reinitializations, such as nrow,ncol,nplot
	oncol=(ststr.window.multi)(1) > 1
	onrow=(ststr.window.multi)(2) > 1
	onplot=n_elements(ststr.window.xmin)
	ncol=oncol & nrow=onrow & nplot=onplot
	zplot=lindgen(nplot)
	nlabel=n_elements(ststr.LOC.X) & zlabel=lindgen(nlabel)
	pnum=0L & lnum=0L & updatelabl=1
    end
    'menu_help': begin				;display help
      widget_control,txts_stts,set_value=''
      if fix(strmid(!version.release,0,1)) lt 5 then $
	for i=0,nhelp-1 do print,k_help(i) else xdisplayfile,$
	'HELP',text=k_help,title='KALPANA HELP',width=70
    end
    'menu_prnt': begin				;show output
      widget_control,menu_file,set_value=outroot
    end
    'prnt_post': begin				;print postscript
      dname=!D.NAME
      if zencap eq 0 then outfile=outroot+'.ps' else outfile=outroot+'.eps'
      c1='Writing to file '+outfile+' ['
      if zcolor eq 1 then c1=c1+'color' else c1=c1+'greyscale'
      if zlandscape eq 1 then c1=c1+' landscape]' else c1=c1+' portrait]'
      widget_control,txts_stts,set_value=c1
      colmax=!d.n_colors
      set_plot,'ps' & device,file=outfile,$
	landscape=zlandscape,encapsulated=zencap,color=zcolor
      darshana,x,y,ststr,colmax=colmax, _extra=e
      device,/close & set_plot,dname
      if !D.N_COLORS gt 256 then device,decomposed=0
    end
    'prnt_gifs': begin				;print GIF
      widget_control,txts_stts,set_value='Writing to file '+outroot+'.gif'
      wset,wshow
      tvlct,r,g,b,/get
      write_gif,outroot+'.gif',tvrd(),r,g,b
      wset,wwork
    end
    'menu_file': begin				;root name of output
      widget_control,menu_file,get_value=oroot & outroot=strtrim(oroot(0),2)
      ;	search for ending .ps,.eps,.gif and remove them
      rlen=strlen(outroot)
      ips=strpos(outroot,'.ps')
      ieps=strpos(outroot,'.eps')
      igif=strpos(outroot,'.gif')
      if ips eq rlen-3 then outroot=strmid(outroot,0,rlen-3) else $
       if ieps eq rlen-4 then outroot=strmid(outroot,0,rlen-4) else $
        if igif eq rlen-4 then outroot=strmid(outroot,0,rlen-4)
      widget_control,menu_file,set_value=outroot
      widget_control,txts_stts,set_value='Output path and file root is '+outroot
    end
    'menu_dpar': begin				;DEVICE keywords
      ;	I am sure this can be tightened up
      widget_control,txts_stts,set_value='Set DEVICE keywords from keyboard'
      dpar_help=['Type E, L, C to toggle',$
	'encapsulated postscript, landscape, and color outputs',$
	'Type Q or X to return to widget']
      c1='type E,L,C, Q, or X' & print,c1
      while c1 ne 'Q' and c1 ne 'X' do begin
        c1=get_kbrd(1) & c1=strupcase(c1)
        case c1 of
	  'E': begin
	    zencap=1-zencap
	    if zencap eq 0 then print,'Set to regular postscript' else $
		print,'Set to encapsulated postscript'
	  end
	  'L': begin
	    zlandscape=1-zlandscape
	    if zlandscape eq 0 then print,'Set to portrait mode' else $
		print,'Set to landscape mode'
	  end
	  'C': begin
	    zcolor=1-zcolor
	    if zcolor eq 0 then print,'greyscale postscript' else $
		print,'color postscript'
	  end
	  'Q': print,'returning to widget control'	;nothing
	  'X': print,'returning to widget control'	;nothing
	  else: begin
	    for i=0,n_elements(dpar_help)-1 do print,dpar_help(i)
	  end
        endcase
      endwhile
    end

    ;1,1}{	second row, first column

    ;2,1}{ third row, first column

    'wind_ncol': begin				;number of columns of plots
      widget_control,wind_ncol,get_value=ncol
      widget_control,wind_nrow,get_value=nrow
      ncol = ncol > 1	;sanity check
      nrow = nrow > 1	;sanity check
      ;nplot = ncol*nrow > 1
      zplot=lindgen(nplot)
      if ncol ne oncol then begin
	rewstr=1
        widget_control,txts_stts,set_value='Changing number of plots in a page'
	oncol=ncol
      endif else widget_control,txts_stts,set_value=$
	'No change in number of plot columns'
    end
    'wind_nrow': begin				;number of rows of plots
      widget_control,wind_ncol,get_value=ncol
      widget_control,wind_nrow,get_value=nrow
      ncol = ncol > 1	;sanity check
      nrow = nrow > 1	;sanity check
      ;nplot = ncol*nrow > 1
      zplot=lindgen(nplot)
      if nrow ne onrow then begin
	rewstr=1
        widget_control,txts_stts,set_value='Changing number of plots in a page'
	onrow=nrow
      endif else widget_control,txts_stts,set_value=$
	'No change in number of plot rows'
    end
    'wind_titl': begin				;overall title
      widget_control,txts_stts,set_value=''
      widget_control,wind_titl,get_value=title
      ststr.window.title=title
      darshana,x,y,ststr,wid=wwork,pid=pnum
      replot=1
    end

    'cols_pltn': begin			;name of plot color
      widget_control,cols_pltn,get_value=pltn
      if strtrim(pltn(0),2) ne '' then setcolor,pltn(0),!D.N_COLORS-2
      widget_control,txts_stts,set_value='Setting color number '+$
	strtrim(!D.N_COLORS-2,2)+' to '+pltn(0)
      widget_control,cols_pltn,set_value=''
    end
    'cols_frmn': begin			;name of frame color
      widget_control,cols_frmn,get_value=frmn
      if strtrim(frmn(0),2) ne '' then setcolor,frmn(0),!D.N_COLORS-1
      widget_control,txts_stts,set_value='Setting color number '+$
	strtrim(!D.N_COLORS-1,2)+' to '+frmn(0)
      widget_control,cols_frmn,set_value=''
    end
    'cols_bkgn': begin			;name of background color
      widget_control,cols_bkgn,get_value=bkgn
      if strtrim(bkgn(0),2) ne '' then setcolor,bkgn(0),0
      widget_control,txts_stts,set_value='Setting color number 0 to '+bkgn(0)
      widget_control,cols_bkgn,set_value=''
    end

    ;3,1}{ fourth row, first column

    'row1_pnum': begin				;display which plot?
      widget_control,row1_pnum,get_value=pnum
      if pnum ge 0 and pnum lt n_elements(ststr.window.xmin) then begin
        widget_control,txts_stts,$
		set_value='Displaying plot window '+strtrim(pnum,2)
	;	update the rest of it
	csiz=(ststr.window.chars)(pnum)
	subt=(ststr.window.subtitle)(pnum)
	xtit=(ststr.window.xtitle)(pnum)
	ytit=(ststr.window.ytitle)(pnum)
	xlog=(ststr.window.xlog)(pnum)
	ylog=(ststr.window.ylog)(pnum)
	xmin=(ststr.window.xmin)(pnum) & xmax=(ststr.window.xmax)(pnum)
	ymin=(ststr.window.ymin)(pnum) & ymax=(ststr.window.ymax)(pnum)
	xsty=(ststr.window.xstyle)(pnum)
	ysty=(ststr.window.ystyle)(pnum)
	cxr=strtrim(xmin,2)+','+strtrim(xmax,2)
	cyr=strtrim(ymin,2)+','+strtrim(ymax,2)
	widget_control,row1_csiz,set_value=csiz
	widget_control,row1_subt,set_value=subt
	widget_control,rowx_xtit,set_value=xtit
	widget_control,rowx_xrng,set_value=cxr
	widget_control,rowx_xsty,set_value=xsty
	widget_control,rowx_xlog,set_droplist_select=xlog
	widget_control,rowy_ytit,set_value=ytit
	widget_control,rowy_yrng,set_value=cyr
	widget_control,rowy_ysty,set_value=ysty
	widget_control,rowy_ylog,set_droplist_select=ylog
	;
	darshana,x,y,ststr,wid=wwork,pid=pnum
	;replot=1
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'row1_subt': begin				;plot (sub)title
      widget_control,row1_pnum,get_value=pnum
      widget_control,row1_subt,get_value=subt
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='changing subtitle for plot '+$
		strtrim(pnum,2)
	ststr.window.SUBTITLE(pnum)=subt(0)
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
	replot=1
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'row1_csiz': begin				;size of axis labels
      widget_control,row1_pnum,get_value=pnum
      widget_control,row1_csiz,get_value=csiz
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='changing CharSize for plot '+$
		strtrim(pnum,2)
	ststr.window.CHARS(pnum)=csiz(0)
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'dstr_plot': begin				;delete this plot
      widget_control,row1_pnum,get_value=pnum
      if pnum ge 0 and pnum lt nplot then begin
	if nplot gt 1 then begin
	  widget_control,txts_stts,set_value='Destroying Window '+$
		strtrim(pnum,2)+': Please confirm at keyboard.'
	  c1='Delete window '+strtrim(pnum,2)+'? [n/y] > ' & c2=''
	  read,prompt=c1,c2
	  c2=strmid(strtrim(strupcase(c2),2),0,1)
	  if c2 eq 'Y' then begin
	    zplot=zplot(where(zplot ne pnum))
	    rewstr=1
            widget_control,txts_stts,set_value='La Plot, she is gone.'
	  endif
	endif else widget_control,txts_stts,$
		set_value='Cannot destroy last remaining window!'
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowx_xtit': begin				;x-axis title
      widget_control,row1_pnum,get_value=pnum
      widget_control,rowx_xtit,get_value=xtit
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='changing X-title for plot '+$
		strtrim(pnum,2)
	ststr.window.XTITLE(pnum)=xtit(0)
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowx_xrng': begin				;x-axis range
      widget_control,row1_pnum,get_value=pnum
      widget_control,rowx_xrng,get_value=cxr
      xlog=widget_info(rowx_xlog,/droplist_select)
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='changing X-range for plot '+$
		strtrim(pnum,2)
	xr=str_2_arr(cxr,/r4,/squish) & nxr=n_elements(xr)
	if nxr eq 1 then begin
	  dlo=abs(xr(0)-(ststr.window.xmin)(pnum))
	  dhi=abs(xr(0)-(ststr.window.xmax)(pnum))
	  if dlo lt dhi then xr=[xr(0),(ststr.window.xmax)(pnum)] else $
		xr=[(ststr.window.xmin)(pnum),xr(0)]
	endif
	if xlog eq 1 then begin
	  if xr(1) le 0 then begin
	    widget_control,txts_stts,set_value=$
		"Whatchew talkin' 'bout, Willis?"
	    xr(0)=(ststr.window.xmin)(pnum) & xr(1)=(ststr.window.xmax)(pnum)
	  endif else begin
	    if xr(0) le 0 then xr(0)=1e-5*xr(1)
	  endelse
	endif
	ststr.window.XMIN(pnum)=xr(0) & ststr.window.XMAX(pnum)=xr(1)
	cxr=strtrim(xr(0),2)+','+strtrim(xr(1),2)
	widget_control,rowx_xrng,set_value=cxr
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowx_xsty': begin				;x-axis style
      widget_control,row1_pnum,get_value=pnum
      widget_control,rowx_xsty,get_value=xsty
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='Setting X-axis style to '+$
		strtrim(xsty,2)+' for plot '+strtrim(pnum,2)
	ststr.window.XSTYLE(pnum)=xsty(0)
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowx_xlog': begin				;x-axis type
      widget_control,row1_pnum,get_value=pnum
      xlog=widget_info(rowx_xlog,/droplist_select)
      if xlog eq 0 then axistype='LINEAR' else axistype='LOG'
      if pnum ge 0 and pnum lt nplot then begin
        widget_control,rowx_xrng,get_value=cxr
	xr=str_2_arr(cxr,/r4,/squish) & nxr=n_elements(xr)
	ok='ok'
	if nxr eq 2 then begin
	  if xr(0) lt 0 or xr(1) lt 0 then ok='not OK for LOG'
	endif else begin
	  if xr(0) lt 0 then ok='not OK for LOG'
	endelse
	if axistype eq 'LINEAR' then ok='ok'
	if ok eq 'ok' then begin
	  widget_control,txts_stts,set_value='Setting X-axis type to '+$
		axistype+' for plot '+strtrim(pnum,2)
	  ststr.window.XLOG(pnum)=xlog
	  ;
          darshana,x,y,ststr,wid=wwork,pid=pnum
	endif else begin
	  widget_control,txts_stts,set_value='Invalid X-range'
	  widget_control,rowx_xlog,set_droplist_select=xlog-1
	endelse
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowy_ytit': begin				;y-axis title
      widget_control,row1_pnum,get_value=pnum
      widget_control,rowy_ytit,get_value=ytit
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='changing Y-title for plot '+$
		strtrim(pnum,2)
	ststr.window.YTITLE(pnum)=ytit(0)
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowy_yrng': begin				;y-axis range
      widget_control,row1_pnum,get_value=pnum
      widget_control,rowy_yrng,get_value=cyr
      ylog=widget_info(rowy_ylog,/droplist_select)
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='changing Y-range for plot '+$
		strtrim(pnum,2)
	yr=str_2_arr(cyr,/r4,/squish) & nyr=n_elements(yr)
	if nyr eq 1 then begin
	  dlo=abs(yr(0)-(ststr.window.ymin)(pnum))
	  dhi=abs(yr(0)-(ststr.window.ymax)(pnum))
	  if dlo lt dhi then yr=[yr(0),(ststr.window.ymax)(pnum)] else $
		yr=[(ststr.window.ymin)(pnum),yr(0)]
	endif
	if ylog eq 1 then begin
	  if yr(1) le 0 then begin
	    widget_control,txts_stts,set_value=$
		"Whatchew talkin' 'bout, Willis?"
	    yr(0)=(ststr.window.ymin)(pnum) & yr(1)=(ststr.window.ymax)(pnum)
	  endif else begin
	    if yr(0) le 0 then yr(0)=1e-5*yr(1)
	  endelse
	endif
	ststr.window.YMIN(pnum)=yr(0) & ststr.window.YMAX(pnum)=yr(1)
	cyr=strtrim(yr(0),2)+','+strtrim(yr(1),2)
	widget_control,rowy_yrng,set_value=cyr
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowy_ysty': begin				;y-axis style
      widget_control,row1_pnum,get_value=pnum
      widget_control,rowy_ysty,get_value=ysty
      if pnum ge 0 and pnum lt nplot then begin
	widget_control,txts_stts,set_value='Setting Y-axis style to '+$
		strtrim(ysty,2)+' for plot '+strtrim(pnum,2)
	ststr.window.YSTYLE(pnum)=ysty(0)
	;
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end
    'rowy_ylog': begin				;y-axis type
      widget_control,row1_pnum,get_value=pnum
      ylog=widget_info(rowy_ylog,/droplist_select)
      if ylog eq 0 then axistype='LINEAR' else axistype='LOG'
      if pnum ge 0 and pnum lt nplot then begin
        widget_control,rowy_yrng,get_value=cyr
	yr=str_2_arr(cyr,/r4,/squish) & nyr=n_elements(yr)
	ok='ok'
	if nyr eq 2 then begin
	  if yr(0) lt 0 or yr(1) lt 0 then ok='not OK for LOG'
	endif else begin
	  if yr(0) lt 0 then ok='not OK for LOG'
	endelse
	if axistype eq 'LINEAR' then ok='ok'
	if ok eq 'ok' then begin
	  widget_control,txts_stts,set_value='Setting Y-axis type to '+$
		axistype+' for plot '+strtrim(pnum,2)
	  ststr.window.YLOG(pnum)=ylog
	  ;
          darshana,x,y,ststr,wid=wwork,pid=pnum
	endif else begin
	  widget_control,txts_stts,set_value='Invalid Y-range'
	  widget_control,rowy_ylog,set_droplist_select=ylog-1
	endelse
      endif else widget_control,txts_stts,$
	set_value='Invalid plot window: '+strtrim(pnum,2)
    end

    ;4,1}{ fifth row, first column

    'row0_sele': begin				;select label
      widget_control,txts_stts,set_value=$
	'MOUSE CONTROL: Select label with mouse button clicks'
      picklabl,x,y,xloc,xlab,ylab,ilab
      widget_control,txts_stts,set_value='Done'
      lnum=ilab(0) & if lnum eq -1 then lnum=0L
      updatelabl=1
    end
    'row0_next': begin				;pick next label
      widget_control,row0_labl,get_value=lnum
      lnum=lnum+1 & if lnum eq nlabel then lnum=0L
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      updatelabl=1
    end
    'row0_prev': begin				;pick previous label
      widget_control,row0_labl,get_value=lnum
      lnum=lnum-1 & if lnum lt 0 then lnum=nlabel-1L
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      updatelabl=1
    end
    'row0_labl': begin				;pick label
      widget_control,row0_labl,get_value=lnum
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      updatelabl=1
    end
    'row1_ltxt': begin				;label text
      widget_control,row0_labl,get_value=lnum
      widget_control,row1_ltxt,get_value=ltxt
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).LABEL=ltxt
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row0_hide': begin				;show/hide label
      widget_control,row0_labl,get_value=lnum
      ihide=widget_info(row0_hide,/droplist_select)
      widget_control,row3_size,get_value=size
      lnum = lnum mod nlabel > 0
      if ihide eq 1 then widget_control,txts_stts,set_value='Label '+$
	strtrim(lnum,2)+' will be hidden from view' else $
	widget_control,txts_stts,set_value='Displaying label '+strtrim(lnum,2)
      if ihide eq 1 then ststr.(lnum+2).SIZE=-abs(size) else $
	ststr.(lnum+2).SIZE=abs(size)
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row0_newl': begin				;new label
      widget_control,txts_stts,set_value='CURSOR CONTROL: '+$
	'Mark reference-location and label position with Click+Drag'
      cursor,x0,y0,/down & cursor,x1,y1,/up
      widget_control,txts_stts,set_value='KEYBOARD CONTROL: '+$
	'Verify reference location at keyboard; <CR> to accept'
      cx0=strtrim(x0,2) & c1='' & read,prompt='Label for ['+cx0+']> ',c1
      if strtrim(c1,2) ne '' then x0=float(c1)
      ;
      mudra,ststr,floc=x0,flab='NEW LABEL'
      nlabel=n_elements(ststr.LOC.X) & zlabel=lindgen(nlabel)
      lnum=nlabel-1 & updatelabl=1
    end
    'dele_labl': begin				;delete label
      widget_control,row0_labl,get_value=lnum
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='KEYBOARD CONTROL: Deleting Label '+$
	strtrim(lnum,2)+'; confirm at keyboard!'
      c1='' & read,prompt='Delete Label '+strtrim(lnum,2)+'? [n/y]> ',c1
      if strlowcase(strtrim(c1,2)) eq 'y' then begin
      mudra,ststr,filter=[lnum]
        nlabel=n_elements(ststr.LOC.X) & lnum = lnum mod nlabel > 0
        zlabel=lindgen(nlabel)
	widget_control,txts_stts,set_value='The deed, he is done.'
        updatelabl=1
      endif else widget_control,txts_stts,set_value='OK, did nothing'
    end
    'row2_styl': begin				;quick menu of styles
      widget_control,txts_stts,set_value=$
	'Short cut to label styles (see RANGOLI.PRO for details)'
      istyl=widget_info(row2_styl,/droplist_select)
      widget_control,row0_labl,get_value=lnum
      lnum = lnum mod nlabel > 0
      widget_control,col2_yfrc,get_value=frac
      frac=abs(frac) & if frac gt 1 then frac=abs(frac)/100. < 0.999
      case istyl of
	1: pstyle=1.+frac
	2: pstyle=2.+frac
	3: pstyle=3.+frac
	4: pstyle=4.+frac
	5: pstyle=5.+frac
	6: pstyle=6.+frac
	7: pstyle=-1
	8: pstyle=-2
	9: pstyle=-3
	10: pstyle=-4-frac
	11: pstyle=-5-frac
	12: pstyle=-6-frac
	else: pstyle=0
      endcase
      ;
      widget_control,col2_ytop,get_value=cyt
      yt=str_2_arr(cyt,/r4,/squish) & nyt=n_elements(yt)
      topy=[0.,0.2]
      case nyt of
	1: begin
	  if yt(0) lt 0.5 then topy(1)=yt(0) else topy=[yt(0),0.]
	end
	else: topy=[yt(0),yt(1)]
      endcase
      widget_control,col2_ybot,get_value=cyb
      yb=str_2_arr(cyb,/r4,/squish) & nyb=n_elements(yb)
      case nyb of
	1: boty=[yb(0),1.1]
	else: boty=[yb(0),yb(1)]
      endcase
      ;
      rangoli,ststr,pstyle,x,y,pnum=pnum,lnum=lnum,ytop=topy,ybot=boty
      ialign=(where(zalign eq ststr.(lnum+2).ALIGN))(0)
      iarrange=(where(zarrange eq ststr.(lnum+2).ARRANGE))(0)
      widget_control,row2_algn,set_droplist_select=ialign
      widget_control,row2_arng,set_droplist_select=iarrange
      updatelabl=1
    end
    'row2_algn': begin				;label alignment
      widget_control,row0_labl,get_value=lnum
      ialign=widget_info(row2_algn,/droplist_select)
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).ALIGN=zalign(ialign)
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row2_arng': begin				;multi-line arrangement
      widget_control,row0_labl,get_value=lnum
      iarrange=widget_info(row2_arng,/droplist_select)
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).ARRANGE=zarrange(iarrange)
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row3_ornt': begin				;label orientation
      widget_control,row0_labl,get_value=lnum
      widget_control,row3_ornt,get_value=orient
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).ORIENT=orient
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row3_size': begin				;label size
      widget_control,row0_labl,get_value=lnum
      widget_control,row3_size,get_value=size
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).SIZE=abs(size)
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row3_thck': begin				;line thickness 
      widget_control,row0_labl,get_value=lnum
      widget_control,row3_thck,get_value=lthick
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).THICK=lthick
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row4_labc': begin				;label color number
      widget_control,row0_labl,get_value=lnum
      widget_control,row4_labc,get_value=labc
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).LABCOLOR=labc
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row4_labn': begin				;label color name
      widget_control,row4_labc,get_value=labc
      widget_control,row4_labn,get_value=labn
      if strtrim(labn(0),2) ne '' then setcolor,labn(0),labc
      widget_control,txts_stts,set_value='Setting color number '+$
	strtrim(labc,2)+' to '+labn(0)
      widget_control,row4_labn,set_value=''
    end
    'row4_linc': begin				;line color number
      widget_control,row0_labl,get_value=lnum
      widget_control,row4_linc,get_value=linc
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).LINCOLOR=linc
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row4_linn': begin				;line color name
      widget_control,row4_linc,get_value=linc
      widget_control,row4_linn,get_value=linn
      if strtrim(linn(0),2) ne '' then setcolor,linn(0),linc
      widget_control,txts_stts,set_value='Setting color number '+$
	strtrim(linc,2)+' to '+linn(0)
      widget_control,row4_linn,set_value=''
    end
    'row5_ulin': begin				;label underline
      widget_control,row0_labl,get_value=lnum
      widget_control,row5_ulin,get_value=ulin
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).UNDERLINE=ulin
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'row5_slin': begin				;label sideline
      widget_control,row0_labl,get_value=lnum
      widget_control,row5_slin,get_value=slin
      lnum = lnum mod nlabel > 0
      widget_control,txts_stts,set_value='Working on label '+strtrim(lnum,2)
      ststr.(lnum+2).SIDELINE=slin
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end

    ;5,1}{ second column

    'col2_move': begin				;adjust label position
      widget_control,txts_stts,set_value=$
	'CURSOR CONTROL: Use cursor to reposition label'
      coarse=1-widget_info(col2_posn,/droplist_select)	;("1-" because 0th
      springy=1-widget_info(col2_ties,/droplist_select)	;element is default)
      widget_control,row1_pnum,get_value=pnum
      widget_control,row0_labl,get_value=lnum
      widget_control,col2_coar,get_value=coar
      pnum = pnum mod nplot > 0
      lnum = lnum mod nlabel > 0
      coar=abs(coar)
      ;
      widget_control,col2_ytop,get_value=cyt
      yt=str_2_arr(cyt,/r4,/squish) & nyt=n_elements(yt)
      topy=[0.,0.2]
      case nyt of
	1: begin
	  if yt(0) lt 0.5 then topy(1)=yt(0) else topy=[yt(0),0.]
	end
	else: topy=[yt(0),yt(1)]
      endcase
      widget_control,col2_ybot,get_value=cyb
      yb=str_2_arr(cyb,/r4,/squish) & nyb=n_elements(yb)
      case nyb of
	1: boty=[yb(0),1.1]
	else: boty=[yb(0),yb(1)]
      endcase
      ;
      savstr=ststr			;for UNDO's sake
      movelabl,x,y,ststr,pnum,lnum,coarse=coarse,dcor=coar,springy=springy,$
	ytop=topy,ybot=boty
      widget_control,txts_stts,set_value='Done.'
      widget_control,col2_posn,set_droplist_select=coarse+1
      widget_control,col2_ties,set_droplist_select=springy+1
      updatelabl=1
    end

    'col2_undo': begin				;undo previous change
      tmpstr=ststr
      ststr=savstr
      savstr=tmpstr	;for that double UNDO's sake!
      widget_control,txts_stts,set_value='Undoing previous change!'
      updatelabl=1
    end

    'grps_alle': begin				;group attributes of all
      widget_control,txts_stts,set_value=$
	'Place all labels into single group; select group leader via cursor'
      picklabl,x,y,xloc,xlab,ylab,ilab
      gllab=ilab(0)
      if gllab ne -1 then begin
	savstr=ststr			;for UNDO's sake
	ststr.LOC.GROUP=gllab
	widget_control,txts_stts,set_value=$
	'Grouping label attributes to properties of Label '+strtrim(gllab,2)
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,set_value=$
	'No group leader selected; ignoring.'
    end
    'grps_sele': begin				;group attributes of selected
      widget_control,txts_stts,set_value='Group attributes of cursor'+$
	' selected labels; first selection is group leader.'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	gllab=ilab(0) & nilab=n_elements(ilab)
	widget_control,txts_stts,set_value=$
	'Grouping label attributes to properties of Label '+strtrim(gllab,2)
	gg=ststr.LOC.GROUP
	for j=0L,nilab-1L do gg(ilab(j))=gllab
	ststr.LOC.GROUP=gg
        darshana,x,y,ststr,wid=wwork,pid=pnum
      endif else widget_control,txts_stts,set_value=$
	'No labels selected; ignoring.'
    end
    'grps_dele': begin				;remove label from group
      widget_control,txts_stts,set_value=$
	'CURSOR CONTROL: Remove cursor-selected labels from existing groupings'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
        ststr.LOC.GROUP(ilab)=ilab
	widget_control,txts_stts,set_value='Done.'
      endif else $
	widget_control,txts_stts,set_value='None selected, naught done'
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end
    'grps_undo': begin				;ungroup all
      widget_control,txts_stts,set_value='Ungroup all labels'
      ststr.LOC.GROUP=lindgen(nlabel)
      darshana,x,y,ststr,wid=wwork,pid=pnum
    end

    'alny_alle': begin				;align all label Y-positions
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'which all others are to be aligned'
      springy=1-widget_info(col2_ties,/droplist_select)	;rigidity
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	ymark=ststr.(ilab(0)+2).POS(1)
	for jlab=0L,nlabel-1L do begin
	  marky=ststr.(jlab+2).POS(1)
	  ypath=ststr.(jlab+2).YPATH
	  if keyword_set(springy) then begin
	    if springy gt 0 then begin
	      yden=marky-ststr.LOC.Y(jlab) & fypath=(ypath-ststr.LOC.Y(jlab))
	      if yden ne 0 then fypath=fypath/yden else fypath=0*fypath
	      ypath=fypath*(ymark-ststr.LOC.Y(jlab))+ststr.LOC.Y(jlab)
	    endif else begin
	      dely=marky-ymark & ypath=ypath-dely
	      ststr.LOC.Y(jlab)=ststr.LOC.Y(jlab)-dely
	    endelse
	    ststr.(jlab+2).YPATH=ypath
	  endif
	  ststr.(jlab+2).POS(1)=ymark
	endfor
	widget_control,txts_stts,set_value='Done'
	lnum=ilab(0) & updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'alny_grps': begin				;align Y-positions of group
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      springy=1-widget_info(col2_ties,/droplist_select)	;rigidity
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for klab=0L,nilab-1L do begin		;{go through all selections
	  jlab=ilab(klab)
	  ymark=(ststr.(jlab+2).POS)(1)
	  oo=where(gg eq gg(jlab),moo)
	  ;if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).POS(1)=ymark
	  for jj=0L,moo-1L do begin		;{groupings exist
	    jjlab=oo(jj)
	    marky=ststr.(jjlab+2).POS(1)
	    ypath=ststr.(jjlab+2).YPATH
	    if keyword_set(springy) then begin
	      if springy gt 0 then begin
		yden=marky-ststr.LOC.Y(jjlab)
		fypath=(ypath-ststr.LOC.Y(jjlab))
	        if yden ne 0 then fypath=fypath/yden else fypath=0*fypath
		ypath=fypath*(ymark-ststr.LOC.Y(jjlab))+ststr.LOC.Y(jjlab)
	      endif else begin
	        dely=marky-ymark & ypath=ypath-dely
		ststr.LOC.Y(jjlab)=ststr.LOC.Y(jjlab)-dely
	      endelse
	      ststr.(jjlab+2).YPATH=ypath
	    endif
	    ststr.(jjlab+2).POS(1)=ymark
	  endfor				;JJ=0,MOO-1}
	endfor					;KLAB=0,NILAB-1}
	widget_control,txts_stts,set_value='Done'
	lnum=ilab(0) & updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'alny_sele': begin				;align selected label Y-pos's
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      springy=1-widget_info(col2_ties,/droplist_select)	;rigidity
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	ymark=(ststr.(ilab(0)+2).POS)(1)
	for klab=1L,nilab-1L do begin
	  jlab=ilab(klab)
	  marky=ststr.(jlab+2).POS(1)
	  ypath=ststr.(jlab+2).YPATH
	  if keyword_set(springy) then begin
	    if springy gt 0 then begin
	      yden=marky-ststr.LOC.Y(jlab) & fypath=(ypath-ststr.LOC.Y(jlab))
	      if yden ne 0 then fypath=fypath/yden else fypath=0*fypath
	      ypath=fypath*(ymark-ststr.LOC.Y(jlab))+ststr.LOC.Y(jlab)
	    endif else begin
	      dely=marky-ymark & ypath=ypath-dely
	      ststr.LOC.Y(jlab)=ststr.LOC.Y(jlab)-dely
	    endelse
	    ststr.(jlab+2).YPATH=ypath
	  endif
	  ststr.(jlab+2).POS(1)=ymark
	endfor
	widget_control,txts_stts,set_value='Done'
	lnum=ilab(0) & updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'ornt_alle': begin				;match all orientations
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'define orientation'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	omark=ststr.(ilab(0)+2).ORIENT
	for jlab=0L,nlabel-1L do ststr.(jlab+2).ORIENT=omark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'ornt_grps': begin				;match grouped orientations
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for jlab=0L,nilab-1L do begin
	  klab=ilab(jlab)
	  omark=ststr.(klab+2).ORIENT
	  oo=where(gg eq gg(klab),moo)
	  if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).ORIENT=omark
	endfor
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'ornt_sele': begin				;match selected orientations
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	omark=ststr.(ilab(0)+2).ORIENT
	for jlab=1L,nilab-1L do ststr.(ilab(jlab)+2).ORIENT=omark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'lsiz_alle': begin				;match all label sizes
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'define label sizes'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	lmark=ststr.(ilab(0)+2).SIZE
	for jlab=0L,nlabel-1L do ststr.(jlab+2).SIZE=lmark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'lsiz_grps': begin				;match grouped label sizes
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for jlab=0L,nilab-1L do begin
	  klab=ilab(jlab)
	  lmark=ststr.(klab+2).SIZE
	  oo=where(gg eq gg(klab),moo)
	  if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).SIZE=lmark
	endfor
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'lsiz_sele': begin				;match selected label sizes
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	lmark=ststr.(ilab(0)+2).SIZE
	for jlab=1L,nilab-1L do ststr.(ilab(jlab)+2).SIZE=lmark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'thck_alle': begin			;match all line thicknesses
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'define line thickness'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	tmark=ststr.(ilab(0)+2).THICK
	for jlab=0L,nlabel-1L do ststr.(jlab+2).THICK=tmark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'thck_grps': begin			;match grouped line thicknesses
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for jlab=0L,nilab-1L do begin
	  klab=ilab(jlab)
	  tmark=ststr.(klab+2).THICK
	  oo=where(gg eq gg(klab),moo)
	  if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).THICK=tmark
	endfor
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'thck_sele': begin			;match selected line thicknesses
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	tmark=ststr.(ilab(0)+2).THICK
	for jlab=1L,nilab-1L do ststr.(ilab(jlab)+2).THICK=tmark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'arng_alle': begin				;match all arrangements
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'define label text arrangement'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	amark=ststr.(ilab(0)+2).ARRANGE
	for jlab=0L,nlabel-1L do ststr.(jlab+2).ARRANGE=amark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'arng_grps': begin				;match grouped arrangements
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for jlab=0L,nilab-1L do begin
	  klab=ilab(jlab)
	  amark=ststr.(klab+2).ARRANGE
	  oo=where(gg eq gg(klab),moo)
	  if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).ARRANGE=amark
	endfor
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'arng_sele': begin				;match selected arrangements
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	amark=ststr.(ilab(0)+2).ARRANGE
	for jlab=1L,nilab-1L do ststr.(ilab(jlab)+2).ARRANGE=amark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'ulin_alle': begin				;match all underlines
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'define underlines'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	umark=ststr.(ilab(0)+2).UNDERLINE
	for jlab=0L,nlabel-1L do ststr.(jlab+2).UNDERLINE=umark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'ulin_grps': begin				;match grouped underlines
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for jlab=0L,nilab-1L do begin
	  klab=ilab(jlab)
	  umark=ststr.(klab+2).UNDERLINE
	  oo=where(gg eq gg(klab),moo)
	  if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).UNDERLINE=umark
	endfor
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'ulin_sele': begin				;match selected underlines
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	umark=ststr.(ilab(0)+2).UNDERLINE
	for jlab=1L,nilab-1L do ststr.(ilab(jlab)+2).UNDERLINE=umark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'slin_alle': begin				;match all sidelines
      widget_control,txts_stts,set_value='CURSOR CONTROL: select label to '+$
	'define sidelines'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	smark=ststr.(ilab(0)+2).SIDELINE
	for jlab=0L,nlabel-1L do ststr.(jlab+2).SIDELINE=smark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected; naught done'
    end
    'slin_grps': begin				;match grouped sidelines
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab) & gg=ststr.LOC.GROUP
	for jlab=0L,nilab-1L do begin
	  klab=ilab(jlab)
	  smark=ststr.(klab+2).SIDELINE
	  oo=where(gg eq gg(klab),moo)
	  if moo gt 0 then for jj=0L,moo-1L do ststr.(oo(jj)+2).SIDELINE=smark
	endfor
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end
    'slin_sele': begin				;match selected sidelines
      widget_control,txts_stts,set_value='CURSOR CONTROL: select labels;'+$
	' first selection defines Y position'
      picklabl,x,y,xloc,xlab,ylab,ilab
      if ilab(0) ne -1 then begin
	savstr=ststr			;for UNDO's sake
	nilab=n_elements(ilab)
	smark=ststr.(ilab(0)+2).SIDELINE
	for jlab=1L,nilab-1L do ststr.(ilab(jlab)+2).SIDELINE=smark
	widget_control,txts_stts,set_value='Done'
	updatelabl=1
      endif else widget_control,txts_stts,set_value='None selected, naught done'
    end

    'col2_posn': begin					;line positioning
      coarse=1-widget_info(col2_posn,/droplist_select)
      widget_control,col2_coar,get_value=coar
      if coarse eq 1 then c1='Line Vertexes will tend to line up' else $
	c1='Lines to labels will be free form'
      widget_control,txts_stts,set_value=c1
    end
    'col2_ties': begin					;relative positioning
      springy=1-widget_info(col2_ties,/droplist_select)
      case springy of
	1: c1='Move label, and the attached line will move like a spring.'
	0: c1='Label and attached lines may be moved independently.'
	-1: c1='Lines attached to label are held as rigid as possible.'
	else: c1='Re RIGIDITY: This should not be happening'
      endcase
      widget_control,txts_stts,set_value=c1
    end
    'col2_coar': begin				;coarseness
      widget_control,col2_coar,get_value=coar
      c1='Coarseness of line vertex positioning, as fraction of label size'
      widget_control,txts_stts,set_value=c1
    end
    'col2_ytop': begin				;YTOP
      widget_control,col2_ytop,get_value=cyt
      yt=str_2_arr(cyt,/r4,/squish) & nyt=n_elements(yt)
      topy=[0.,0.2]
      case nyt of
	1: begin
	  if yt(0) lt 0.5 then topy(1)=yt(0) else topy=[yt(0),0.]
	end
	else: topy=[yt(0),yt(1)]
      endcase
      widget_control,txts_stts,set_value='Setting Label Y position to: '+$
        'YMAX-('+strtrim(topy(0),2)+'+'+strtrim(topy(1),2)+'*(YMAX-YMIN))'
    end
    'col2_ybot': begin				;YBOT
      widget_control,col2_ybot,get_value=cyb
      yb=str_2_arr(cyb,/r4,/squish) & nyb=n_elements(yb)
      case nyb of
	1: boty=[yb(0),1.1]
	else: boty=[yb(0),yb(1)]
      endcase
      widget_control,txts_stts,set_value='Setting line end Y position to: '+$
        strtrim(boty(0),2)+'+'+strtrim(boty(1),2)+'*Y(XLOC)'
    end
    'col2_yfrc': begin
      widget_control,col2_yfrc,get_value=frac
      c1="length of label's line as fraction of window height"
      widget_control,txts_stts,set_value=c1
    end

    ;1,2}

    else: ;nothing
  endcase					;handle EV}

  if keyword_set(rewstr) then begin
    ststr.WINDOW.MULTI(1)=ncol & ststr.WINDOW.MULTI(2)=nrow
    subt=ststr.WINDOW.SUBTITLE & char=ststr.WINDOW.CHARS
    xt=ststr.WINDOW.XTITLE & yt=ststr.WINDOW.YTITLE
    xsty=ststr.WINDOW.XSTYLE & ysty=ststr.WINDOW.YSTYLE
    xlog=ststr.WINDOW.XLOG & ylog=ststr.WINDOW.YLOG
    xmin=ststr.WINDOW.XMIN & ymin=ststr.WINDOW.YMIN
    xmax=ststr.WINDOW.XMAX & ymax=ststr.WINDOW.YMAX
    if ncol*nrow gt nplot then nplot=ncol*nrow
    if nplot gt onplot then begin	;(add window(s)
      nn=nplot-onplot
      vv=replicate(subt(0),nn) & subt=[subt,vv]
      vv=replicate(xt(0),nn) & xt=[xt,vv]
      vv=replicate(yt(0),nn) & yt=[yt,vv]
      vv=replicate(xsty(0),nn) & xsty=[xsty,vv]
      vv=replicate(ysty(0),nn) & ysty=[ysty,vv]
      vv=replicate(xlog(0),nn) & xlog=[xlog,vv]
      vv=replicate(ylog(0),nn) & ylog=[ylog,vv]
      vv=replicate(char(0),nn) & char=[char,vv]
      vv=replicate(xmin(0),nn) & xmin=[xmin,vv]
      vv=replicate(xmax(0),nn) & xmax=[xmax,vv]
      vv=replicate(ymin(0),nn) & ymin=[ymin,vv]
      vv=replicate(ymax(0),nn) & ymax=[ymax,vv]
      onplot=nplot & zplot=lindgen(nplot)
    endif else begin			;)(delete a window
      nplot=n_elements(zplot) & onplot=nplot
      subt=subt(zplot) & xt=xt(zplot) & yt=yt(zplot) & char=char(zplot)
      xsty=xsty(zplot) & ysty=ysty(zplot) & xlog=xlog(zplot) & ylog=ylog(zplot)
      xmin=xmin(zplot) & xmax=xmax(zplot) & ymin=ymin(zplot) & ymax=ymax(zplot)
    endelse				;NPLOT?ONPLOT)
    wstr=create_struct('MULTI',ststr.WINDOW.MULTI,'TITLE',ststr.WINDOW.TITLE,$
	'SUBTITLE',subt,'XTITLE',xt,'YTITLE',yt,'XSTYLE',xsty,$
	'YSTYLE',ysty,'XLOG',xlog,'YLOG',ylog,'CHARS',char,$
	'XMIN',xmin,'XMAX',xmax,'YMIN',ymin,'YMAX',ymax)
    tstr=ststr & ststr=create_struct('WINDOW',wstr,'LOC',tstr.LOC)
    for k=0L,nlabel-1L do ststr=create_struct(ststr,$
	'L'+strtrim(k+1L,2),tstr.(k+2L))
    replot=1
  endif

  ;	update label info in widget
  if keyword_set(updatelabl) then begin
    ialign=widget_info(row2_algn,/droplist_select)
    iarrange=widget_info(row2_arng,/droplist_select)
    ihide=widget_info(row0_hide,/droplist_select)
    widget_control,row0_labl,set_value=lnum
    widget_control,row0_hide,set_droplist_select=ihide
    widget_control,row1_ltxt,set_value=ststr.(lnum+2).LABEL
    widget_control,row2_styl,set_droplist_select=0
    widget_control,row2_algn,set_droplist_select=ialign
    widget_control,row2_arng,set_droplist_select=iarrange
    widget_control,row3_ornt,set_value=ststr.(lnum+2).ORIENT
    widget_control,row3_size,set_value=ststr.(lnum+2).SIZE
    widget_control,row3_thck,set_value=ststr.(lnum+2).THICK
    widget_control,row4_labc,set_value=ststr.(lnum+2).LABCOLOR
    widget_control,row4_linc,set_value=ststr.(lnum+2).LINCOLOR
    widget_control,row5_ulin,set_value=ststr.(lnum+2).UNDERLINE
    widget_control,row5_slin,set_value=ststr.(lnum+2).SIDELINE
    darshana,x,y,ststr,wid=wwork,pid=pnum
    xmark=ststr.LOC.X(lnum) & ymark=max(!y.crange)
    xyouts,xmark,ymark,'!95!3',charsize=3,color=!d.n_colors-1,align=0.5
    updatelabl=0
  endif

  ;	replot
  if keyword_set(replot) then begin
    darshana,x,y,ststr,wid=wshow,_extra=e
    wset,wwork
    replot=0
  endif

endwhile				;end of endless loop}

return
end
