function pickrange,x,y,sigmay=sigmay,xsize=xsize,ysize=ysize,wid=wid,$
	help=help,dynrng=dynrng,markx=markx,marky=marky,markp=markp,$
	markc=markc,marks=marks,marko=marko,markl=markl,legg=legg,$
	verbose=verbose,xo=xo,yo=yo,oxr=oxr,oyr=oyr,xtitle=xtitle,$
	ytitle=ytitle,title=title,opx=opx,opy=opy,opc=opc, _ref_extra=ex
;+
;function	pickrange
;	interactively select intervals in data -- returns a long-integer
;	array of position indices of selected bins
;
;syntax
;	oo=pickrange(x,y,sigmay=sigmay,xsize=xsize,ysize=ysize,wid=wid,$
;	/help,dynrng=dynrng,markx=markx,marky=marky,markp=markp,markc=markc,$
;	marks=marks,marko=marko,markl=markl,/legg,verbose=verbose,$
;	xo=xo,yo=yo,oxr=oxr,oyr=oyr,xtitle=xtitle,ytitle=ytitle,$
;	title=title,opx=opx,opy=opy,opc=opc,$
;	/nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol,$
;	stasiz=stasiz,stathk=stathk)
;		
;parameters
;	x	[INPUT; required] points where "reference curve" is defined
;	y	[INPUT] the reference curve (usually data)
;		* if not given, X is taken to be Y, and the indices of X are
;		  taken to be X!
;
;keywords
;	sigmay	[INPUT] if given, plots error bars on Y
;	xsize	[INPUT] width of window to display curve in
;		* ignored if WID is not set
;	ysize	[INPUT] height of window to display curve in
;		* ignored if WID is not set
;	wid	[INPUT] window ID
;		* use current/default window if -ve/0
;	help	[INPUT] prints out usage summary
;	dynrng	[INPUT; default: 3 dex] maximum dynamic range in log scale
;	markx	[INPUT] mark these x-coordinate points
;	marky	[INPUT] mark at these y-coordinate positions
;		* ignored if incompatible with MARKX
;		* if not given, labels and marks will be placed 2/3 of
;		  the way up the plot.
;	markp	[INPUT: default: 1] PSYMs for each (MARKX,MARKY) pair
;	markc	[INPUT: default: 222] COLOR for each (MARKX,MARKY) pair
;	marks	[INPUT: default: 1] SYMSIZE for each (MARKX,MARKY) pair
;	marko	[INPUT; default: 90 deg] ORIENTATION for MARKL
;		* MARK[P,C,S,O]($) overrides the hardcoded defaults
;	markl	[INPUT] label for each (MARKX,MARKY) pair
;		* marked only if given, unless there's only one element,
;		  in which case it is accepted as default
;	legg	[INPUT] if set, drops a leg to the X-axis
;	xo	[I/O] returns the x-position of last cursor click
;	yo	[I/O] returns the y-position of last cursor click
;	oxr	[I/O] 2-element array containing the final XRANGE
;	oyr	[I/O] 2-element array containing the final YRANGE
;	xtitle	[INPUT] passed to PLOT w/o checking
;	ytitle	[INPUT] passed to PLOT w/o checking
;	title	[INPUT] passed to PLOT w/o checking
;	opx	[INPUT] X-coords of a curve to be overplotted
;	opy	[INPUT] the curve to be overplotted
;		* OPX and OPY may be 2-D, in which case the 1st dimension
;		  defines the number of curves to be overplotted
;		* sizes of OPX and OPY >must< match
;	opc	[INPUT] color for the overplotting
;		* default is 222
;	verbose	[INPUT] controls chatter
;		* replaces QUIET
;
;	_ref_extra	[INPUT ONLY] pass defined keywords to
;		STAMPLE: NUTHIN,NODAY,NOTIME,NOUSER,NOPACK,STACOL,STASIZ,STATHK
;
;restrictions
;	* works only under X
;
;usage summary
;	see variable USAGE below
;
;history
;	vinay kashyap (Oct96)
;	added keyword SIGMA; changed to return only position indices of
;	  selected bins; modified click/drag behavior (VK; Dec96)
;	added keyword MARKX (VK; Feb97)
;	added keywords MARKY, MARKP, MARKC, MARKS, MARKL, MARKO, LEGG, QUIET,
;	  XO, YO, OXR, OYR; modified behavior of LEFT CLICK; added keyboard
;	  control via LEFT CLICK (VK; Aug98)
;	allowed XO,YO,OXR,OYR to be I/O (VK; Sep98)
;	now MARKL(0) does not act as default, XYOUTS forced to CLIP (VK; Nov98)
;	now zoom is remembered if OXR or OYR are set on input; made XR,YR
;	  double precision; MARKY bug corrected; renamed keyword SIGMA as
;	  SIGMAY; changed _EXTRA to _REF_EXTRA (VK; Dec98)
;	bug correction: reversed X-axis was preventing zoom (VK; Apr99)
;	added keywords OPX,OPY,OPC (VK; Oct99)
;	cursor changes shape when mouse button is pressed (VK; FebMM)
;	undid cursor shape changes (VK; AprMM)
;	MARKC can now be I*4; added extra keyboard options; added call to
;	  STAMPLE (VK; DecMM)
;	changed MARKY default from 0.5 to 0.67 (VK; JanMMI)
;	added keyword TITLE; improved color-scale setting for 24-bit
;	  consoles (VK; FebMMI)
;	streamlined X,Y input behavior (VK; May02)
;	handle color tables in 24-bit displays (VK; Jun02)
;	replaced keyword QUIET with VERBOSE, changed color defaults
;	  (VK; Jul02)
;	speeded up marker plots if no variations exist (VK; Dec02)
;	button press status now stored in !MOUSE, not !ERR (VK; Apr09)
;-

;	help
usage=[	'',$
	'		USAGE SUMMARY',$
	' * call as function with at least one parameter',$
	' * returns position indices of selected array elements',$
	' * hold down left button and drag to select range',$
	' * hold down middle button and drag (any amount) to delete selection',$
	'   -- if cursor is not on any selection, delete previous',$
	'   -- if cursor on some selection, delete selection',$
	' * hold down right button and drag to select area to zoom',$
	' * click left button to:',$
	'   -- zoom to area selected with right button drag',$
	'   -- if no zoom possible, give program control to keyboard',$
	'      ?,q,x,z,h,r,<,>,2,^,_ have meaning',$
	'      any other key to list available commands',$
	' * click middle button to:',$
	'   -- if zooomed, then unzoom',$
	'      (if click is outside plot window, then toggle log scale)',$
	'   -- if not zoomed, then toggle log scaling of axes',$
	'      (inside plot window ... toggle both axes)',$
	'      (below/above x-axis ... toggle in X only)',$
	'      (left/right of y-axis ... toggle in Y only)',$
	'      (if above/right, toggle even if some elements are .LE. 0)',$
	' * click right button to exit',$
	'' ]
nuse=n_elements(usage)
cmsg=[	'	CLICK & DRAG mouse buttons',$
	'  LEFT: select interval; MIDDLE: deselect; RIGHT: select zoom area',$
  	'	CLICK mouse buttons',$
	'LEFT: zoom/keyboard; MIDDLE: unzoom/toggle-logscale; RIGHT: QUIT']
nmsg=n_elements(cmsg)

;	verbosity
vv=0 & if keyword_set(verbose) then vv=long(verbose(0))>1

;	calling sequence
ok='ok' & np=n_params(0) & nx=n_elements(x) & ny=n_elements(y)
if np eq 0 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if nx eq 1 then ok='Cannot plot scalars'
if ok ne 'ok' then begin
  print,'Usage: oo=pickrange(x,y,sigmay=s,xsize=xs,ysize=ys,wid=window_id,$'
  print,'       dynrng=dynamic_range,/help,markx=mx,marky=my,markp=psym,$'
  print,'       markc=color,marks=symsize,markl=label,marko=orientation,$'
  print,'	  /legg,verbose=verbose,xo=xout,yo=yout,oxr=outxr,oyr=outyr,$'
  print,'       opx=opx,opy=opy,opc=opc)
  print,'  interactively select intervals in data'
  if vv gt 0 or keyword_set(help) then $
    for i=0,nuse-1 do print,usage(i)
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	check variables
xx=x
if ny ne nx then begin
  if ny eq nx-1L then begin
    xx=0.5*(xx(1:*)+xx)	;assuming that what came in were bin-boundaries
    yy=y
  endif else begin
    yy=xx & xx=lindgen(nx)
  endelse
endif else yy=y & ny=n_elements(yy)
;
if keyword_set(sigmay) then begin		;for error bars
  sig=sigmay & szs=size(sig) & su=0*xx & sd=su
  case szs(0) of
  1: begin					;symmetric error bars
    if szs(1) lt nx then su(0)=sig else su=sig(0:nx-1)
    sd=su
  end
  2: begin					;asymmetric error bars
    if szs(1) eq 2 then begin
      if szs(2) lt nx then su(0)=sig(0,*) else su=sig(0,0:nx-1)
      if szs(2) lt nx then sd(0)=sig(1,*) else sd=sig(1,0:nx-1)
    endif else begin
      su(*)=sig(0) & sd(*)=sig(1)
    endelse
  end
  else: begin					;constant, or whatever
    su(*)=sig(0) & sd=su
  end
  endcase
endif

;	marker stuff
mx=n_elements(markx) & my=n_elements(marky)
mp=n_elements(markp) & mc=n_elements(markc)
ms=n_elements(marks) & ml=n_elements(markl) & mo=n_elements(marko)
markall=1
if mx gt 0 then begin
  psim=intarr(mx) & kol=lonarr(mx)
  ssiz=fltarr(mx) & lbl=strarr(mx) & ori=fltarr(mx)
  if mp eq 0 then begin & markp=1 & mp=1 & endif
  if mc eq 0 then begin & markc=200 & mc=1 & endif
  if ms eq 0 then begin & marks=1. & ms=1 & endif
  if mo eq 0 then begin & marko=90. & mo=1 & endif
  if ml ne 0 then pltlbl=1 else pltlbl=0
  if mp gt 1 or mc gt 1 or ms gt 1 or mo gt 1 then markall=0
  for i=0L,mx-1L do psim(i)=markp([i])
  for i=0L,mx-1L do kol(i)=markc([i])
  for i=0L,mx-1L do ssiz(i)=marks([i])
  for i=0L,mx-1L do ori(i)=marko([i])
  if pltlbl then begin
    if ml le mx then lbl(0L:ml-1L)=string(markl) else $
      lbl(0L:mx-1L)=string(markl(0L:mx-1L))
  endif
endif

;	overplot stuff
sox=size(opx) & soy=size(opy) & loc=n_elements(opc)
nsox=n_elements(sox) & nsoy=n_elements(soy)
oplok='ok'	;by default, we'll overplot
if sox(0) eq 0 then oplok='OPX is scalar or undefined' else $
 if sox(0) gt 2 then oplok='cannot handle 3D data' else $
  if sox(0) ne soy(0) then oplok='OPX and OPY have different Ds' else $
   if sox(1) ne soy(1) then oplok='OPX and OPY incompatible' else $
    if sox(nsox-2) gt 5 then oplok='OPX is a perfectly horrid' else $
     if soy(nsoy-2) gt 5 then oplok='OPY is a perfectly horrid'
if oplok ne 'ok' and oplok ne 'OPX is scalar or undefined' then message,$
	oplok,/info
if oplok eq 'ok' then begin
  if sox(0) eq 2 then nopl=sox(1) else nopl=1	;number of curves to overplot
  opcol=lonarr(nopl)+222
  if loc gt 0 and loc le nopl then opcol(0L:loc-1L)=opc(*)
  if loc gt 0 and loc gt nopl then opcol(*)=opc(0L:nopl-1L)
endif

;initialize
minx=min(xx,max=maxx,/nan) & miny=min(yy,max=maxy,/nan)
xr=[minx,maxx] & yr=[miny,maxy]
if xr(0) eq xr(1) then xr(1)=xr(0)+1. & dx=xr(1)-xr(0)
if yr(0) eq yr(1) then yr(1)=yr(0)+1. & dy=(yr(1)-yr(0))
xr=double(xr) & yr=double(yr)
oo=lonarr(nx)-1L & hh=-1L & go_on=1 & xlog=0 & ylog=0 & izm=0 & jzm=0
storange=[xr,yr]
if n_elements(oxr) eq 2 then begin
  jzm=1 & xr=double(oxr)
  if xr(1) lt xr(0) then xr=reverse(xr)
endif else oxr=xr
if n_elements(oyr) eq 2 then begin
  jzm=1 & yr=double(oyr)
  if yr(1) lt yr(0) then yr=reverse(yr)
endif else oyr=yr
currange=storange
ikbrd=0 & k1='' & cret=string("15b)

;	window descriptors
if not keyword_set(xsize) then xsize=(nx<900)>300
if not keyword_set(ysize) then ysize=(ny<800)>256
if n_elements(wid) eq 0 then begin
  if !d.window lt 0 then wid=0 else wid = -(!d.window)
endif

;	dynamic range in log scale
if keyword_set(dynrng) then derange=10.^(-abs(dynrng)) else derange=1e-3

;	usage
if keyword_set(help) then for i=0,nuse-1 do print,usage(i)

;	open window
if wid ge 0 and (!d.window ne wid) then window,wid,xsize=xsize,ysize=ysize
;	plot
plot,xx,yy,psym=10,xr=xr,yr=yr,/xs,/ys,col=200,xtitle=xtitle,ytitle=ytitle,$
	title=title
;	define plotting symbol "-"
usersym,[-0.5,0.5],[0.,0.]
if keyword_set(sigmay) then begin
  oplot,xx,yy+su,psym=8 & oplot,xx,yy-sd,psym=8
  for i=0,nx-1 do oplot,xx(i)*[1,1],yy(i)+[su(i),-sd(i)]
endif
;
;	extra markers?
if mx gt 0 then begin
  xm=[markx]
  recompy=mx-my
  if my eq 0 then marky=0.67*(yr(0)+yr(1)) & ym=0*xm+marky(0)
  ;if recompy ne 0 then begin
  ;  ymark=0.5*(yr(0)+yr(1)) & ym=0*xm+ymark
  ;  if my gt 0 and my le mx then ym(0:my-1)=marky else ym(0:mx-1)=marky(0:mx-1)
  ;endif else ym=marky
  my=n_elements(marky)
  if my le mx then ym(0:my-1)=marky else ym(0:mx-1)=marky(0:mx-1)
  dncolors=256. > !D.N_COLORS	;24-bit color screen temporary fix
  ;kol=long((float(kol)/dncolors)*!d.n_colors+0.5)
  for i=0L,mx-1L do begin
    if keyword_set(markall) then begin
      if i eq 0 then oplot,[xm],[ym],psym=psim[0],symsize=ssiz[0],col=kol[0]
    endif else oplot,[xm(i)],[ym(i)],psym=psim(i),symsize=ssiz(i),col=kol(i)
    if keyword_set(legg) then oplot,xm(i)*[1,1],[yr(0),ym(i)],col=kol(i)
    if pltlbl gt 0 then xyouts,xm(i),ym(i),' '+lbl(i),orientation=ori(i),$
	col=kol(i),noclip=0
  endfor
endif
;
;	extra curves?
if oplok eq 'ok' then begin
  for iopl=0,nopl-1 do begin
    xxo=reform(opx(iopl,*)) & yyo=reform(opy(iopl,*))
    if nopl eq 1 then begin
      xxo=opx & yyo=opy
    endif
    oplot,xxo,yyo,psym=10,color=opcol(iopl)
  endfor
endif
if n_elements(xo) eq 1 and n_elements(yo) eq 1 then tvcrs,xo,yo,/data

;	go ahead and select intervals
if vv gt 0 then for i=0,nmsg-1 do print,cmsg(i) else begin
  ;szq=size(quiet) & nszq=n_elements(szq)
  ;if szq(nszq-2) eq 7 then for i=0,szq(nszq-1)-1 do print,quiet(i)
endelse

idlvers=float(!VERSION.release) & if idlvers le 0 then idlvers=7.

while go_on eq 1 do begin			;{begin selection
  ;	click and/or drag
  if ikbrd eq 0 then begin
    cursor,x0,y0,/down,/data
    if idlvers lt 5 then mbutton=!err else mbutton=!MOUSE.button
    if !D.NAME eq 'X' and vv ge 5 then begin
      if mbutton eq 1 then device,cursor_standard=74
      if mbutton eq 2 then device,cursor_standard=82
      if mbutton eq 4 then device,cursor_standard=100
    endif
    cursor,x1,y1,/up,/data
    if !D.NAME eq 'X' and vv ge 5 then device,/cursor_original
    ;	drag?
    drag=abs(x1-x0)+abs(y1-y0)
    draglim=xr(1)-xr(0) < (yr(1)-yr(0))	;shaky mice
    ;	toggle log?
    wreg=1		;1=>outside plot area; 0=>inside plot area
    if y0 ge yr(0) and y0 le yr(1) and $
       x0 ge xr(0) and x0 le xr(1) then wreg=0
    if x0 ge xr(0) and x0 le xr(1) then begin	;toggle xlog?
      xreg=1		;1=>below plot area; -1=>above; 0=>inside
      if y0 gt yr(1) then xreg=-1
      ;if y0 ge yr(0) and y0 le yr(1) then wreg=0
    endif else xreg=0 ;& xreg=wreg*xreg
    if y0 ge yr(0) and y0 le yr(1) then begin	;toggle ylog?
      yreg=1		;1=>left of plot area; -1=>right; 0=>inside
      if x0 gt xr(1) then yreg=-1
      ;if x0 ge xr(0) and x0 le xr(1) then wreg=0
    endif else yreg=0 ;& yreg=wreg*yreg
  endif 				;else go directly to keyboard control

  if drag lt 1e-4*draglim then begin		;(CLICK!

    if mbutton eq 1 then begin		;(	left click
      if izm gt 0 then begin		;(		zoom
	xr(0)=currange(0) > (minx-0.1*dx)
	xr(1)=currange(1) < (maxx+0.1*dx)
	yr(0)=currange(2) > (miny-0.1*dy)
	yr(1)=currange(3) < (maxy+0.1*dy)
	if xr(0) gt xr(1) then xr=reverse(xr)
	if yr(0) gt yr(1) then yr=reverse(yr)
	izm=0 & jzm=1
      endif else begin			;)(		KEYBOARD CONTROL
	ikbrd=1 & if not keyword_set(k0) then k0=''
	if k0 ne 'h' then begin
          for ik=1,10 do k1=get_kbrd(0)			;flush keyboard
	  k1='KEYBOARD: h for hardcopy, ? for help, x for cursor control'
	  print,form='(a,$)',k1+cret
	  k1=get_kbrd(1) & k11=strlowcase(k1)		;read keyboard input
	endif
	case k11 of					;{keyboard input

	  '?': begin					;help
	    if fix(strmid(!version.release,0,1)) lt 5 then $
		for i=0,nuse-1 do print,usage(i) else xdisplayfile,$
		'help',text=usage,title='PICKRANGE HELP',width=70
	  end
	  'x': begin					;return cursor control
	    ikbrd=0
	    k1='returning control to cursor                                 '
	    if vv eq 0 then print,k1 else $
		for i=0,nmsg-1 do print,cmsg(i)
	  end
	  'z': begin					;force halt
	    print,''
	    message,'HALTING!',/info
	    print,''
	    print,' why would one use this?  one could, for example, manually'
	    print,' set all the DEVICE keywords one wishes here, the fonts,'
	    print,' the bits/pixel, etc., then type .CON, then halt again and'
	    print,' do DEVICE,/CL -- just a bit more freedom to screw up.'
	    print,''
	    print,'	caveat usuator.'
	    print,''
	    stop,'type .CON to continue'
	  end
	  'q': begin					;quit
	    print,'' & go_on=0
	  end
	  ;
	  'h': begin			;(		hardcopy
	    if !d.name eq 'X' then begin	;(if not PS..
	      nops=0		;do make postscript, please.
	      ;	set some useful DEVICE keywords
	      c1='device'
	      if not keyword_set(landscape) then begin
		landscape=0 & c1=c1+',Landscape?'
	      endif else c1=c1+',/landscape'
	      if not keyword_set(color) then begin
		color=0 & c1=c1+',Color?'
	      endif else c1=c1+',/color'
	      if not keyword_set(encapsulated) then begin
		encapsulated=0 & c1=c1+',Encapsulated?'
	      endif else c1=c1+',/encapsulated'
	      if not keyword_set(psdir) then psdir='.'
	      if not keyword_set(psfil) then psfil='idl.ps'
	      ;
	      c0='	= to print, q to abort, or Dir,Fil,l,c,e to set DEVICE kwrds'
	      print,'' & print,c0
	      c1=c1+',file='+psdir+'/'+psfil & print,form='("'+c1+'",a,$)',cret
	      ;
	      while k0 ne '=' do begin			;{set DEVICE keywords
		k0=get_kbrd(1) & k00=strlowcase(k0)
		case k0 of				;(set DEV kwrds
		  'f': begin
		    fil='' & print,fil & read,prompt='hardcopy to file> ',fil
	            fil=strtrim(fil,2) & if fil ne '' then psfil=fil
		  end
		  'd': begin
		    dir='' & print,dir & read,prompt='output directory> ',dir
		    dir=strtrim(dir,2) & if dir ne '' then psdir=dir
		  end
		  'l': if landscape eq 1 then landscape=0 else landscape=1
		  'c': if color eq 1 then color=0 else color=1
		  'e': if encapsulated eq 1 then encapsulated=0 else $
			encapsulated=1
		  'q': begin
		    k0='=' & nops=1
		  end
		  else: begin				;do nothing
		    print,'' & print,c0
		  end
		endcase					;DEVICE kwrds)
		c1='device'
		if keyword_set(landscape) then c1=c1+',/landscape' else $
			c1=c1+',Landscape?'
		if keyword_set(color) then c1=c1+',/color' else $
			c1=c1+',Color?'
		if keyword_set(encapsulated) then c1=c1+',/encapsulated' else $
		  c1=c1+',Encapsulated?'
		c1=c1+',file='+psdir+'/'+psfil
		print,form='("'+c1+'",a,$)',cret
	      endwhile					;K0="="}

	      if nops eq 0 then begin
		tvlct,red,green,blue,/get
	        owin=!d.name & set_plot,'ps'
	        device,filename=psdir+'/'+psfil,landscape=landscape,$
			color=color,encapsulated=encapsulated
		help,landscape,color,encapsulated
		tvlct,red,green,blue
	        k0='h'
	      endif else begin
		print,'' & print,'changed mind about making a hardcopy?'
		k0=''
	      endelse
	    endif else begin				;X)(not X .. PS?
	      k0=''
	      device,/close & if !D.N_COLORS gt 256 then device,decomposed=0
	      print,'' & print,'hardcopy in file: ',psdir+'/'+psfil
	      if not keyword_set(owin) then begin
		help,!D,/str & return,-1L		;aut X, aut nihil
	      endif
	      set_plot,owin
	    endelse					;!X)
	  end				;		end hardcopy)
	  'r': begin					;reset range
	    print,'' & print,form='("range in ",$)'
	    k2=get_kbrd(1) & k22=strlowcase(k2)
	    case k22 of
	      'x': begin				;X-range
		print,form='("X",$)'
		read,rx0,rx1
		if rx0 lt rx1 then xr=[rx0,rx1] else xr=[rx1,rx0]
	      end
	      'y': begin				;Y-range
		print,form='("Y ",$)'
		read,ry0,ry1
		if ry0 lt ry1 then yr=[ry0,ry1] else yr=[ry1,ry0]
	      end
	      else: begin				;unzoom
		xr=[storange(0),storange(1)]
		yr=[storange(2),storange(3)]
		print,form='(a,$)',cret
	      end
	    endcase
	  end
	  '<': begin					;scroll left
	    rx0=xr(0) & rx1=xr(1) & drx=abs(rx1-rx0)/5. & xr=[rx0-4*drx,rx1-drx]
	  end
	  '>': begin					;scroll right
	    rx0=xr(0) & rx1=xr(1) & drx=abs(rx1-rx0)/5. & xr=[rx1-drx,rx1+4*drx]
	  end
	  '2': begin					;zoom out 2x
	    rx0=xr(0) & rx1=xr(1) & drx=abs(rx1-rx0)/2. & xr=[rx1-drx,rx1+drx]
	  end
	  '^': begin					;extend upwards
	    ry0=yr(0) & ry1=yr(1) & dry=abs(ry1-ry0)/10. & yr=[ry0,ry1+dry]
	  end
	  '_': begin					;extend downwards
	    ry0=yr(0) & ry1=yr(1) & dry=abs(ry1-ry0)/10. & yr=[ry0-dry,ry1]
	  end
	  else: begin
	    print,'	?: general help,'
	    print,'	x: return to cursor control'
	    print,'	z: stop'
	    print,'	q: quit program'
	    print,'	h: hardcopy'
	    print,'	rx: set x-range'
	    print,'	ry: set y-range'
	    print,'	r*: reset (x,y)-range'
	    print,'	<: scroll left'
	    print,'	>: scroll right'
	    print,'	2: zoom out in X by 2X'
	    print,'	^: extend upwards'
	    print,'	_: extend downwards'
	  end
	endcase						;keyboard input}
      endelse				;	END KEYBOARD CONTROL)
    endif				;	end left click)
    ;
    if mbutton eq 2 then begin		;(	middle click
      if jzm gt 0 and wreg eq 0 then begin ;		unzoom
	currange=storange & jzm=0
	xr=[storange(0),storange(1)]
	yr=[storange(2),storange(3)]
      endif else begin			;		toggle log scales
	if keyword_set(xreg) then begin	;			xlog
	  if xlog eq 0 then xlog=1 else xlog=0	;==(xlog=1-xlog)
	  ox=where(xx ge xr(0) and xx le xr(1),mox)
	  if mox ne 0 then ox=where(xx(ox) le 0,mox)
	  if mox ne 0 and xreg gt 0 then xlog=0
	endif
	if keyword_set(yreg) then begin	;			ylog
	  if ylog eq 0 then ylog=1 else ylog=0	;==(ylog=1-ylog)
	  oy=where(yy ge yr(0) and yy le yr(1),moy)
	  if moy ne 0 then oy=where(yy(oy) le 0,moy)
	  if moy ne 0 and yreg gt 0 then ylog=0
	endif
      endelse
    endif				;	end mid click)
    ;
    if mbutton eq 4 then go_on=0	;	right click

  endif else begin			;CLICK)(CLICK+DRAG!

    if mbutton eq 1 then begin		;	left drag (select)
      hh=where(xx ge x0 and xx le x1)
      if hh(0) ne -1 then oo(hh)=1
    endif				;	end left drag
    ;
    if mbutton eq 2 then begin		;	middle drag (deselect)
      z=abs(xx-x0) & ih=(where(z eq min(z)))(0)
      if oo(ih) eq 1 then begin		;		deselect this interval
        o0=ih & o1=ih & do0=0 & do1=-1
        hh=[ih]
        while do0 ne do1 do begin
	  do0=o1-o0
	  oi=o0-1 > 0
	  if oo(oi) eq 1 then begin & o0=oi & hh=[oi,hh] & endif
	  oi=o1+1 < nx
	  if oo(oi) eq 1 then begin & o1=oi & hh=[hh,oi] & endif
	  do1=o1-o0
        endwhile
        hh=hh(uniq(hh,sort(hh)))
        oo(hh)=-1L
        hh=hh-1L
      endif else begin			;		deselect previous
        if hh(0) ne -1 then oo(hh)=-1L
      endelse
    endif				;	end middle drag
    ;
    if mbutton eq 4 then begin		;	right drag (select zoom area)
      izm=1
      if x0 eq x1 then begin & x0=xr(0) & x1=xr(1) & endif	;default
      if y0 eq y1 then begin & y0=yr(0) & y1=yr(1) & endif	;default
      if x0 gt x1 then currange=[x1,x0] else currange=[x0,x1]
      if y0 gt y1 then currange=[currange,y1,y0] else currange=[currange,y0,y1]
      ;
      ;	if drag extends past plot area, extend by half the viewing area
      if currange(0) lt xr(0) then currange(0)=$
	xr(0)-0.5*(xr(1)-xr(0)) > (minx-0.1*dx)
      if currange(1) gt xr(1) then currange(1)=$
	xr(1)+0.5*(xr(1)-xr(0)) < (maxx+0.1*dx)
      if currange(2) lt yr(0) then currange(2)=$
	yr(0)-0.5*(yr(1)-yr(0)) > (miny-0.1*dy)
      if currange(3) gt yr(1) then currange(3)=$
	yr(1)+0.5*(yr(1)-yr(0)) < (maxy+0.1*dy)
      ;
      oplot,[x0,x1,x1,x0,x0],[y0,y0,y1,y1,y1],linestyle=1
    endif					;	end right drag

  endelse				;CLICK+DRAG)

  ;	make sure the plot does not go off to -infinity!
  if xlog eq 1 then begin
    ;ox=where(xx ge xr(0) and xx le xr(1) and yy ge yr(0) and yy le yr(1))
    ox=where(xx ge xr(0) and xx le xr(1),mox)
    if mox ne 0 then begin
      xmax=max(xx(ox))
      if xmax le 0 then begin
	message,'Cannot convert X-axis to log scaling',/info & xlog=0
      endif
      if xr(1) le 0 then xr(1)=1.
      if xr(0) le 0 then xr(0)=derange*xr(1)
    endif else begin
      message,'	no data points in selected window!',/info
      xlog=0
    endelse
  endif
  if ylog eq 1 then begin
    ;oy=where(xx ge xr(0) and xx le xr(1) and yy ge yr(0) and yy le yr(1),moy)
    oy=where(yy ge yr(0) and yy le yr(1),moy)
    if moy ne 0 then begin
      ymax=max(yy(oy))
      if ymax le 0 then begin
	message,'Cannot convert Y-axis to log scaling',/info & ylog=0
      endif
      if yr(1) le 0 then yr(1)=1.
      if yr(0) le 0 then yr(0)=derange*yr(1)
    endif else begin
      message,'	no data points in selected YRANGE!',/info
      ylog=0
    endelse
  endif

  ;	plot
  plot,xx,yy,psym=10,xr=xr,yr=yr,/xs,/ys,xlog=xlog,ylog=ylog,$
	xtitle=xtitle,ytitle=ytitle,title=title
  stample, _extra=e
  oh=where(oo gt 0)
  if oh(0) ne -1 then oplot,xx(oh),0*yy(oh)+0.5*(yr(0)+yr(1)),psym=8
  if mx gt 0 then begin
    if recompy ne 0 then begin
      ymark=0.67*(yr(0)+yr(1)) & ym=0*xm+ymark
      if my gt 0 and my le mx then ym(0:my-1)=marky else $
	ym(0:mx-1)=marky(0:mx-1)
    endif
    for i=0L,mx-1L do begin
      if keyword_set(markall) then begin
        if i eq 0 then oplot,[xm],[ym],psym=psim[0],symsize=ssiz[0],col=kol[0]
      endif else oplot,[xm(i)],[ym(i)],psym=psim(i),symsize=ssiz(i),col=kol(i)
      if keyword_set(legg) then oplot,xm(i)*[1,1],[yr(0),ym(i)],col=kol(i)
      if pltlbl gt 0 then xyouts,xm(i),ym(i),' '+lbl(i),orientation=ori(i),$
	col=kol(i),noclip=0
    endfor
  endif
  if keyword_set(sigmay) then begin
    oplot,xx,yy+su,psym=8 & oplot,xx,yy-sd,psym=8
    for i=0,nx-1 do oplot,xx(i)*[1,1],yy(i)+[su(i),-sd(i)]
  endif
  if oplok eq 'ok' then begin
    for iopl=0,nopl-1 do begin
      xxo=reform(opx(iopl,*)) & yyo=reform(opy(iopl,*))
      if nopl eq 1 then begin
	xxo=opx & yyo=opy
      endif
      oplot,xxo,yyo,psym=10,color=opcol(iopl)
    endfor
  endif
  if izm gt 0 then begin
    x0=currange(0)>xr(0) & x1=currange(1)<xr(1)
    y0=currange(2)>yr(0) & y1=currange(3)<yr(1)
    oplot,[x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],linestyle=1
  endif
endwhile					;end selection}

;	shorten the array
oo=where(oo ge 0)

;	optional outputs
xo=x1		;last click x-pos
yo=y1		;last click y-pos
oxr(*)=xr		;final x-range
oyr(*)=yr		;final y-range

return,oo
end
