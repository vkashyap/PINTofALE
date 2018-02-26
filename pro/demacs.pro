function demacs,logt,dem0=dem0,logt0=logt0,norm=norm,slope=slope,$
	group=group,igroup=igroup, _extra=e
;+
;function	demacs
;	a DEM editor that allows interactively "tweaking" a specified DEM
;
;syntax
;	dem=demacs(logt,dem0=dem0,logt0=logt0,norm=norm,slope=slope,$
;	group=group,igroup=igroup, PLOT_KEYWORDS)
;
;parameters
;	logt	[INPUT; required] log10(Temperature [K]) at which DEM
;		is to be determined
;
;keywords
;	dem0	[INPUT] initial DEM (just for guidance).
;		* if size does not match LOGT0, gets stretched to suit.
;		* overrides SLOPE, but not NORM
;		* if max(DEM0)<100, assumed to be in log form -- i.e.,
;		  plots are in "linear" form and +vity is not enforced.
;	logt0	[INPUT] logT at which DEM0 are defined.  if not given,
;		assumed to be LOGT.
;	norm	[INPUT] if set, first element of initial DEM is set to NORM
;	slope	[INPUT; default=1] if given, generates initial
;		DEM=NORM*T0^(SLOPE)
;		* if DEM0 is specified, SLOPE is ignored
;	group	[OUTPUT] long-int array of same size as LOGT containing
;		grouping information
;	igroup	[OUTPUT] array of same size as GROUP, containing index
;		of "representative element" of each group
;	_extra	[INPUT] allows passing defined keywords (e.g., YRANGE)
;		to PLOT.
;
;restrictions
;	requires X-windows, display, and mouse capability.
;
;history
;	vinay kashyap (Feb97)
;	corrected log behavior, added _extra to PLOT (VK; Apr97)
;	cleaned up log behavior, added keystrokes *,/,v,^,z (VK; Aug98)
;	button press status now stored in !MOUSE, not !ERR (VK; Apr09)
;-

;	usage
nt=n_elements(logt)
if nt lt 3 then begin
  print,'Usage: dem=demacs(logt,dem0=dem0,logt0=logt0,norm=norm,slope=slope,$'
  print,'       group=group,igroup=igroup,PLOT_KEYWORDS)'
  print,'  DEM editor' & return,0.
endif

;	check keywords
;
nt0=n_elements(logt0)			;is LOGT0 defined?
if nt0 eq 0 then begin			;nope.
  lt0=logt & nt0=nt			;use LOGT
endif else lt0=[logt0]			;make sure 'tis an array
;
if keyword_set(dem0) then d0=[dem0] else begin	;DEM0 as an array
  if n_elements(slope) eq 0 then slope=1.	;no DEM0, so what slope?
  d0=10.D^(slope*lt0+10)			;compute DEM0
endelse
if n_elements(d0) eq 1 then d0=dblarr(nt0)+d0(0)
nd0=n_elements(d0)
;if max(d0) lt 100 then d0=10.^(d0)		;in log form?
if max(d0) lt 100 then ylog=0 else ylog=1	;in log form?
;	NOTE: ylog=0 says that the *input* is in LOG form!
if nd0 ne nt0 then begin
  t0=lt0(0) & t1=lt0(nt0-1) & ti=findgen(nd0)*(t1-t0)/(nd0-1.)+t0
  d0=interpol(d0,ti,lt0)			;stretch to fit
endif
;
if keyword_set(norm) then d0=d0*norm/d0(0)	;input normalization
;
;	groups of points
if n_elements(group) eq nt then grp=group else grp=lindgen(nt)
;	representative point in each group
if n_elements(igroup) eq nt then igrp=igroup else igrp=grp

;	initialize
mbt=0					;buttons pressed
kbrd=0					;input from keyboard
go_on=1					;active mouse
act='?'					;action to take
factor=1.1				;multiplicative "nudge" factor
togint=1				;spline/poly/linear interpolation
;dem=interpol(d0,lt0,logt)>1		;output
dem=(spline(lt0,d0,logt,1.)>min(d0))<max(d0)		;output
dd=dem & olddd=dd & ig0=0 & ig1=0 & jg0=ig0 & jg1=ig0 & ip=0
hilfe=[	'?: help;	q: quit;	x: cursor;	z: halt',$
	'h,l: move left/right by 1 element;',$
	'b,w: by 4;	B,W: by 16',$
	'0: move to first element;		$: move to last element',$
	'j,k: divide/multiply by FACTOR;	J,K: by 2;	-,+: by 10',$
	'f: set FACTOR (default=1.1);		F: reset VALUE of FACTOR',$
	'n: set normalization',$
	'(,): start/end grouping;		{,}: start/end ungroup',$
	'.: select "representative" element of given group',$
	'*,/: move block up/down by FACTOR;	v,^: ask for factor',$
	't: toggle spline/polynomial/linear interpolation',$
	'u: undo last change;			U: remove point from group']

;	edit DEM
print,'		CLICK & DRAG mouse buttons'
print,'	LEFT: group points; MIDDLE: ungroup points; RIGHT: control via keyboard'
print,'		CLICK mouse buttons'
print,'	LEFT: position; MIDDLE: undo last; RIGHT: quit'
print,'	(click on right margin to toggle interpolation style)'

idlvers=float(!VERSION.release) & if idlvers le 0 then idlvers=7

while go_on eq 1 do begin				;{interactive edit

  ;	plot
  plot,logt,dem,/xs,/ys,ylog=ylog,xtitle='log!d10!n(T [K])',ytitle='DEM',$
	_extra=e
  oplot,logt,olddd,linestyle=1
  oplot,[logt(ip)],[dem(ip)],psym=4
  xyouts,logt(ip),1.15*max(dem),'logT='+string(logt(ip),'(f5.2)'),align=0.5
  ;	plot the grouped point pointers
  ugrp=grp(uniq(grp,sort(grp))) & ngrp=n_elements(ugrp)
  if ngrp lt nt then begin
    for i=0,ngrp-1 do begin
      oo=where(grp eq ugrp(i),moo)
      if moo gt 1 then begin
	xyouts,logt(oo(0)),dem(oo(0)),'(',align=0.1,chars=2,chart=2
	xyouts,logt(oo(moo-1)),dem(oo(moo-1)),')',align=-0.1,chars=2,chart=2
	oplot,[logt(igrp(oo(0)))],[dem(igrp(oo(0)))],psym=5
      endif
    endfor
  endif
  recomp=0

  if kbrd eq 0 then begin		;{click &/or drag
    cursor,x0,y0,/change,/data
    if idlvers lt 5 then mbt=!err else $	;OLD IDL
    	mbt=!MOUSE.BUTTON
    while mbt eq 0 do begin	;display cursor position until button up
      cursor,x0,y0,/change,/data
      if idlvers lt 5 then mbt=!err	;OLD IDL
      	mbt=!MOUSE.BUTTON
      tmpx=min(abs(logt-x0),ii) & tmpx=logt(ii)
      c1='logT='+string(tmpx,'(f6.3)')+' DEM='+string(y0,'(g10.4)')
      print,form='($,"'+c1+'",a)',string("15b)
    endwhile
    cursor,x1,y1,/up,/data
    drag=abs(x1-x0)+abs(y1-y0)
    xmin=x0 & if x1 lt xmin then xmin=x1
    xmax=x1 & if x0 gt xmax then xmax=x0

    if drag gt 0 then begin		;(DRAG

      if mbt eq 1 then begin		;	left drag
	tmp=min(abs(logt-xmin),ig0) & tmp=min(abs(logt-xmax),ig1)
	act='>'
      endif
      if mbt eq 2 then begin		;	middle drag
	tmp=min(abs(logt-xmin),jg0) & tmp=min(abs(logt-xmax),jg1)
	act=']'
      endif
      if mbt eq 4 then begin		;	right drag
	kbrd=1 & print,''
	print,'?: help;	q: quit;	x: cursor control'
        act=get_kbrd(1)
      endif

    endif else begin			;DRAG)(CLICK

      if mbt eq 1 then begin		;	left click
	if xmin ge logt(0) and xmin le logt(nt-1) then begin
	  tmp=min(abs(logt-xmin),ip)
	  act='g'
	endif else begin
	  act='l'
	  recomp=1
	endelse
	gg=ip
      endif
      if mbt eq 2 then begin		;	middle click
	if xmin ge logt(0) or xmin le logt(nt-1) then begin
	  tmp=min(abs(logt-xmin),ip)
	  act='u'
	endif else begin
	  act='h'
	  recomp=1
	endelse
      endif
      if mbt eq 4 then begin 		;	right click
	act='q' & print,''
      endif
      ;
      if xmin gt logt(nt-1) then togint=togint+1
      ;if xmin lt logt(0) then togint=togint-1

    endelse				;CLICK)

  endif else begin			;}{get input from keyboard
    act=get_kbrd(1)
  endelse				;kbrd=1}

  case act of				;{action!
    'h': ip=(ip-1)>0			;	hop left
    'l': ip=(ip+1)<(nt-1)		;	hop right
    ' ': ip=(ip+1)<(nt-1)		;	hop right
    'b': ip=(ip-4)>0			;	skip left
    'w': ip=(ip+4)<(nt-1)		;	skip right
    'B': ip=(ip-16)>0			;	jump left
    'W': ip=(ip+16)<(nt-1)		;	jump right
    '0': ip=0				;	move to 1st element
    '$': ip=nt-1			;	move to last element
    'j': begin				;	divide
      if ylog eq 1 then dd(igrp(ip))=dd(igrp(ip))/factor else $
	dd(igrp(ip))=dd(igrp(ip))-alog10(factor)
      recomp=1
    end
    'J': begin				;	divide
      if ylog eq 1 then dd(igrp(ip))=dd(igrp(ip))/2 else $
	dd(igrp(ip))=dd(igrp(ip))-alog10(2.)
      recomp=1
    end
    '-': begin				;	divide
      if ylog eq 1 then dd(igrp(ip))=dd(igrp(ip))/10 else $
	dd(igrp(ip))=dd(igrp(ip))-1.
      recomp=1
    end
    'k': begin				;	multiply
      if ylog eq 1 then dd(igrp(ip))=dd(igrp(ip))*factor else $
	dd(igrp(ip))=dd(igrp(ip))+alog10(factor)
      recomp=1
    end
    'K': begin				;	multiply
      if ylog eq 1 then dd(igrp(ip))=dd(igrp(ip))*2 else $
	dd(igrp(ip))=dd(igrp(ip))+alog10(2.)
      recomp=1
    end
    '+': begin				;	multiply
      if ylog eq 1 then dd(igrp(ip))=dd(igrp(ip))*10 else $
	dd(igrp(ip))=dd(igrp(ip))+1.
      recomp=1
    end
    '*': begin				;	multiply whole group
      og=where(grp eq grp(ip),mog)
      if mog eq 0 then message,'bug!'
      if ylog eq 1 then dd(og)=dd(og)*factor else $
	dd(og)=dd(og)+alog10(factor)
      olddd=dem & dem=dd
    end
    '/': begin				;	divide whole group
      og=where(grp eq grp(ip),mog)
      if mog eq 0 then message,'bug!'
      if ylog eq 1 then dd(og)=dd(og)/factor else $
	dd(og)=dd(og)-alog10(factor)
      olddd=dem & dem=dd
    end
    'v': begin				;	divide whole group
      og=where(grp eq grp(ip),mog) & tmp=0.
      if mog eq 0 then message,'bug!'
      while tmp le 0 do $
        read,prompt='divide this group by what factor? ',tmp
      if ylog eq 1 then dd(og)=dd(og)/tmp else $
	dd(og)=dd(og)-alog10(tmp)
      olddd=dem & dem=dd
    end
    '^': begin				;	multiply whole group
      og=where(grp eq grp(ip),mog) & tmp=0.
      if mog eq 0 then message,'bug!'
      while tmp le 0 do $
        read,prompt='multiply this group by what factor? ',tmp
      if ylog eq 1 then dd(og)=dd(og)*tmp else $
	dd(og)=dd(og)+alog10(tmp)
      olddd=dem & dem=dd
    end
    'f': begin				;	x by...
      read,prompt='nudge factor: ',factor
      if factor eq 0 then factor=1.1
      if factor lt 1 then factor=abs(1./factor)
    end
    'F': begin				;	ask what value
      tmp=0.
      while tmp le 0 do $
        read,prompt='DEM(logT='+strtrim(logt(igrp(ip)),2)+') =? ',tmp
      dd(igrp(ip))=tmp
      recomp=1
    end
    'g': begin				;	fix value
      dd(ip)=y0
      recomp=1
    end
    'n': begin				;	normalization
      tmp=0.
      while tmp le 0 do $
        read,prompt='set normalization @ logT='+strtrim(logt(ip),2)+': ',tmp
      if ylog eq 1 then dd=dd*tmp/dd(ip) else dd=dd+tmp-dd(ip)
      recomp=1
    end
    '(': begin				;	start grouping
      ig0=ip & act2=''
      while act2 ne ')' do begin		;(group
	act2=get_kbrd(1)
	case act2 of
          'h': ip=(ip-1)>0			;	move left
          'l': ip=(ip+1)<(nt-1)			;	move right
          ' ': ip=(ip+1)<(nt-1)			;	move right(
	  else: print,'h: left; l: right; ): exit'
	endcase
	oplot,[logt(ip)],[dd(ip)],psym=6 & wait,0.02
	oplot,[logt(ip)],[dd(ip)],psym=6,col=0
	oplot,[logt(ip)],[dd(ip)],psym=7
      endwhile					;end group)
      act='>' & ig1=ip
    end
    '{': begin					;	start ungrouping
      jg0=ip & act2=''
      while act2 ne '}' do begin		;(ungroup
	act2=get_kbrd(1)
	case act2 of
          'h': ip=(ip-1)>0			;	move left
          'l': ip=(ip+1)<(nt-1)			;	move right
          ' ': ip=(ip+1)<(nt-1)			;	move right{
	  else: print,'h: left; l: right; }: exit'
	endcase
	oplot,[logt(ip)],[dd(ip)],psym=6 & wait,0.02
	oplot,[logt(ip)],[dd(ip)],psym=6,col=0
	oplot,[logt(ip)],[dd(ip)],psym=7
      endwhile					;end ungroup)
      act=']' & jg1=ip
    end
    'u': begin				;	undo last change
      dem=olddd & olddd=dd
      ;recomp=1
    end
    'U': begin				;	ungroup just one point
      jg0=ip & jg1=ip & act=']'
    end
    't': begin				;	toggle spline/poly/lin
      print,'1:cubic spline, 2:polynomial, 3:linear'
      act2='' & act2=get_kbrd(1)
      if act2 eq ' ' or act2 eq '' then togint=togint+1 else togint=fix(act2)
      recomp=1
    end
    '.': gg=ip				;	representative point
    '?': print,hilfe			;	help
    'x': kbrd=0				;	cursor
    'q': go_on=0			;	quit
    'z': stop,'halting'			;	stop
    '>':				;	same as group
    ']':				;	same as ungroup
    else: print,'?: help; x: cursor; q: quit'
  endcase				;action}

  oplot,[logt(ip)],[dem(ip)],psym=4

  if act eq '>' then begin		;group
    if ig1 gt ig0 then begin
      grp(ig0:ig1)=grp(ig0)
      igrp(ig0:ig1)=igrp(ig0)
      oplot,logt(ig0:ig1),dem(ig0:ig1),psym=1
    endif
  endif

  if act eq '.' or act eq 'g' then begin	;representative point
    oo=where(grp eq grp(ip))
    igrp(oo)=gg
    oplot,logt(ig0:ig1),dem(ig0:ig1),psym=4,col=0
    oplot,logt(ig0:ig1),dem(ig0:ig1),psym=1
    oplot,[logt(ip)],[dem(ip)],psym=4
  endif

  if act eq ']' then begin		;ungroup
    if jg1 ge jg0 then begin
      grp(jg0:jg1)=lindgen(jg1-jg0+1)+jg0
      igrp(jg0:jg1)=grp(jg0:jg1)
      oplot,logt(jg0:jg1),dem(jg0:jg1),psym=1,col=0
      oplot,logt(jg0:jg1),dem(jg0:jg1),psym=4,col=0
      oplot,logt(jg0:jg1),dem(jg0:jg1)
      if jg1 lt nt-1 and jg0 gt 1 then begin
	if grp(jg0-1) eq grp(jg1+1) then begin
	  o1=where(grp eq grp(jg0-1) and lindgen(nt) lt jg0)
	  o2=where(grp eq grp(jg1+1) and lindgen(nt) gt jg1)
	  grp(o2)=jg1+1
	  if igrp(jg0-1) ge jg0 then igrp(o1)=jg0-1
	  if igrp(jg1+1) le jg1 then igrp(o2)=jg1+1
	endif
      endif
    endif
  endif

  ilogt=igrp(uniq(igrp,sort(igrp)))
  if n_elements(ilogt) lt 3 then recomp=0

  if recomp eq 1 then begin			;(recompute
    lt0=logt(ilogt) & d0=dd(ilogt)

    if total(abs(dem(ilogt)-dd(ilogt)))/total(dd(ilogt)) gt 1e-7 then olddd=dem

    if togint gt 4 then togint=1
    case togint of
      1: begin				;cubic spline
        ;dem=spline(lt0,d0,logt,0.1)
        dem=spline(lt0,d0,logt,1.)
      end
      2: begin				;polynomial
        if ylog eq 1 then dem=10.D^(spline(lt0,alog10(d0),logt,5)) else $
	  dem=spline(lt0,d0,logt,5)
      end
      3: begin				;linear
        dem=interpol(d0,lt0,logt)
      end
      else: begin			;spline interpolate
        y2=spl_init(lt0,d0,/double)
        dem=spl_interp(lt0,d0,y2,logt,/double)
      end
    endcase

    oo=where(dem gt 0,moo)
    if moo gt 0 then begin
      if ylog eq 1 then dem=dem>(min(dem(oo)))
      dd=dem
    endif else begin
      dd=olddd & dem=dd
    endelse
  endif						;recomputed)

endwhile						;go_on=1}

group=grp & igroup=igrp & return,dem			;outputs
end
