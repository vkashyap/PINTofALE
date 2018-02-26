pro poaintsym, symb, psize=psize,pthick=pthick,pfill=pfill,$
	pie=pie,halftri=halftri,npoint=npoint,aspect=aspect,angle=angle,$
	adent=adent,phase=phase,verbose=verbose, _extra=e
;+
;procedure	poaintsym
;	Define useful plotting symbols not in the standard !PSYM
;	definitions.  After symbol has been defined with this
;	routine, a plotting command should follow with either
;	PSYM = 8 or !P.PSYM = 8 (see USERSYM).
;
;syntax
;	poaintsym,symb,psize=psize,pthick=pthick,/pfill,color=color,$
;	pie=pie,/halftri,npoint=npoint,aspect=aspect,angle=angle,adent=adent,$
;	verbose=verbose
;
;parameters
;	symb	[INPUT; required] string describing the type of
;		symbol to define.  The following are accepted:
;		  'circle'
;		  'dash'
;		  'star' 'starfish' 'asterisk'
;		  'mercury' 'venus' 'earth' 'moon' 'mars'
;		  'jupiter' 'saturn' 'uranus' 'neptune' 'pluto'
;		  'square' 'diamond' 'box' 'spoke'
;		  'triangle' 'down triangle' 'right triangle' 'left triangle'
;		  'down arrow' 'up arrow' 'left arrow' 'right arrow'
;		  'upper limit' 'lower limit'
;		  'tilde ANGLE OVERFLOW'
;		* if something unknown is specified, the program
;		  prints a list of accepted options and quits
;		  without making a change
;		* special cases:
;		  -- STARFISH: is 5-pointed
;		  -- ASTERISK: is (NPOINT/2)-pointed and uses keyword
;		     ADENT to figure out indentation
;		  -- EARTH: forces PFILL=0
;		  -- MOON: use ADENT to pull the horns out, and PHASE to
;		     control thickness
;		  -- URANUS: forces PFILL=0
;		  -- PLUTO: forces PFILL=0
;		  -- SPOKE: essentially DIAMOND with an aspect x:y::1:10
;		  -- UPPER/LOWER LIMIT: synonyms for DOWN/UP ARROW respec.
;		  -- TILDE: set ANGLE (default=0) to rotate the symbol,
;		     and OVERFLOW (default=1) to pull the ends out a bit
;
;keywords
;	psize	[INPUT] size of the plotting symbol in
;		multiples of the standard size (see SYMSIZE); does
;		not need to be an integer
;	pthick	[INPUT] thickness of symbol being drawn
;	pfill	[INPUT] if set, fills the symbol
;		* does not affect arrows or character symbols
;	pie	[INPUT] range of angles in [degrees] to include while
;		making a circle
;		* if scalar, includes from PIE[0] to 360.
;		* if vector, includes from PIE[0] to PIE[1]
;	halftri	[INPUT] if set, lops off the triangle such that only the
;		part that protrudes past the axes (above x-axis for up
;		triangle, below x-axis for down triangle, etc.) is plotted
;	npoint	[INPUT] number of vertices in the symbol
;		* note that usersym maximum is 49
;		* currently, only matters to TILDE
;	aspect	[INPUT] changes the aspect ratio of the symbols
;		by scaling the axes prior to rotation
;		* if +ve, scales X-axis
;		* if -ve, scales Y-axis
;		* if array, only the first element is used
;		* if 0, nothing is done
;	rotate	[INPUT] angle by which to rotate the symbol [deg]
;	adent	[INPUT] the magnitude of the indentation for STARFISH
;		and ASTERISK and the extent of the horns for MOON
;		* default is 0.4 for STARFISH and ASTERISK,
;		  and 10 for MOON
;	phase	[INPUT] controls thickness of moon
;		* default is 0
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		USERSYM: COLOR
;
;method
;	Appropiate X,Y vectors are used to define the symbol and passed
;	to USERSYM
;
;examples
;	Plot with filled stars as the symbol, twice the default size
;	  POAINTSYM,'star',psize=2,/pfill & PLOT,findgen(10),PSYM=8
;	Plot with open circle as the symbol
;	  POAINTSYM,'circle' & PLOT,findgen(10),PSYM=8
;	Plot filled color circles
;	  POAINTSYM,'circle',/pfil,psize=3,col=2 & PLOT,findgen(10),PSYM=8
;	Plot arrows all around
;	  POAINTSYM,'left arrow',pthick=2,psize=2 & plot,findgen(11),psym=8
;	  POAINTSYM,'up arrow',pthick=3,psize=3 & oplot,findgen(11),psym=8
;	  POAINTSYM,'right arrow',pthick=4,psize=4 & oplot,findgen(11),psym=8
;	  POAINTSYM,'down arrow',pthick=5,psize=5 & oplot,findgen(11),psym=8
;	Plot tildes
;	  POAINTSYM,'tilde',pthick=5,psize=5 & plot,findgen(11),psym=8
;	  POAINTSYM,'tilda 90',pthick=5,psize=5 & plot,findgen(11),psym=8
;	  POAINTSYM,'tilda 90 3',pthick=5,psize=5 & plot,findgen(11),psym=8
;	  POAINTSYM,'tilda 0 20',pthick=2,psize=5,npoint=49 & plot,findgen(11),psym=8
;
;history
;	modified from JDSYM (Jeremy Drake 2001; JJD added square and
;	  diamond symbols to PLOTSYM [ASTROLib; Wayne Landsman 1992])
;	  by Vinay Kashyap (Oct03)
;	added option TILDE and keyword NPOINT (VK; Feb04)
;	added keywords ASPECT, ROTATE, and ADENT; added options SPOKE,
;	  MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS,
;	  NEPTUNE, PLUTO (VK; Jun05)
;	added keyword PIE (VK; Aug12)
;	added option DASH (VK; Jul15)
;-

;On_error,2

;	usage
ok='ok' & np=n_params() & ns=n_elements(symb) & ss=size(symb,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SYMB: undefined' else $
  if ns gt 1 then ok='SYMB cannot be a vector' else $
   if ss ne 7 then ok='SYMB must be a character string'
if ok ne 'ok' then begin
  print,'Usage: poaintsym,symb,psize=psize,pthick=pthick,/pfill,color=color,$'
  print,'       pie=pie,/halftri,npoint=npoint,aspect=aspect,angle=angle,adent=adent,$
  print,'       phase=phase,verbose=verbose'
  print,'  define useful plotting symbols not in standard !PSYM definitions'
  if np ne 0 then message,ok,/informational
  return
endif

;	figure out keywords
psiz=1.0 & if keyword_set(psize) then psiz=psize[0]
pthk=1.0 & if keyword_set(pthick) then pthk=pthick[0]
pfil=0 & if keyword_set(pfill) then pfil=1
npt=24L & if keyword_set(npoint) then npt=long(npoint[0])
if npt lt 2 or npt gt 49 then begin
  message,'Number of points must be between 2 and MAX<50',/informational
  if npt lt 2 then npt=2L
  if npt ge 50 then npt=49L
endif
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
aspx=1. & aspy=1.
if keyword_set(aspect) then begin
  if aspect[0] gt 0 then aspx=aspect[0]
  if aspect[0] lt 0 then aspy=abs(aspect[0])
endif
ang=0. & if keyword_set(angle) then ang=angle[0]*!pi/180.

;	figure out what SYMB contains
ss=strtrim(strlowcase(symb[0]),2) & s1=strmid(ss,0,1) & csymb=ss
if s1 eq 'c' then begin				;(C*
  if strpos(ss,'ci',0) ge 0 then csymb='CIRCLE'	;expand to CIR.. if other shapes are coded
endif						;CIRCLE)
if s1 eq 'b' then begin				;(B*
  if strpos(ss,'bo',0) ge 0 then csymb='SQUARE'
endif						;BOX)
if s1 eq 'd' then begin				;(D*
  if strpos(ss,'da',0) ge 0 then csymb='DASH'
  if strpos(ss,'di',0) ge 0 then csymb='DIAMOND'
endif						;DASH,DIAMOND)
if s1 eq 'm' then begin				;(M*
  if strpos(ss,'me',0) ge 0 then csymb='MERCURY'
  if strpos(ss,'ma',0) ge 0 then csymb='MARS'
  if strpos(ss,'mo',0) ge 0 then csymb='MOON'
endif						;MERCURY, MARS, MOON)
if s1 eq 'v' then begin				;(V*
  if strpos(ss,'ve',0) ge 0 then csymb='VENUS'
endif						;VENUS)
if s1 eq 'e' then begin				;(E*
  if strpos(ss,'ea',0) ge 0 then csymb='EARTH'
endif						;EARTH)
if s1 eq 'j' then begin				;(J*
  if strpos(ss,'ju',0) ge 0 then csymb='JUPITER'
endif						;JUPITER)
if s1 eq 's' then begin				;(S*
  if strpos(ss,'sa',0) ge 0 then csymb='SATURN'
  if strpos(ss,'st',0) ge 0 then csymb='STAR'
  if strpos(ss,'starf',0) ge 0 then csymb='STARFISH'
  if strpos(ss,'sp',0) ge 0 then csymb='SPOKE'
  if strpos(ss,'sq',0) ge 0 then csymb='SQUARE'
endif						;SATURN, STAR, STARFISH SPOKE, SQUARE)
if s1 eq 'a' then begin				;(A*
  if strpos(ss,'as',0) ge 0 then csymb='ASTERISK'
endif						;ASTERISK)
if s1 eq 'u' then begin				;(U*
  if strpos(ss,'ur',0) ge 0 then csymb='URANUS'
endif						;URANUS)
if s1 eq 'n' then begin				;(N*
  if strpos(ss,'ne',0) ge 0 then csymb='NEPTUNE'
endif						;NEPTUNE)
if s1 eq 'p' then begin				;(P*
  if strpos(ss,'pl',0) ge 0 then csymb='PLUTO'
endif						;PLUTO)
if strpos(ss,'tr',0) ge 0 then begin		;(TR*
  csymb='UP TRIANGLE'
  if strpos(ss,'do',0) ge 0 then csymb='DOWN TRIANGLE'
  if strpos(ss,'rig',0) ge 0 then csymb='RIGHT TRIANGLE'
  if strpos(ss,'lef',0) ge 0 then csymb='LEFT TRIANGLE'
endif						;TRIANGLE)
if strpos(ss,'arr',0) ge 0 then begin		;(ARR*
  if strpos(ss,'do',0) ge 0 then csymb='DOWN ARROW'
  if strpos(ss,'up',0) ge 0 then csymb='UP ARROW'
  if strpos(ss,'rig',0) ge 0 then csymb='RIGHT ARROW'
  if strpos(ss,'lef',0) ge 0 then csymb='LEFT ARROW'
endif						;ARROW)
if strpos(ss,'lim',0) ge 0 then begin		;(LIM*
  if strpos(ss,'upp',0) ge 0 then csymb='DOWN ARROW'	;upper limit == down arrow
  if strpos(ss,'low',0) ge 0 then csymb='UP ARROW'	;lower limit == up arrow
endif						;LIMIT==ARROW)
if strpos(ss,'til',0) ge 0 then begin		;(TIL*
  csymb='TILDE'
  cc=strsplit(ss,'[ '+STRING(9B)+']+',/extract) & ncc=n_elements(cc)
  tilde_ang=0. & if ncc gt 1 then tilde_ang=float(cc[1])
  tilde_ovr=1. & if ncc gt 2 then tilde_ovr=float(cc[2])
  if vv gt 100 then stop
endif						;TILDE)
if vv gt 0 then message,'Setting point symbol to: '+csymb,/informational

case csymb of					;{CSYMB
  'CIRCLE': begin
    thtdeg = 360.*findgen(49)/48.
    ii=lindgen(49)
    if keyword_set(pie) then begin
      npie=n_elements(pie)
      piemin=(pie[0] mod 360.)
      if npie gt 1 then piemax=(pie[1] mod 360.) else piemax=360.
      if piemin gt piemax then begin
        tmp=piemin & piemin=piemax & piemax=tmp
      endif
      oo=where(thtdeg ge piemin and thtdeg le piemax,moo)
      thtdeg=[piemin,thtdeg[oo],piemax]
      thtdeg=thtdeg[uniq(sort(thtdeg))]
    endif
    tht = thtdeg*!pi/180.
    xarr = aspx*psiz*cos(tht)  &  yarr = aspy*psiz*sin(tht)
  end
  'DASH': begin
    xarr=psiz*0.5*[-1,1] & yarr=0.*[1,1]
  end
  'STAR': begin
    r = psiz   
    tht = (720. / 5*findgen(6) + 45) / !RADEG  ;Define star angles every 144 deg
    xarr = aspx*r*cos(tht)  & yarr = aspy*r*sin(tht)
  end
  'STARFISH': begin
    indent=0.4 & if keyword_set(adent) then indent=adent[0]
    tht=2.*!pi*findgen(11)/10.+!pi/2.
    xarr=cos(tht) & yarr=sin(tht)
    ii=2*lindgen(5)+1 & xarr[ii]=indent*xarr[ii] & yarr[ii]=indent*yarr[ii]
    xarr = xarr * psiz*aspx
    yarr = yarr * psiz*aspy
  end
  'ASTERISK': begin
    indent=0.4 & if keyword_set(adent) then indent=adent[0]
    mpt = 2*(npt/2)+1
    tht=2.*!pi*findgen(mpt)/(mpt-1.)+!pi/2.
    xarr=cos(tht) & yarr=sin(tht)
    ii=2*lindgen(mpt/2)+1 & xarr[ii]=indent*xarr[ii] & yarr[ii]=indent*yarr[ii]
    xarr = xarr * psiz*aspx
    yarr = yarr * psiz*aspy
  end
  'MERCURY': begin
    tht = 2*!PI*findgen(31)/30.-(!pi/2.) & xarr1=cos(tht) & yarr1=sin(tht)
    jnk=min(abs(tht-70.*!pi/180.),i1) & jnk=min(abs(tht-110.*!pi/180.),i2)
    xarrh=[0.,0.2,0.4,0.2,0] & yarrh=[0.,0.2,0.5,0.2,0]
    xarrc=[0.0, 0.0, -0.7,0.0, 0.7, 0.0] +xarr1[0]
    yarrc=[-1.2,-0.6,-0.6,-0.6,-0.6,-0.6]+yarr1[0]
    xarr = [xarrc, xarr1[0:i1],xarr1[i1]+xarrh,xarr1[i1:i2],xarr1[i2]-xarrh,xarr1[i2:*]] * psiz*aspx
    yarr = [yarrc, yarr1[0:i1],yarr1[i1]+yarrh,yarr1[i1:i2],yarr1[i2]+yarrh,yarr1[i2:*]] * psiz*aspy
  end
  'VENUS': begin
    tht = 2*!PI*findgen(31)/30.-(!pi/2.) & xarr1=cos(tht) & yarr1=sin(tht)
    xarrc=[0.0, 0.0, -0.7,0.0, 0.7, 0.0] +xarr1[0]
    yarrc=[-1.2,-0.6,-0.6,-0.6,-0.6,-0.6]+yarr1[0]
    xarr = [xarrc, xarr1] * psiz*aspx
    yarr = [yarrc, yarr1] * psiz*aspy
  end
  'EARTH': begin
    tht1 = (!PI/2.)*findgen(11)/10. & xarr1=cos(tht1) & yarr1=sin(tht1)
    tht2 = tht1 + !pi/2. & xarr2=cos(tht2) & yarr2=sin(tht2)
    tht3 = tht2 + !pi/2. & xarr3=cos(tht3) & yarr3=sin(tht3)
    tht4 = tht3 + !pi/2. & xarr4=cos(tht4) & yarr4=sin(tht4)
    xarr = [xarr2[10],xarr1,xarr3[10],xarr1[10],xarr2,xarr3,xarr4] * psiz*aspx
    yarr = [yarr2[10],yarr1,yarr3[10],yarr1[10],yarr2,yarr3,yarr4] * psiz*aspy
    pfil = 0
  end
  'MOON': begin
    ff=0.8 & if keyword_set(phase) then ff=phase[0]
    if ff lt 0 then ff=0.9 & if ff gt 1 then ff=1./ff
    dx=0.1 & if keyword_set(adent) then dx=adent[0]
    dd=(1-ff)+dx
    ;	R*cos(theta)=r*cos(alpha)+delta
    ;	R*sin(theta)=r*sin(alpha)
    ;	R = delta*cos(pi-theta) + r*cos(theta-alpha)
    ;	==> sin(theta)=f*sin(alpha), f=r/R
    ;	==> sqrt(1-f^2*sin(alpha)^2)=f*cos(alpha)+delta/R
    ;	==> 1=f^2+2*f*(delta/R)*cos(alpha)+(delta/R)^2
    ;	==> cos(alpha) ==
    cc=(1.-ff^2-dd^2)/(2.*ff*dd)
    alpha=acos(cc) & theta=asin(ff*sin(alpha))
    thout=findgen(24)*((2*!pi-theta)-theta)/23.+theta
    thinn=reverse(findgen(24)*((2*!pi-alpha)-alpha)/23.+alpha)
    xo=cos(thout) & yo=sin(thout)
    xi=ff*cos(thinn)+dd & yi=ff*sin(thinn)
    xarr=[xo,xi] * psiz*aspx
    yarr=[yo,yi] * psiz*aspy
  end
  'MARS': begin
    tht = 2*!PI*findgen(41)/40. & xarr1=cos(tht) & yarr1=sin(tht)
    jnk=min(abs(tht-!pi/3.),imn)
    xarr2 = [ 0, 2, 1.4, 2, 1.4, 2]
    yarr2 = [ 0, 0, 0.5, 0, -.5, 0]
    xarr3 = xarr2*cos(-!pi/3.) + yarr2*sin(-!pi/3.)
    yarr3 = -xarr2*sin(-!pi/3.) + yarr2*cos(-!pi/3.)
    xarr = [xarr1[0:imn], xarr3+xarr1[imn], xarr1[imn:*]] * psiz*aspx
    yarr = [yarr1[0:imn], yarr3+yarr1[imn], yarr1[imn:*]] * psiz*aspy
  end
  'JUPITER': begin
    tht = (!PI/2.)*findgen(41)/40.-!pi/4. & xarr1=cos(tht) & yarr1=sin(tht)
    xarrc=[1.0,0.7,0.7,0.7,0.7,0.]+xarr1[0]
    yarrc=[0.0,0.0,-0.6,0.7,0.0,0.]+yarr1[0]
    xarr = [xarrc,xarr1] * psiz*aspx
    yarr = [yarrc,yarr1] * psiz*aspy
  end
  'SATURN': begin
    tht = (3.*!PI/2.)*findgen(37)/36.-3*!pi/4. & xarr1=0.5*cos(tht) & yarr1=0.5*sin(tht)
    xarr2 = xarr1*cos(-!pi/4.)+yarr1*sin(-!pi/4.)
    yarr2 = -xarr1*sin(-!pi/4.)+yarr1*cos(-!pi/4.)
    xarrc = [0.0,0.0,-0.5,0.0,0.5,0.0,0.0] + xarr2[36]
    yarrc = [1.5,1.1, 1.1,1.1,1.1,1.1,0.0] + yarr2[36]
    xarrt = [0.0,0.1,0.2] + xarr2[0]
    yarrt = [0.0,-0.1,-0.2] + yarr2[0]
    xarr = [xarrc,reverse(xarr2),xarrt] * psiz*aspx
    yarr = [yarrc,reverse(yarr2),yarrt] * psiz*aspy
  end
  'URANUS': begin
    tht = 2*!PI*findgen(31)/30.+!pi/2. & xarr1=0.3*cos(tht) & yarr1=0.3*sin(tht)-1.3
    xarrl=[-1.0,-1.0,-1.01,-1.01,-1.0,-1.0,0.0]
    yarrl=[ 0.0, 1.0, 1.0, -1.0, -1.0, 0.0,0.0]
    xarrc=[0.0,0.0]
    yarrc=[0.4,-1.0]
    xarrr=[0.0,1.0,1.0,1.01,1.01,1.0,1.0]
    yarrr=[0.0,0.0,-1.0,-1.0,1.0,1.0,0.0]
    xarr = [xarrl,xarrc,xarr1,xarrr] * psiz*aspx
    yarr = [yarrl,yarrc,yarr1,yarrr] * psiz*aspy
    pfil = 0
  end
  'NEPTUNE': begin
    xarrm = [0.0, 0.0, -0.6, 0.0, 0.6, 0.0, 0.0]
    yarrm = [-1.0,-0.5,-0.5,-0.5,-0.5,-0.5, 0.0]
    thtr = (!pi/2.)*findgen(10)/9.-(!pi/2.)
    xarrr = cos(thtr)
    yarrr = sin(thtr)+1.0
    xarrra = [0.2,0.0,-0.2,0.0] + xarrr[9]
    yarrra = [-0.2,0.0,-0.2,0.0] + yarrr[9]
    xarrmu = [0.0, 0.2, 0.0, -0.2, 0.0]
    yarrmu = [0.8, 0.6, 0.8,  0.6, 0.8]
    thtl = reverse((!pi/2.)*findgen(10)/9.+!pi)
    xarrl = cos(thtl)
    yarrl = sin(thtl)+1.0
    xarrla = [0.2,0.0,-0.2] + xarrl[9]
    yarrla = [-0.2,0.0,-0.2] + yarrl[9]
    xarr = [xarrm, xarrr, xarrra, reverse(xarrr), xarrmu, xarrl, xarrla] * psiz*aspx
    yarr = [yarrm, yarrr, yarrra, reverse(yarrr), yarrmu, yarrl, yarrla] * psiz*aspy
  end
  'PLUTO': begin
    tht = !PI*findgen(31)/30.-!pi/2. & xarrP=0.6*cos(tht) & yarrP=0.4*sin(tht)
    xarrL = [0.6,0.0,0.0]
    yarrL = [0.0,0.0,1.5]
    xarr = [xarrL, xarrP] * psiz*aspx
    yarr = [yarrL, yarrP+1.1] * psiz*aspy
    pfil = 0
  end
  'SQUARE': begin
    xarr = [1,-1,-1,1,1]*psiz*aspx
    yarr = [1,1,-1,-1,1]*psiz*aspy
  end
  'DIAMOND': begin
    xarr = [0,1.41,0,-1.41,0]*psiz*aspx
    yarr = [1.41,0,-1.41,0,1.41]*psiz*aspy
  end
  'SPOKE': begin
    xarr = [0.0,0.1,0.0,-0.1,0.0]*psiz*aspx
    yarr = (1+[1.0,0.0,-1.0,0.0,1.0])*psiz*aspy
  end
  'UP TRIANGLE': begin
    xarr = [-1,0,1,-1]*psiz*aspx
    yarr = [-1,1,-1,-1]*psiz*aspy
    if keyword_set(halftri) then yarr=yarr+psiz*aspy
  end
  'DOWN TRIANGLE': begin
    xarr = [-1, 0, 1, -1]*psiz*aspx
    yarr = [ 1,-1, 1, 1]*psiz*aspy
    if keyword_set(halftri) then yarr=yarr-psiz*aspy
  end
  'RIGHT TRIANGLE': begin
    xarr = [-1,1,-1,-1]*psiz*aspx
    yarr = [1,0,-1,1]*psiz*aspy
    if keyword_set(halftri) then xarr=xarr+psiz*aspx
  end
  'LEFT TRIANGLE': begin
    xarr = [1,1,-1,1]*psiz*aspx
    yarr = [1,-1,0,1]*psiz*aspy
    if keyword_set(halftri) then xarr=xarr-psiz*aspx
  end
  'DOWN ARROW': begin
    xarr = [0,0,.5,0,-.5]*psiz*aspx
    yarr = [0,-2,-1.4,-2,-1.4]*psiz*aspy
    pfil = 0
  end
  'UP ARROW': begin
    xarr = [0,0,.5,0,-.5]*psiz*aspx
    yarr = [0,2,1.4,2,1.4]*psiz*aspy
    pfil = 0
  end
  'LEFT ARROW': begin
    yarr = [0, 0, 0.5, 0, -.5]*psiz*aspy
    xarr = [0,-2,-1.4,-2,-1.4]*psiz*aspx
    pfil = 0
  end
  'RIGHT ARROW': begin
    yarr = [ 0, 0, 0.5, 0, -.5] * psiz * aspy
    xarr = [ 0, 2, 1.4, 2, 1.4] * psiz * aspx
    pfil = 0
  end
  'TILDE': begin
    narr=npt/2L & marr=(2L*narr+1L < 49 ) > 2
    xarr = (findgen(marr)-marr/2)/(marr/2)
    if narr ne tilde_ovr then xarr = (findgen(marr)-narr)/(narr-tilde_ovr)
    yarr = (-sin(xarr*!pi))
    if tilde_ang ne 0 then begin
      xarr2 = xarr*cos(tilde_ang*!pi/180.) + yarr*sin(tilde_ang*!pi/180.)
      yarr2 = -xarr*sin(tilde_ang*!pi/180.) + yarr*cos(tilde_ang*!pi/180.)
      xarr = xarr2 & yarr=yarr2
    endif
    xarr = xarr * psiz * aspx
    yarr = yarr * psiz * aspy
    if vv gt 100 then stop
  end
  else: begin
    message,CSYMB+': plotting symbol not supported',/informational
    print,'Available options are:'
    print,'CIrcle'
    print,'STar / STARFish / ASterisk'
    print,'MErcury/VEnus/EArth/MOon/MArs/JUpiter/SAturn/URanus/NEptune/PLuto'
    print,'SQuare / Diamond / Box / SPoke'
    print,'up TRIangle / DOWn TRIangle / RIGht TRIangle /LEFt TRIangle'
    print,'UP ARRow / DOWn ARRow / LEFt ARRow / RIGht ARRow'
    print,'UPPer LIMit / LOWer LIMit'
    print,'TILde <angle> <overshoot>'
    return
  end
endcase						;CSYMB}

if ang ne 0 then begin
  xarr2 = xarr*cos(ang) + yarr*sin(ang)
  yarr2 = -xarr*sin(ang) + yarr*cos(ang)
  xarr = xarr2 & yarr=yarr2
endif

usersym, xarr, yarr, FILL = pfil, THICK=pthk, _extra=e

return
end          
