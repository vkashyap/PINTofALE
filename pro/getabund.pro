function getabund,hint,elem=elem,source=source,norm=norm,metals=metals,$
	fipbias=fipbias,abunde=abunde,_extra=e
;+
;function	getabund
;		returns elemental abundances in an array
;
;syntax
;	abund=getabund(hint,elem=elem,source=source,norm=norm,$
;	metals=metals,fipbias=fipbias,abunde=abunde)
;
;parameters
;	hint	[INPUT] some hint as to what output is required.  examples:
;		'help' -- prints out usage
;		'angr', 'anders', 'anders & grevesse', '1989', etc. -- Anders & Grevesse (1989)
;		'meyer', '1985', etc. -- Meyer (1985)
;		'allen', 'cosmic', '1973' -- Allen (1973)
;		'ross', 'aller', 'ross & aller', '1976' etc. -- Ross & Aller (1976)
;		'grevesse', 'grevesse et al.', etc. -- Grevesse et al (1992)
;			modifications to Anders & Grevesse (1989)
;		'feld', 'feldman', '1992', etc. -- Feldman et al (1992)
;		'waljeski', '1994', etc. -- Waljeski et al. (1994)
;		'murphy', 'chromospheric', etc. -- Murphy (1985) chromospheric
;		'grevesse & anders', '1991' -- Grevesse & Anders (1991)
;		'grsa', 'grevesse & sauval', '1998' -- Grevesse & Sauval (1998) "standard
;			abundances"
;		'young', 'photo', '1997' -- photospheric abundances from Young et
;			al. (1997)
;		'fludra', 'fludra & schmelz', 'hybrid', '1999' -- Fludra & Schmelz (1999) "hybrid coronal"
;		'asplund', '2005' -- Asplund, Grevesse, & Sauval (2005) 3D hydro
;			(NOTE: even though it has a different reference, is identical to aspl)
;		'schmelz', 'schmelz et al', 'coronal', '2012' -- Schmelz et al (2012) "recommended coronal"
;		'aspl', 'asplund 2009', 'asp ... & scott', '2009', etc. -- Asplund, Grevesse, Sauval, & Scott (2009)
;			model of solar atmosphere to constrain convection
;		'aneb', 'anders & ebihara', 'ebihara', '1982', etc. - Anders & Ebihara (1982)
;		'wilm', 'wilms, allen, mccray', '2000', etc. -- Wilms, Allen, & McCray (2000)
;		'lodd', 'lodders', '2003', etc. -- Lodders (2003)
;		'path/file', etc. -- if string contains a "/" in it, then
;			read from disk file which contains one entry per
;			line, with atomic symbol first, and abundance next
;
;keywords
;	elem	[INPUT] element -- if given, return abundances for only
;		the specified element
;		* if string, converts to atomic number using SYMB2ZION
;		* if byte, integer, or float, assumed to be atomic number
;		* may be an array
;	source	[INPUT] another way to say "hint"
;		1: Anders & Grevesse (1989)
;		2: Meyer (1985)
;		3: Allen (1973)
;		4: Ross & Aller (1976)
;		5: Grevesse et al. (1992)
;		6: Feldman et al. (1992)
;		7: Waljeski et al. (1994)
;		8: Murphy (1985)
;		9: Grevesse & Anders (1991)
;		10: Grevesse & Sauval (1998)
;		11: Young et al. (1997)
;		12: Fludra & Schmelz (1999)
;		13: Asplund, Grevesse, & Sauval (2005)
;		14: Schmelz et al. (2012)
;		15: Asplund, Grevesse, Sauval, & Scott (2009)
;		16: Anders & Ebihara (1982)
;		17: Wilms, Allen, & McCray (2000)
;		18: Lodders (2003)
;		* overrides HINT iff HINT is not understood
;	norm	[INPUT] normalization
;		* use this to set output format
;		  NORM=1: output will be n(Z)/n(H) (this is the default)
;		  NORM>1: output will be NORM(0)+log10(n(Z)/n(H))
;		  	(e.g., NORM=12 is the standard log form)
;		  0<NORM<1: output will be NORM(0)*n(Z)/n(H)
;		  	(this can be useful when ELEM is set)
;		  -1<=NORM<0: output will be n(Z)/DEFAULT_n(Z)
;		  	(if you want to compare different abundance tables)
;		  NORM<-1: output will be log10(n(Z)/DEFAULT_n(Z))
;		  	(I have no idea why anyone would use this)
;	metals	[INPUT] metallicity.  if set, modifies the abundances of
;		all elements past He by the appropriate amount
;		* if not set, makes no changes
;		* assumed to be in log10 form if
;		  -- abs(NORM) > 1, or
;		  -- < 0
;		* assumed to be multiplicative factor otherwise
;		* NOTE: for example, if norm=12 then metals=0 means
;		  Solar abundances, while if norm=1 then metals=0
;		  means complete metal depletion!
;	fipbias	[INPUT] if set to a non-zero scalar, multiplies the
;		abundances of low-FIP elements by this factor.  be
;		careful not to use it for abundance lists that already
;		include the FIP effect.
;		* if not set, or set to 0, does nothing.
;	abunde	[OUTPUT] errors on abundances
;	
;	_extra	[INPUT] junk -- here only to prevent program from crashing
;
;subroutines
;	RDABUND
;	SYZE
;	SYMB2ZION
;
;history
;	vinay kashyap (Dec96)
;	avoid analyzing HINT if SOURCE is set (VK; Jan97)
;	added 'angr' and tightened 'Anders&Grevesse' v/s 'Grevesse' (VK; Jul97)
;	added option of reading from disk file (VK; Jun98)
;	added Murphy (1985) data, and allowed ELEM to be array (VK; Jul98)
;	corrected bug associated with non-array ELEMs (VK; Sep98)
;	added keyword METALS, deleted keywords ALTER, FORM, removed call to
;	  SETABUND, added call to RDABUND (VK; Dec98)
;	added Grevesse & Anders (1991), Grevesse & Sauval (1998), Young et
;	  al (1998) (VK; Apr99)
;	added keyword FIPBIAS (VK; Jul01)
;	added Fludra & Schmelz (1999); now HINT overrides SOURCE and program
;	  warns if there is a conflict (VK; Jun02)
;	corrected Meyer abundances for O and Na (VK; Apr04)
;	corrected Fludra & Schmelz abundances (VK; Jun04)
;	added Asplund, Grevesse, & Sauval abundances; added keyword ABUNDE
;	  (VK; Apr05)
;	B was off by 2 orders of magnitude (thx Marc Audard; VK Jul06)
;	streamlined input HINTs; added Schmelz et al. (2012); updated
;	  bracket notation to IDL5+ (VK Jul12)
;	added 'aspl', 'aneb', 'wilm', and 'lodd' from the elemental abundance
;	  table used in XSPEC plasma code models, as given in ahelp get_xsabund;
;	  changed convection model from Asplund et al (2005) to aspl; streamlined
;	  selection keywords (VK Feb13)
;	changed the behavior of default to accept !ABREF if defined (VK May28)
;-

;	check input
np=n_elements(hint)
if np eq 0 then begin
  todo='help'
  defsysv,'!ABREF',exists=ivar
  if ivar eq 1 then begin
    message,'getting abundances for !ABREF='+!ABREF,/informational
    jnk=execute("todo=strlowcase(!ABREF)")
  endif
endif else begin
  if np gt 1 then todo='help' else todo=strlowcase(string(hint))
endelse
;if np eq 0 or np gt 1 then todo='help' else todo=strlowcase(string(hint))
if not keyword_set(norm) then norm=1
;
if n_elements(metals) gt 0 then begin
  met=metals[0]
  if abs(norm) gt 1 or met lt 0 then lmet=1 else lmet=0
endif

codemax=19				;so far this many have been coded.
					;CODEMAX would be file input

code_s=0 & code_h=0

if keyword_set(source) then begin	;(SOURCE has been set
  szs=size(source) & nszs=n_elements(szs)
  code_s=0
  if szs[nszs-2] lt 7 then code_s=fix(source[0])
  if szs[nszs-2] eq 7 then begin
    code_s=codemax	;must be a filename
    hint=strtrim(source[0],2)	;the name of the file
  endif
endif					;SOURCE)

;	now analyze the HINT
  code_h=0					;default
  			;ANDERS & GREVESSE?
  i0=strpos(todo,'ander') & i1=strpos(todo,'gre') & i2=strpos(todo,'1989') & i3=strpos(todo,'angr')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 or i3 ge 0 then code_h=1
  			;MEYER?
  i0=strpos(todo,'mey') & i1=strpos(todo,'mayer') & i2=strpos(todo,'1985')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=2
  			;ALLEN?
  i0=strpos(todo,'allen') & i1=strpos(todo,'allan') & i2=strpos(todo,'cos') & i3=strpos(todo,'1973')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 or i3 ge 0 then code_h=3
  			;ROSS & ALLER?
  i0=strpos(todo,'ros') & i1=strpos(todo,'aller') & i2=strpos(todo,'1976')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=4
  			;GREVESSE, NOELS, SAUVAL?
  i0=strpos(todo,'gre') & i1=strpos(todo,'et. al') & i2=strpos(todo,'et al')
  i3=strpos(todo,'noe') & i4=strpos(todo,'sau')
  if i0 ge 0 and (i1 ge 0 or i2 ge 0) then code_h=5
  i5=strpos(todo,'ander')
  if i0 ge 0 and i5 ge 0 then code_h=1		;<-- Anders & Grevesse
  if i0 ge 0 and i5 lt 0 then code_h=5		;<-- no Anders w. Grevesse!
  if i5 gt i0 and i0 ge 0 then code_h=9		;<-- Grevesse & Anders
  if i4 ge 0 then begin
    if i0 ge 0 and i3 lt 0 then code_h=10	;<-- G&S, no N
    if i0 ge 0 and i3 ge 0 then code_h=5	;<-- G, N *and* S!
  endif
  i6=strpos(todo,'1992') & i7=strpos(todo,'feld')
  if i6 ge 0 and i7 lt 0 then code_h=5		;<-- not Feldman et al 1992
  			;FELDMAN, MANDELBAUM, SEELY, DOSCHEK, GURSKY?
  i0=strpos(todo,'fel') & i3=strpos(todo,'mand')
  i4=strpos(todo,'see') & i5=strpos(todo,'dos') & i6=strpos(todo,'gur')
  if i0 ge 0 and (i1 ge 0 or i2 ge 0) then code_h=6
  if i0 ge 0 or i3 ge 0 or i4 ge 0 or i5 ge 0 or i6 ge 0 then code_h=6
  i7=strpos(todo,'grev') & i8=strpos(todo,'1992')
  if i7 lt 0 and i8 ge 0 then code_h=6		;<-- Feldman et al, not Grevesse et al
  			;WALJESKI, MOSES, DERE, SABA, STRONG, WEBB, ZARRO?
  i0=strpos(todo,'wal') & i3=strpos(todo,'mos') & i4=strpos(todo,'dere')
  i5=strpos(todo,'sab') & i6=strpos(todo,'stro') & i7=strpos(todo,'we')
  i8=strpos(todo,'zar') & i9=strpos(todo,'1994')
  if i0 ge 0 and (i1 ge 0 or i2 ge 0) then code_h=7
  if i0 ge 0 or i3 ge 0 or i4 ge 0 or i5 ge 0 or i6 ge 0 $
	or i7 ge 0 or i8 ge 0 or i9 ge 0 then code_h=7
			;MURPHY (1985)?
  i0=strpos(todo,'mur') & i1=strpos(todo,'chromo') & i2=strpos(todo,'1985')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=8
			;GREVESSE & ANDERS (1991)?
  i0=strpos(todo,'gre') & i1=strpos(todo,'ander') & i2=strpos(todo,'1991')
  if (i0 ge 0 and i1 ge 0 and i1 gt i0) or (i2 ge 0) then code_h=9
			;GREVESSE & SAUVAL (1998)?
  i0=strpos(todo,'gre') & i1=strpos(todo,'sau') & i2=strpos(todo,'1998') & i3=strpos(todo,'grsa')
  if (code_h ne 5 and i0 ge 0 and i1 ge 0) or i2 ge 0 or i3 ge 0 then code_h=10
			;YOUNG ET AL (1997)/PHOTOSPHERIC?
  i0=strpos(todo,'you') & i1=strpos(todo,'phot') & i2=strpos(todo,'1997')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=11
			;FLUDRA & SCHMELZ (1999)
  i0=strpos(todo,'flu') & i1=strpos(todo,'hyb') & i2=strpos(todo,'1999')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=12
			;ASPLUND, GREVESSE, & SAUVAL (2005)
  i0=strpos(todo,'asplund') & i1=strpos(todo,'scott') & i2=strpos(todo,'2005')
  i3=strpos(todo,'et al') & i4=strpos(todo,'& sauval') & i5=strpos(todo,'and sauval')
  if (i0 ge 0 and i3 ge 0) or (i0 ge 0 and i4 ge 0) or (i0 ge 0 and i5 ge 0) or $
  	(i0 ge 0 and i1 lt 0) or i2 ge 0 then code_h=13
  			;SCHMELZ, REAMES, VON STEIGER, BASU (2012)?
  i0=strpos(todo,'sch') & i1=strpos(todo,'rea') & i2=strpos(todo,'stei')
  i3=strpos(todo,'von') & i4=strpos(todo,'basu') & i5=strpos(todo,'2012')
  i6=strpos(todo,'flu') & i7=strpos(todo,'coro')
  if (i0 ge 0 and i6 lt 0) or i1 ge 0 or i2 ge 0 or i3 ge 0 or $
  	i4 ge 0 or i5 ge 0 or i7 ge 0 then code_h=14
  			;Asplund, Grevesse, Sauval, & Scott (2009)?
  i0=strpos(todo,'aspl') & i1=strpos(todo,'scott') & i2=strpos(todo,'2009')
  i3=strpos(todo,'convection')
  if (i0 ge 0) or (i0 ge 0 and i1 ge 0) or (i2 ge 0) or (i3 ge 0) then code_h=15
			;Anders & Ebihara (1982)?
  i0=strpos(todo,'aneb') & i1=strpos(todo,'ebihara') & i2=strpos(todo,'1982')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=16
			;Wilms, Allen, & McCray (2000)?
  i0=strpos(todo,'wilm') & i1=strpos(todo,'mccray') & i2=strpos(todo,'2000')
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=17
			;Lodders (2003)?
  i0=strpos(todo,'lodd') & i1=strpos(todo,'2003')
  if i0 ge 0 or i1 ge 0 then code_h=18
			;disk file?
  i0=strpos(todo,'/') & i1=strpos(todo,'\') & i2=strpos(todo,path_sep())
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=codemax
  ;
  			;AIUTO!
  i0=strpos(todo,'help') & if i0 ge 0 then code_h=0

;	the above block used to be inside the first if block below,
;	until the possibility of a SOURCE/HINT conflict was addressed
;if not keyword_set(source) then begin			;{analyze the hint
;endif else begin				;}{hint, shint..
;  szs=size(source) & nszs=n_elements(szs)
;  code=0
;  if szs[nszs-2] lt 7 then code=fix(source[0])
;  if szs[nszs-2] eq 7 then begin
;    code=codemax	;must be a filename
;    hint=strtrim(source[0],2)	;the name of the file
;  endif
;endelse						;loud and clear}

;	catch some trivial errors

;	what if there is a conflict between SOURCE and HINT?
if not keyword_set(source) then code_s = code_h else $
 if np ne 1 then code_h = code_s
if code_s ne code_h then begin
  message,'keyword value incompatible with parameter value',/info
  if code_h ne 0 then begin
    message,'choosing to go with HINT="'+todo+'"',/info
    code=code_h
  endif else begin
    message,'choosing to go with the SOURCE='+strtrim(code_s,2),/info
    code=code_s
  endelse
endif else code=code_h
;
if code lt 0 then code=0 & if code gt codemax then code=0

;	usage
if code eq 0 then begin
  print,'Usage: abund=getabund(hint,elem=elem,source=source,norm=norm,$'
  print,'       metals=metals,fipbias=fipbias,abunde=abunde)'
  print,'  returns elemental abundances'
  print,''
  print,'  SOURCE=1: Anders & Grevesse (1989) **[angr]** <<<DEFAULT>>>'
  print,'  SOURCE=2: Meyer (1985)'
  print,'  SOURCE=3: Allen (1973) **[cosmic]**'
  print,'  SOURCE=4: Ross & Aller (1976)'
  print,'  SOURCE=5: Grevesse et al. (1992)'
  print,'  SOURCE=6: Feldman et al. (1992) **[feld]**'
  print,'  SOURCE=7: Waljeski et al. (1994)'
  print,'  SOURCE=8: Murphy (1985) **(chromospheric)**'
  print,'  SOURCE=9: Grevesse & Anders (1991)'
  print,'  SOURCE=10: Grevesse & Sauval (1998) **[grsa]**'
  print,'  SOURCE=11: Young et al. (1997) **(photospheric)**'
  print,'  SOURCE=12: Fludra & Schmelz (1999) **(hybrid)**'
  print,'  SOURCE=13: Asplund, Grevesse, & Sauval (2005)'
  print,'  SOURCE=14: Schmelz, Reames, von Steiger, & Basu (2012) **(coronal)**'
  print,'  SOURCE=15: Asplund, Grevesse, Sauval & Scott (2009) **[aspl]**'
  print,'  SOURCE=16: Anders & Ebihara (1982) **[aneb]**'
  print,'  SOURCE=17: Wilms, Allen, & McCray (2000) **[wilm]**'
  print,'  SOURCE=18: Lodders (2003) **[lodd]**'
  print,''
  print,'  0<NORM<=1:	ABUND(Z)=NORM*n(Z)/n(H) {DEFAULT: NORM=1}'
  print,'  NORM>1:	ABUND(Z)=log10(n(Z)/n(H))+NORM'
  print,'  -1<=NORM<0:	ABUND(Z)=ABUND(Z)/DEFAULT_ABUND(Z)'
  print,'  NORM<-1:	ABUND(Z)=log10(ABUND(Z)/DEFAULT_ABUND(Z))'
  print,''
endif

;	default abundances
defabu=[1., 0.0977, 1.45e-11, 1.15e-11, 3.98e-10, 3.63e-4, 1.12e-4,$
	8.51e-4, 3.63e-8, 1.23e-4, 2.14e-6, 3.8e-5, 2.95e-6, 3.55e-5,$
	2.82e-7, 1.62e-5, 3.16e-7, 3.63e-6, 1.32e-7, 2.29e-6, 1.26e-9,$
	9.77e-8, 1e-8, 4.68e-7, 2.45e-7, 4.68e-5, 8.32e-8, 1.78e-6,$
	1.62e-8, 3.98e-8] & nelem=n_elements(defabu)
;Anders, E. \& Grevesse, N.\ 1989, {\it Geochimica et Cosmochimica Acta},
;53, 197
abunde=0.*defabu	;errors on abundances

;	database of abundances
case code of

  2: begin	;Meyer, J.P.\ 1985, ApJS, 57, 173
		;NOTE: CHIANTI lists abund(7,10)=[8.39,6.44]
    abund=fltarr(nelem)
    abund[0]=12.00 & abund[1]=10.99 & abund[5]=8.37 & abund[6]=7.59
    ;abund[7]=8.33 & abund[9]=7.55 & abund[11]=7.57 & abund[12]=6.44
    abund[7]=8.39 & abund[9]=7.55 & abund[10]=6.44 & abund[11]=7.57 & abund[12]=6.44
    abund[13]=7.59 & abund[15]=6.94 & abund[17]=6.33 & abund[19]=6.47
    abund[25]=7.59 & abund[27]=6.33
    oo=where(abund ne 0) & abund[oo]=10.^(abund[oo]-abund[0])
    oo=where(abund eq 0) & abund[oo]=defabu[oo]
  end

  3: begin	;Allen, C.W.\ 1973, {it Astrophysical Quantities}, 3$^{rd}$
		;Ed., Athlone Press, London
    frac=[1., 0.87, 0.35, 0.89, 2.51, 0.91, 0.81, 0.78, 1.10, 0.68, 0.83,$
	0.69, 0.83, 0.93, 1.17, 0.98, 1.26, 1.74, 0.68, 0.87, 1.32, 1.38,$
	2.51, 1.51, 1.02, 0.85, 1.51, 1.12, 1.95, 0.4]
    abund=frac*defabu
  end

  4: begin	;Ross, J.E., \& Aller, L.H.\ 1976, {\it Science}, 191, 1223
    frac=[1., 0.65, 0.69, 1., 0.32, 1.15, 0.78, 0.81, 1., 0.3, 0.89, 1.02,$
	  1.12, 1.26, 1.12, 0.98, 1., 0.28, 1.1, 0.98, 0.87, 1.15, 1.05,$
	  1.1, 1.07, 0.68, 0.95, 1.07, 0.71, 0.71]
    abund=frac*defabu
  end

  5: begin	;Grevesse, N., Noels, A., \& Sauval, A.J.\ 1992, in {\it
		;Coronal Streamers, Coronal Loops and Coronal and Solar
		;Wind Composition}, Proc.\ 1st SOHO Workshop, ESA SP-348, ESA
		;Publ.\ Div., ESTEC, Noordwijk, p.305
    abund=defabu
    abund[1]=0.095 & abund[5]=3.55e-4 & abund[6]=9.33e-5 & abund[7]=7.41e-4
    abund[8]=3.63e-8 & abund[9]=1.2e-4 & abund[16]=3.16e-7 & abund[17]=3.31e-6
    abund[20]=1.58e-9 & abund[21]=1.1e-7 & abund[25]=3.24e-5
  end

  6: begin	;Feldman, U., Mandelbaum, P., Seely, J.L., Doschek, G.A., \&
		;Gursky H.\ 1992, ApJS, 81, 387
    abund=[-3.41, -4.00, -3.11, -7.4, -3.92, -5.07, -3.85, -4.96,$
	-3.90, -6.48, -4.73, -6.4, -5.42, -7.05, -5.07, -8.78, -6.87, -7.6,$
	-6.15, -6.60, -3.90, -6.9, -5.16, -7.5, -7.8]
    abund=10.^(abund) & abund=[1.,10.^(-1.1),defabu[2:4],abund]
  end

  7: begin	;Waljeski, K., Moses, D., Dere, K.P., Saba, J.L., Strong, K.T.,
		;Webb, D.F., \& Zarro, D.M.\ 1994, ApJ, 429, 909
    abund=fltarr(nelem)
    abund[0]=12. & abund[1]=10.90 & abund[5]=9.28 & abund[6]=8.50
    abund[7]=9.30 & abund[9]=8.45 & abund[11]=8.48 & abund[12]=7.35
    abund[13]=8.50 & abund[15]=7.84 & abund[17]=7.24 & abund[19]=7.38
    abund[25]=8.50 & abund[27]=7.24
    oo=where(abund ne 0) & abund[oo]=10.^(abund[oo]-abund[0])
    oo=where(abund eq 0) & abund[oo]=defabu[oo]
  end

  8: begin	;Murphy, R.J.\ 1985, Ph.D.\ Thesis, Univ. of Maryland
    abund=fltarr(nelem)
    abund[1-1]=1. & abund[2-1]=7.57e-2 & abund[6-1]=1.14e-4
    abund[7-1]=4.82e-5 & abund[8-1]=1.72e-4 & abund[10-1]=7.84e-5
    abund[12-1]=2.63e-5 & abund[14-1]=3.92e-5 & abund[16-1]=1.88e-5
    abund[20-1]=6.27e-6 & abund[26-1]=2.47e-5
    oo=where(abund eq 0) & abund[oo]=defabu[oo]
  end

  9: begin	;Grevesse, N. \& Anders, E.\ 1991, in Solar Interior and
		;Atmosphere, ed. A.N.Cox, W.C.Livingston, & M.S.Matthews,
		;(Tucson: Univ. Arizona Press), 1227.
    abund=fltarr(nelem)
    abund[1-1]=12. & abund[2-1]=10.90
    abund[6-1:30-1]=[8.60,8.00,8.93,4.56,8.09,6.33,7.58,6.47,7.55,5.45,$
	7.21,5.5,6.56,5.12,6.36,3.10,4.99,4.00,5.67,5.39,7.67,4.92,$
	6.25,4.21,4.60]
    oo=where(abund eq 0,moo)
    abund=10.^(abund-12.) & if moo gt 0 then abund[oo]=defabu[oo]
  end

  10: begin	;Grevesse, N. \& Sauval, A.J.\ 1998, in Space Science Reviews,
		;85, 161-174
		;CHIANTI says they are same as Grevesse, Noels, & Sauval 1996.
    abund=[12.00,10.93,1.10,1.40,2.55,8.52,7.92,8.83,4.56,8.08,6.33,7.58,$
	6.47,7.55,5.45,7.33,5.50,6.40,5.12,6.36,3.17,5.02,4.00,5.67,5.39,$
	7.50,4.92,6.25,4.21,4.60]
    abund=10.^(abund-12.)
  end

  11: begin	;Young, P.R., Mason, H.E., Keenan, F.P., \& Widing K.G.
		;1997, A&A, 323, 243
		;also see CHIANTI's photospheric_may97.abund
    abund=fltarr(nelem)
    abund[1-1]=12. & abund[2-1]=10.99
    abund[6-1:30-1]=[8.55,7.97,8.87,4.56,8.08,6.33,7.58,6.47,7.55,5.45,$
	7.21,5.5,6.47,5.12,6.36,3.20,5.02,4.00,5.67,5.39,7.51,4.92,6.25,$
	4.21,4.60]
    oo=where(abund eq 0,moo)
    abund=10.^(abund-12.) & if moo gt 0 then abund[oo]=defabu[oo]
  end

  12: begin	;Fludra, A. & Schmelz, J. T., 1999, A&A, 348, 286
		;CHIANTI's fludra_schmelz_hybrid.abund
    abund=fltarr(nelem)
    ;{VK/Jun04: CHIANTI had this incorrectly, so redid, direct from source
    ;iz=[1,    2,    6,   7,   8,   10,  11,  12,  13,  14,  15,  16,$
    ;	17,  18,  19,  20,  22,  24,  26,  28,  30]
    ;xz=[12.00,10.80,8.41,7.81,8.74,7.95,6.63,7.90,6.80,7.87,5.44,7.32,$
    ;	5.08,6.36,5.46,6.66,5.25,6.00,7.83,6.56,4.09]
    ;abund[iz]=xz
    ;oo=where(abund eq 0,moo)
    ;abund=10.^(abund-12.) & if moo gt 0 then abund[oo]=defabu[oo]
    ;}
    ;
    iz=[19,       11,       13,       20,       24,       22,$
	12,       28,       26,       14,       6,        17,$
	8,        7,        18,       10,       2]
    xz=[2.79e-07, 4.25e-06, 6.30e-06, 4.54e-06, 1.00e-06, 1.78e-07,$
	7.94e-05, 3.65e-06, 6.74e-05, 7.42e-05, 0.000257, 1.20e-07,$
	0.000548, 6.47e-05, 2.31e-06, 8.83e-05, 0.0631]
    abund[iz-1]=xz
    oo=where(abund eq 0,moo) & if moo gt 0 then abund[oo]=defabu[oo]
  end

  13: begin	;Asplund, M., Grevesse, N., & Sauval, A.J., 2005,
		;Invited review presented at "Cosmic abundances as
		;records of stellar evolution and nucleosynthesis",
		;F.N. Bash & T.G Barnes(editors). ASP conf. series,
		;in press.
		;http://arxiv.org/abs/astro-ph/0410214
    xz=[12.00, 10.93, 1.05, 1.38, 2.70, 8.39, 7.78, 8.66, 4.56, 7.84,$
	 6.17,  7.53, 6.37, 7.51, 5.36, 7.14, 5.50, 6.18, 5.08, 6.31,$
	 3.05,  4.90, 4.00, 5.64, 5.39, 7.45, 4.92, 6.23, 4.21, 4.60]
    xze=[0.00, 0.01, 0.10, 0.09, 0.20, 0.05, 0.06, 0.05, 4.43, 0.06,$
	 0.04, 0.09, 0.06, 0.04, 0.04, 0.05, 0.30, 0.08, 0.07, 0.04,$
	 0.08, 0.06, 0.02, 0.10, 0.03, 0.05, 0.08, 0.04, 0.04, 0.03]
    abund=10.^(xz-12.)
    abunde=abund*alog(10.)*xze
  end

  14: begin	;Schmelz, J.T., Reames, D.V., von Steiger, R., & Basu, S.
  		;2012, ApJ, in press(?)
		;from Table 1
    abund=defabu
    iz=[2,    6,   7,   8,   10,  11,   12,  13,   14,  15,    16,  17,    18,   19,    20,   22,    24,    26,  28,   30]
    xz=[6080.,22.6,51.7,41.1,8.01,0.417,7.41,0.616,7.24,0.0309,1.69,0.0127,0.226,0.0275,0.436,0.0182,0.0955,7.08,0.355,0.00826]/1e5
    abund[iz-1]=xz
  end

  15: begin	;Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)
    abund=[1., 8.51e-02, 1.12e-11, 2.40e-11, 5.01e-10, 2.69e-04,$
    	6.76e-05, 4.90e-04, 3.63e-08, 8.51e-05, 1.74e-06, 3.98e-05,$
    	2.82e-06, 3.24e-05, 2.57e-07, 1.32e-05, 3.16e-07, 2.51e-06,$
    	1.07e-07, 2.19e-06, 1.41e-09, 8.91e-08, 8.51e-09, 4.37e-07,$
    	2.69e-07, 3.16e-05, 9.77e-08, 1.66e-06, 1.55e-08, 3.63e-08]
  end

  16: begin	;Anders E. & Ebihara (1982, Geochimica et Cosmochimica Acta 46, 2363)
    abund=[1., 8.01e-2, 2.19e-9, 2.87e-11, 8.82e-10, 4.45e-4,$
    	9.12e-5, 7.39e-4, 3.10e-8, 1.38e-4, 2.10e-6, 3.95e-5,$
    	3.12e-6, 3.68e-5, 3.82e-7, 1.89e-5, 1.93e-7, 3.82e-6,$
    	1.39e-7, 2.25e-6, 1.24e-9, 8.82e-8, 1.08e-8, 4.93e-7,$
    	3.50e-7, 3.31e-5, 8.27e-8, 1.81e-6, 1.89e-8, 4.63e-8]
  end

  17: begin	;Wilms, Allen & McCray (2000, ApJ 542, 914)
  		;unknown values filled in with default
    abund=[1, 9.77e-2, -1, -1, -1, 2.40e-4,$
    	7.59e-5, 4.90e-4, -1, 8.71e-5, 1.45e-6, 2.51e-5,$
    	2.14e-6, 1.86e-5, 2.63e-7, 1.23e-5, 1.32e-7, 2.57e-6,$
    	-1, 1.58e-6, -1, 6.46e-8, -1, 3.24e-7,$
    	2.19e-7, 2.69e-5, 8.32e-8, 1.12e-6, -1, -1]
    o0=where(abund lt 0,mo0) & if mo0 gt 0 then abund[o0]=defabu[o0]
  end

  18: begin	;Lodders, K (2003, ApJ 591, 1220)
    abund=[1., 7.92e-2, 1.90e-9, 2.57e-11, 6.03e-10, 2.45e-4,$
    	6.76e-5, 4.90e-4, 2.88e-8, 7.41e-5, 1.99e-6, 3.55e-5,$
    	2.88e-6, 3.47e-5, 2.88e-7, 1.55e-5, 1.82e-7, 3.55e-6,$
    	1.29e-7, 2.19e-6, 1.17e-9, 8.32e-8, 1.00e-8, 4.47e-7,$
    	3.16e-7, 2.95e-5, 8.13e-8, 1.66e-6, 1.82e-8, 4.27e-8]
  end

  codemax: begin	;read from disk file
    abund=rdabund(hint,norm=1, _extra=e)
  end

  else: abund=defabu

endcase

;	size check
nab=n_elements(abund)
if nab lt nelem then abund=[abund,defabu[nab:*]]

;	include metallicity
if n_elements(metals) gt 0 then begin
  if keyword_set(lmet) then met=10.^(met)
  abund[2:*]=met*abund[2:*]
endif

;	include FIP bias
if keyword_set(fipbias) then begin
  fip=abs(float(fipbias[0]))
  ;	the following are elements with FIP < 10 eV
  zlofip=[3,4,5,11,12,13,14,19,20,21,22,23,23,25,26,27,28,29,30]-1L
  if fip gt 0 then abund[zlofip]=abund[zlofip]*fip else message,$
	'FIPBIAS of '+strtrim(fip,2)+' not understood; ignoring',/info
endif

;	extract element(s)
sze=size(elem) & nsze=n_elements(sze) & icode=sze[nsze-2] & nn=sze[nsze-1]
case icode of
  1: atno=((elem > 1) < nelem)
  2: atno=((elem > 1) < nelem)
  3: atno=((elem > 1) < nelem)
  4: atno=((elem > 1) < nelem)
  7: begin
     atno=intarr(nn)
     for i=0,nn-1 do begin
       symb2zion,elem[i],x,y & atno[i]=x
     endfor
  end
  else: atno=lindgen(nelem)+1
endcase
;
if icode gt 0 then atno=atno[uniq(atno,sort(atno))]
abund=[abund[atno-1]]

;	renorm
if norm[0] gt 0 and norm[0] le 1 then abund=norm[0]*abund
if norm[0] gt 1 then abund=alog10(abund)+norm[0]
if norm[0] lt 0 and norm[0] ge -1 then abund=abund/defabu
if norm[0] lt -1 then abund=alog10(abund/defabu)

return,abund
end
