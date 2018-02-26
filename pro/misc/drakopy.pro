pro drakopy,milieu,help=help,undo=undo,magnify=magnify,verbose=verbose,$
	_extra=e
;+
;procedure	drakopy
;	set up the environment to make plots per individual preferences
;
;syntax
;	drakopy,milieu,/HELP,/UNDO,MAGNIFY=MAGNIFY,VERBOSE=VERBOSE
;
;parameters
;	milieu	[INPUT] string specifying what sort of environment to set to.
;		accepted values are:
;		-- 'jeremy'
;		-- 'olivia'
;		-- 'jjd' (jeremy old preferences)
;		-- 'drake' (jeremy older preferences)
;		-- 'vinay' (different depending device at call time)
;		-- 'poster'
;		-- 'original'
;
;keywords
;	help	[INPUT] if set, prints out a somewhat verbose help
;	undo	[INPUT] if set, restores original state
;	verbose	[INPUT] controls chatter level
;	magnify	[INPUT] scale all the variables by this factor
;	_extra	[JUNK] here only to prevent crashing the program
;
;commons
;	drakopy	{old_P, old_X, old_Y, old_Z}
;
;description
;	sets the !P, !X, !Y, and !Z variables to special and specific values.
;	stores the old values in a common block for later access.
;
;example
;	set_plot,'ps' & device,file='/tmp/test_drakopy.ps'
;	  plot,findgen(10),title='MILIEU=Original'
;	  drakopy,'jeremy' & plot,findgen(10),title='MILIEU=Jeremy'
;	  drakopy,'olivia' & plot,findgen(10),title='MILIEU=Olivia'
;	  drakopy,'original' & plot,findgen(10),title='MILIEU=Original'
;	  drakopy,/undo & plot,findgen(10),title='MILIEU=previous (Olivia)'
;	device,/close & set_plot,'x'
;
;history
;	vinay kashyap (SepMM)
;	added option "drake" (VK; Aug01)
;	added option "jjd"; added keyword MAGNIFY (VK; Feb02)
;	added option "vinay" (VK; May03)
;	added option "poster" (VK; Jul03)
;	told the code that "poster" is a valid option (VK; Sep03)
;	now prints out help and exits if no input is given (VK; Jul13)
;	changed "jeremy" defaults (JD; Nov15)
;-

;	usage
if n_params() eq 0 or keyword_set(help) then begin
  print,'Usage: drakopy,milieu,/HELP,/UNDO,MAGNIFY=MAGNIFY,VERBOSE=VERBOSE'
  print,"  set up the environment to make good-lookin' hardcopy plots"
  print,"  acceptable values for MILIEU:"
  print,"  'jeremy','olivia','drake','jjd','poster','vinay','original'"
  message,'Returning without making any changes',/info
  return
endif

;	populate the defaults
common drakopy,old_P,old_X,old_Y,old_Z
if n_tags(old_P) eq 0 then old_P=!P
if n_tags(old_X) eq 0 then old_X=!X
if n_tags(old_Y) eq 0 then old_Y=!Y
if n_tags(old_Z) eq 0 then old_Z=!Z
;
vv=0 & if keyword_set(verbose) then vv=fix(verbose[0]) > 1
;
zm=1.0 & if keyword_set(magnify) then zm=float(magnify[0])
if zm lt 0 then begin
  if vv gt 0 then message,$
	'MAGNIFY < 0 makes no sense; assuming abs(MAGNIFY)',/info
  zm=abs(zm)
endif

;	undo previous change
if keyword_set(undo) then begin
  if vv gt 0 then message,'Reverting to old environment',/info
  !P=old_P & !X=old_X & !Y=old_Y & !Z=old_Z
  old_P=!P & old_X=!X & old_Y=!Y & old_Z=!Z
  return		;All done, quit.
endif

;	save old values
old_P=!P & old_X=!X & old_Y=!Y & old_Z=!Z

;	set new values
nm=n_elements(milieu) & setto='splurg'
if nm gt 0 then begin
  mm=strupcase(milieu[0])
  if strpos(mm,'JE',0) ge 0 then setto='jeremy' else $
   if strpos(mm,'OL',0) ge 0 then setto='olivia' else $
    if strpos(mm,'OR',0) ge 0 then setto='original' else $
     if strpos(mm,'D',0) ge 0 then setto='drake' else $
      if strpos(mm,'JJ',0) ge 0 then setto='drake' else $
       if strpos(mm,'PO',0) ge 0 then setto='poster' else $
	if strpos(mm,'VI',0) ge 0 then setto='vinay'
endif
if vv gt 0 then message,'Setting milieu to: '+setto,/info

case setto of
  'original': begin
    pcharsize=1. & pcharthick=1. & psymsize=1. & pthick=1.
    xthick=1. & xcharsize=1.
    ythick=1. & ycharsize=1.
    zthick=1. & zcharsize=1.
  end
  'olivia': begin		;== to avoid 'spidery' plots
    pcharsize=2. & pcharthick=2. & psymsize=2. & pthick=2.
    xthick=2. & xcharsize=2.
    ythick=2. & ycharsize=2.
    zthick=2. & zcharsize=2.
  end
  'drake': begin		;== older 'jeremy'
    pcharsize=1.7 & pcharthick=10. & psymsize=1.7 & pthick=10.
    xthick=10. & xcharsize=1.7
    ythick=10. & ycharsize=1.7
    zthick=10. & zcharsize=1.7
  end
  'jjd': begin			;== old 'jeremy'
    pcharsize=1.3 & pcharthick=5. & psymsize=1.3 & pthick=7.
    xthick=4. & xcharsize=1.2
    ythick=4. & ycharsize=1.2
    zthick=4. & zcharsize=1.2
  end
  'poster': begin			;== 'large posters'
    pcharsize=1.1 & pcharthick=3. & psymsize=1.3 & pthick=3.
    xthick=3. & xcharsize=1.1
    ythick=3. & ycharsize=1.1
    zthick=3. & zcharsize=1.1
  end
  'jeremy': begin			;== 'jeremy'
    pcharsize=1.2 & pcharthick=4. & psymsize=1.2 & pthick=6.
    xthick=5. & xcharsize=1.1
    ythick=5. & ycharsize=1.1
    zthick=5. & zcharsize=1.1
  end
  'vinay': begin		;== 'vinay'
    if !d.NAME eq 'PS' or !d.NAME eq 'PRINTER' then begin
      pcharsize=1.4 & pcharthick=5. & psymsize=1.3 & pthick=5.
      xthick=5. & xcharsize=1.5
      ythick=5. & ycharsize=1.5
      zthick=5. & zcharsize=1.5
    endif else begin
      pcharsize=1.3 & pcharthick=2. & psymsize=1.2 & pthick=2.
      xthick=2. & xcharsize=1.4
      ythick=2. & ycharsize=1.4
      zthick=2. & zcharsize=1.4
    endelse
  end
  else: begin
    message,'Unrecognized MILIEU; not doing anything!',/info
    return
  end

endcase

!P.CHARSIZE=pcharsize*zm
!P.CHARTHICK=pcharthick*zm
!P.SYMSIZE=psymsize*zm
!P.THICK=pthick*zm
;
!X.THICK=xthick*zm
!X.CHARSIZE=xcharsize*zm
;
!Y.THICK=ythick*zm
!Y.CHARSIZE=ycharsize*zm
;
!Z.THICK=zthick*zm
!Z.CHARSIZE=zcharsize*zm

return
end
