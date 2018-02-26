pro chandra_deroll,skyx,skyy,roll,xout,yout,xcen=xcen,ycen=ycen,$
	hrci=hrci,hrcs=hrcs,acisi=acisi,aciss=aciss,verbose=verbose,$
	_extra=e
;+
;procedure	chandra_deroll
;	convert the input SKY coordinates to how it would look with
;	ROLL_NOM=0
;
;syntax
;	chandra_deroll,skyx,skyy,roll,xout,yout,xcen=xcen,ycen=ycen,$
;	/hrci,/hrcs,/acisi,/aciss,verbose=verbose
;
;parameters
;	skyx	[INPUT; required] SKYX in [pix]
;	skyy	[INPUT; required] SKYY in [pix]
;	roll	[INPUT; required] ROLL_NOM [degree]
;		* get it from the event list using
;		  dmkeypar file ROLL_NOM echo+
;		* for HRC-I, an additional non-negotiable 45 deg correction
;		  is made to account for the orientation of the detector
;	xout	[OUTPUT; required] derolled X-coordinate, oriented
;		along the LETG dispersion direction
;	yout	[OUTPUT; required] derolled Y-coordinate, oriented
;		along the cross-dispersion direction
;
;keywords
;	xcen	[INPUT] the X center point around which to roll
;	ycen	[INPUT] the Y center point around which to roll
;		* if given, overrides the defaults set up using
;		  /HRCI,/HRCS,/ACISI,/ACISS
;		* if no keywords are given, assumed to be (0,0)
;	hrci	[INPUT] if set, uses default of XCEN=YCEN=16384
;	hrcs	[INPUT] if set, uses default of XCEN=YCEN=32768
;	acisi	[INPUT] if set, uses default of XCEN=YCEN=4096.5
;	aciss	[INPUT] if set, uses default of XCEN=YCEN=4096.5
;		* precedence: HRCI >> HRCS >> ACISI >> ACISS
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (MMXII.VI)
;	bug fix: deroll was not lining up along the SIM axis
;	also, changed default XCEN,YCEN to 4096.5 for ACIS
;	  (VK; MMXIII.VIII)
;	added HRC-I offset correction; allow XCEN and YCEN to be
;	  set to zero explicitly (VK; MMXVI.V)
;-

;usage
ok='ok' & np=n_params() & nx=n_elements(skyx) & ny=n_elements(skyy)
nr=n_elements(roll)
if np lt 5 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='SKYX is not defined' else $
  if ny eq 0 then ok='SKYY is not defined' else $
   if nr eq 0 then ok='ROLL is not defined' else $
    if nx ne ny then ok='SKYX and SKYY are incompatible' else $
     if nr gt 1 then ok='ROLL cannot be an array'
if ok ne 'ok' then begin
  print,'Usage: chandra_deroll,skyx,skyy,roll,xout,yout,xcen=xcen,ycen=ycen,$'
  print,'       /hrci,/hrcs,/acisi,/aciss,verbose=verbose'
  print,'  deroll SKY coordinates'
  if np ne 0 then message,ok,/informational
  return
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
xc=0. & yc=0.
if keyword_set(aciss) then begin & xc=4096.5 & yc=xc & endif
if keyword_set(acisi) then begin & xc=4096.5 & yc=xc & endif
if keyword_set(hrcs) then begin & xc=32768. & yc=xc & endif
if keyword_set(hrci) then begin & xc=16384. & yc=xc & endif
if n_elements(xcen) gt 0 then xc=xcen[0]
if n_elements(ycen) gt 0 then yc=ycen[0]
zroll=roll[0] & if keyword_set(hrci) then zroll=roll[0]-45.	;be sure to change this HRCI offset (currently 45) to whatever Mike J. and Ping Z. converge to
theta=zroll*!pi/180.	;convert to radians
if vv gt 10 then message,'Rolling '+strtrim(roll[0],2)+' deg around (Xc,Yc)='+$
	strtrim(xc,2)+','+strtrim(yc,2),/informational

;	deroll
;	(this originally had +theta, and the + and - of second term interchanged)
sinth=sin(-theta) & costh=cos(-theta)
xout=(skyx-xc)*costh + (skyy-yc)*sinth
yout=(skyx-xc)*sinth - (skyy-yc)*costh

;	plot if asked
if vv gt 50 then begin
  nn=nx < 10000L
  ii=long(randomu(seed,nn)*nx)
  plot,skyx[ii]-xc,skyy[ii]-yc,psym=3,/xs,/ys,title='Roll = '+strtrim(roll[0],2)+' [deg]'
  oplot,xout[ii],yout[ii],psym=3,col=200
endif

if vv gt 1000 then stop,'halting; type .CON to continue'

return
end
