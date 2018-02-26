pro splac,x,y,ytol,u,v,w,noband=noband,discon=discon,verbose=verbose, _extra=e
;+
;procedure	splac
;	compute Segmented Piecewise Linear Approximation to Curve.
;	generates piecewise linear approximations with fewest line
;	segments within given tolerances.
;
;about the name
;	have you heard the one about the Stanford Linear Accelerator (SLAC)
;	and what would happen after the Big One hits?  Obviously, it would
;	have to be renamed SPLAC, the Stanford Piecewise Linear Accelerator..
;
;syntax
;	splac,x,y,ytol,u,v,w,/noband,/discon,verbose=verbose
;
;parameters
;	x	[INPUT; required] points at which data are defined
;		* must be in ascending order
;	y	[INPUT; required] values of data array, Y(X)
;		* must match size of X
;	ytol	[INPUT; required] tolerance of approximation: allowed
;		>deviation< from the true curve
;		* if scalar, single absolute-value tolerance for entire curve.
;		* if vector, tolerance at each X -- size must match X
;		* if -ve, ABS(YTOL) is taken to be the fractional factor
;		  relative to the true value.
;	u	[OUTPUT] a partition of the interval [X(0),X(N[X]-1)]
;		with U(0)=X(0) and U(N[U]-1)=X(N[X]-1)
;	v	[OUTPUT] ordinates, V(U)
;		* the approximating segment is (u[i],v[i]) -> (u[i+1],v[i+1])
;	w	[OUTPUT] ordinates if the approximation may be discontinuous
;		* the approximating segment is (u[i],w[i]) -> (u[i+1],v[i+1])
;
;keywords
;	noband	[INPUT] indicates whether the approximation is to be
;		restricted to the 'tolerance band' about the data
;		(the 'tolerance band' is a piecewise linear band
;		centered at the data whose width is determined by the
;		tolerances at the data points.)
;		* if set, indicates no band restriction
;		* if not set, or set to 0, indicates apply this restriction
;	discon	[INPUT] indicates whether or not the approximation must be
;		continuous
;		* if set, indicates continuity not required
;		* if not set, or set to 0, indicates continuity is required
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;examples
;	x=findgen(100) & y=sqrt(x)
;	splac,x,y,0.1,u,v,w,/noband,verbose=10
;	splac,x,y,0.3,u,v,w,/discon,verbose=10
;	splac,x,y,-0.1,u,v,w & plot,x,y & oplot,u,v,psym=-1,col=100
;
;subroutines
;	NONE
;
;restrictions
;	Will not work with IDL versions < 5
;
;these are all the smoothing tools in PINTofALE
;	ALALOESS() makes a loess curve
;	CLTSMOOTH() removes -ves from background subtracted spectra or light curves
;	CONV_RMF convolves with specified RMF
;	HAARTRAN() smoothes by filtering on Haar wavelet coefficients
;	LINEREM() removes lines from spectrum by iterative S/N filtering
;	NOISMOOTH() does boxcar accumulation a la Ebeling's asmooth
;	REGROUP() accumulates from one end
;	SCRMF does fast convolution using optimized RMF
;	SMOOTHIE does peak-descent accumulation to build up S/N
;	SPLAC computes a segmented piecewise linear approximations
;	UNKINK() removes sharp kinks from a curve
;	VARSMOOTH() does local smoothing at specified scales
;	VOORSMOOTH() does fast local smoothing at specified scales
;
;history
;	From D.G. Wilson, 1976, ACM Transactions on Mathematical Software,2,388
;	http://www.netlib.org/toms/510
;	  converted to IDL by Vinay Kashyap, who is very unhappy about all
;	  those GOTOs, but doesn't have the patience to rewrite from scratch
;	  (FebMM)
;	bug correction for when YTOL was -ve; made into IDL5 (VK; JunMM)
;	bug correction in plotting with /DISCON (VK; AprMMVII)
;-

;	usage
ok='ok'
np=n_params() & nx=n_elements(x) & ny=n_elements(y) & ntol=n_elements(ytol)
if np lt 4 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='Input data points undefined' else $
  if ny eq 0 then ok='Input data values undefined' else $
   if ntol eq 0 then ok='Tolerance undefined' else $
    if nx ne ny then ok='Input data mismatch' else $
     if ntol gt 1 and ntol ne nx then ok=$
	'Tolerance is a virtue sadly lacking here'
if ok ne 'ok' then begin
  print,'Usage: splac,x,y,ytol,u,v,w,/noband,/discon,verbose=verbose'
  print,'  compute a Segmented Piecewise Linear Approximation to Curve'
  if np ne 0 then message,ok,/info
  return
endif

;	check inputs
; there must be at least 3 points in curve
if nx lt 3 then begin
  message,'Nothing to approximate here, is there?',/info
  u=x & v=y & w=y & return
endif
; X must be in ascending order
dx=x[1:*]-x & oo=where(dx lt 0,moo)
if moo gt 0 then begin
  message,'Input data points not in ascending order; doing nothing.',/info
  u=x & v=y & w=y & return
endif
; YTOL must be +ve, and if not is taken to be fractional deviations
toly=0.*y
if ntol eq 1 then begin
  if ytol[0] lt 0 then toly=abs(ytol[0]*y) else toly[*]=ytol[0]
endif else begin
  oo=where(ytol lt 0,moo)
  if moo gt 0 then toly[oo]=abs(ytol[oo]*y[oo])
  oo=where(ytol ge 0,moo)
  if moo gt 0 then toly[oo]=ytol[oo]
endelse

;	keywords
if not keyword_set(noband) then noband=0	;noband=(0,1)::i2=(2,1)
if not keyword_set(discon) then discon=0	;discon=(0,1)::i3=(3,1)
i2=2 & if keyword_set(noband) then i2=1
i3=3 & if keyword_set(discon) then i3=1
itch=i2*i3
vv=0 & if keyword_set(verbose) then vv=fix(verbose) > 1

;	outputs
u=x & v=y & w=y

;	INITIALIZATION FOR ENTIRE PROGRAM
form="($,a)"
epsln=toly[1-1l] & sgn=1.0 & keep=1L & i=1L & u[1-1l]=x[1-1l] & j=2L
init=1L & indc=0L
goto, L30					;(30->

;	INITIALIZATION FOR EACH SEGMENT
L20: if vv ge 5 then print,form=form,'-20'		;(20<-

j=j+1L & init=i & indc=0L
if itch lt 3 then keep=i	;DISCON=1
if abs(itch-4) ne 2 then $	;NOBAND=1
	goto, L30				;(30->

;	RESTRICTED TO TOLERANCE BAND
xeye=u[j-1L-1l] & yeye=v[j-1L-1l]
temp1=epsln+(sgn*toly[i-1L-1l]-epsln)*(x[i-1l]-u[j-1L-1l])/(x[i-1l]-x[i-1L-1l])
yinit=yeye-temp1-temp1
goto, L40					;(40->

L30: if vv ge 5 then print,form=form,'-30'		;->30))
if vv ge 1 and vv lt 5 then print,form=form,'.'

;	NOT RESTRICTED TO TOLERANCE BAND
xeye=x[i-1l] & yeye=y[i-1l]+epsln & yinit=y[i-1l]-epsln
if itch eq 1 or i eq 1 then $	;NOBAND=DISCON=1
	goto, L40				;(40->
temp1 = sgn*toly[i+1L-1l]
smin = (y[i+1L-1l]-yeye-temp1)/(x[i+1L-1l]-xeye)
temp1 = sgn*toly[i-1L-1l]
smax = (yeye-y[i-1L-1l]+temp1)/(xeye-x[i-1L-1l])
if keep eq i-1L then goto, L50			;(50->
it=i-2L & xinit=xeye & ipiv=i & igraze=i & i=i+1L
goto, L150					;(150->

L40: if vv ge 5 then print,form=form,'-40'		;->40))

if xeye ge x[i-1l] then i=i+1L
epsln = sgn*toly[i-1l]
dx = x[i-1l]-xeye
smax = (y[i-1l]+epsln-yeye)/dx
smin = (y[i-1l]-epsln-yeye)/dx

L50: if vv ge 5 then print,form=form,'-50'		;->50)

xinit=xeye & ipiv=i & igraze=i

;	DETERMINATION OF INDIVIDUAL SEGMENT
L60: if vv ge 5 then print,form=form,'-60'		;((60<-
if vv eq 4 then print,form=form,'.'

if i eq nx then goto, L260			;(260->
i = i+1L

L70: if vv ge 5 then print,form=form,'-70'		;((70<-

;	TEST FOR NEW *MAX* SLOPE
dx=x[i-1l]-xeye & epsln=sgn*toly[i-1l]
temp1=(y[i-1l]+epsln-yeye)/dx
test = temp1-smax
if sgn le 0 then test=-test
if test eq 0 then goto, L90			;(90->
if test gt 0 then goto, L100			;(100->

;	TEST FOR END OF CANDIDATE SEGMENT
test = temp1-smin
if sgn le 0 then test=-test
if test lt 0 then goto, L210			;(210->
smax = temp1

L90: if vv ge 5 then print,form=form,'-90'		;->90)

;	TEST FOR NEW *MIN* SLOPE
ipiv = i

L100: if vv ge 5 then print,form=form,'-100'		;->100)

temp2 = (y[i-1l]-epsln-yeye)/dx
test = temp2 - smax
if sgn le 0 then test=-test
if test lt 0 then goto, L110			;(110->
if test eq 0 then goto, L120			;(120->
if test gt 0 then goto, L140			;(140->
L110: if vv ge 5 then print,form=form,'-110'		;->110)
test=smin-temp2
if sgn le 0 then test=-test
if test lt 0 then goto, L120			;(120->
if test eq 0 then goto, L130			;(130->
if test gt 0 then goto, L60			;->60)
L120: if vv ge 5 then print,form=form,'-120'		;->120))
smin=temp2
L130: if vv ge 5 then print,form=form,'-130'		;->130))
igraze=i
goto, L60					;->60)

;	CHECK FOR PIVOT AT NEW EYE POINT
L140: if vv ge 5 then print,form=form,'-140'		;->140)

if xeye eq x[ipiv-1l] then goto, L220		;(220->
epsln = sgn*toly[ipiv-1l]
indc=1
svx=xeye & svy=yeye & svmn=smin & svmx=smax
xeye=x[ipiv-1l] & yeye=y[ipiv-1l]+epsln
smin=smax & smax=(yinit-yeye)/(xinit-xeye)
if keep ge ipiv then goto, L170			;(170->
it = ipiv-1L

L150: if vv ge 5 then print,form=form,'-150'		;->150)

temp2=yeye+epsln
for l=keep,it do begin				;{160->
  temp2=yeye+sgn*toly[l-1l]
  temp1=(y[l-1l]-temp2)/(x[l-1l]-xeye)
  test=temp1-smax
  if sgn le 0 then test=-test
  if test lt 0 then smax=temp1
endfor						;L=KEEP,IT ->160}

L170: if vv ge 5 then print,form=form,'-170'		;->170)

if ipiv ge i-1L then goto, L70			;70->)
it=i-2L
temp2=yeye-epsln
idiot=ipiv
for l=idiot,it do begin				;{200->
  dx=x[l+1L-1l]-xeye
  temp2=yeye-sgn*toly[l+1L-1l]
  temp1=(y[l+1L-1l]-temp2)/dx
  test=temp1-smax
  if sgn le 0 then test=-test
  if test lt 0 then smax=temp1			;(->180<-)
  if test le 0 then ipiv=l+1L			;(->190<-)
endfor						;L=IDIOT,IT ->200}

goto, L70					;70->)

;	END OF CURRENT SEGMENT
L210: if vv ge 5 then print,form=form,'-210'		;->210)
temp2=smin
if i eq nx then goto, L240			;(240->
keep=igraze
goto, L250					;(250->

L220: if vv ge 5 then print,form=form,'-220'		;->220)

temp2=smax
if i eq nx then goto, L230			;(230->
sgn=-sgn & epsln=-epsln & keep=ipiv
goto, L250					;(250->

L230: if vv ge 5 then print,form=form,'-230'		;->230)

if indc eq 0 or xeye ne x[nx-2L-1l] then goto, L240;(240->
xeye=svx & yeye=svy & smin=svmn & smax=svmx

L240: if vv ge 5 then print,form=form,'-240'		;->240))

u[j-1l]=x[nx-1L-1l] & yinit=y[nx-1L-1L]
goto, L270					;(270->

L250: if vv ge 5 then print,form=form,'-250'		;->250))

if abs(itch-4) ne 2 then $	;NOBAND=1
	goto, L300				;(300->

;	DETERMINE KNOT ON EDGE OF TOLERANCE BAND
temp1=0.0
temp1=epsln-sgn*toly[i-1L-1l]
temp1=(y[i-1l]-y[i-1L-1l]+temp1)/(x[i-1l]-x[i-1L-1l])
u[j-1l] = (y[i-1l]+epsln-yeye-temp1*x[i-1l]+temp2*xeye)/(temp2-temp1)
goto, L310					;(310->

L260: if vv ge 5 then print,form=form,'-260'		;->260)

u[j-1l]=x[nx-1l] & yinit=y[nx-1l]

L270: if vv ge 5 then print,form=form,'-270'		;->270)

;	CONTINUITY CHECK FOR LAST SEGMENT
if itch ge 3 or init eq 0 then $	;DISCON=0
	goto, L290				;(290->
it=init-1L & svmx=smax+sgn & temp2=yeye+epsln
for l=kp,it do begin				;{280->
  temp2=yeye+sgn*toly[l-1l]
  temp1=(y[l-1l]-temp2)/(x[l-1l]-xeye)
  test=temp1-svmx
  if sgn le 0 then test=-test
  if test lt 0 then svmx=temp1
endfor						;L=KP,IT ->280}

if abs(svmx-smax+svmx-smin) le abs(smax-smin) then smax=svmx

L290: if vv ge 5 then print,form=form,'-290'		;->290)

;	NEARNESS CHECK FOR LAST SEGMENT
temp2=smax & temp1=yeye+smax*(u[j-1l]-xeye) & test=yinit-temp1
if sgn lt 0 then test=-test
if test gt 0 then goto, L310			;(310->
temp2=smin & temp1=yeye+smin*(u[j-1l]-xeye) & test=yinit-temp1
if sgn lt 0 then test=-test
if test gt 0 then goto, L310			;(310->
temp2=(yinit-yeye)/(u[j-1l]-xeye) & v[j-1l]=yinit
goto, L320					;(320->

L300: if vv ge 5 then print,form=form,'-300'		;->300)

if itch ge 3 then $	;DISCON=0
	goto, L330				;(330->
u[j-1l]=0.5*(x[i-1l]+x[i-1L-1l])
L310: if vv ge 5 then print,form=form,'-310'		;->310)))
v[j-1l]=yeye+temp2*(u[j-1l]-xeye)
L320: if vv ge 5 then print,form=form,'-320'		;->320)

if xeye ne xinit then goto, L330		;(330->
if itch eq 2 then $	;NOBAND=0 && DISCON=1
	goto, L360				;(360->
if itch ne 6 then $	;NOBAND=1 || DISCON=1
	goto, L330				;(330->
if j le 2 then goto, L380			;(380->
goto, L390					;(390->

L330: if vv ge 5 then print,form=form,'-330'		;->330)))

;	RECOMPUTATION OF KNOT FOR CONTINUITY
if j le 2 then goto, L370			;(370->
if slope eq temp2 then goto, L360		;(360->
yinit = v[j-2L-1l]
if itch lt 3 then yinit=w[j-2L-1l]	;DISCON=1
temp1=(xeye*temp2-u[j-2L-1l]*slope+yinit-yeye)/(temp2-slope)
if itch ge 3 then $	;DISCON=0
	goto, L350				;(350->
if temp1 gt xinit then goto, L360		;(360->
test=abs(epsln) & idiot=init-kp
for l=1,idiot do begin				;{340->
  it=init-l
  if temp1 ge x[it-1l] then goto, L350		;(350->
  dx=y[it-1l]-yeye-temp2*(x[it-1l]-xeye)
  test=toly[it-1l]
  if abs(dx) gt test then goto, L360		;(360->

endfor						;L=1,IDIOT ->340}

L350: if vv ge 5 then print,form=form,'-350'		;->350))

u[j-1L-1l]=temp1 & v[j-1L-1l]=yeye+temp2*(u[j-1L-1l]-xeye)
if itch lt 3 then w[j-1L-1l]=v[j-1L-1l]	;DISCON=1
goto, L390					;(390->

L360: if vv ge 5 then print,form=form,'-360'		;(->360))))

w[j-1L-1l]=yeye+temp2*(u[j-1L-1l]-xeye)
goto, L390					;(390->

L370: if vv ge 5 then print,form=form,'-370'		;->370)

if itch lt 3 then $	;DISCON=1
	goto, L360				;<-360)

L380: if vv ge 5 then print,form=form,'-380'		;->380)

v[1-1l]=yeye+temp2*(u[1-1l]-xeye)

L390: if vv ge 5 then print,form=form,'-390'		;->390)))

slope=temp2 & kp=keep
if i lt nx then goto, L20			;<-20)
if x[nx-1l] eq u[j-1l] then goto, L400		;(400->
if itch lt 3 then w[j-1l]=v[j-1l]	;DISCON=1
j=j+1L & u[j-1l]=x[nx-1l] & v[j-1l]=y[nx-1l]

L400: if vv ge 5 then print,form=form,'-400'		;->400)
if vv ge 1 and vv lt 5 then print,form=form,'.'

if j ge 2 and itch lt 3 then v[1-1l]=w[1-1l]	;DISCON=1
k=j-1L

u=u[0L:k] & v=v[0L:k] & w=w[0L:k-1l]

tt='Compression: '+strtrim(nx,2)+' -> '+strtrim(n_elements(u),2)
if vv ge 1 then message,tt,/info
if vv gt 0 then print,''
if vv ge 3 then begin
  if vv ge 10 then begin
    print,"plot,x,y,/xs,/ys,xtitle='X',ytitle='Y',title='"+tt+"',line=1,thick=3"
    if not keyword_set(discon) then print,'oplot,u,v,psym=-1' else $
    	print,'for i=0L,k do oplot,[u[i],u[i+1L]],[w[i],v[i+1L]],psym=-1'
  endif
  plot,x,y,/xs,/ys,xtitle='X',ytitle='Y',title=tt,line=1,thick=3
  if not keyword_set(discon) then oplot,u,v,psym=-1 else $
  	for i=0L,k-1L do oplot,[u[i],u[i+1L]],[w[i],v[i+1L]],psym=-1
endif
if vv ge 100 then stop,'HALTING.  type .CON to continue'

return
end
