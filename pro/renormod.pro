function renormod,data,model,err,uchisq=uchisq,uabsdev=uabsdev,$
	umaxdev=umaxdev,ureldev=ureldev,statval=statval,$
	dstep=dstep,eps=eps,maxstep=maxstep,verbose=verbose, _extra=e
;+
;function	renormod
;	given model values and data, find a normalization factor for
;	the model that makes it a better fit to the data and return
;	this new value.
;
;syntax
;	norm=renormod(data,model,err,/uchisq,/uabsdev,/umaxdev,/ureldev,$
;	statval=statval,dstep=dstep,eps=eps,verbose=verbose)
;
;parameters
;	data	[INPUT; required] data values
;	model	[INPUT; required] model values
;		* size must match DATA
;	err	[INPUT; optional] errors to use to determine the goodness
;		of fit statistic.  if not given, assumed to be the Gehrels
;		approximation to Poisson, sqrt(abs(DATA+0.75)+1)
;		* size must match DATA, or must be 1-element
;		* if 1-element and +ve, constant value errors are assumed
;		* if 1-element and -ve, a constant fractional value for
;		  the errors, DATA*abs(ERR) are assumed
;		* WARNING: data points with ERR=0 are ignored when it
;		  is necessary to divide by the errors
;
;keywords
;	uchisq	[INPUT] if set, use the chi-square statistic as the
;		measure of goodness of fit
;		* this is the default
;	uabsdev	[INPUT] if set, use abs(DATA-MODEL) as statistic
;	umaxdev	[INPUT] if set, use min(abs(DATA-MODEL)) as statistic
;	ureldev	[INPUT] if set, use total(abs(DATA-MODEL)/ERR) as statistic
;		* priority: UCHISQ > UABSDEV > UMAXDEV > URELDEV
;	statval	[OUTPUT] "best-fit" value of the statistic
;	dstep	[INPUT] value of the first step
;		* default is 10% of initial value
;	eps	[INPUT] a small number
;		* default is 1e-6
;	maxstep	[INPUT] maximum number of iterations
;		* default is 100
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
;
;history
;	vinay kashyap (Jun2001)
;-

;	usage
ok='ok' & np=n_params()
nd=n_elements(data) & nm=n_elements(model) & ns=n_elements(err)
if np lt 2 then ok='Insufficient input' else $
 if nd eq 0 then ok='Data are missing' else $
  if nm eq 0 then ok='Model is not defined' else $
   if nm ne nd then ok='Data and model size mismatch'
if ok ne 'ok' then begin
  print,'Usage: norm=renormod(data,model,err,/uchisq,/uabsdev,/umaxdev,/ureldev,$
  print,'       statval=statval,dstep=dstep,eps=eps,verbose=verbose)'
  print,'  return renormalization factor for model'
  if np ne 0 then message,ok,/info
  return,0.
endif

;	check error
sigdat=sqrt(abs(DATA+0.75)+1.)
if ns eq nd then sigdat[*]=err[*] else begin
  if ns eq 1 then begin
    if err[0] ge 0 then begin
      message,'assuming constant error-bars',/info
      sigdat[*]=err[0]
    endif else begin
      message,'assuming constant fractional error',/info
      sigdat[*]=data[*]*abs(err[0])
    endelse
  endif else message,'assuming Poisson errors',/info
endelse

;	check keywords
if keyword_set(ureldev) then statqty='RELDEV' else $
 if keyword_set(umaxdev) then statqty='MAXDEV' else $
  if keyword_set(uabsdev) then statqty='ABSDEV' else $
   if keyword_set(uchisq) then statqty='CHISQ' else statqty='CHISQ'
vv=0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1
if not keyword_set(eps) then eps=1e-6
stepmax=100L & if keyword_set(maxstep) then stepmax=long(maxstep) > 1

;	the first value that will be tried out is total(DATA)/total(MODEL)
;	and from there onwards, 3 points are chosen, 2 on either side that
;	are equidistant from the current point, and an extra one on the
;	downstream side at half the step size.  the point with the smallest
;	value of the statistic becomes the central point for the next
;	iteration, and if it happens to be the inner downstream point,
;	or the current point, the step size is halved.  the step size is
;	doubled if the current point turns out to be the maximum.

;	first point
tmpn=total(data) & tmpd=total(model)
z=1. & z0=1. & if tmpd ne 0 then z=tmpn/tmpd
stepsz=0.1*z & if keyword_set(dstep) then stepsz=dstep
;	compute the statistic
case statqty of
  'ABSDEV': statval=total(abs(data-z*model))
  'MAXDEV': statval=max(abs(data-z*model))
  'RELDEV': begin
    ox0=where(sigdat ne 0,mox0)
    if mox0 eq 0 then begin
      message,'No data points found with non-zero errors',/info
      message,'assuming all errors are equal to 1',/info
      sigdat[*]=1.
      ox0=lindgen(nd)
    endif
    statval=total(abs(data[ox0]-z*model[ox0])/sigdat[ox0])
  end
  else: begin
    ox0=where(sigdat ne 0,mox0)
    if mox0 eq 0 then begin
      message,'No data points found with non-zero errors',/info
      message,'assuming all errors are equal to 1',/info
      sigdat[*]=1.
      ox0=lindgen(nd)
    endif
    statval=total((data[ox0]-z*model[ox0])^2/sigdat[ox0]^2)
  end
endcase
oldstatval=statval+100.	;just something to start with

;	iterate
go_on=1 & k=0L
while go_on eq 1 do begin		;{go wander
  ;	pick the downstream points
  if z lt z0 then downsgn=-1 else downsgn=1
  zdownfar=z+downsgn*stepsz & zdownear=z+downsgn*0.5*stepsz
  ;	pick the upstream point
  zup=z-downsgn*stepsz
  ;	get the value of the statistic at all the points
  zz=[zdownfar,zdownear,zup] & ss=fltarr(3)
  for j=0,2 do begin
    case statqty of
      'ABSDEV': ss[j]=total(abs(data-zz[j]*model))
      'MAXDEV': ss[j]=max(abs(data-zz[j]*model))
      'RELDEV': ss[j]=total(abs(data[ox0]-zz[j]*model[ox0])/sigdat[ox0])
      else: ss[j]=total((data[ox0]-zz[j]*model[ox0])^2/sigdat[ox0]^2)
    endcase
  endfor
  ;	which is the new minimum?
  smin=min([statval,ss],imin)
  smax=max([statval,ss],imax)
  ;	update the controls
  if imax eq 0 then stepsz=2.*stepsz
  case imin of
    0: begin	;Z
      if imax ne 0 then stepsz=0.5*stepsz
    end
    1: begin	;ZDOWNFAR
      z0=z & z=zdownfar & oldstatval=statval & statval=smin
    end
    2: begin	;ZDOWNEAR
      z0=z & z=zdownear & oldstatval=statval & statval=smin
      if imax ne 0 then stepsz=0.5*stepsz
    end
    3: begin	;ZUP
      z0=z & z=zup & oldstatval=statval & statval=smin
    end
  endcase
  if vv gt 2 then kilroy
  if vv ge 5 then print,z,statval
  k=k+1L
  ;	stopping rules: stop if
  ;	-- step size gets too small
  if abs(stepsz) lt eps then go_on=0
  ;	-- statval changes by very little
  if abs(statval-oldstatval) le 10.*eps then go_on=0
  ;	-- too many steps
  if k ge stepmax then begin
    message,'too many steps!  quitting here',/info
    go_on=0
    if vv gt 10 then stop,'HALTing! type .CON to continue'
  endif
endwhile				;GO_ON}

return,z
end
