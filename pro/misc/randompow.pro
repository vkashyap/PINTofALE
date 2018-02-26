function randompow,seed,num,alpha=alpha,cumul=cumul,outpoi=outpoi,$
	Nerr=Nerr,Xmin=Xmin,Xmax=Xmax,help=help,$
	verbose=verbose, _extra=e
;+
;function	randompow
;	return a sample drawn from a power-law distribution
;
;	let the power-law be defined as
;	  f(x)=N*x^(-A)
;	the cumulative function,
;	  F(X)= (N/(1-A))*(X^(1-A)-Xmin^(1-A))	if A.ne.1
;	      = N*alog(X)-N*alog(Xmin) 		if A.eq.1
;	which may be normalized to go from 0 to 1:
;	  G(X)=(X^(1-A)-Xmin^(1-A))/(Xmax^(1-A)-Xmin^(1-A))	if A.ne.1
;	      =(alog(X)-alog(Xmin))/(alog(Xmax)-alog(Xmin))	if A.eq.1
;	given a random number 0<r==G(X)<1, solve for X:
;	  X=(r*(Xmax^(1-A)-Xmin^(1-A))+Xmin^(1-A))^(1/(1-A))	if A.ne.1
;	   =Xmin*exp(r*alog(Xmax/Xmin))				if A.eq.1
;
;syntax
;	x=randompow(seed,num,alpha=alpha,/cumul,/outpoi,Nerr=Nerr,$
;	Xmin=Xmin,Xmax=Xmax,/help,verbose=verbose)
;
;parameters
;	seed	[INPUT; necessary] seed for random number generator
;		* this parameter must be present, but there is no
;		  need for it to be defined on input.
;		* will be passed w/o modification to RANDOMU()
;	num	[INPUT] this has a different meaning depending on
;		whether the keyword CUMUL is set or not.
;		* if CUMUL is not set, this is the number of random
;		  samples that are returned, with a minimum of 1
;		* if CUMUL is set, this is taken to be the sum of the
;		  numbers sampled, with a default of max{1,XMIN}
;
;keywords
;	alpha	[INPUT] index of power-law
;		* default is 1.8
;	cumul	[INPUT] if set, generates random samples until the
;		sum of the values matches NUM.  otherwise, the sum
;		is irrelevant and NUM samples are extracted
;	outpoi	[INPUT] if set, a Poisson deviate is obtained from
;		the sample generated
;		* if CUMUL is set, the Poisson deviates are what are
;		  summed to match against NUM
;	Nerr	[INPUT] error on NUM.  The actual SUM{X} matches a random
;		realization of NCT
;		* if given, assumed to be 1-sigma Gaussian error
;		* if not given, is set to sqrt(NUM+0.75)+1
;	Xmin	[I/O] minimum value in sample
;		* default: 2
;	Xmax	[I/O] maximum value in sample
;		* default: max{Normal(NUM,NERR),XMIN}
;	help	[INPUT] if set, prints out usage and exits
;	verbose	[INPUT] controls chatter
;
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run randompow
;
;history
;	vinay kashyap (Nov09; based loosely on samplpl.pro)
;	allowed XMIN,XMAX to be defined on return, modified example
;	  calls to be self-consistent (VK; JanMMX)
;	allowed large photon counts by using N(mu,sqrt(mu)) to take
;	  over as Poisson (VK; JanMMXIII)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(seed) & nn=n_elements(num)
if np eq 0 then ok='Insufficient parameters' else $
 if keyword_set(help) then ok='exiting' else $
  if nn gt 1 then ok='NUM cannot be an array'
if ok ne 'ok' then begin
  print,'Usage: x=randompow(seed,num,alpha=alpha,/cumul,/outpoi,Nerr=Nerr,$'
  print,'       Xmin=Xmin,Xmax=Xmax,/help,verbose=verbose)'
  if np gt 0 then message,ok,/informational
  print,'  return a sample drawn from a power-law distribution'
  return,-1L
endif

;	check input
Nct=1L & if nn ne 0 then Nct=abs(num[0])

;	check keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
aa=1.8 & if keyword_set(alpha) then aa[0]=alpha[0]
;
Nsig=sqrt(abs(Nct)+0.75)+1. & if keyword_set(Nerr) then Nsig=abs(Nerr[0])
;
Zmin=2. & if keyword_set(Xmin) then Zmin=Xmin[0]
;
Zmax=randomu(seed)*Nsig+Nct & if keyword_set(Xmax) then Zmax=Xmax[0]
Zmax=Zmax>Zmin
;
Xmin=Zmin & Xmax=Zmax

;...............................................................................
;	compute
if not keyword_set(cumul) then begin	;(sample "with replacement"
  rr=randomu(seed,Nct)
  if aa eq 1 then X=Zmin*exp(rr*alog(Zmax/Zmin)) else $
    X=(rr*(Zmax^(1.-aa)-Zmin^(1.-aa))+Zmin^(1.-aa))^(1./(1.-aa))
  if keyword_set(outpoi) then begin
    Z=lonarr(Nct)
    for i=0L,Nct-1L do if X[i] gt 0 then Z[i]=randomu(seed,poisson=X[i])
    X=Z
  endif
endif else begin			;CUMUL=0)(sample to accumulate
  oldway=0
  go_on=1 & storc=0 & lamN=Nct & kk=0L
  while go_on do begin
    if lamN gt 1e5 then begin
      ir=(randomn(seed)*sqrt(lamN)+lamN) > 1
    endif else ir=randomu(seed,poisson=lamN)>1
    xx=randompow(seed,ir,alpha=aa,Xmin=Zmin,Xmax=Zmax,outpoi=outpoi)
    zc=total(xx) & pz=-(Nct-zc)^2/Nsig^2/2. & rz=alog(randomu(seed))
    if pz gt rz then begin
      go_on=0
      if vv gt 5 then print,vv
      if vv gt 3 then print,ir,lamN,Nct,zc,pz,rz
    endif else begin
      if zc gt 0 then lamN=lamN*(Nct/zc) else lamN=lamN*2.
    endelse
    if keyword_set(storc) eq 0 then begin
      storc=zc
      storp=pz
      storr=rz
      pmin=pz & xbest=xx
    endif else begin
      storc=[storc,zc]
      storp=[storp,pz]
      storr=[storr,rz]
      if pz lt pmin then xbest=xx
    endelse
    if kk eq (100L/((vv>1L)<100))*fix(kk*(vv<100)/100L) then kilroy
    if kk gt 10000L then begin
      if vv gt 0 then message,'Not converging; exiting with closest match',/informational
      if vv gt 1000 then stop,'HALTing; type .CON to continue'
      xx=xbest & go_on=0
      if vv gt 3 then print,'**'
      if vv gt 1 then print,ir,lamN,Nct,zc,pz,rz
    endif
    X=xx & kk=kk+1L
  endwhile

endelse					;CUMUL=1)

return,X
end

;	usage
jnk=randompow()

;	example
if not keyword_set(v_alpha) then v_alpha=1.5	;set ALPHA
if n_elements(verbose) eq 0 then verbose=1	;set verbosity

print,''
print,';	get a set of random deviates from the power-law distribution'
print,';	such that they add up to a specific quantity'
print,'X=randompow(seed,5e4,alpha=v_alpha,/cumul,outpoi=outpoi)
minx=0 & maxx=0
x=randompow(seed,5e4,alpha=v_alpha,/cumul,verbose=verbose,outpoi=outpoi,xmin=minx,xmax=maxx)
nx=n_elements(x)	;there are these many in the sample
xmin=min(X,max=xmax)	;need this range for direct draws

print,';	plot up the sample'
print,'plot,X[sort(X)],findgen(nx)/(nx-1),/xl'
plot,x[sort(x)],findgen(nx)/(nx-1),/xl,xtitle='X',ytitle='F(<X)',title='!4a!X='+strtrim(v_alpha,2),xr=[1,xmax]

print,''
print,';	you can also draw the same number directly'
print,'Z=randompow(seed,NX,alpha=v_alpha,Xmin=min(X),Xmax=total(X),outpoi=outpoi)'
z=randompow(seed,nx,alpha=v_alpha,Xmin=minx,Xmax=maxx,verbose=verbose,outpoi=outpoi)
nz=n_elements(z)	;this should be identical to NX
zmin=min(z,max=zmax)	;this won't necessarily match [XMIN,XMAX]

print,';	overplot the second sample'
peasecolr
print,'oplot,Z[sort(Z)],findgen(nz)/(nz-1),col=2'
oplot,z[sort(z)],findgen(nz)/(nz-1),col=2

print,''
print,';	summary'
print,'X : total='+strtrim(total(X),2)+'; number='+strtrim(nx,2)+'; minmax='+strtrim(xmin,2)+','+strtrim(xmax,2)
print,'Z : total='+strtrim(total(Z),2)+'; number='+strtrim(nZ,2)+'; minmax='+strtrim(zmin,2)+','+strtrim(zmax,2)

end
