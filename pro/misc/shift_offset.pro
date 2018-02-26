function shift_offset,arr1,arr2,barr1,barr2,bkgscal=bkgscal,nboot=nboot,$
	arange=arange,deltax=deltax,sigx=sigx,maxstep=maxstep,$
	verbose=verbose, _extra=e
;+
;function	shift_offset
;	computes and returns the best offset match between two arrays
;
;syntax
;	dx=shift_offset(arr1,arr2,barr1,barr2,bkgscal=bkgscal,nboot=nboot,$
;	arange=arange,deltax=deltax,sigx=sigx,verbose=verbose)
;
;parameters
;	arr1	[INPUT; required] list of values
;	arr2	[INPUT; required] list of values to compare to ARR1
;		* ARR2 is swept across ARR1 from when their maxima match
;		  till their minima match.  So make sure there is some
;		  extra buffer in the range of ARR2 relative to ARR1,
;		  otherwise this routine does nothing
;	barr1	[INPUT] list of values to use as background for ARR1
;		* currently ignored
;	barr2	[INPUT] list of values to use as background for ARR2
;		* if not given, uses BARR1
;		NOTE: background arrays are currently not implemented
;		and are silently ignored
;		NOTE: all arrays must be of size >2
;
;keywords
;	bkgscal	[INPUT] the background area scaling
;		* default is 1
;		* if 2-element array, first element is for ARR1/BARR1 and
;		  second element is for ARR2/BARR2
;		* currently ignored because all background is ignored
;	nboot	[INPUT] number of bootstrap samples to pick to get the error bar
;		* default is 50
;	arange	[INPUT] the range of ARR1 to actually consider for purposes
;		of matching
;		* ignored unless a 2-element array
;		* if range falls outside of actual bounds of ARR1, the latter
;		  are used
;	deltax	[OUTPUT] the best offset found, same as the main return value,
;		here only for redundancy and clarity
;		* _add_ DELTAX to ARR2 to match to ARR1
;	sigx	[OUTPUT] the error on DELTAX, computed via bootstrap
;	maxstep	[INPUT] maximum number of steps to step through to find offset
;		* default is 1000
;		* cannot set this to a number smaller than 10
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2016-aug; based on call_findoffsets)
;-

;	usage
ok='ok' & np=n_params() & n1=n_elements(arr1) & n2=n_elements(arr2)
nb1=n_elements(barr1) & nb2=n_elements(barr2)
if np lt 2 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='ARR1 is not defined' else $
  if n2 eq 0 then ok='ARR2 is not defined' else $
   if n1 lt 3 then ok='ARR1 is too small' else $
    if n2 lt 3 then ok='ARR2 is too small'
if ok ne 'ok' then begin
  print,'Usage: dx=shift_offset(arr1,arr2,barr1,barr2,bkgscal=bkgscal,nboot=nboot,$'
  print,'       arange=arange,deltax=deltax,sigx=sigx,verbose=verbose)'
  print,'  computes and returns shift offset between ARR1 and ARR2'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
ok='ok'
if nb2 eq 0 then ok='BARR2 not defined' else if nb2 lt 3 then ok='BARR2 is too small'
if ok ne 'ok' then begin
  nb2=nb1
  if nb1 ge 3 then BARR2=BARR1
endif
;
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
bkgscl=1. & if keyword_set(backscal) then bkgscl=float(backscal[0])
;
mboot=50L & if keyword_set(nboot) then mboot=long(nboot[0])>0L
;
maxnumstep=1000L & if keyword_set(maxstep) then maxnumstep=long(maxstep[0])>10

;	figure out the range to test
a1min=min(arr1,max=a1max) & a2min=min(arr2,max=a2max)
xrange=[a1min,a1max]
if n_elements(arange) eq 2 then begin
  mina1=min(arange,max=maxa1)
  if mina1 gt a1min then xrange[0]=mina1
  if maxa1 lt a1max then xrange[1]=maxa1
endif
dA=xrange[1]-a2max & dB=xrange[0]-a2min & dmin=min([dA,dB],max=dmax)

;	don't bother if ARR2 range is smaller than ARR1
if a2max-a2min le xrange[1]-xrange[0] or dmax-dmin eq 0 then begin
  message,'ARR2 must have a buffer zone to sweep through ARR1',/informational
  message,'Perhaps you should consider switching the order of the inputs',/informational
  message,'As is, no point in continuing; returning with an error',/informational
  sigx=!values.F_NAN
  return,-1L
endif

o1=where(arr1 ge xrange[0] and arr1 le xrange[1],mo1)
if mo1 lt 4 then begin
  message,'nothing to compare to in ARR1 in selected range!',/informational
  return,-1L
endif

;	how many steps in this range?
ndx=(2L*(n1+n2))<maxnumstep
ddx=(dmax-dmin)/ndx & dxarr=findgen(ndx)*ddx+dmin & dksarr=fltarr(ndx,mboot+1L)
for i=0L,ndx-1L do begin	;{step through the offsets
  tmp=arr2+dxarr[i]
  o2=where(tmp ge xrange[0] and tmp le xrange[1],mo2)
  if mo2 gt 3 then begin
    kstwo,arr1[o1],tmp[o2],dks,pks & dksarr[i,mboot]=dks
  endif else dksarr[i,mboot]=!values.F_NAN
;  for j=0L,mboot-1L do begin	;{bootstrap on arr2
;    ii=long(randomu(seed,n2)*n2)
;    kstwo,arr1,arr2[ii]+dxarr[i],dks,pks & dksarr[i,j]=dks
;  endfor			;J=0,MBOOT-1}
endfor				;I=0,NDX-1}

;	bootstrap the error
for j=0L,mboot-1L do begin	;{bootstrap on ARR2
  ii=long(randomu(seed,n2)*n2) & altarr2=arr2[ii] & os=sort(altarr2) & altarr2=altarr2[os]
  for i=0L,ndx-1L do begin	;{step through the offsets
    tmp=altarr2+dxarr[i]
    o2=where(tmp ge xrange[0] and tmp le xrange[1],mo2)
    if mo2 gt 3 then begin
      kstwo,arr1[o1],tmp[o2],dks,pks & dksarr[i,j]=dks
    endif else dksarr[i,j]=!values.F_NAN
  endfor			;I=0,NDX-1}
endfor				;J=0,MBOOT-1}

;	where's the minimum?
jnk=min(dksarr[*,mboot],imn,/nan) & deltax=dxarr[imn]

;	what's the error?
dtmp=fltarr(mboot+1L) & xtmp=dtmp
for j=0L,mboot do begin & dtmp[j]=min(dksarr[*,j],imn,/nan) & xtmp[j]=dxarr[imn] & endfor
avx=total(xtmp/dtmp)/total(1./dtmp)
sdx=sqrt(total(xtmp^2/dtmp)/total(1./dtmp)-avx^2)
deltax=avx & sigx=sdx
;sigx=stddev(tmp,/nan)
if vv gt 5 then print,deltax,mean(tmp,/nan),sigx

;DEBUG
;cdf1=findgen(mo1)/(mo1-1.)
;tmp=arr2+deltax & o2=where(tmp ge xrange[0] and tmp le xrange[1],mo2) & cdf2=findgen(mo2)/(mo2-1.)
;plot,arr1[o1],cdf1,/xs & oplot,tmp[o2],cdf2,col=2

if vv gt 1000 then stop,'halting; type .CON to continue'

return,deltax
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;example
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print,'Calling Sequence' & print,''
jnk=shift_offset()
print,''

;	example -- make two Gaussians and find their offset
if not keyword_set(n1) then n1=1000
if not keyword_set(n2) then n2=100
if not keyword_set(nboot) then nboot=50L
if not keyword_set(dx) then dx=0.5
if not keyword_set(maxstep) then maxstep=100L

;	make simulated datasets
arr1=randomn(seed,n1)		;<- centered on 0
arr2=randomn(seed,n2)+dx	;<- centered on DX, should add -DX to match to ARR1
;print,mean(arr1),mean(arr2),mean(arr2)-mean(arr1)

;	find shift
arange=[-1.,1.]
o1=where(arr1 ge arange[0] and arr1 le arange[1],mo1)
o2=where(arr2 ge arange[0] and arr2 le arange[1],mo2)
delx=shift_offset(arr1,arr2,arange=[-1,1],nboot=50,deltax=deltax,sigx=sigx,verbose=verbose,maxstep=maxstep)
print,'measured shift is '+strtrim(deltax,2)+'+-'+strtrim(sigx,2)+' (cf. actual='+strtrim(-dx,2)+'; approx expected='+strtrim(mean(arr1)-mean(arr2),2)+')'

end
