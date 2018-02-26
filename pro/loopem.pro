function loopem,logTmax,EM,logT=logT,sloop=sloop,verbose=verbose, _extra=e
;+
;function	loopem
;	convert a set of integrated emission measures to a differential
;	emission measure distribution assuming a hydrostatic-loop type
;	distribution and return computed DEM
;
;	DEM(T) is generally n^2*dX/dlnT, and
;	EM(T) == \int_0^T dlnT' DEM(T') = \int_0^T dT' DEM(T')/T'
;	In particular, DEM is parameterized as DEM(T)=a*T^b
;	So for a given EM(T), we obtain
;		EM(T)=(a/b)*(Tmax^(b)-Tmin^(b)), Tmin=0, b>0
;	with Tmin=0 for b>0
;
;syntax
;	dem=loopem(logTmax,EM,logT=logT,sloop=sloop,verbose=verbose)
;
;parameters
;	logTmax	[INPUT; required] loop top temperature(s) in log10(deg K)
;	EM	[INPUT] if given, normalize the output DEM to give a
;		total Emission Measure equal to this.
;		* size must match that of logTmax -- if fewer, fills out
;		  with 1st element; if more, ignored
;
;keywords
;	logT	[INPUT] temperature grid over which output DEM is returned
;		* default is 4..8 in steps of 0.05
;	sloop	[INPUT] slope of DEM in log(DEM)-log(T) space
;		* default is 1.5
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Nov01)
;	changed NORM calc to allow for large values of SLOOP (VK; Feb05)
;-

;	usage
ok='ok' & np=n_params() & ntmx=n_elements(logTmax)
if np eq 0 then ok='Insufficient parameters' else $
 if ntmx eq 0 then ok='logTmax is undefined'
if ok ne 'ok' then begin
  print,'Usage: dem=loopem(logTmax,EM,logT=logT,sloop=sloop,verbose=verbose)'
  print,'  return DEM(T) [cm^-5/logK] of a static RTV loop'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	check inputs
vv=0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;
logTm=logTmax
;
EMtot=dblarr(ntmx)+1.0 & nem=n_elements(EM)
if nem gt 0 then begin
  EMtot=dblarr(ntmx)+double(EM[0])
  if nem le ntmx then EMtot[0L:nem-1L]=EM else EMtot=EM[0L:ntmx-1L]
endif
;
tlog=findgen(81)*0.05+4. & nT=n_elements(logT)
if nT ne 0 then tlog=logT & nT=n_elements(tlog)
DEM=dblarr(nT)
;
ss=1.5 & if keyword_set(sloop) then ss=sloop[0]
ok='ok'
if ss le 0 then begin
  message,'slope cannot be -ve',/info
  if vv gt 4 then message,'Assuming delta-fn T distributions',/info
  for i=0L,ntmx-1L do begin
    jnk=max(abs(logTm[i]-tlog),imx)
    dlogT=(logT[imx]-logT[[imx-1L]]) > (logT[[imx+1L]]-logT[imx])
    if dlogT eq 0 then dlogT=1.
    DEM[imx]=DEM[imx]+EMtot[i]/dlogT
  endfor
  return,DEM
endif

;	initialize
tmin=0. & if ss eq 0 then tmin=10.^(min(logT))

;	compute the normalization
if tmin eq 0 then norm=10.D^(alog10(EMtot)+alog10(ss)-ss*logTm) else $
	norm = EMtot*ss/((10.D^(logTm))^(ss)-(tmin)^(ss))

;	accumulate
for i=0L,ntmx-1L do begin
  oo=where(tlog le logTm[i],moo)
  if moo gt 0 then DEM[oo]=DEM[oo]+norm[i]*(10.D^(tlog[oo]))^(ss)
endfor

if vv gt 1000 then message,'HALT FOR DEBUGGING'

return,DEM
end
