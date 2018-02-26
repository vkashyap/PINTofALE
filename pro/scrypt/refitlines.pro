;+
;REFITLINES.PRO
;	so you ran through FITLINES and got your fits.  and you are
;	ready to publish.
;	but wait, did you use the right model?  sigh.
;	ok, let's have a go at it again.
;
;	need the spectrum, of course.  should be in variables
;	LAMBDA and SPECT[LAMBDA] and SPECTSIG[LAMBDA]
;
;	need the old fit structure.  should be in variable
;	OLDFITSTR
;
;	new fit structure will be in variable, what else, NEWFITSTR
;
;vinay kashyap (MMJul)
;	included HISTERIX (VK; MMIVFeb)
;-

;	check if spectrum exists
nx=n_elements(LAMBDA) & ny=n_elements(SPECT) & nye=n_elements(SPECTSIG)
if nx eq 0 then message,'LAMBDA is undefined'
if ny eq 0 then message,'SPECT[LAMBDA] is undefined'
if nx ne ny and nx ne (ny+1L) then message,$
	'LAMBDA and SPECT[LAMBDA] are incompatible'
if nx eq ny then xLAMBDA=LAMBDA
if nx eq ny+1L then xLAMBDA=0.5*(LAMBDA[1:*]+LAMBDA)	;if bin-boundaries
if nye ne ny then begin
  message,'assuming poisson errors on SPECT',/info
  ySPECTSIG=sqrt(abs(SPECT)+0.75)+1.
endif else ySPECTSIG=SPECTSIG

;	check if fit structure exists
nf=n_tags(OLDFITSTR)
if nf eq 0 then message,'OLDFITSTR does not exist'
if nf lt 19 then message,'OLDFITSTR in unknown format'

;	set some obvious keywords
if n_elements(dchi) eq 0 then dchi=1.0
if n_elements(dumb) eq 0 then dumb=1
if n_elements(funcs) eq 0 then funcs='x3model'
if n_elements(intens) eq 0 then intens=0
if n_elements(verbose) eq 0 then verbose=5

;	extract continuum info from OLDFITSTR
conlev=OLDFITSTR.CONLEV & conlevx=OLDFITSTR.CONLEVX
consig=OLDFITSTR.CONSIG & consigx=OLDFITSTR.CONSIGX
yCONLEV=interpol(conlev,conlevx,xLAMBDA)
yCONSIG=interpol(consig,consigx,xLAMBDA)

;	extract ties, if any
namf=tag_names(OLDFITSTR)
onamf=where(namf eq 'TIES',monamf)
if monamf ne 0 then OLDTIES=OLDFITSTR.TIES else OLDTIES=''
onamf=where(namf eq 'EPITHET',monamf)
if monamf ne 0 then OLDEPITHET=OLDFITSTR.EPITHET else OLDEPITHET=''

;	call FITLINES
message,'calling FITLINES',/info
NEWFITSTR=fitlines(xLAMBDA,SPECT,$
	ysig=ySPECTSIG,funcs=funcs,intens=intens,dchi=dchi,$
	pos=OLDFITSTR.POS,perrp=OLDFITSTR.PERRP,perrm=OLDFITSTR.PERRM,$
	wdt=OLDFITSTR.WDT,werrp=OLDFITSTR.WERRP,werrm=OLDFITSTR.WERRM,$
	flx=OLDFITSTR.FLX,ferrp=OLDFITSTR.FERRP,ferrm=OLDFITSTR.FERRM,$
	perrc=OLDFITSTR.PERRC,werrc=OLDFITSTR.WERRC,ferrc=OLDFITSTR.FERRC,$
	thaw=OLDFITSTR.THAW,type=OLDFITSTR.TYPE,ties=OLDTIES,$
	epithet=OLDEPITHET,conlev=yCONLEV,consig=yCONSIG,comment=comment,$
	histerix=histerix,$
	dumb=dumb,verbose=verbose,$
	xsize=xsize,ysize=ysize,wid=wid,dynrng=dynrng,/posve,/negve,$
	itmax=itmax,chithr=chithr,jumpup=jumpup,jumpdn=jumpdn,$
	svdthr=svdthr,missing=missing)

message,'Fit output in: NEWFITSTR',/info

end
