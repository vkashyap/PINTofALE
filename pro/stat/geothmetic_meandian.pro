function geothmetic_meandian,YY,epsilon=epsilon,verbose=verbose, _extra=e
;+
;function	geothmetic_meandian
;	iteratively compute the arithmetic and geometric means and the median until convergence
;
;syntax
;	gmdn=geothmetic_meandian(Y,epsilon=epsilon,verbose=verbose)
;
;parameters
;	Y	[INPUT; required] sample for which to compute the geothmetic meandian
;		* must have at least 3 elements
;
;keywords
;	epsilon	[INPUT; default=1e-6] a small number
;	verbose	[INPUT; default=0] controls chatter
;	
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	https://xkcd.com/2435/
;
;WARNING
;	DO NOT USE THIS FOR ACTUAL SCIENCE!
;	The Geothmetic Meandian is biased and non-optimal -- see example run
;
;example
;	.run geothmetic_meandian
;	(makes plots)
;
;history
;	Vinay Kashyap (2021-Mar-11)
;-

;	check input
ok='ok' & np=n_params() & nY=n_elements(YY)
if np eq 0 then ok='Insufficient parameters' else $
 if nY eq 0 then ok='Input sample is not defined' else $
  if nY lt 3 then ok='Y is too small'
if ok ne 'ok' then begin
  print,'Usage: gmdn=geothmetic_meandian(Y,epsilon=epsilon,verbose=verbose)'
  print,' compute the Geothmetic Meandian of a sample'
  print,' see https://xkcd.com/2435/'
  if np ne 0 then message,ok,/informational else print,''
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
zeps=1d-6 & if keyword_set(epsilon) then zeps=abs(epsilon[0])

;	check for negatives
oneg=where(YY lt 0,moneg) & sgn=(-1L)^(moneg)

;	exclude zero from geometric mean
ok=where(YY ne 0,mok)
if mok eq 0 then begin
  if vv gt 0 then message,'all elements are zero; returning zero',/informational
  return,YY[0]
endif

kk=0L & ff=[mean(YY,_extra=e,/nan),sgn*exp(mean(alog(abs(YY[ok])),_extra=e,/nan)),median(YY)] & df=stddev(ff) & fk=ff
go_on=1
while go_on do begin
  kk=kk+1L
  fp=fk[3*kk-3L:3*kk-1L]
  oneg2=where(fp lt 0,moneg2) & sgn2=(-1L)^(moneg2)
  ok2=where(fp ne 0,mok2)
  fp2=[mean(fp,_extra=e,/nan),sgn2*exp(mean(alog(abs(fp[ok])),_extra=e,/nan)),median(fp)] & df2=stddev(fp)
  gmdn=mean(fp2,/nan)
  fk=reform([fk[*],fp2],3,kk+1L)
  if vv ge 1000 then begin
    plot,/nodata,[0,kk],minmax(fk),xtitle='iteration',ytitle='gmdn',/xs,/ys
    oplot,lindgen(kk),fk[0,*],psym=-1,thick=2,symsize=2
    oplot,lindgen(kk),fk[1,*],psym=-2,thick=2,symsize=2
    oplot,lindgen(kk),fk[2,*],psym=-4,thick=2,symsize=2
  endif
  if df2 le zeps then go_on=0
endwhile

if vv gt 1000 then stop,'halting; type .CON to continue'

return,gmdn
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example usage
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	display usage
jnk=geothmetic_meandian()

if not keyword_set(verbose) then verbose=1000L

;	try the XKCD originating example
window,0
print,'geothmetic_meandian([1,1,2,3,5]) = ',geothmetic_meandian([1,1,2,3,5],verbose=verbose)

window,1
;	how does it do for a normal distribution compared to the regular mean?
nsim=500L & nct=[50,80,130,210,340,550,890L] & mct=n_elements(nct) & Ngmdn=fltarr(mct,nsim) & Navg=Ngmdn
for i=0L,mct-1L do begin
  for j=0L,nsim-1L do begin
    rr=randomn(seed,nct[i]) & Navg[i,j]=mean(rr)
    Ngmdn[i,j]=geothmetic_meandian(rr)
  endfor
endfor
xmin=min(Nct,max=xmax)
ymin=min([Navg[*],Ngmdn[*]],max=ymax)
plot,nct,0*nct,xrange=[xmin*0.5,xmax*2],/xl,/xs,yrange=[ymin,ymax],/nodata,xtitle='SAMPLE SIZE',ytitle='CENTRAL MEASURE',title='E[Normal] vs GMDN[Normal]'
for i=0L,mct-1L do oplot,nct[i]*0.9+fltarr(nsim),Navg[i,*],psym=1
for i=0L,mct-1L do oplot,nct[i]*1.1+fltarr(nsim),Ngmdn[i,*],psym=4

window,2
;	how does it do for a Poisson distribution compared to the regular mean?
nsim=100L & nsamp=500L & Pct=[10,20,30,50,80,130,210,340,550,890L] & mPct=n_elements(Pct) & Pgmdn=fltarr(mPct,nsim) & Pavg=Pgmdn
for i=0L,mPct-1L do begin
  for j=0L,nsim-1L do begin
    rr=randomn(seed,nsamp,poisson=Pct[i]) & Pavg[i,j]=mean(rr)-Pct[i]
    Pgmdn[i,j]=geothmetic_meandian(rr,epsilon=0.01)-Pct[i]
  endfor
endfor
xmin=min(Pct,max=xmax)
ymin=min([Pavg[*],Pgmdn[*]],max=ymax)
plot,Pct,0*Pct,xrange=[xmin*0.5,xmax*2],/xl,/xs,yrange=[ymin,ymax],/nodata,xtitle='COUNTS',ytitle='CENTRAL MEASURE',title='E[Poisson] vs GMDN[Poisson]',subtitle='sample size = '+strtrim(nsamp,2)
for i=0L,mPct-1L do oplot,Pct[i]*0.9+fltarr(nsim),Pavg[i,*],psym=1
for i=0L,mPct-1L do oplot,Pct[i]*1.1+fltarr(nsim),Pgmdn[i,*],psym=4

end
