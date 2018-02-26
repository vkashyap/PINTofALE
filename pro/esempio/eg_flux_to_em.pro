;+
;EG_FLUX_TO_EM.PRO
;	example calling program for FLUX_TO_EM.
;
;	start from a set of ID'd lines whose fluxes have
;	been measured, and make EMs appropriate for the
;	given IDs.
;
;	1. read in ID structure from $SCARDIR/pro/esempio/eg_flux_to_em.sav
;	2. extract the line emissivities for each ID'd feature
;	   (and sum them up if multiply ID'd)
;	3. extract the fluxes too
;	4. call FLUX_TO_EM
;	5. make plots showing the output
;
;vinay kashyap (MMJul)
;	added call to GETPOADEF (VK; MMXVAug)
;-

;	read in ID structure
if not keyword_set(dir_SCAR) then dir_SCAR=getpoadef()
if not keyword_set(scardir) then scardir=getenv('SCAR')
if not keyword_set(scardir) then scardir=getenv('SCARDIR')
if not keyword_set(scardir) then scardir=dir_SCAR
idsav=scardir+'/pro/esempio/eg_flux_to_em.sav'
restore,idsav,/verbose
if n_tags(all_fit_ids) eq 0 then message,idsav+$
	': does not contain id structure named ALL_FIT_IDS'

;	squish the IDs
if n_elements(abund) lt 30 then begin
  message,'using Anders & Grevesse abundances',/info
  abund=getabund('anders & grevesse')
endif
sq_allids=squishem(ALL_FIT_IDS,abund=abund)

;	initialize
wvl=sq_allids.WVL & logT=sq_allids.(1).LOGT
nW=n_elements(wvl) & nT=n_elements(logT)

;	extract the emissivities, fluxes, and flux errors
allemis=dblarr(nT,nW) & allflux=fltarr(nW)
allferrp=allflux & allferrm=allflux	;for forwards compatibility
allZ=intarr(nW) & allION=allZ & allsymb=strarr(nW)
for iw=0L,nW-1L do begin
  tmp=sq_allids.(iw+1L).EMIS
  if n_elements(tmp) eq nT then allemis[*,iw]=tmp
  allflux[iw]=(sq_allids.(iw+1L).FLUX)[0]
  ferr=(sq_allids.(iw+1L).FLUXERR)[0]
  allferrp[iw]=ferr & allferrm[iw]=ferr
  allZ[iw]=(sq_allids.(iw+1L).Z)[0]
  allION[iw]=(sq_allids.(iw+1L).ION)[0]
endfor
zion2symb,allZ,allION,allsymb,ziform='ZION'

;	warn if effective areas are not specified
ok='ok' & nea=n_elements(effar) & nwa=n_elements(wvlar)
if nea eq 0 then ok='EFFAR not defined' else $
 if nwa eq 0 then ok='WVLAR not defined' else $
  if nwa ne nea then ok='WVLAR and EFFAR(WVLAR) incompatible'
if ok ne 'ok' then begin
  message,ok,/info
  cc='' & read,prompt='continue? [y/n]',cc
  cc=strtrim(strlowcase(cc),2)
  if cc ne '' and cc ne 'y' then stop,'HALTing.  type .CON to continue'
endif

;	call FLUX_TO_EM
em=flux_to_em(allemis,flux=allflux,logT=logT,wvl=wvl,Z=allZ,NH=NH,$
	defEM=defEM,noph=noph,thresh=thresh, abund=abund,effar=effar,$
	wvlar=wvlar, fh2=fh2,he1=he1,heii=heii,fano=fano)
emm=flux_to_em(allemis,flux=allflux+allferrp,logT=logT,wvl=wvl,Z=allZ,NH=NH,$
	defEM=defEM,noph=noph,thresh=thresh, abund=abund,effar=effar,$
	wvlar=wvlar, fh2=fh2,he1=he1,heii=heii,fano=fano)
emp=flux_to_em(allemis,flux=allflux-allferrm,logT=logT,wvl=wvl,Z=allZ,NH=NH,$
	defEM=defEM,noph=noph,thresh=thresh, abund=abund,effar=effar,$
	wvlar=wvlar, fh2=fh2,he1=he1,heii=heii,fano=fano)

;	plot
if not keyword_set(ymax) then ymax=max(emp)*2.
if not keyword_set(ymin) then begin
  ymin=ymax/100. & oo=where(emm gt 0,moo)
  if moo gt 0 then ymin=min(emm[oo])/2.
endif
if not keyword_set(xmin) then xmin=min(logT)
if not keyword_set(xmax) then xmax=max(logT)
if not keyword_set(labelcolor) then labelcolor='yellow'
if not keyword_set(boundcolor) then boundcolor='red'
if not keyword_set(labelindex) then labelindex=2
if not keyword_set(boundindex) then boundindex=3
;
plot,logT,em,/nodata,xrange=[xmin,xmax],yrange=[ymin,ymax],/xs,/ys,/ylog,$
	xtitle='log!d10!n(T [K])',ytitle='EM [cm!u-5!n]'
setcolor,labelcolor,2
setcolor,boundcolor,3
for iw=0L,nw-1L do begin
  y=reform(em[*,iw]) & oo=where(y gt 0,moo)
  if moo gt 0 then begin
    oplot,logT[oo],y[oo]
    y0=min(y[oo],i0) & x0=logT[oo[i0]]
    c='!3'+allsymb[iw]+' !4k!3'+strtrim(string(wvl[iw],'(f10.2)'),2)
    xyouts,x0,y0,c,color=labelindex
  endif
  ;
  y=reform(emm[*,iw]) & oo=where(y gt 0,moo)
  if allferrp[iw] eq 0 then moo=0
  if moo gt 0 then oplot,logT[oo],y[oo],line=1,color=boundindex
  ;
  y=reform(emp[*,iw]) & oo=where(y gt 0,moo)
  if allferrm[iw] eq 0 then moo=0
  if moo gt 0 then oplot,logT[oo],y[oo],line=1,color=3
endfor

end
