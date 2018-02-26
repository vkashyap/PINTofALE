pro oplot_bowtie,tgr,tgd,caldb=caldb,bowtie=bowtie,foracis=foracis,$
	hetg=hetg,srccol=srccol,bkgcol=bkgcol,npoint=npoint, _extra=e
;+
;procedure	oplot_bowtie
;	display a scatter plot of LETG dispersed photons in (TG_R,TG_D)
;	space and overplot the source and background extraction regions
;	on it.
;
;syntax
;	oplot_bowtie,tgr,tgd,caldb=caldb,bowtie=bowtie,/foracis,/hetg,$
;	srccol=srccol,bkgcol=bkgcol,npoint=npoint,$
;	and PLOT/OPLOT keywords
;
;parameters
;	tgr	[INPUT; required] TG_R coordinates from the evt2 file
;	tgd	[INPUT; required] TG_D coordinates from the evt2 file
;
;keywords
;	caldb	[INPUT] location of the Chandra CALDB
;		* default is /soft/ciao/CALDB
;		* the bowtie region file is assumed to in the directory
;		  $CALDB/data/chandra/hrc/bcf/tgmask2/
;		  and named letgD1999-07-22regN0002.fits
;	bowtie	[I/O] a structure that contains the bowtie read in
;		from the bowtie region file using MRDFITS()
;		* if defined on input, will not use CALDB and will
;		  not read it in
;	foracis	[INPUT] if set, ignores both the CALDB and BOWTIE
;		keywords and instead just plots the source and
;		background regions defined in tgextract, scaled
;		by FORACIS (e.g.,
;		FORACIS*tgextract.min_tg_d, FORACIS*tgextract.max_tg_d
;		for the source region, etc.)
;	hetg	[INPUT] if set, the default min_tg_d and max_tg_d are
;		set to values appropriate for the HETG.
;	srccol	[INPUT] color index to plot the source region in
;		* default is 3 (corresponds to green in PEASECOLR)
;	bkgcol	[INPUT] color index to plot the source region in
;		* default is 2 (corresponds to red in PEASECOLR)
;	npoint	[INPUT] plots only NPOINT points chosen randomly
;		* if -ve, plots the _first_ abs(NPOINT) points
;	_extra	[INPUT] use this to pass defined keywords to PLOT
;
;subroutines
;	STAMPLE
;
;history
;	vinay kashyap (Oct04)
;	added keyword HETG (VK; Nov04)
;-

;	usage
ok='ok' & np=n_params() & nr=n_elements(tgr) & nd=n_elements(tgd)
if np lt 2 then ok='Insufficient parametes' else $
 if nr eq 0 then ok='TG_R is undefined' else $
  if nd eq 0 then ok='TG_D is undefined' else $
   if nr ne nd then ok='TG_R and TG_D are incompatible'
if ok ne 'ok' then begin
  print,'Usage: oplot_bowtie,tgr,tgd,caldb=caldb,bowtie=bowtie,/foracis,/hetg,$'
  print,'       srccol=srccol,bkgcol=bkgcol,npoint=npoint,$'
  print,'       and PLOT keywords'
  print,'  display a scatter plot of LETG dispersed photons and'
  print,'  overplot the bowtie extraction region on it'
  if np ne 0 then message,ok,/informational
  return
endif

;	read in bowtie region if necessary
if keyword_set(hetg) then begin
  ;	assume HETG is never used together with HRC
  if not keyword_set(foracis) then foracis=1.
endif
if not keyword_set(foracis) then begin
  nb=n_tags(bowtie)
  if nb eq 0 then begin
    tgcal='/soft/ciao/CALDB'
    if keyword_set(caldb) then tgcal=strtrim(caldb,2)
    tgcal=tgcal+'/data/chandra/hrc/bcf/tgmask2'
    bowtiefil=tgcal+'/letgD1999-07-22regN0002.fits'
    bowtie=mrdfits(bowtiefil,1,hbowtie)
  endif
endif

;	make the plot
xmin=min(tgr,max=xmax,/nan) & ymin=min(tgd,max=ymax,/nan)
oo=lindgen(nr)
if keyword_set(npoint) then begin
  if abs(npoint[0]) lt nr then begin
    if npoint[0] lt 0 then oo=oo[0:-npoint[0]] else $
      oo=long(randomu(seed,npoint[0])*nr)
  endif
endif
plot,tgr[oo],tgd[oo],psym=3,xrange=[xmin,xmax],yrange=[ymin,ymax],$
	xtitle='TG!dR!n',ytitle='TG!dD!n', _extra=e

;	plot the source region
if not keyword_set(foracis) then begin
  ;	with bowtie
  scol=3 & if keyword_set(srccol) then scol=fix(srccol[0])
  oplot,(bowtie.tg_r)[*,3],(bowtie.tg_d)[*,3],col=scol, _extra=e
  scol=3 & if keyword_set(srccol) then scol=fix(srccol[[1]])
  oplot,(bowtie.tg_r)[*,0],(bowtie.tg_d)[*,0],col=scol, _extra=e
endif else begin
  ;	without bowtie
  max_tg_d=1.33e-3 & if keyword_set(hetg) then max_tg_d=6.6E-04
  scol=3 & if keyword_set(srccol) then scol=fix(srccol[0])
  oplot,[-1000,1000],max_tg_d*[1,1],col=scol, _extra=e
  min_tg_d=-1.33e-3 & if keyword_set(hetg) then min_tg_d=-6.6E-04
  scol=3 & if keyword_set(srccol) then scol=fix(srccol[[1]])
  oplot,[-1000,1000],min_tg_d*[1,1],col=scol, _extra=e
endelse

;	plot the background region
if not keyword_set(foracis) then begin
  ;	with bowtie
  bcol=2 & if keyword_set(bkgcol) then bcol=fix(bkgcol[0])
  oplot,(bowtie.tg_r)[*,4],(bowtie.tg_d)[*,4],col=bcol, _extra=e
  oplot,(bowtie.tg_r)[*,5],(bowtie.tg_d)[*,5],col=bcol, _extra=e
  bcol=2 & if keyword_set(bkgcol) then bcol=fix(bkgcol[[1]])
  oplot,(bowtie.tg_r)[*,1],(bowtie.tg_d)[*,1],col=bcol, _extra=e
  oplot,(bowtie.tg_r)[*,2],(bowtie.tg_d)[*,2],col=bcol, _extra=e
endif else begin
  ;	without bowtie
  max_upbkg_tg_d=1.33e-2 & if keyword_set(hetg) then max_upbkg_tg_d=6e-3
  bcol=2 & if keyword_set(bkgcol) then bcol=fix(bkgcol[0])
  oplot,[-1000,1000],max_tg_d*[1,1],col=bcol,line=2, _extra=e
  oplot,[-1000,1000],max_upbkg_tg_d*[1,1],col=bcol, _extra=e
  min_downbkg_tg_d=-1.33e-2 & if keyword_set(hetg) then min_downbkg_tg_d=-6e-3
  bcol=2 & if keyword_set(bkgcol) then bcol=fix(bkgcol[[1]])
  oplot,[-1000,1000],min_tg_d*[1,1],col=bcol,line=2, _extra=e
  oplot,[-1000,1000],min_downbkg_tg_d*[1,1],col=bcol, _extra=e
endelse

;	stamp with authority
stample

return
end
