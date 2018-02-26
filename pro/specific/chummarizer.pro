function chummarizer,counts,wave,effar=effar,wvlar=wvlar,lsffwhm=lsffwhm,$
	type=type,hawastr=hawastr,kkcont=kkcont,ctline=ctline,wscale=wscale,$
	wcont=wcont,lpos=lpos,lerrpos=lerrpos,lflx=lflx,lerrflx=lerrflx,$
	lwdt=lwdt,lerrwdt=lerrwdt,verbose=verbose, _extra=e
;+
;function	chummarizer
;	analyzes a high-resolution Chandra grating spectrum and
;	extracts some meaningful summarizing quantities from it,
;	such as the ratio of lines-to-continuum fluxes, etc., and
;	returns the results in a structure
;
;syntax
;	chstr=chummarizer(counts,wave,effar=effar,wvlar=wvlar,$
;	lsffwhm=lsffwhm,type=type,hawastr=hawastr,kkcont=kkcont,$
;	ctline=ctline,wscale=wscale,wcont=wcont,lpos=lpos,lerrpos=lerrpos,$
;	lflx=lflx,lerrflx=lerrflx,lwdt=lwdt,lerrwdt=lerrwdt,$
;	verbose=verbose,maxkern=maxkern,clev=clev,maxiter=maxiter)
;
;parameters
;	counts	[INPUT; required] counts spectrum
;	wave	[INPUT; required] wavelength grid for the spectrum
;		* if size matches that of COUNTS, taken to be mid-bin values
;		* if size exceeds COUNTS by 1, taken to be grid boundaries
;		* if size is double that of COUNTS, taken to be two
;		  concatenated arrays of lower and upper boundaries of grid
;		* it is assumed that WAVE is on a regular grid, i.e., that
;		  WAVE[1:*]-WAVE is constant
;
;keywords
;	effar	[INPUT] effective area [cm^2]
;		* if set, the counts are converted to flux prior to
;		  summarizing
;	wvlar	[INPUT] wavelengths at which effective area is defined
;		* if not set, assumed to be WAVE
;		* size MUST match EFFAR, else both are ignored
;		** currently, both EFFAR and WVLAR are ignored
;	lsffwhm	[I/O] passed w/o check to HAWALINER()
;	type	[INPUT] string denoting which model function to use
;		in fitting the detected lines (see LIBMODEL)
;		* default: 'beta=2.5' (like Chandra grating LSFs)
;	hawastr	[OUTPUT] output from HAWALINER()
;	kkcont	[OUTPUT] output from HAWALINER()
;	ctline	[OUTPUT] output from HAWALINER()
;	wcont	[OUTPUT] output from HAWALINER()
;	wscale	[I/O] pass defined wavelet scales to HAWALINER()
;	lpos	[OUTPUT] output from HAWALINER()
;	lerrpos	[OUTPUT] output from HAWALINER()
;	lflx	[OUTPUT] output from HAWALINER()
;	lerrflx	[OUTPUT] output from HAWALINER()
;	lwdt	[OUTPUT] output from HAWALINER()
;	lerrwdt	[OUTPUT] output from HAWALINER()
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		* HAWALINER: MAXKERN, CLEV, MAXITER
;
;restrictions
;	backgrounds are ignored
;	apply only to MEG
;	no error bars yet
;
;subroutines
;	HAWALINER()
;	HIPD_INTERVAL()
;	MID2BOUND()
;	GETLOCMAX()
;
;history
;	vinay kashyap (May07)
;	multiple bug corrections (VK; Jun07)
;	modified calling sequence to HAWALINER, changed behavior of
;	  how counts are collected (VK; Jul08)
;-

;	usage
ok='ok' & & np=n_params() & nc=n_elements(counts) & nw=n_elements(wave)
if np lt 2 then ok='Insufficient parameters' else $
 if nc eq 0 then ok='COUNTS array undefined' else $
  if nw eq 0 then ok='WAVE array undefined' else $
   if nw ne nc and nw ne nc+1L and nw ne 2L*nc then $
   	ok='COUNTS and WAVE arrays incompatible' else $
    if nc lt 9 then ok='arrays too small to bother'
if ok ne 'ok' then begin
  print,'Usage: chstr=chummarizer(counts,wave,effar=effar,wvlar=wvlar,$'
  print,'       lsffwhm=lsffwhm,type=type,hawastr=hawastr,kkcont=kkcont,$'
  print,'       ctline=ctline,wscale=wscale,wcont=wcont,lpos=lpos,lerrpos=lerrpos,$'
  print,'       lflx=lflx,lerrflx=lerrflx,lwdt=lwdt,lerrwdt=lerrwdt,$'
  print,'       verbose=verbose,maxkern=maxkern,clev=clev,maxiter=maxiter)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
ct=float(counts) & wvl=findgen(nc) & wgrid=findgen(nc+1L)-0.5
if nw eq nc then begin	;(WAVE are mid-bin values
  wvl=float(wave)
  wgrid=mid2bound(wvl, _extra=e)
endif			;NW==NC)
if nw eq nc+1L then begin	;(WAVE are grid boundaries
  wgrid=float(wave)
  wvl=0.5*(wgrid[1:*]+wgrid)
endif				;NW==NC+1)
if nw eq 2L*nc then begin	;(WAVE are BIN_LO,BIN_HI arrays
  w1=wave[0L:nc-1L] & w2=wave[nc:*]
  if w1[0] lt w2[0] then wgrid=[w1[0],w2] else wgrid=[w1,w2[0]]
  wvl=0.5*(w1[0L:nc-1L]+w2[nc:*])
endif				;NW==2*NC)

;	keywords
;VERBOSE
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;EFFAR(WVLAR)
arwvl=wvl & areff=0.*wvl+1.
ok='ok' & nea=n_elements(effar) & nwa=n_elements(wvlar)
if nea eq 0 then ok='EFFAR not defined' else $
 if nwa eq 0 then ok='WVLAR not defined' else $
  if nea ne nwa then ok='EFFAR and WVLAR are incompatible'
if ok eq 'ok' then areff=(interpol(effar,wvlar,arwvl)>0)<(max(effar)) else $
 if nea gt 0 then message,ok,/informational
lsftype_def='beta=2.5'
if not keyword_set(type) then lsftype=lsftype_def else $
 if size(type,/type) ne 7 then lsftype=lsftype_def else $
  if n_elements(type) gt 1 then lsftype=type[0] else lsftype=type[0]

;	call HAWALINER
linewvl=hawaliner(ct,wvl,lsffwhm=lsffwhm,type=lsftype,hawastr=hawastr,$
	kkcont=kkcont,ctline=ctline,wscale=wscale,$
	wcont=wcont,ewcont=ewcont,lpos=lpos,lerrpos=lerrpos,$
	lflx=lflx,lerrflx=lerrflx,$
	lwdt=lwdt,lerrwdt=lerrwdt,$
	wrange=[min(wgrid),max(wgrid)],/fullerr,verbose=vv, _extra=e)
xgrid=hawastr.X
yct=hawastr.Y
hwcont=hawastr.YCONT
ehwcont=hawastr.YCONTERR

;	we are interested in these:
;	Ne X 12.13
;	O VIII 18.97
;	O VIIrif 21.6,21.8,22.1
;	Mg XII 8.419
;	Si XIV 6.18
;	S XV 5.1
;	Fe XVII 15.01,17.1
;	all lines shortward and longward of Ne X
;	and the continnum in these ranges:
;	3-6, 6-12, 13-25

;	get the counts in all the locations
wNe10=12.13 & wO8=18.97 & wO7r=21.60 & wO7i=21.80 & wO7f=22.10
wMg12=8.419 & wSi14=6.18 & wS15=5.1 & wFe17_15=15.01 & wFe17_17=17.1

;
oNe10=where(abs(wNe10-hawastr.LINEPOS) lt lsffwhm,moNe10) & fNe10=0. & efNe10=1e10
if moNe10 ne 0 then begin
  fNe10=(hawastr.LINEFLX)[oNe10[0]] & efNe10=(hawastr.LINEFLXERR)[oNe10[0]]
endif
if fNe10 eq 0 then begin
  oNe10=where(xgrid ge wNe10-2*lsffwhm and xgrid le wNe10+2*lsffwhm,moNe10)
  if moNe10 gt 0 then begin
    fNe10=(total(yct[oNe10])-total(hwcont[oNe10]))>0
    efNe10=sqrt(total(yct[oNe10])+total(ehwcont[oNe10]^2))
  endif
endif

;
;oO8=where(xgrid ge wO8-2*lsffwhm and xgrid le wO8+2*lsffwhm,moO8)
;
oO8=where(abs(wO8-hawastr.LINEPOS) lt lsffwhm,moO8) & fO8=0. & efO8=1e10
if moO8 ne 0 then begin
  fO8=(hawastr.LINEFLX)[oO8[0]] & efO8=(hawastr.LINEFLXERR)[oO8[0]]
endif
if fO8 le 0 then begin
  oO8=where(xgrid ge wO8-2*lsffwhm and xgrid le wO8+2*lsffwhm,moO8)
  if moO8 gt 0 then begin
    fO8=(total(yct[oO8])-total(hwcont[oO8]))>0
    efO8=sqrt(total(yct[oO8])+total(ehwcont[oO8]^2))
  endif
endif

;
oO7r=where(abs(wO7r-hawastr.LINEPOS) lt lsffwhm,moO7r) & fO7r=0. & efO7r=1e10
if moO7r ne 0 then begin
  fO7r=(hawastr.LINEFLX)[oO7r[0]] & efO7r=(hawastr.LINEFLXERR)[oO7r[0]]
endif
if fO7r le 0 then begin
  oO7r=where(xgrid ge wO7r-2*lsffwhm and xgrid le wO7r+2*lsffwhm,moO7r)
  if moO7r gt 0 then begin
    fO7r=(total(yct[oO7r])-total(hwcont[oO7r]))>0
    efO7r=sqrt(total(yct[oO7r])+total(ehwcont[oO7r]^2))
  endif
endif

;
oO7i=where(abs(wO7i-hawastr.LINEPOS) lt lsffwhm,moO7i) & fO7i=0. & efO7i=1e10
if moO7i ne 0 then begin
  fO7i=(hawastr.LINEFLX)[oO7i[0]] & efO7i=(hawastr.LINEFLXERR)[oO7i[0]]
endif
if fO7i le 0 then begin
  oO7i=where(xgrid ge wO7i-2*lsffwhm and xgrid le wO7i+2*lsffwhm,moO7i)
  if moO7i gt 0 then begin
    fO7i=(total(yct[oO7i])-total(hwcont[oO7i]))>0
    efO7i=sqrt(total(yct[oO7i])+total(ehwcont[oO7i]^2))
  endif
endif

;
oO7f=where(abs(wO7f-hawastr.LINEPOS) lt lsffwhm,moO7f) & fO7f=0. & efO7f=1e10
if moO7f ne 0 then begin
  fO7f=(hawastr.LINEFLX)[oO7f[0]] & efO7f=(hawastr.LINEFLXERR)[oO7f[0]]
endif
if fO7f le 0 then begin
  oO7f=where(xgrid ge wO7f-2*lsffwhm and xgrid le wO7f+2*lsffwhm,moO7f)
  if moO7f gt 0 then begin
    fO7f=(total(yct[oO7f])-total(hwcont[oO7f]))>0
    efO7f=sqrt(total(yct[oO7f])+total(ehwcont[oO7f]^2))
  endif
endif

;
oMg12=where(abs(wMg12-hawastr.LINEPOS) lt lsffwhm,moMg12) & fMg12=0. & efMg12=1e10
if moMg12 ne 0 then begin
  fMg12=(hawastr.LINEFLX)[oMg12[0]] & efMg12=(hawastr.LINEFLXERR)[oMg12[0]]
endif
if fMg12 le 0 then begin
  oMg12=where(xgrid ge wMg12-2*lsffwhm and xgrid le wMg12+2*lsffwhm,moMg12)
  if moMg12 gt 0 then begin
    fMg12=(total(yct[oMg12])-total(hwcont[oMg12]))>0
    efMg12=sqrt(total(yct[oMg12])+total(ehwcont[oMg12]^2))
  endif
endif

;
oSi14=where(abs(wSi14-hawastr.LINEPOS) lt lsffwhm,moSi14) & fSi14=0. & efSi14=1e10
if moSi14 ne 0 then begin
  fSi14=(hawastr.LINEFLX)[oSi14[0]] & efSi14=(hawastr.LINEFLXERR)[oSi14[0]]
endif
if fSi14 le 0 then begin
  oSi14=where(xgrid ge wSi14-2*lsffwhm and xgrid le wSi14+2*lsffwhm,moSi14)
  if moSi14 gt 0 then begin
    fSi14=(total(yct[oSi14])-total(hwcont[oSi14]))>0
    efSi14=sqrt(total(yct[oSi14])+total(ehwcont[oSi14]^2))
  endif
endif

;
oS15=where(abs(wS15-hawastr.LINEPOS) lt lsffwhm,moS15) & fS15=0. & efS15=1e10
if moS15 ne 0 then begin
  fS15=(hawastr.LINEFLX)[oS15[0]] & efS15=(hawastr.LINEFLXERR)[oS15[0]]
endif
if fS15 le 0 then begin
  oS15=where(xgrid ge wS15-2*lsffwhm and xgrid le wS15+2*lsffwhm,moS15)
  if moS15 gt 0 then begin
    fS15=(total(yct[oS15])-total(hwcont[oS15]))>0
    efS15=sqrt(total(yct[oS15])+total(ehwcont[oS15]^2))
  endif
endif

;
oFe17_15=where(abs(wFe17_15-hawastr.LINEPOS) lt lsffwhm,moFe17_15) & fFe17_15=0. & efFe17_15=1e10
if moFe17_15 ne 0 then begin
  fFe17_15=(hawastr.LINEFLX)[oFe17_15[0]] & efFe17_15=(hawastr.LINEFLXERR)[oFe17_15[0]]
endif
if fFe17_15 le 0 then begin
  oFe17_15=where(xgrid ge wFe17_15-2*lsffwhm and xgrid le wFe17_15+2*lsffwhm,moFe17_15)
  if moFe17_15 gt 0 then begin
    fFe17_15=(total(yct[oFe17_15])-total(hwcont[oFe17_15]))>0
    efFe17_15=sqrt(total(yct[oFe17_15])+total(ehwcont[oFe17_15]^2))
  endif
endif

;	for Fe 17 doublet, don't bother with trying to _match_ the wavelength
fFe17_17=0. & efFe17_17=1e10
oFe17_17=where(xgrid ge wFe17_17-5*lsffwhm and xgrid le wFe17_17+5*lsffwhm,moFe17_17)
if moFe17_17 gt 0 then begin
  fFe17_17=(total(yct[oFe17_17])-total(hwcont[oFe17_17]))>0
  efFe17_17=sqrt(total(yct[oFe17_17])+total(ehwcont[oFe17_17]^2))
endif

oc36=where(xgrid ge 3. and xgrid le 6,moc36)
if moc36 gt 0 then fc36=total(hwcont[oc36]) else fc36=0.
if moc36 gt 0 then efc36=sqrt(total(ehwcont[oc36]^2)) else efc36=1e10

oc612=where(xgrid ge 6. and xgrid le 12,moc612)
if moc612 gt 0 then fc612=total(hwcont[oc612]) else fc612=0.
if moc612 gt 0 then efc612=sqrt(total(ehwcont[oc612]^2)) else efc612=1e10

oc1325=where(xgrid ge 13. and xgrid le 25,moc1325)
if moc1325 gt 0 then fc1325=total(hwcont[oc1325]) else fc1325=0.
if moc1325 gt 0 then efc1325=sqrt(total(ehwcont[oc1325]^2)) else efc1325=1e10

oNes=where(hawastr.LINEPOS lt wNe10-2*lsffwhm,moNes)
fNes=0. & if moNes gt 0 then fNes=total((hawastr.LINEFLX)[oNes])
if moNes gt 0 then efNes=sqrt(total(((hawastr.LINEFLXERR)[oNes])^2)) else efNes=1e10

oNel=where(hawastr.LINEPOS gt wNe10+2*lsffwhm,moNel)
fNel=0. & if moNel gt 0 then fNel=total((hawastr.LINEFLX)[oNel])
if moNel gt 0 then efNel=sqrt(total(((hawastr.LINEFLXERR)[oNel])^2)) else efNel=1e10

;	now get some useful summarizing ratios
;behr=behr_hug(fNe10,fO8,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fNe10d=reform(var[0,*]) & fO8d=reform(var[1,*])
;behr=behr_hug(fO7i,fO7f,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fO7id=reform(var[0,*]) & fO7fd=reform(var[1,*])
;behr=behr_hug(fO7r,fMg12,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fO7rd=reform(var[0,*]) & fMg12d=reform(var[1,*])
;behr=behr_hug(fSi14,fS15,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fSi14d=reform(var[0,*]) & fS15d=reform(var[1,*])
;behr=behr_hug(fFe17_15,fFe17_17,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fFe17_15d=reform(var[0,*]) & fFe17_17d=reform(var[1,*])
;behr=behr_hug(fNes,fNel,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fNesd=reform(var[0,*]) & fNeld=reform(var[1,*])
;behr=behr_hug(fc36,fc612,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fc36d=reform(var[0,*]) & fc612d=reform(var[1,*])
;behr=behr_hug(fc1325,fNe10,0,0,1.,1.,1.,1.,algo='gibbs',nsim=10000L,nburnin=5000L,/outputMC,verbose=vv) & var=fltarr(2,5000)
;openr,umc,'BEHR_draws.txt',/get_lun & readf,umc,var & close,umc & free_lun,umc & fc1325d=reform(var[0,*]) & fNe10d=reform(var[1,*])

Ne2O=fNe10/fO8 & eNe2O=Ne2O*sqrt((efNe10/fNe10)^2+(efO8/fO8)^2)
Oif=fO7i/fO7f & eOif=Oif*sqrt((efO7i/fO7i)^2+(efO7f/fO7f)^2)
O78=fO7r/fO8 & eO78=O78*sqrt((efO7r/fO7r)^2+(efO8/fO8)^2)
O8Fe=fO8/(fFe17_15+fFe17_17) & efnum=sqrt(efFe17_15^2+efFe17_17^2) & eO8Fe=O8Fe*sqrt((efO8/fO8)^2+(efnum/(fFe17_15+fFe17_17))^2)
Mg2Fe=fMg12/(fFe17_15+fFe17_17) & efnum=sqrt(efFe17_15^2+efFe17_17^2) & eMg2Fe=Mg2Fe*sqrt((efMg12/fMg12)^2+(efnum/(fFe17_15+fFe17_17))^2)
Si2Fe=fSi14/(fFe17_15+fFe17_17) & efnum=sqrt(efFe17_15^2+efFe17_17^2) & eSi2Fe=Si2Fe*sqrt((efSi14/fSi14)^2+(efnum/(fFe17_15+fFe17_17))^2)
S2Fe=fS15/(fFe17_15+fFe17_17) & efnum=sqrt(efFe17_15^2+efFe17_17^2) & eS2Fe=S2Fe*sqrt((efS15/fS15)^2+(efnum/(fFe17_15+fFe17_17))^2)

sl2c=fNes/(fc36+fc612+fc1325) & efnum=sqrt(efc36^2+efc612^2+efc1325^2) & esl2c=sl2c*sqrt((efNes/fNes)^2+(efnum/(fc36+fc612+fc1325))^2)
ll2c=fNel/(fc36+fc612+fc1325) & efnum=sqrt(efc36^2+efc612^2+efc1325^2) & ell2c=ll2c*sqrt((efNel/fNel)^2+(efnum/(fc36+fc612+fc1325))^2)
l2c=(fNes+fNel+fNe10)/(fc36+fc612+fc1325) & efnum1=sqrt(efNes^2+efNel^2+efNe10^2) & efnum2=sqrt(efc36^2+efc612^2+efc1325^2) & el2c=l2c*sqrt((efnum1/(fNes+fNel+fNe10))^2+(efnum2/(fc36+fc612+fc1325))^2)
cl2cs=fc1325/(fc36+fc612) & efnum=sqrt(efc36^2+efc612^2) & ecl2cs=cl2cs*sqrt((efc1325/fc1325)^2+(efnum/(fc36+fc612))^2)
cl2cm=fc1325/fc612 & ecl2cm=cl2cm*sqrt((efc1325/fc1325)^2+(efc612/fc612)^2)
cm2cs=fc612/fc36 & ecm2cs=cm2cs*sqrt((efc612/fc612)^2+(efc36/fc36)^2)

;	make the output structure
chhlp=create_struct('Ne2O','Ne10/O8: tracks Ne/O',$
	'Oif','O7i/O7f: tracks density',$
	'O78','O7r/O8: tracks T',$
	'O8Fe','O8/Fe: tracks O abun',$
	'Mg2Fe','Mg12/Fe17: tracks Mg abun',$
	'Si2Fe','Si14/Fe17: tracks Si abun',$
	'S2Fe','S15/Fe17: tracks S abun',$
	'sl2c','lines(w<Ne10)/cont: line-to-continuum for lines shortward of Ne 10',$
	'll2c','lines(w>Ne10)/cont: line-to-continuum for lines longward of Ne 10',$
	'l2c','lines/cont: line-to-continuum ratio, tracks metallicity',$
	'cl2cs','continnum at short/long wvl: tracks T',$
	'cl2cm','continnum at short/medium wvl: tracks T',$
	'cm2cs','continnum at medium/long wvl: tracks T',$
	'HAWAstr','HAWALINER output')
chstr=create_struct('help',chhlp,$
	'fNe10',fNe10,'efNe10',efNe10,'fO8',fO8,'efO8',efO8,$
	'fO7r',fO7r,'efO7r',efO7r,'fO7i',fO7i,'efO7i',efO7i,$
	'fO7f',fO7f,'efO7f',efO7f,'fMg12',fMg12,'efMg12',efMg12,$
	'fSi14',fSi14,'efSi14',efSi14,'fS15',fS15,'efS15',efS15,$
	'fFe17_15',fFe17_15,'efFe17_15',efFe17_15,'fFe17_17',fFe17_17,'efFe17_17',efFe17_17,$
	'fc36',fc36,'efc36',efc36,'fc612',fc612,'efc612',efc612,$
	'fc1325',fc1325,'efc1325',efc1325,'fNes',fNes,'efNes',efNes,$
	'fNel',fNel,'efNel',efNel,'LSFFWHM',lsffwhm,$
	'Ne2O',Ne2O,'Oif',Oif,'O78',O78,'O8Fe',O8Fe,'Mg2Fe',Mg2Fe,'Si2Fe',Si2Fe,'S2Fe',S2Fe,$
	'sl2c',sl2c,'ll2c',ll2c,'l2c',l2c,'cl2cs',cl2cs,'cl2cm',cl2cm,'cm2cs',cm2cs,$
	'eNe2O',eNe2O,'eOif',eOif,'eO78',eO78,'eO8Fe',eO8Fe,'eMg2Fe',eMg2Fe,'eSi2Fe',eSi2Fe,'eS2Fe',eS2Fe,$
	'esl2c',esl2c,'ell2c',ell2c,'el2c',el2c,'ecl2cs',ecl2cs,'ecl2cm',ecl2cm,'ecm2cs',ecm2cs,$
	'HAWAstr',hawastr)

if vv gt 500 then stop,'HALTing; type .CON to continue'

return,chstr
end
