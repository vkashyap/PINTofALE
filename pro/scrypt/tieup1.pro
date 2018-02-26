;+
;script 	tieup1
;
;given an ID structure (IDSTR)
;	-- sets the line width of the first component to LSFWDT,
;	-- ties all the other widths to the first one,
;	-- freezes all but the first width parameter,
;	-- ties all line locations to that of the first one, but
;	   with a leeway defined by LINERR, and
;	-- if variable POSITIVF is set, then constrains all fluxes
;	   to be +ve
;
;outputs are
;	TIES	a string array containing all the ties
;	THAW	integer array signifying frozen(0) and thawed(1) parameters
;	POS	a float array containing the initial positions
;	WDT	a float array containing the initial widths
;	HGT	a float array containing the initial heights
;
;then, if
;	SPECT(LAMBDA) is well defined, goes ahead and calls FITLINES
;	and makes a new copy of IDSTR -- IDSTR2 -- with the new fluxes.
;
;the fit parameters are in
;	FPAR	a structure containing everything of interest
;	POS,WDT,FLX	best-fit values
;	[PWF]ERR[UL]	upper and lower deviations on best-fit values
;	[PWF]ERRC	delta-chisq threshold for which deviations are computed
;
;vinay kashyap (MarMM)
;	added EPITHET to call (VK; MMaug)
;	now LINERR can be 0 (VK; MVIMaug)
;-

;	absolutely need the following
if n_tags(idstr) eq 0 then message,'IDSTR not available'

;	this is where we start from
tmp=cat_id(idstr)

;	initialize some useful parameters
if not keyword_set(lsfwdt) then lsfwdt=0.015	;width of line
if n_elements(linerr) eq 0 then linerr=0.01	;uncertainty in line positions
if not keyword_set(positivf) then positivf=0	;flux positivity?
print,'' & help,lsfwdt,linerr,positivf

nf=n_elements(idstr.WVL)	;# features
nc=0L & for i=0L,nf-1L do nc=nc+n_elements(idstr.(i+1).WVL)	;# components

;	extract the parameters
pos=fltarr(nc) & wdt=pos+float(lsfwdt[0]) & hgt=pos
kc=0L
for kf=0L,nf-1L do begin
  idf=idstr.(kf+1L)
  wvlc=abs(idf.WVL) & flxc=idf.FLUX & mc=n_elements(wvlc)
  pos[kc]=wvlc[*] & hgt[kc]=flxc[*] & kc=kc+mc
endfor

;	freeze all widths except the first one
thaw=intarr(3L*nc)+1 & ii=indgen(nc)
thaw[3L*ii+1]=0	& thaw[1]=1

;	tie up all the widths to the first one
for i=1L,nc-1L do begin
  ic=strtrim(3L*i+1,2) & ac='a'+ic
  ;	aI=a1
  cc=ac+' = a1'
  if i eq 1 then ties=[cc] else ties=[ties,cc]
endfor

;	(loosely) tie up the relative positions
for i=1L,nc-1L do begin
  ic=strtrim(3L*i,2) & ac='a'+ic & dpos=pos[i]-pos[0]
  if linerr gt 0 then begin
    ;	ac = (ac > (a0+(dpos-linerr))) < (a0+(dpos+linerr))
    cc=ac+' = ('+ac+' > (a0+('+strtrim(dpos-linerr[0],2)+'))) < (a0+('+$
	strtrim(dpos+linerr[0],2)+'))'
  endif else begin
    ;	ac = a0 + dpos
    cc=ac+' = a0+'+strtrim(dpos,2)
  endelse
  if not keyword_set(ties) then ties=[cc] else ties=[ties,cc]
endfor

;	flux positivity constraint
if keyword_set(POSITIVF) then begin
  for i=0L,nc-1L do begin
    ic=strtrim(3L*i+2,2) & ac='a'+ic
    ;	aI=aI > 0
    cc=ac+' = '+ac+' > 0'
    if not keyword_set(ties) then ties=[cc] else ties=[ties,cc]
  endfor
endif

;	make labels for the IDs
epithet=idlabel(idstr,Zbeda=Zbeda,Ibeda=Ibeda,Wbeda=Wbeda,Lbeku=Lbeku,$
	wform=wform,wstyle=wstyle, ziform=ziform)

;	report
help,ties,thaw,pos,wdt,hgt

;	ready to call FITLINES?
nl=n_elements(lambda) & ns=n_elements(spect)
if nl eq 0 then stop,'LAMBDA undefined'
if ns eq 0 then stop,'SPECT(LAMBDA) undefined'
if nl ne ns then stop,'LAMBDA and SPECT(LAMBDA) incompatible'

;	so some initialization for FITLINES
nse=n_elements(sperr) & if nse ne ns then sperr=sqrt(abs(SPECT)+0.75)+1.
if not keyword_set(funcnam) then funcnam='x3model'
if not keyword_set(functyp) then functyp='lorentz'
if n_elements(dumb) eq 0 then dumb=1
if n_elements(verbose) eq 0 then verbose=5

;	call FITLINES
print,'*****************************************'
print,'calling FITLINES with following as input:'
print,'*****************************************'
help,lambda,spect,sperr,funcnam,functyp,pos,wdt,hgt,thaw,ties,$
	conlev,consig,dumb,verbose
print,'*****************************************'
print,string("7b) & print,'OK?  type any key to continue, n, q, or z to stop'
print,'*****************************************'
c=get_kbrd(1) & c=strlowcase(c)
if c eq 'n' or c eq 'q' or c eq 'z' then stop,'HALTING: type .CON to continue'
print,'*****************************************'
;
fpar=fitlines(lambda,spect,ysig=sperr,funcs=funcnam,pos=pos,wdt=wdt,flx=hgt,$
	perrp=perrp,perrm=perrm,perrc=perrc,$
	werrp=werrp,werrm=werrm,werrc=werrc,$
	ferrp=ferrp,ferrm=ferrm,ferrc=ferrc,$
	ties=ties,thaw=thaw,type=functyp,epithet=epithet,$
	conlev=conlev,consig=consig,dumb=dumb,verbose=verbose)

;	overwrite the fluxes in IDSTR
flx=fpar.FLX & flxp=fpar.FERRP & flxm=fpar.FERRM & flxc=fpar.FERRC
flxe=0.5*(flxp+flxm)
oo=where(flxc eq 2.7,moo) & if moo gt 0 then flxe[oo]=0.5*flxe[oo]
idstr2=updatid(idstr,flx,flxerr=flxe)

;	report
help,cat_id(idstr2),idstr2

end
