;+
;WRAP_ID2FIT
;	program that calls ID_TO_FITPAR and FITLINES in sequence
;
;history
;	vinay kashyap (Nov98)
;	added call to UPDATID (VK; Dec98)
;	improved PICK handling, included keyword LSFWDT
;	  (VK; Mar99)
;-

;	input for ID_TO_FITPAR
nid=n_tags(idstr) & if nid eq 0 then stop,'IDSTR (for ID_TO_FITPAR) undefined'
npk=n_elements(pick) & if npk eq 0 or npk ge nid then pick=lindgen(nid-1L)
id1str=cat_id(idstr,pick=pick)
if not keyword_set(shackle) then begin
  shackle=0 & message,'Multiple IDs will not be shackled',/info
endif else help,shackle
if not keyword_set(leeway) then begin
  leeway=0 & message,'Parameters are completely unconstrained',/info
endif else help,leeway
if not keyword_set(parPerr) then parPerr=0.1
if not keyword_set(parWerr) then parWerr=0.01
if not keyword_set(parFerr) then parFerr=0.1
if not keyword_set(lsfwdt) then lsfwdt=parPerr

;	call ID_TO_FITPAR
id_to_fitpar,id1str,pars,ties,thaw,pos=parpos,wdt=parwdt,flx=parflx,$
	lsfwdt=lsfwdt,shackle=shackle,leeway=leeway,$
	perr=parPerr,werr=parWerr,ferr=parFerr

;parflx=parflx*randomu(seed)
;parwdt(*)=instr_wdt

;	input for FITLINES
nlam=n_elements(lambda) & if nlam eq 0 then stop,'LAMBDA (for FITLINES) not given'
nspc=n_elements(spect) & if nspc ne nlam then stop,'SPECT(LAMBDA) not given'
nspe=n_elements(sperr) & if nspe ne nspc then message,'Default errors on SPECT assumed',/info
if not keyword_set(funcnam) then funcnam='x3model'	;use X3MODEL
if not keyword_set(functyp) then functyp='gauss'	;Lorentz also possible
if not keyword_set(betha) then betha=1			;for MK_LORENTZ
if not keyword_set(dumb) then dumb=1			;for FIT_LEVMAR
poserr=0 & wdterr=0 & flxerr=0

;	call FITLINES
fitpar=fitlines(lambda,spect,ysig=sperr,funcs=funcnam,pos=parpos,perrp=poserr,$
	wdt=parwdt,werrp=wdterr,flx=parflx,ferrp=flxerr,thaw=thaw,ties=ties,$
	type=functyp,conlev=conlev,consig=consig, betha=betha,dumb=dumb)

;	print out result
ncomp=n_elements(fitpar.POS)
print,'Component	Position +-Err	Width +-Err		Flux +-Err'
form='(i4,"	",f10.4," +- ",f8.5,f8.4," +- ",f7.4,g12.4," +- ",g12.4)'
for i=0,ncomp-1 do print,i,fitpar.POS(i),fitpar.PERRP(i),$
	fitpar.WDT(i),fitpar.WERRP(i),fitpar.FLX(i),fitpar.FERRP(i),form=form

;	overwrite the ID fluxes
c1='Overwrite the fluxes into the ID list? [y/n]' & print,c1
c1=get_kbrd(1)
if strlowcase(c1) ne 'n' then begin
  if ncomp eq n_elements(parpos) then begin
    if not keyword_set(outfil) then outfil='/tmp/id2fit.lst'
    c2='' & read,prompt='filename to dump output in ['+outfil+'] ',c2
    if strtrim(c2,2) ne '' and strlowcase(c2) ne 'n' then outfil=strtrim(c2,2)
    id2str=updatid(id1str,fitpar.FLX)		;overwrite fluxes
    ;
    if strlowcase(c2) ne 'n' then help,$	;make hardcopy
	cat_id(id2str,flst=outfil,fluxerr=fitpar.FERRP)
    ;
    id3str=cat_id(idstr,id2str,ask='r')		;replace old with new
  endif else print,'Number of components does not match number of IDs'
endif

end
