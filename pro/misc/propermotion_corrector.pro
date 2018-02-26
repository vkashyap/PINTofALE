pro propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,$
	catyr=catyr,newyr=newyr,$
	i24=i24,osxg=osxg,$
	verbose=verbose, _extra=e
;+
;procedure	propermotion_corrector
;	computes proper motion and new coordinates
;
;syntax
;	propermotion_corrector.catra,catdec,pmra,pmdec,newra,newdec,$
;	catyr=catyr,newyr=newyr,verbose=verbose,$
;	/i24,/osxg,/hrs,/idnum,/idpre,/idform,/pm,/chop,/osep
;
;parameters (all REQUIRED)
;	catra	[INPUT] catalog RA
;	catdec	[INPUT] catalog Dec
;		* CATRA and CATDEC must have the same sizes
;		* may be either float or string array
;		* if float, assumed to be decimal degrees unless /I24 is set
;		* if string, assumed to be decimal degrees unless ':' or ' '
;		  are found
;	pmra	[INPUT] proper motion in RA [mas/yr]
;	pmdec	[INPUT] proper motion in Dec [mas/yr]
;		* PMRA,PMDEC must be scalars,
;		  or have the same size as CATRA,CATDEC
;	newra	[OUTPUT] new RA
;	newdec	[OUTPUT] new Dec
;		* newRA, newDEC will have the same size as CATRA,CATDEC
;		* default is to return decimal degrees
;		* set keyword /OSXG to return the output in sexagesimal
;
;keywords
;	catyr	[INPUT] epoch of CATRA,CATDEC
;		* default is 2000.0
;	newyr	[INPUT] epoch for which to compute new coordinates
;		* default is today
;	verbose	[INPUT] controls chatter
;	i24	[INPUT] caught and passed on to SKYPOS's first run
;		* use iff input are _not_ in decimal degrees
;	osxg	[INPUT] caught and passed on to SKYPOS's second run
;		* use iff you want the output in sexagesimal
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		SKYPOS: HRS, IDNUM, IDPRE, IDFORM, PM, CHOP, OSEP
;
;requires
;	SKYPOS [DEG2HRS [STR_2_ARR, SYZE]]
;
;history
;	vinay kashyap (2014sep05)
;-

;	usage
ok='ok' & np=n_params() & nra=n_elements(catra) & ndec=n_elements(catdec)
mra=n_elements(pmra) & mdec=n_elements(pmdec)
if np lt 6 then ok='Insufficient parameters' else $
 if nra eq 0 then ok='CATRA is undefined' else $
  if ndec eq 0 then ok='CATDEC is undefined' else $
   if nra ne ndec then ok='CATRA and CATDEC are not compatible' else $
    if mra eq 0 then ok='PMRA is undefined' else $
     if mdec eq 0 then ok='PMDEC is undefined' else $
      if mra ne ndec then ok='PMRA and PMDEC are not compatible' else $
       if mra ne 1 or mra ne nra then ok='PM* should be scalar or match CAT*'
if ok ne 'ok' then begin
  print,'Usage: propermotion_corrector.catra,catdec,pmra_masperyr,pmdec_masperyr,newra,newdec,$'
  print,'       catyr=catyr,newyr=newyr,verbose=verbose,$'
  print,'       /i24,/osxg,/hrs,/idnum,/idpre,/idform,/pm,/chop,/osep'
  print,'  computes new coordinates including proper motion'
  if np ne 0 then message,ok,/informational
  return
endif

;	regularize inputs
rapm=dblarr(nra) & decpm=rapm
if nra eq 1 then begin & rapm[0]=pmra[0] & decpm[0]=pmdec[0] & endif
if nra gt 1 and mra eq 1 then begin & rapm[*]=pmra[0] & decpm[*]=pmdec[0] & endif
if nra gt 1 and mra gt 1 then begin & rapm=pmra & decpm=pmdec & endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
catepoch=fltarr(nra)+2000.0
if keyword_set(catyr) then begin
  ncyr=n_elements(catyr)
  if ncyr eq 1 then catepoch[*]=float(catyr[0])
  if ncyr gt 1 then catepoch[0:(nra<ncyr)-1]=float(catyr[0:(ncyr<nra)-1])
endif
caldat,julday(),mon,day,yr,hr,min,sec
curepoch=2000.+(julday()-julday(1,1,2000))/365.256360417
newepoch=fltarr(nra)+curepoch
if keyword_set(newyr) then begin
  nnyr=n_elements(newyr)
  if nnyr eq 1 then newepoch[*]=float(newyr[0])
  if nnyr gt 1 then newepoch[0:(nra<nnyr)-1]=float(newyr[0:(nnyr<nra)-1])
endif

;	do each row separately
for i=0L,nra-1L do begin	;{I=0,NRA-1

  zira=catra[i] & zidec=catdec[i]
  skypos,zira,zidec,zora,zodec,i24=i24
  if vv gt 0 then print,zira,' ',zidec,zora[0],zodec[0]
  ;
  dyr=newepoch[i]-catepoch[i]
  if vv gt 0 then print,catepoch[i],newepoch[i],dyr
  cosdec=cos(zodec[0]*!pi/180.)
  zpra=rapm[i] & zpdec=decpm[i]
  dra=(zpra*1e-3/3600.)*dyr	;[deg/yr]
  ddec=(zpdec*1e-3/3600.)*dyr	;[deg/yr]
  yra=(zora[0]*cosdec+dra)/cosdec
  ydec=zodec[0]+ddec
  ;
  skypos,yra,ydec,outra,outdec,osxg=osxg, _extra=e
  if vv gt 0 then print,yra,ydec,' ',outra[0],' ',outdec[0]

  if n_elements(xra) eq 0 then begin
    xra=[outra[0]] & xdec=[outdec[0]]
  endif else begin
    xra=[xra,outra[0]] & xdec=[xdec,outdec[0]]
  endelse

endfor				;I=0,NRA-1}

if vv gt 1000 then stop,'PROPERMOTION_CORRECTOR: type .CON to continue'

;	outputs
newra=xra & newdec=xdec

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	example
;	Prox Cen has ICRS2000.0 coordinates 14:29:42.94853,-62:40:46.1631
;	and proper motion -3775.75,765.34 mas/yr
;
;	this example calls propermotion_corrector for 4 different epochs,
;	2000.0, 2012.5, now, and 2015.5

catra='14:29:42.94853' & catdec='-62:40:46.1631'
pmra=-3775.75 & pmdec=765.34

print,'' & help,catra,catdec,pmra,pmdec & print,''

;	this is for the epoch 2000.0, should produce no corrections
propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,newyr=2000.0,/i24,/osxg
print,'propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,newyr=2000.0,/i24,/osxg'
print,'OLD 2000.0: ',catra,' ',catdec
print,'NEW 2000.0: ',newra[0],' ',newdec[0],' (should show no change)'
print,''

;	this is for the epoch 2012.5
propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,newyr=2012.5,/i24,/osxg
print,'propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,newyr=2012.5,/i24,/osxg'
print,'OLD 2000.0: ',catra,' ',catdec
print,'NEW 2012.5: ',newra[0],' ',newdec[0]
print,'    2012.5: 14:29:36.09 -62:40:36.59 (expected)'
print,''

;	starting from 2012.5 and going back to 2000.0
propermotion_corrector,newra,newdec,pmra,pmdec,altra,altdec,catyr=2012.5,newyr=2000.0,/i24,/osxg
print,'propermotion_corrector,newra,newdec,pmra,pmdec,altra,altdec,catyr=2012.5,newyr=2000.0,/i24,/osxg'
print,'OLD 2012.5: ',newra[0],' ',newdec[0]
print,'NEW 2000.0: ',altra[0],' ',altdec[0]
print,'    2000.0: ',catra[0],' ',catdec[0],' (actual)'
print,''

;	this is for *now*, output in decimal degrees
propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,/i24
print,'propermotion_corrector,catra,catdec,pmra,pmdec,newra,newdec,/i24'
print,'OLD 2000.0: ',catra,' ',catdec
print,'NEW NOW: ',newra[0],' ',newdec[0],form='(a,f12.7,a1,f12.7)'
print,''

;	now show the calling sequence
print,'propermotion_corrector'
propermotion_corrector

end
