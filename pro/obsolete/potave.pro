function potave,atom,fx,wvl,mwvl,midx,sigma=sigma,ewt=ewt,tmax=tmax,$
	emavg=emavg,emsig=emsig,unpot=unpot, _extra=e
;+
;function	potave
;	returns Pottasch point EM [cm^-5] ([cm^-3] if WVLAR and
;	EFFAR are supplied) and corresponding standard deviation
;	and peak temperature by averaging the Pottasch corrected
;	EMs derived for each of the supplied line fluxes.
;
;	for each line, compute the pottasch corrected EM using FLX2EM
;	and POTTASCH.  average them all together and return the result.
;
;parameters
;	atom	[INPUT; required] specifies the origin of the lines
;		* if scalar, all lines are assumed to be from the
;		  same source
;		* if vector, must match the size of MWVL
;		* e.g., 'Fe XIII'
;	fx	[INPUT; required] observed line fluxes [ct/s]
;	wvl	[INPUT; required] wavelengths at which FX is observed [Ang]
;	mwvl	[INPUT; required] matched wavelengths in the database [Ang]
;		* N(MWVL) must be .GE. N(WVL)
;		* if N(MWVL) .GT. N(WVL), look in MIDX for how they match up
;	midx	[INPUT] array of position indices -- must be of same size
;		as MWVL, and points to matching WVL.
;		* required if N(MWVL) .GT. N(WVL)
;
;keywords
;	sigma	[INPUT; default: 1+sqrt(abs(FLX)+0.75)] error on FLX
;	ewt	[INPUT] if set, weights the averages by the inverse of SIGMA
;	tmax	[OUTPUT] returns the position of the peak in the line
;		contribution for each line.
;	emavg	[OUTPUT] the full averaged, Pottasch corrected emissivity
;		curve
;	emsig	[OUTPUT] the error on EMAVG
;	unpot	[INPUT] if set, does not do Pottasch corrections
;	_extra	[INPUT] allows setting defined keywords to
;		FLX2EM [ABUND, NH, SKIPRD, DEFEM]
;		RD_LINE [LOGP, N_E, DBDIR, ALLAH]
;		FOLD_IONEQ [CHIFIL [EQFILE, CHIDIR]]
;		LINEFLX [EFFAR, WVLAR]
;		ISMTAU [FH2, HE1, HEII, FANO]
;		POTTASCH [LEVEL]
;
;restrictions
;	* LINEID must've been run first, else where'll you get MWVL and MIDX?
;	* requires subroutines FLX2EM, RD_LINE, FOLD_IONEQ, GETABUND,
;	  LINEFLX, ISMTAU, POTTASCH, RD_IONEQ, READ_IONEQ, SYMB2ZION,
;	  LAT2ARAB, KILROY, WHEE
;
;history
;	vinay kashyap (Jun97)
;	added keywords EMAVG and EMSIG (VK; Jul97)
;	corrected normalization factors (bug; VK; Sep97)
;	added message saying routine is obsolete (VK; Oct98)
;-

message,'This routine is obsolete; use ID_TO_EMIS, FLUX_TO_EM, and POTTASCH!',/info

;	check input for problems
ok=''			;assume everything is OK
na=n_elements(atom) & nf=n_elements(fx) & nw=n_elements(wvl)
mw=n_elements(mwvl) & mi=n_elements(midx)
if (na eq 0) then ok='(Z,ION) not given' else $
 if (nf eq 0) then ok='Observed flux(es) missing' else $
  if (nw eq 0) then ok='Observed wavelength(s) missing' else $
   if (mw eq 0) then ok='Matching wavelengths missing' else $
    if (na gt 1) and (na ne mw) then ok='mismatched MWVL(Z,ION)' else $
     if (nf ne nw) then ok='mismatched FLX(WVL)' else $
      if (mw lt nw) then ok='mismatched matching' else $
       if (nw gt 1) and (mw gt nw) and (mw ne mi) then ok='confused matching'

;	usage
if ok ne '' then begin
  message,ok,/info & print,''
  print,'Usage: emavg=potave(atom,fx,wvl,mwvl,midx,sigma=sigma,ewt=ewt,$'
  print,'       tmax=tmax,emavg=emavg,emsig=emsig,unpot=unpot)'
  print,'  return [Pottasch corrected ion-averaged EMs, std. deviation]'
  print,'  also accepts defined keywords ABUND,NH,SKIPRD,DEFEM (FLX2EM);'
  print,'  LOGP,N_E,DBDIR,ALLAH (RD_LINE); CHIFIL,EQFILE,CHIDIR (FOLD_IONEQ);'
  print,'  EFFAR,WVLAR (LINEFLX); FH2,HE1,HEII,FANO (ISMTAU); LEVEL (POTTASCH)'
  return,-1L
endif

;	initialize
Z=intarr(mw) & ion=Z & mtch=Z
for iz=0,na-1 do begin		;{decode ATOM
  symb2zion,atom(iz),zz,ii
  if ii eq 0 then begin
    message,'Ionic state not given -- setting to I!',/info
    ii=1
  endif
  if na eq 1 then begin
    Z(*)=zz & ion(*)=ii		;all are of same type
  endif else begin
    Z(iz)=zz & ion(iz)=ii	;each may be different
  endelse
endfor				;IZ=0,NA-1}
if mi eq 0 then mtch(*)=lindgen(nw) else mtch=[midx]	;position indices
tmax=fltarr(nw)			;store location of peak contribution

;	weights and measures
sig=1.+sqrt(abs(fx)+0.75) & ns=n_elements(sigma)
if ns le nw then sig(0)=sigma(*) else sig(*)=sigma(0:nw-1)
wts=fltarr(nw)+1.

for iw=0,nw-1 do begin			;{Pottasch corrected EMs for each line
  lfx=0
  flux=fx(iw)			;line flux
  ow=where(mtch eq iw,mow)	;any matching wavelengths?
  if mow gt 0 then begin	;(if any matches exist for this line...
    wav=mwvl(ow)		;matching wavelengths
    zz=Z(ow)			;atomic numbers
    ii=ion(ow)			;ionic states
    em=flx2em(flux,wav,zz,ii,lfx,logT=logT, _extra=e)	;get EM(T)
    pot=pottasch(lfx,logT,Z=zz,ion=ii, _extra=e)	;get correction factors
    if iw eq 0 then begin
      nt=n_elements(logT)
      avgem=dblarr(nt,nw)	;store Pottasch corrected EMs
    endif
    emnorm=intarr(nt)		;normalization for avgem
    for jw=0,mow-1 do begin	;{do the correction
      oo=where(em(*,jw) gt 0,moo) & count=0
      if moo gt 0 then begin
	if not keyword_set(unpot) then em(oo,jw)=em(oo,jw)/pot(jw)
							;corrected EM
	avgem(oo,iw)=avgem(oo,iw)+em(oo,jw)
	emnorm(oo)=emnorm(oo)+1 	;for the normalization
      endif
    endfor			;JW=0,MOW-1}
    oo=where(emnorm gt 0,moo)
    if moo gt 0 then begin
      avgem(oo,iw)=avgem(oo,iw)/emnorm(oo)
      y0=min(avgem(oo,iw),imn)
      tmax(iw)=(logT(oo))(imn)
    endif
  endif				;MOW>0)
endfor					;IW=0,NW-1}

;	average
emavg=dblarr(nt) & sdev=dblarr(nt) & nwt=fltarr(nt)
if keyword_set(ewt) then wts=1./sig
for it=0,nt-1 do begin
  tmp=reform(avgem(it,*))
  oo=where(tmp gt 0,moo)
  if moo gt 0 then begin
    wt=wts & wt(oo)=wt(oo)/total(wt(oo))
    emavg(it)=total(tmp(oo)*wt(oo))
    sdev(it)=sqrt(total(wt(oo)*tmp(oo)^2)-emavg(it)^2)
  endif
endfor
oo=where(emavg gt 0,moo)
if moo gt 0 then begin
  y=min(emavg(oo),imn)
  s=(sdev(oo))(imn)
  t=(logT(oo))(imn)
endif

return,[y,s,t]
end
