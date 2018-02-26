function fold_ioneq,emis,z,ion,logt=logt,tmax=tmax,trng=trng,level=level,$
	userarr=userarr,eqfile=eqfile,chifil=chifil,verbose=verbose, _extra=e
;+
;function	fold_ioneq
;	folds in ionization equilibrium fractions of the various ions
;	of a given element with line emissivities (which include only
;	level population info; see RD_LINE()).  returns
;	emis(T;Z,ion)*ioneq(T;Z,ion)
;
;syntax
;	lflx=fold_ioneq(emis,Z,ion,logt=logt,tmax=tmax,trng=trng,$
;	level=level,userarr=userarr,chifil=chifil,chidir=chidir,$
;	eqfile=eqfile,verbose=v)
;
;parameters
;	emis	[INPUT; required] 1- or 2-D array of line emissivities as
;		a function of temperature
;		* if 2D, assumed to be EMIS(Temperature,Wavelength)
;		* if 1D, assumed to be EMIS(Temperature), unless the
;		  number of elements matches that of Z
;	Z	[INPUT; required] atomic number(s) corresponding to
;		each wavelength
;		* if size does not match the 2nd dimension of EMIS, then
;		  if >, ignore extra elements
;		  if <, but >0, use available ones, and use Z[0] for rest
;	ion	[INPUT; default: Z+1] ionic state(s) corresponding to
;		each wavelength
;		* if size does not match the 2nd dimension of EMIS, then
;		  if >, ignore extra elements
;		  if <, but >0, use available ones, and use ION[0] for rest
;		* NOTE: this is the ionic state >>that matters<<, i.e.,
;		  the ionic state of the species that populates the
;		  upper level.  See keyword JON of RD_LINE, which is
;		  what is needed here.
;
;keywords
;	logt	[INPUT] 1D array of log10(temperature [K]) at which EMIS
;		is given. (default: findgen(81)*0.05+4)
;		* if size does not match that of EMIS, use default
;	level	[INPUT; default: 0.5] level at which to determine the
;		temperature range on line flux (in other words, default
;		will be to return the FWHM)
;		* if <1, TRNG=LOGT(where(FLUX>LEVEL*MAX(FLUX))
;		* if >=1, TRNG=LOGT(where(FLUX>(LEVEL/100.)*MAX(FLUX))
;		* if <=0 or >=100, set to 0.5
;	tmax	[OUTPUT] log10(temperature) at which contribution is maximum
;	trng	[OUTPUT] range in LOGT (LEVEL of MAX)
;		-- TRNG(0,*) is lower bound, TRNG(1,*) is upper bound
;	userarr	[INPUT] an array of ionization fractions in the same
;		format as is returned by RD_IONEQ(), and has the size
;		[n(T),max(Z)+1,n(Z)]
;		* if this is given and is legit (i.e., matches the
;		  supplied EMIS and LOGT), then EQFILE and CHIFIL are
;		  ignored
;	eqfile	[I/O] file from which to read in ion-balance data
;		* default is to read in the CHIANTI file, !IONEQF
;		* hard-coded default is
;		  ioneq/chianti.ioneq
;		* if not given, but CHIFIL is given, uses that.
;	chifil	[I/O] if set, reads in from CHIANTI-type database; must be
;		set to the name of the file containing the ionization
;		equilibrium info (i.e., CHIDIR/CHIFIL)
;		* may be overridden with EQFILE
;		* default is ioneq/chianti.ioneq
;		* currently, set automatically if not given.  
;	verbose	[INPUT] higher the number, greater the chatter
;	_extra	[INPUT] use to transmit defined keywords to called functions
;		* RD_IONEQ [CHIDIR]
;
;subroutines
;	-- WHEE
;	-- SETSYSVAL
;	-- RD_IONEQ [READ_IONEQ (%)]
;	(%) CHIANTI subroutine, used as is
;
;history
;	vinay kashyap (Dec96)
;	removed call to KILROY and added call to WHEE, hacked default
;	  use of READ_IONEQ until more options become available; corrected
;	  bug that was reading in ionstate-1 rather than ionstate; added
;	  code to handle Z=0 (VK; Feb97)
;	bug: if Z is scalar take em all to be same element (VK; Apr97)
;	added keyword VERBOSE (VK; MarMM)
;	streamlined meanings of EQFILE, CHIFIL (VK; DecMM)
;	changed default for EQFILE; now pass VERBOSE to RD_IONEQ; added
;	  call to SETSYSVAL (VK; Jul01)
;	allowed EQFILE to be "none" (VK; Feb04)
;	changed default behavior of Z and ION when sizes don't match EMIS;
;	  now uses first element of input if given (VK; Nov04)
;	bug correction, ION is offset by 1 from array index (LL/VK; Dec04)
;	added keyword USERARR, slightly modified behavior when EMIS
;	  was input as 1D (VK; Aug08)
;	changed default of EQFILE and CHIFIL to ioneq/chianti.ioneq (VK; Aug13)
;-

;	usage
ok='ok' & np=n_params() & nem=n_elements(emis) & szl=size(emis)
nz=n_elements(z) & ni=n_elements(ion)
if np lt 2 then ok='Insufficient parameters' else $
 if nem eq 0 then ok='EMIS is not defined' else $
  if nz eq 0 then ok='Z is not defined' else $
   if szl(0) gt 2 then ok='EMIS is not in recognizable format'
if ok ne 'ok' then begin
  print,'Usage: lflx=fold_ioneq(emissivity,Z,ion,logt=logt,tmax=tmax,trng=trng,$'
  print,'       level=level,userarr=userarr,eqfile=eqfile,chifil=chifil,$'
  print,'       chidir=chidir,verbose=v)'
  print,'  fold line emissivities read in with RD_LINE with ion balance'
  print,'  also accepts defined keywords CHIDIR and CHIFIL for RD_IONEQ'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check input
if nz eq 0 then begin & z=1 & nz=1 & endif	;default element is H
if ni eq 0 then begin & ion=z+1 & ni=nz & endif	;default ionic state is Z+1
case szl(0) of
  1: begin				;emis=emis(LOGT) or emis(WVL)
    if szl(1) eq nz then begin
      nt=1 & nw=nz & atno=[z(*)] & jon=[ion(*)]-1	;the "-1" coz of IDL
    endif else begin
      nt=szl(1) & nw=1 & atno=[z(0)] & jon=[ion(0)]-1	;the "-1" coz of IDL
    endelse
  end
  2: begin				;emis=emis(LOGT,WVL)
    nt=szl(1) & nw=szl(2) & atno=intarr(nw)+1 & jon=atno
    if nz gt 0 then atno(*)=z(0)
    if ni gt 0 then jon(*)=ion(0)
    if nz gt 0 and nz le nw then atno(0L:nz-1L)=z
    if ni gt 0 and ni le nw then jon(0L:nz-1L)=ion-1L
    ;if nz le nw then atno(0:nz-1)=([z])(*) else atno(*)=z(0:nw-1)
    ;if ni le nw then jon(0:ni-1)=([ion])(*)-1 else jon(*)=ion(0:nw-1)-1
  end
  else: return,-1L			;nolle comprendo
endcase
;
if not keyword_set(logt) then begin
  tlog=findgen(nt)*(4./((nt>2)-1.))+4.
endif else tlog=logt
if n_elements(tlog) ne nt then begin
  message,'Line emissivity does not match temperature grid',/info
  tlog=findgen(nt)*(4./((nt>2)-1.))+4.
endif
;
if not keyword_set(level) then lev=0.5 else lev=level
if lev le 0. or lev ge 100. then lev=0.5
if lev ge 1. then lev=lev/100.
;
;	HACK ALERT!!!  until there are other ways to compute the
;	ion balance, use CHIANTI!!!  to avoid using CHIANTI, explicitly
;	set CHIFIL=0
if n_elements(chifil) eq 0 then chifil=1
if keyword_set(userarr) then begin
  ok='ok' & szu=size(userarr)
  if szu(0) lt 2 then ok='USERARR must be 3D (Temp,Z,ION)' else $
   if szu(1) ne nt then ok='USERARR does not match LOGT' else $
    if szu(2) ne nz then ok='USERARR does not match Z' else $
     if szu(3) ne max(Z)+1 then ok='USERARR does not match ION'
  if ok eq 'ok' then useuserarr=1 else $	;no need to call CHIANTI!
    useuserarr=0
endif
;
zEQFILE='ioneq/chianti.ioneq'
ivar=0 & defsysv,'!IONEQF',exists=ivar  ;if !IONEQF exists
if ivar ne 0 then setsysval,'IONEQF',zEQFILE,/getval
if not keyword_set(eqfile) then eqfile=zEQFILE
;def_eqfil='ioneq/arnaud_rothenflug.ioneq'
sze=size(eqfile) & nsze=n_elements(sze)
szc=size(chifil) & nszc=n_elements(szc)
if not keyword_set(eqfile) or sze(nsze-2) ne 7 then begin	;no EQFILE
  if szc(nszc-2) eq 7 then eqfile=chifil else eqfile=zEQFILE ;def_eqfil
endif

;	chatter
vv=0 & if keyword_set(verbose) then vv=long(verbose(0)) > 1

;	initialize
lflx=0*emis & tmax=fltarr(nw) & trng=fltarr(2,nw)

;	find all the unique elements
zz=atno(uniq(atno,sort(atno))) & mz=n_elements(zz)

;	for each element, read in ionization equilibrium and multiply
for iz=0L,mz-1L do begin
  oz=where(atno eq zz(iz),noz)
  if zz(iz) le 0 then begin		;{error! error!
    ioneq=fltarr(nt)+1.
    if noz gt 0 then jon(oz)=0
  endif else begin			;}{all OK
    ioneq=fltarr(nt,zz(iz)+1)+1.	;default -- multiply by 1
    ;
    whee,spoke=spoke,/moveit	;we're working!
    if keyword_set(useuserarr) then ioneq=userarr(*,*,iz) else begin
      if keyword_set(chifil) then begin	;use CHIANTI
        szc=size(chifil) & nszc=n_elements(szc)
        if szc(nszc-2) eq 7 then eqfile=chifil else begin
	  if not keyword_set(eqfile) then eqfile=zEQFILE ;'ioneq/chianti.ioneq'
        endelse
        ;eqfile=chifil
        ;if szc(nszc-2) ne 7 then eqfile=zEQFILE ;'ioneq/chianti.ioneq'
        if strlowcase(eqfile) ne 'none' then begin
          if vv gt 0 then message,'Using ion balances from: '+eqfile,/info
          ioneq=rd_ioneq(zz(iz),tlog,eqfile=eqfile,verbose=vv,_extra=e)
        endif
      endif
    endelse
  endelse				;}
  ;
  for ii=0L,noz-1L do begin
    if vv gt 1 then whee,spoke=spoke,/moveit	;we're working!
    ;	fold in ionization equilibrium fractions
    lflx(*,oz(ii))=emis(*,oz(ii))*ioneq(*,jon(oz(ii)))
    ;	and here figure out TMAX ...
    tmp=reform(lflx(*,oz(ii))) ;& imax=where(tmp eq max(tmp),moo)
    tmpmax=max(tmp,imax)
    tmax(oz(ii))=total(tlog(imax));/moo
    ;	... and TRNG
    oo=where(tmp gt lev*tmpmax,moo)
    ;oo=where(tmp gt lev*tmp(imax(0)),moo)
    if moo ne 0 then begin
      trng(0,oz(ii))=tlog(oo(0)) & trng(1,oz(ii))=tlog(oo(moo-1))
    endif
  endfor
  if vv gt 1 then whee,spoke=spoke,/moveit	;we're working!
  ;
endfor

return,lflx
end
