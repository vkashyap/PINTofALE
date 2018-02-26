function rd_ioneq,z,logt,eqfile=eqfile,chidir=chidir,zlost=zlost,$
	verbose=verbose, _extra=e
;+
;function	rd_ioneq
;	return the values of the ionization equilibrium for ions
;	of element Z at specified temperatures as ARRAY(NLOGT,NIONS,NZ).
;
;	right now this is simply a wrapper to the CHIANTI routine
;	READ_IONEQ.  someday other formats will be supported.
;
;syntax
;	ioneq=rd_ioneq(z,logt,eqfile=eqfile,chidir=chidir,zlost=zlost,$
;	verbose=verbose)
;
;parameters
;	z	[INPUT; required] atomic number
;		* must be integers
;	logt	[INPUT/OUTPUT] 1D array of log10(temperatures [K]) at which to
;		provide ionization equilibrium values.  if not defined, returns
;		the temperatures at which the values in EQFILE are given.
;
;keywords
;	eqfile	[INPUT] pathname, relative to CHIDIR, of file containing
;		ionization equilibrium values
;		* default: !IONEQF
;		* hard-coded default: ioneq/chianti.ioneq
;	chidir	[INPUT] path name to the CHIANTI installation
;		* default: !CHIDIR
;		* hard-coded default: /data/fubar/CHIANTI/dbase
;	zlost	[OUTPUT] average number of lost electrons at each temperature
;		* will be of size (NLOGT,NZ) or (NLOGT)
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT] junk -- here only to prevent program from crashing
;
;description
;	reads in ionization equilibrium values a la CHIANTI
;	and interpolates to appropriate temperatures
;
;	for example, the ionization fraction of Fe XVII would be
;		(rd_ioneq(26,logT))[*,17-1]
;	or for that matter,
;		(rd_ioneq([8,26,10],logT))[*,17-1,1]
;
;requires
;	SETSYSVAL
;	READ_IONEQ (from CHIANTI)
;	GETPOADEF
;
;history
;	vinay kashyap (Nov96; based on SYNTHETIC.PRO)
;	added keyword _EXTRA (VK; Dec96)
;	interpolation in log to avoid stepped look (VK; Feb97)
;	smooth after interpolation to avoid segmented look (VK; Apr97)
;	bug: force maximum IONEQ to be 1, not 10^(1)! (VK; 99Sep)
;	added keyword VERBOSE, now takes CHIDIR and EQFILE defaults
;	  from !CHIDIR and !IONEQF if they are defined; added call to
;	  SETSYSVAL (VK; Jul01)
;	cleaned up the interface, and allowed Z to be an array
;	  (VK; Jun02)
;	added keyword ZLOST, converted "()" to "[]" (VK: Apr06)
;	allowed EQFILE to be a plain filename; changed hardcoded
;	  defaults of EQFILE and CHIDIR (VK; Jul13)
;	changed default of EQFILE to ioneq/chianti.ioneq (VK; Aug13)
;	added call to GETPOADEF (VK; Aug15)
;-

;	usage
ok='ok' & np=n_params(0) & nz=n_elements(Z)
if np lt 1 then ok='Insufficient parameters' else $
 if nz eq 0 then ok='Z is undefined'
if ok ne 'ok' then begin
  print,'Usage: ioneq=rd_ioneq(Z,logT,chidir=!CHIDIR,eqfile=!IONEQF,verbose=v)'
  print,'  return 2D (or 3D) array of ionization equilibrium values for given'
  print,'  element at specified temperatures'
  if np ne 0 then message,ok,/info
  return,fltarr(1,1)
endif

;	initialize
zslash='/'
case !version.OS_FAMILY of
  'unix': zslash='/'
  'windows': zslash='\'
  'macos': zslash=':'
  'vms': zslash='\'
  else: zslash='/'	;unknown OS, assume UNIX
endcase

;	check keywords
;	first figure out where CHIANTI is installed
;zCHIDIR='/data/fubar/CHIANTI/dbase'
zCHIDIR=getpoadef('CHIDIR')
ivar=0 & defsysv,'!CHIDIR',exists=ivar	;if !CHIDIR exists
if ivar ne 0 then setsysval,'CHIDIR',zCHIDIR,/getval
if not keyword_set(chidir) then chidir=zCHIDIR
;	do the same for the ion-balance file
;zEQFILE='ioneq/chianti.ioneq'
zEQFILE=getpoadef('IONEQF')
ivar=0 & defsysv,'!IONEQF',exists=ivar  ;if !IONEQF exists
if ivar ne 0 then setsysval,'IONEQF',zEQFILE,/getval
if not keyword_set(eqfile) then eqfile=zEQFILE
;	verbosity
vv=0
ivar=0 & defsysv,'!VERBOSE',exists=ivar	;if !VERBOSE exists
if ivar ne 0 then setsysval,'VERBOSE',vv,/getval
if n_elements(verbose) gt 0 then vv=long(verbose[0])
if vv gt 5 then message,'Assuming CHIDIR='+chidir,/info
if vv gt 5 then message,'Using EQFILE='+eqfile,/info

;	call CHIANTI routine, which returns the ion balance as
;	an array of the form (NLOGT,NZ,NION)
;	NOTE: this is currently the only option.  some day we will
;	have other formats available which will be accessible via
;	appropriate keywords

;	figure out the name of the input file
if strpos(eqfile,zslash) ge 0 then begin	;(EQFILE is fragment of a path
  infile=chidir+zslash+eqfile
endif else begin				;path fragment)(EQFILE is filename
  infile=(file_search(chidir,eqfile))[0]	;pick the first option
endelse						;isolated filename)
read_ioneq,infile,tt,ieq,ref
sieq=size(ieq) & ntt=n_elements(tt)
if max(tt)/min(tt) gt 100 then begin
  if vv gt 1 then message,$
	'converting temperatures from ion balance file to log10',/info
  ltt=alog10(tt)
endif else ltt=tt

;	check inputs and define output array
nT=n_elements(logT)
if nT eq 0 then begin
  logT=ltt & nT=n_elements(logT)
endif
if nz eq 1 then begin
  nion=Z[0]+1L & ioneq=fltarr(nT,nION)
  zlost=fltarr(nT)
endif else begin
  nion=max(z)+1L > 1
  ioneq=fltarr(nT,nION,nZ)
  zlost=fltarr(nT,nZ)
endelse

;	if file temperature grid is different from requested grid
k=0	;implies same grid
if ntt ne nt then k=1	;number of elements in grid differ
if k eq 0 then for i=0,nt-1 do k=k+(logt[i] ne ltt[i])	;even so, values differ

;	and now put the data that were read in into the output array
for iz=0L,nZ-1L do begin	;{for each Z
  kZ=Z[iZ]
  mion=kZ+1L > 1L
  if sieq[2] ge kZ and kZ gt 0 then begin	;(data exist in IEQ
    jeq=reform(ieq[*,kZ-1L,0L:mion-1L])
    joneq=fltarr(nT,mION)
    if k ne 0 then begin	;(grids differ, so interpolate
      for i=0L,mion-1L do begin		;{for each ion
	tmp=reform(jeq[*,i]) > 1e-17
	tmp=alog10(tmp)
	ineq=interpol(tmp,ltt,logt) < 0
	;NO	ineq=spl_interp(ltt,tmp,spl_init(ltt,tmp,/d),logt,/d)<1
	;NO	ineq=spline(ltt,tmp,logt,2)<1
	ineq=10.^(ineq)
	oo=where(ineq lt 1e-9,moo) & if moo gt 0 then ineq[oo]=0.
	joneq[*,i]=ineq[*]
      endfor				;I=0,MION-1}
    endif else joneq[*,0L:mion-1L]=reform(ieq[*,kZ-1L,0:mion-1L])	;k.NE.0)
    if nz eq 1 then begin
      ioneq=joneq
      cion=ioneq
      for ii=0L,kZ do cion[*,ii]=total(reform(ioneq[*,ii]),/cumul)
      for ii=0L,kZ do cion[*,ii]=cion[*,ii]/cion[nT-1L,ii]
      for ii=0L,kZ do zlost=zlost+cion[*,ii]
    endif else begin
      ioneq[*,0L:mion-1L,iz]=joneq[*,*]
      cion=ioneq
      for ii=0L,kZ do cion[*,ii,iz]=total(reform(ioneq[*,ii,iz]),/cumul)
      for ii=0L,kZ do cion[*,ii,iz]=cion[*,ii,iz]/cion[nT-1L,ii,iz]
      for ii=0L,kZ do zlost[*,iz]=zlost[*,iz]+cion[*,ii,iz]
    endelse
  endif else begin			;SIEQ[2].GE.Z[IZ])(no ion balance data
    if vv gt 1 then message,$
	'Ion balance data do not exist for Z='+strtrim(kZ,2),/info
    if nZ eq 1 then begin
      ioneq[*,0L:mion-1L]=1./float(mion)
      cion=ioneq
      for ii=0L,mion-1L do cion[*,ii]=total(reform(ioneq[*,ii]),/cumul)
      for ii=0L,mion-1L do cion[*,ii]=cion[*,ii]/cion[mion-1L,ii]
      for ii=0L,mion-1L do zlost=zlost+cion[*,ii]
    endif else begin
      ioneq[*,0L:mion-1L,iz]=1./float(mion)
      cion=ioneq
      for ii=0L,mion-1L do cion[*,ii,iz]=total(reform(ioneq[*,ii,iz]),/cumul)
      for ii=0L,mion-1L do cion[*,ii,iz]=cion[*,ii,iz]/cion[mion-1L,ii,iz]
      for ii=0L,mion-1L do zlost[*,iz]=zlost[*,iz]+cion[*,ii,iz]
    endelse
  endelse				;SIEQ[2]<Z[IZ])
endfor			;IZ=0,NZ-1}

;;catch trivial error
;sieq=size[ieq] & nt=n_elements(logt) & ntt=n_elements(tt) & nion=z+1
;if nt eq 0 then begin & nt=ntt & logt=tt & endif
;if sieq[2] lt Z then begin
;  message,"er...doesn't look like we got enough elements in "+infile,/info
;  return,fltarr(nt,z+1)+(1./float(z+1.))
;endif
;
;;now, IEQ is a 3D array (T,Z,ION); from which take a slice of Z
;ieq=reform(ieq[*,Z-1,0:nion-1])
;
;;if tt is different from t...
;k=0 & if ntt ne nt then k=1			;k=0: no difference
;if k eq 0 then for i=0,nt-1 do k=k+(logt[i] ne tt[i])
;;
;;...linear interpolate
;if k ne 0 then begin
;  ioneq=fltarr(nt,nion)
;  for i=0,nion-1 do begin
;    tmp=reform(ieq[*,i])>1e-17
;    tmp=alog10(tmp)
;    ineq=interpol(tmp,tt,logt)<0
;    	;NO	ineq=spl_interp(tt,tmp,spl_init(tt,tmp,/d),logt,/d)<1
;    	;NO	ineq=spline(tt,tmp,logt,2)<1
;    ;if ntt lt nt then ineq=smooth(ineq,3)
;    ineq=10.^(ineq)
;    oo=where(ineq lt 1e-9) & if oo[0] ne -1 then ineq[oo]=0.
;    ioneq[*,i]=ineq[*]
;    x=tt & y=tmp & x2=logt
;  endfor
;endif else ioneq=ieq

return,ioneq
end
