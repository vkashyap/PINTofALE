function nenh,logT,abund=abund,eqfile=eqfile,chidir=chidir,$
	elem=elem,xelem=xelem,nproton=nproton, _extra=e
;+
;function 	nenH
;	calculate the ratio N(e)/N(H) (or N(e)/N(p) -- see keyword
;	NPROTON below) in a plasma of specified temperature
;
;syntax
;	r=nenh(logT,abund=abund,eqfile=eqfile,chidir=chidir,elem=elem,$
;	xelem=xelem,/nproton)
;
;parameters
;	logT	[INPUT; required] temperature(s) (log_10(T [K])) at which
;		to compute ratio.
;
;keywords
;	abund	[I/O] abundances of elements.  default is Allen abundances.
;	eqfile	[INPUT] pathname, relative to CHIDIR, of file containing
;		ionization equilibrium values, passed w/o comment to RD_IONEQ
;	chidir	[INPUT] path name to the CHIANTI installation, passed w/o
;		comment to RD_IONEQ
;	elem	[INPUT] return data for only these elements (default: all)
;	xelem	[INPUT] exclude these elements (overrides ELEM)
;		* default: none
;		  (note that older CHIANTI ion balance files do not have
;		  Cu, Zn, so originally XELEM was set by default to [29,30])
;	nproton	[INPUT] if set, returns the ratio N(e)/N(p)
;	_extra	[JUNK] ignore -- here only to prevent crashing program
;
;subroutines
;	GETABUND [RDABUND]
;	RD_IONEQ [READ_IONEQ (a CHIANTI routine)]
;	SYMB2ZION [LAT2ARAB]
;	KILROY
;
;history
;	vinay kashyap (Sep97)
;	modified to have just one call to RD_IONEQ instead of one for
;	  each Z, changed default of XELEM to none (VK; Jun02)
;-

;	usage
if n_elements(logT) eq 0 then begin
  print,'Usage: r=nenh(logT,abund=abund,eqfile=eqfile,chidir=chidir,$'
  print,'       elem=elem,xelem=xelem,/nproton)'
  print,'  compute ratio of e/H number density at specified temperatures'
  return,0.
endif

;	initialize
tlog=[logT(*)] & nt=n_elements(tlog)		;at which temperatures..
atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al','Si',$
	'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',$
	'Cu','Zn']				;elements from 1-30
nz=n_elements(atom)				;{H..Zn}
if n_elements(abund) ne nz then abund=getabund('allen')
reH=fltarr(nt)					;OUTPUT

;	if ELEM and/or XELEM are set...
if not keyword_set(xelem) then xelem=[-1]
sze=size(elem) & szx=size(xelem) & nsze=n_elements(sze) & nszx=n_elements(szx)
etyp=sze(nsze-2) & xtyp=szx(nszx-2) & elm=intarr(nz)
if etyp eq 0 then elm=indgen(nz)+1	;default -- all elements
if etyp gt 0 and etyp le 7 then begin
  if etyp eq 7 then begin qut=!quiet & !quiet=1 & endif
  nelm=sze(nsze-1)
  for i=0,nelm-1 do begin
    if etyp eq 7 then symb2zion,elem(i),z,ion else z=fix(elem(i))
    if z gt 0 then elm([z-1])=z
  endfor
  if etyp eq 7 then !quiet=qut
endif
if xtyp gt 0 and xtyp le 7 then begin
  if xtyp eq 7 then begin qut=!quiet & !quiet=1 & endif
  nxlm=szx(nszx-1)
  for i=0,nxlm-1 do begin
    if xtyp eq 7 then symb2zion,xelem(i),z,ion else z=fix(xelem(i))
    if z gt 0 then elm([z-1])=0
  endfor
  if xtyp eq 7 then !quiet=qut
endif
;
oz=where(elm gt 0,mz)		;select elements

;	compute electron number
ieq=rd_ioneq(elm(oz),tlog,eqfile=eqfile,chidir=chidir)
nprot=fltarr(nt)
for iz=0,mz-1 do begin			;{for each element in list
  zz=elm(oz(iz)) & nion=zz+1L > 1
  kilroy; was here
  if mz eq 1 then ioneq=ieq else ioneq=reform(ieq(*,0L:nion-1L,iz),nt,nion)
  numelec=findgen(zz+1)		;# of electrons for each ionization stage
  for it=0,nt-1 do begin		;{for each temperature
    reH(it)=reH(it)+abund(zz-1)*total(numelec*ioneq(it,*))
    if zz eq 1 then nprot(it)=abund[oz[iz]]*ioneq[it,1]
  endfor				;IT=0,NT-1}
endfor					;IZ=0,MZ-1}

;	If the proton density were required

;	this part because of CHIANTI routine PROTON_DENS()
;	-- note that the default operation of NENH() is equivalent
;	to PROTON_DENS(/hydrogen) and the default operation of
;	PROTON_DENS() is equivalent to NENH(/nproton))
;
;	example:
;		r=nenh(!logt,/nproton) & rc=proton_dens(!logT)
;		plot,!logT,rc,/ynoz & oplot,!logT,1./r,col=2
;	example:
;		r=nenh(!logt) & rc=proton_dens(!logT,/hydrogen)
;		plot,!logT,1./rc & oplot,!logT,r,col=2

oH=where(nprot gt 0,moH)
if moH eq 0 then begin
  ;	was H included in XELEM?  tch tch.
  ieq=rd_ioneq(1,tlog,eqfile=eqfile,chidir=chidir)
  nprot(*)=abund(0)*ioneq(*,1)
endif
oH=where(nprot gt 0,moH)
if moH eq 0 then message,'BUG?!'
if keyword_set(nproton) then begin
  tmp=0.*reH+1.
  tmp[oH]=reH[oH]/nprot[oH] & reH=tmp
endif

return,reH
end
