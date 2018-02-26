function pottasch,emis,logt,z=z,ion=ion,level=level,wdth=wdth,trange=trange,$
	_extra=e
;+
;function	pottasch
;	returns the Pottasch factor that is a measure of the width of
;	the ion emissivity in temperature space over which most of
;	the line emission occurs.
;
;	the Pottasch factor is the ratio
;		{\int [EMIS*ION_BALANCE] dlogT}/{WDTH*max([EMIS*ION_BALANCE])}
;	where WDTH is the range in logT where EMIS*ION_BALANCE drops to a
;	set fraction of the maximum.
;
;	output has as many elements as the 2nd dimension of EMIS, one point
;	for each line.
;
;syntax
;	pot=pottasch(emis,logt,z=z,ion=ion,level=level,wdth=wdth,$
;	trange=trange,chifil=chifil)
;
;parameters
;	emis	[INPUT; required] 1- or 2-D array of line intensities as
;		a function of temperature
;		* if 1D, assumed to be EMIS(Temperature)
;		* if 2D, assumed to be EMIS(Temperature,{Z,ION,Wavelength})
;		* if Z is not specified, assumed to be the product of
;		  line emissivity and ion balance
;	logt	[INPUT] 1D array of log10(T[K])
;		* WARNING: if not given, or if size does not match 1st dim
;		  of EMIS, makes a dumb guess
;
;keywords
;	Z	[INPUT] atomic numbers corresponding to the 2nd dimension
;		* if not set, ION is ignored and FOLD_IONEQ is not called
;		* if size does not match 2nd dim of EMIS,
;		  >: ignore extra elements
;		  <,>1: use as many as specified, take rest to be H (Z=1)
;		  1: all are of same element
;	ion	[INPUT; default: Z+1] ionic state
;	level	[INPUT; default: 0.1] defines what WDTH means
;	wdth	[OUTPUT] full width in logT at LEVEL*maximum
;		* minimum value at any point is the resolution in LOGT at
;		  that point.
;	trange	[OUTPUT] 2-column array of the boundaries of logT that
;		define WDTH
;		* trange(0,*)=lower bound, trange(1,*)=upper bound
;	_extra	[INPUT] use to pass defined keywords
;		* currently, only CHIFIL to FOLD_IONEQ
;
;subroutines
;	FOLD_IONEQ [WHEE, RD_IONEQ [READ_IONEQ]]
;
;history
;	vinay kashyap (Apr97)
;	changed temperature width calculation to account for bumpy
;	  emissivities (VK; Oct98)
;-

;	usage
sze=size(emis) & & nsze=n_elements(sze) & nt=n_elements(logt)
if sze(0) eq 0 then begin
  print,'Usage: b=pottasch(emissivity,logT,Z=Z,ion=ion,wdth=wdth)'
  print,'  return Pottasch factor describing "width" of line emissivity'
  return,-1L
endif

;	check input
if sze(0) gt 2 then begin
  message,'Input not understandable',/info & return,-1L
endif
if sze(1) ne nt then begin		;guess the temperature grid
  nt=sze(1)
  logt=findgen(nt)*(4./(nt-1))+4.
  if nt eq 21 then logt=findgen(nt)*0.1+5
  if nt eq 41 then logt=findgen(nt)*0.05+5
endif
;
mz=sze(2) & if sze(0) eq 1 then mz=1L & nz=n_elements(z) & zz=intarr(mz)+1
ieq=nz
if mz ne nz then begin			;atomic numbers
  if nz eq 0 then z=[1] else begin
    if nz lt mz then zz(0:nz-1)=([z])(*) else zz=z(0:mz-1)
  endelse
endif else zz=z
;
ni=n_elements(ion) & jon=zz+1
if ni ne mz then begin			;ionic states
  if ni eq 0 then ion=jon else begin
    if ni lt mz then jon(0:ni-1)=([ion])(*) else jon=ion(0:mz-1)
  endelse
endif else jon=ion

;	fold ion balance if needed
if ieq gt 0 then line=fold_ioneq(emis,zz,jon,logt=logt,_extra=e) else line=emis

;	integrate over logT
numer=dblarr(mz)
for i=0,mz-1 do numer(i)=int_tabulated(logt,reform(line(*,i)))

;	find maxima & widths
lmax=dblarr(mz) & wdth=fltarr(mz) & trange=fltarr(2,mz) & dlogT=logT(1:*)-logT
if n_elements(level) eq 0 then dw=0.1 else dw=float(level)
for i=0,mz-1 do begin
  tmp=reform(line(*,i))
  tmx=max(tmp,imx) & lmax(i)=tmx
  oo=where(tmp ge dw*tmx,moo)
  wdth_min=max([abs(logt([imx-1])-logt(imx)),abs(logt([imx+1])-logt(imx))])
    ;{VK: the following used to be the original width calculation, now
    ;corrected to include the possibility that EMIS may have multiple bumps
    ;
    ;wdth(i)=abs(logt(moo-1)-logt(0)) > wdth_min
    ;
  if moo gt 0 then wdth(i)=total(dlogT(oo)) else message,'bug!'
    ;
    ;VK}
  trange(0,i)=logt(oo(0)) & trange(1,i)=logt(oo(moo-1))
endfor

;	compute Pottasch factor
bet=dblarr(mz)+1.
denom=wdth*lmax & oo=where(denom gt 0,moo)
if moo gt 0 then bet(oo)=numer(oo)/denom(oo)

return,bet
end
