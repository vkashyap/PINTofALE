;+
;MK_ROSAT_HR.PRO
;	program to compute hardness ratios of counts in the ROSAT/PSPC
;	obtained for 1-T plasmas
;
;	read in the appropriate effective area curve, read in the line
;	and continuum spectral databases, compute predicted count rates
;	at the various temperatures for a specified emission measure in
;	the various passbands of interest, save results to an IDL save
;	file for later use with ROSAT_HR.PRO
;
;	requires:
;		subroutines ABUND, GET_PIMMS_FILE, WC, RD_LINE, RD_CONT,
;			FOLD_IONEQ, ISMTAU
;		Line and continuum SCAR databases
;
;vinay kashyap (Mar98)
;changed defaults and added call to getpoadef (VK; Aug15)
;-

;	set up keywords and various variables
mission='ROSAT' & instrument='PSPC'
if not keyword_set(configuration) then configuration='OPEN'
if not keyword_set(pdir) then pdir='/soft/prop_cli/config/pimms/data/'
if not keyword_set(arfile) then begin		;(read in effective areas
  arfile=get_pimms_file(mission,instrument,configuration,pdir=pdir)
  nlin=wc(arfile) & nwrd=wc(arfile,/word) & ncol=nwrd/nlin
  openr,ua,arfile,/get_lun & var=fltarr(ncol,nlin)
  readf,ua,var & close,ua & free_lun,ua
endif						;ARFILE)
if not keyword_set(ldbdir) then ldbdir=getpoadef('LDBDIR')
if not keyword_set(cdbdir) then cdbdir=getpoadef('CDBDIR')
if not keyword_set(chifil) then chifil=getpoadef('CHIFIL')
if not keyword_set(eden) then eden=1e9
if not keyword_set(abund) then abund=getabund(getpoadef('ABREF'))
;
if not keyword_set(flstr) then flstr=1
if not keyword_set(fcstr) then fcstr=1
if not keyword_set(ecol) then ecol=0L
if not keyword_set(acol) then acol=1L
if not keyword_set(EM) then EM=1d10	;[cm^-5]
if not keyword_set(e1l) then e1l=0.1
if not keyword_set(e1u) then e1u=0.4
if not keyword_set(e2l) then e2l=0.5
if not keyword_set(e2u) then e2u=2.0
if not keyword_set(e3l) then e3l=0.5
if not keyword_set(e3u) then e3u=0.9
if not keyword_set(e4l) then e4l=0.9
if not keyword_set(e4u) then e4u=2.0

;	define backup keywords
if not keyword_set(oarfil) then oarfil=arfile
if not keyword_set(opdir) then opdir=pdir
if not keyword_set(oldbdir) then oldbdir=ldbdir
if not keyword_set(ocdbdir) then ocdbdir=cdbdir
if not keyword_set(ochifil) then ochifil=chifil
if not keyword_set(oeden) then oeden=eden
if not keyword_set(oabnd) then oabnd=abund

;	read in effective areas
if oarfil ne arfile then begin
  arfile=get_pimms_file(mission,instrument,configuration,pdir=pdir)
  nlin=wc(arfile) & nwrd=wc(arfile,/word) & ncol=nwrd/nlin
  openr,ua,arfile,/get_lun & var=fltarr(ncol,nlin)
  readf,ua,var & close,ua & free_lun,ua
endif
arf_e=reform(var(ecol,*)) & arf_a=reform(var(acol,*))
arf_w=12.3985/arf_e

;	set the wavelength range
if not keyword_set(emin) then emin=min(arf_e)
if not keyword_set(emax) then emax=max(arf_e)
wmin=12.3985/emax & wmax=12.3985/emin
if not keyword_set(owmn) then owmn=wmin
if not keyword_set(owmx) then owmx=wmax
w1u=12.3985/e1l & w1l=12.3985/e1u & w2u=12.3985/e2l & w2l=12.3985/e2u
w3u=12.3985/e3l & w3l=12.3985/e3u & w4u=12.3985/e4l & w4l=12.3985/e4u

;	read in databases IF..
if owmn ne wmin or owmx ne wmax then begin flstr=1 & fcstr=1 & endif
if oldbdir ne ldbdir then flstr=1
if ocdbdir ne cdbdir then fcstr=1
if total(abund-oabnd) gt min([abund,oabnd]) then fcstr=1

;	read in the line database
if n_tags(flstr) ne 8 then begin
  lem=rd_line(atom,n_e=eden,wrange=[wmin,wmax],dbdir=ldbdir,fstr=flstr)
  lem=fold_ioneq(lem,flstr.Z,flstr.jon,logT=flstr.logT,chifil=chifil)
  flstr.line_int=lem
endif else lem=flstr.line_int
nrgl=6.626176e-27*2.9979e10/(abs(flstr.wvl)/1e8)
w=abs(flstr.wvl)

;	read in continuum database
if n_tags(fcstr) ne 6 then begin
  cem=rd_cont('cie',n_e=eden,wrange=[wmin,wmax],dbdir=cdbdir,abund=abund,$
	wavl=wavl,logT=logT,fcstr=fcstr)
endif else cem=fcstr.cont_int
wave=fcstr.midwvl
nrgc=6.626176e-27*2.9979e10/(wave/1e8)

;	get interstellar absorption
etaul=1.+0.*flstr.wvl & etauc=0.*wave+1.
if keyword_set(NH) then begin
  etaul=exp(-ismtau(abs(flstr.wvl),NH=NH,fH2=fH2,He1=He1,HeII=HeII,Fano=Fano))
  etauc=exp(-ismtau(wave,NH=NH,fH2=fH2,He1=He1,HeII=HeII,Fano=Fano))
endif

;	interpolate the effective area arrays
arl=interpol(arf_a,arf_w,abs(flstr.wvl))
arc=interpol(arf_a,arf_w,wave)

;	compute total line and continuum counts at each temperature
logT=flstr.logT & nT=n_elements(logT) & fxl=fltarr(nT) & fxc=fxl
fxl1=fxl & fxl2=fxl & fxl3=fxl & fxl4=fxl
fxc1=fxc & fxc2=fxc & fxc3=fxc & fxc4=fxc
hrA=fxl & hrB=fxl
nZ=30
for iT=0,nT-1 do begin
  dW=wavl(1:*)-wavl
  fxc(iT)=total(EM*cem(iT,*)*etauc(*)*arc(*)*dW(*)/nrgc/1e23)
  o1=where(wavl ge w1l and wavl lt w1u,mo1)
  o2=where(wavl ge w2l and wavl lt w2u,mo2)
  o3=where(wavl ge w3l and wavl lt w3u,mo3)
  o4=where(wavl ge w4l and wavl lt w4u,mo4)
  if mo1 ne 0 then fxc1(iT)=total(EM*cem(iT,o1)*etauc(o1)*arc(o1)*dW(o1)/nrgc/1e23)
  if mo2 ne 0 then fxc2(iT)=total(EM*cem(iT,o2)*etauc(o2)*arc(o2)*dW(o2)/nrgc/1e23)
  if mo3 ne 0 then fxc3(iT)=total(EM*cem(iT,o3)*etauc(o3)*arc(o3)*dW(o3)/nrgc/1e23)
  if mo4 ne 0 then fxc4(iT)=total(EM*cem(iT,o4)*etauc(o4)*arc(o4)*dW(o4)/nrgc/1e23)
  ;	[cm^-5]*[1e-23 ergs cm^3/s/A]*[]*[cm^2]*[A]/[ergs/ph]*[1e23]
  for iZ=1,nZ do begin
    oZ=where(flstr.Z eq iZ,moZ)
    if moZ gt 0 then fxl(iT)=fxl(iT)+$
	abund(iZ-1)*total(EM*lem(iT,oZ)*etaul(oZ)*arl(oZ)/nrgl(oZ)/1e23)
  	;	[]*[cm^-5]*[1e-23 ergs cm^3/s]*[]*[cm^2]/[ergs/ph]*[1e23]
    oZ1=where(flstr.Z eq iZ and (w ge w1l and w lt w1u),moZ1)
    oZ2=where(flstr.Z eq iZ and (w ge w2l and w lt w2u),moZ2)
    oZ3=where(flstr.Z eq iZ and (w ge w3l and w lt w3u),moZ3)
    oZ4=where(flstr.Z eq iZ and (w ge w4l and w lt w4u),moZ4)
    if moZ1 gt 0 then fxl1(iT)=fxl1(iT)+$
	abund(iZ-1)*total(EM*lem(iT,oZ1)*etaul(oZ1)*arl(oZ1)/nrgl(oZ1)/1e23)
    if moZ2 gt 0 then fxl2(iT)=fxl2(iT)+$
	abund(iZ-1)*total(EM*lem(iT,oZ2)*etaul(oZ2)*arl(oZ2)/nrgl(oZ2)/1e23)
    if moZ3 gt 0 then fxl3(iT)=fxl3(iT)+$
	abund(iZ-1)*total(EM*lem(iT,oZ3)*etaul(oZ3)*arl(oZ3)/nrgl(oZ3)/1e23)
    if moZ4 gt 0 then fxl4(iT)=fxl4(iT)+$
	abund(iZ-1)*total(EM*lem(iT,oZ4)*etaul(oZ4)*arl(oZ4)/nrgl(oZ4)/1e23)
  endfor
endfor

fx=fxl+fxc	;[ph/s]
fx1=fxl1+fxc1 & fx2=fxl2+fxc2 & fx3=fxl3+fxc3 & fx4=fxl4+fxc4
;oA=where(fx2 gt 1e-6*max(fx2),moA) & oB=where(fx4 gt 1e-6*max(fx4),moB)
;hrA(oA)=fx1(oA)/fx2(oA) & hrB(oB)=fx3(oB)/fx4(oB)
oA=where(fx1+fx2 gt 0,moA) & oB=where(fx3+fx4 gt 0,moB)
hrA(oA)=(fx2(oA)-fx1(oA))/(fx1(oA)+fx2(oA))
hrB(oB)=(fx4(oB)-fx3(oB))/(fx3(oB)+fx4(oB))

save,file='ROSAT_HRT.sav',logT,hrA,hrB,fx,fx1,fx2,fx3,fx3,oA,oB

oarfil=arfile & opdir=pdir & oldbdir=ldbdir & ocdbdir=cdbdir & ochifil=chifil
oeden=eden

end
