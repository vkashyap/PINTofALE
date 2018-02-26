function translabel,chianti,wmn=wmn,wmx=wmx,tex=tex,config=config, _extra=e
;+
;function	translabel
;	returns 2D array describing the levels involved in each transition
;	array(0,*) are lower levels, and array(1,*) are upper levels
;	
;	if DR transitions, then a "(DR)" tag is appended to the lower level
;	designations.
;
;usage
;	trans=translabel(chianti,wmn=wmn,wmx=wmx,/tex,config=config)
;
;parameters
;	chianti	[INPUT; required] structure that contains all the relevant
;		information gleaned from the CHIANTI database (type
;			DBHELP=RD_CHIANTI() & HELP,DBHELP,/STRUCTURE
;		for a description of what this structure contains)
;
;keywords
;	wmn	[INPUT; default=0 Ang] only include wavelengths > WMN
;	wmx	[INPUT; default=900 Ang] maximum wavelength to include
;	tex	[INPUT] if set, returns the levels formatted for TeX
;	config	[OUTPUT] returns the electron configurations of the lower
;		and upper levels in 2D array
;		* config(0,*) are for lower levels, and config(1,*) are
;		  for upper levels
;	_extra	[INPUT] junk -- ignore
;
;history
;	vinay kashyap (Nov96; based on LATEX_WVL_DEM.PRO of CHIANTI)
;	made loops go over long integers (VK; Jun10)
;-

np=n_params(0)
if np eq 0 then begin
  print,'Usage: trans=translabel(chianti,config=config,wmn=wmn,wmx=wmx,/tex)'
  print,'  returns 2D array describing levels involved in transitions'
  return,-1L
endif

;check keywords
if not keyword_set(wmn) then wmn=0.
if not keyword_set(wmx) then wmx=900.
trans=strarr(2,1) & config=trans

;extract from structure
wvl=chianti.wvl					;wavelengths
llo=chianti.lev1 & lup=chianti.lev2		;lower and upper levels
jj=chianti.j & ll=chianti.l & ss=chianti.ss	;J,L,2S+1
term=chianti.term				;e configurations
delec=chianti.DR				;DR transition?

;initialize
spd=['S','P','D','F','G','H','I','K','?']	;designation
intj=where(fix(jj) eq jj,nintj)			;integer Js
fltj=where(fix(jj) ne jj,nfltj)			;half-integer Js
strj=strarr(n_elements(jj))
if nintj gt 0 then strj(intj)=strtrim(string(jj(intj),'(i2)'),2)
if nfltj gt 0 then begin
  strj(fltj)=strtrim(string(2*jj(fltj),'(i2)'),2)
  strj(fltj)=strj(fltj)+'/2'
endif

;	any transitions in this range?
oo=where(abs(wvl) gt wmn and abs(wvl) le wmx,moo)
if moo eq 0 then begin				;...no
  message,'none in this wavelength range',/info & return,trans
endif
trans=strarr(2,moo) & config=trans		;...yes

for i=0L,moo-1L do begin			;for each transition...
  j=oo(i) & jl=llo(j)-1 & ju=lup(j)-1
  spin_lo=strtrim(string(ss(jl),'(i2)'),2)	;lower level...
  spin_up=strtrim(string(ss(ju),'(i2)'),2)	;upper level...
  ispd_lo=ll(jl) & ispd_up=ll(ju)
  spd_lo=spd([ispd_lo]) & spd_up=spd([ispd_up])	;S/P/D/...
  j_lo=strj(jl) & j_up=strj(ju)			;J value
  t_lo='' & t_up=''				;e configuration
  cl=strtrim(term(jl),2) & cu=strtrim(term(ju),2)
  bl=rstrpos(cl,' ') & if bl eq -1 then t_lo='' else t_lo=strmid(cl,0,bl)
  bu=rstrpos(cu,' ') & if bu eq -1 then t_up='' else t_up=strmid(cu,0,bu)
  ;
  drstr='' & if keyword_set(delec) then drstr=' (DR) '
  if keyword_set(tex) then begin
    trans(0,i)='$\^{'+spin_lo+'} '+spd_lo+'_{'+j_lo+'}'+drstr+'$'
    trans(1,i)='$\^{'+spin_up+'} '+spd_up+'_{'+j_up+'}$'
  endif else begin
    trans(0,i)=spin_lo+spd_lo+'_'+j_lo+drstr
    trans(1,i)=spin_up+spd_up+'_'+j_up
  endelse
  config(0,i)=t_lo & config(1,i)=t_up
endfor

return,trans
end
