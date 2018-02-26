;+
;WRAP_SHOW_LINE
;	wrapper for SHOW_LINE
;
;vinay kashyap (Nov98)
;-

;	verify that input wavelengths and fluxes exist
if n_tags(linstr) eq 0 then begin
  ok='ok'
  if n_elements(linwvl) eq 0 then ok='Wavelengths missing' else $
   if n_elements(linZ) eq 0 then ok='Atomic numbers missing' else $
    if n_elements(ionion) eq 0 then ok='Ionic states missing'
  if ok ne 'ok' then stop,ok
endif else begin
  linwvl=linstr.WVL
  linZ=linstr.Z
  linion=linstr.ION
endelse
;
if n_elements(linflx) eq 0 then begin
  print,''
  print,'setting all fluxes in LINFLX to 0'
  print,''
  linflx=0.*linwvl
endif

;	which to show?
if not keyword_set(elems) then elems='all'
c1='Lines of which ion? (type "," separated list, e.g, "Fe10, C , Si V")'
print,c1 & read,prompt='<CR> to accept ['+elems+'] > ',c1
if strtrim(c1,2) eq '' then c1=elems else elems=c1
cc=str_sep(c1,',') & ncc=n_elements(cc)
;
i=0L & olin=-1L & molin=0L
while i lt ncc do begin
  c1=strtrim(cc(i),2)
  if strlowcase(c1) eq 'all' then begin
    molin=n_elements(linwvl) & olin=lindgen(molin) & i=ncc
  endif else begin
    if c1 ne '' then begin
      symb2zion,c1,zz,jon
      moo=0
      if zz ne 0 and jon ne 0 then oo=where(linZ eq zz and linion eq jon,moo)
      if zz ne 0 and jon eq 0 then oo=where(linZ eq zz,moo)
      if moo gt 0 then begin
	if molin gt 0 then olin=[olin,oo] else olin=oo
	molin=molin+moo
	print,strtrim(moo,2)+' lines in '+c1
      endif
    endif
  endelse
  i=i+1L
endwhile
;
if molin eq 0 then stop,'No lines of '+elem+' in current dataset!'
;
if n_elements(olin) ne molin then stop,'BUG!'
;
olin=olin(uniq(olin,sort(olin))) & molin=n_elements(olin)

if n_elements(lambda) lt 2 or n_elements(spec) lt 2 then begin
  print,'	SPEC(LAMBDA) not defined'
  c1='Continue? [y/n]' & print,c1 & c1=get_kbrd(1)
  if strlowcase(c1) eq 'n' then stop,'halting; type .CON to continue'
endif

print,''
c1='There are '+strtrim(molin,2)+' lines'
linflxmin=min(linflx(olin),max=linflxmax) & hlin=lindgen(molin)
if linflxmax gt linflxmin then begin
  c1=c1+', with fluxes ranging from '+strtrim(linflxmin,2)+' to '+$
	strtrim(linflxmax,2)
  print,c1
  if not keyword_set(linflxthr) then linflxthr=0 & c1=''
  if not keyword_set(thrlinflx) then thrlinflx=linflxthr
  c1='Set a threshold cut? (type fraction of MAX, actual threshold, or'
  print,c1
  read,prompt='-ve of actual threshold, or <CR> to accept '+$
	strtrim(linflxthr,2)+') ',c1
  if strtrim(c1,2) ne '' then linflxthr=float(c1)
  thrlinflx=linflxthr
  if thrlinflx lt 0 or thrlinflx ge 1 then thrlinflx=abs(linflxthr)/linflxmax
  hlin=where(linflx(olin) ge thrlinflx*linflxmax,mhlin)
  if mhlin eq 0 then stop,'No lines selected!'
endif else print,c1

show_line,linwvl(olin(hlin)),linflx(olin(hlin)),Z=linZ(olin(hlin)),$
  ion=linion(olin(hlin)),markw=markw,markf=markf,lambda=lambda,spec=spec,$
  order=order,oxr=oxr,oyr=oyr,markp=markp,markc=markc,marks=marks,$
  marko=marko,xtitle=xtitle,ytitle=ytitle

end
