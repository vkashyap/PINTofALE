;+
;poakeit.pro
;
;	given an ascii table of energies and line strengths,
;	and filenames of the RMF and ARF, generate a spectrum
;	and simulate counts from the spectrum
;
;usage
;	specfil=specfil
;	rmf=rmf & arf=arf
;	exptim=exptim
;	xlog=xlog & ylog=ylog
;	xrange=xrange & yrange=yrange
;	title=title
;	.run poakeit
;
;vinay k (Apr 2005)
;-

;	requires PINTofALE
defsysv,'!PoA',exists=iPoA
if iPoA eq 0 then stop,'PINTofALE not initialized; .run /data/fubar/SCAR/pro/scrypt/initale'

;	input list of energies
if not keyword_set(specfil) then begin
  specfil='' & read,prompt='filename containing line energies and strengths> ',specfil
  print,'NOTE: '+specfil+' is assumed to have two columns --'
  print,'       energy [keV] and line strength [ph/s/cm^2]'
  ok='ok' & nnrg=n_elements(nrg) & nspc=n_elements(spc)
  if nnrg eq 0 then ok='(NRG) line energies are undefined' else $
   if nspc eq 0 then ok='(SPC) line strengths are undefined' else $
    if nnrg ne nspc then ok='line energies (NRG) and strengths (SPC) incompatible'
  if ok ne 'ok' then message,ok
endif else begin
  nlin=wc(specfil) & var=fltarr(2,nlin)
  openr,usp,specfil,/get_lun & readf,usp,var & close,usp & free_lun,usp
  nrg=reform(var[0,*]) & spc=reform(var[1,*])
endelse

;	ARF
if n_tags(arstr) eq 0 then begin
  if not keyword_set(arf) then arf='NONE' & cc=''
  read,prompt='ARF ['+arf+']> ',cc
  if strtrim(cc,2) ne '' then arf=cc
  if strupcase(arf) ne 'NONE' then begin
    areff=rdarf(arf,arstr)
    arnrg=0.5*(arstr.ELO+arstr.EHI)
  endif
endif
nareff=n_elements(areff) & narnrg=n_elements(arnrg)
if nareff gt 0 and nareff ne narnrg then message,$
   'effective area (AREFF) and energies (ARNRG) are incompatible'
if narnrg eq 0 then begin
  arnrg=nrg & areff=0.*nrg+1.
endif

;	RMF
if n_tags(rmstr) eq 0 then begin
  if not keyword_set(rmf) then rmf='NONE' & cc=''
  read,prompt='RMF ['+rmf+']> ',cc
  if strtrim(cc,2) ne '' then rmf=cc
  if strupcase(rmf) ne 'NONE' then begin
    rmstr=rd_ogip_rmf(rmf,effar=effar)
    elo=rmstr.elo & ehi=rmstr.ehi & os=sort(elo)
    egrid=[elo[os],max(ehi)]
    emid=0.5*(egrid[1:*]+egrid) & dE=egrid[1:*]-egrid
  endif
endif
nelo=n_elements(elo) & nehi=n_elements(ehi)
if nelo gt 0 and nelo ne nehi then message,$
   'energy grid (ELO and EHI) are incompatible'
if nelo eq 0 then begin
  if not keyword_set(nbin) then nbin=1024L
  emin=min(nrg,max=emax)
  delE=(emin-emax)/(nbin+1.)
  elo=findgen(nbin)*delE+emin & ehi=elo+delE
  egrid=[elo,max(ehi)]
  emid=0.5*(egrid[1:*]+egrid) & dE=fltarr(nbin)+delE
endif

;	exptime
if not keyword_set(exptim) then exptim=1e-3 & cc=''
read,prompt='exposure time ['+strtrim(exptim,2)+' ks]> ',cc
if strtrim(cc,2) ne '' then exptim=double(cc)

;	define spectrum grid

;	interpolate effar to the spectrum grid
areff2=(interpol(areff,arnrg,emid)>0)<(max(areff))

;	rebin the input spectrum into the spectrum grid
spec=hastogram(nrg,egrid,wts=spc)	;[ph/s/cm^2]

;	convert
spec=spec*areff2	;[ph/s]
spec=spec*exptim*1e3	;[ph]

;	push through RMF
conv_rmf,egrid,spec,chan,spec,rmstr

;	make simulation
nspec=n_elements(spec) & simspec=lonarr(nspec)
for i=0L,nspec-1L do if spec[i] gt 0 then simspec[i]=randomu(seed,poisson=spec[i])

;	outputs
help,chan,spec,simspec
plot,chan,simspec,psym=10,xtitle='E [keV]',ytitle='counts',subtitle=specfil,$
	title=title,xrange=xrange,yrange=yrange,xlog=xlog,ylog=ylog
oplot,chan,spec,color=2
stample

print,"" & print,""
print,"plot,chan,simspec,psym=10,xtitle='E [keV]',ytitle='counts',subtitle=specfil,$"
print,"	title=title,xrange=xrange,yrange=yrange,xlog=xlog,ylog=ylog"
print,"oplot,chan,spec,color=2"
print,""
print,"To force rereading of '+specfil+', type specfil='' and rerun"
print,"To force rereading of RMF, type rmstr=0 and rerun"
print,"To force rereading of ARF, type arstr=0 and rerun"

end
