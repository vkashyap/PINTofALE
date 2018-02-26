;+
;WRAP_STARFLUX.PRO
;	wrap routine for starflux.pro
;
;	basically makes sure that the line and continuum emissivities
;	are read in (and reads them in if not) and match the desired
;	wavelength range.
;
;vinay kashyap (1999May)
;-

;	make sure all necessary variables are defined
if not keyword_set(lindir) then begin
  lindir='$SPEX' & print,'setting LINDIR to: '+lindir
endif
if not keyword_set(chidir) then begin
  chidir='/data/fubar/SCAR/CHIANTI/'
  print,'CHIANTI installation is at: '+chidir
endif
if not keyword_set(chifil) then begin
  chifil='ioneq/arnaud_raymond_lmf.ioneq'
  print,'ion balance data in file: '+chidir+'/'+chifil
endif
if not keyword_set(condir) then begin
  condir='$CONT' & print,'setting CONDIR to: '+condir
endif
if not keyword_set(conroot) then begin
  conroot='cie' & print,'setting CONROOT to: '+conroot
endif
if not keyword_set(n_e) then begin
  n_e=1e10 & print,'setting electron density, N_E='+strtrim(n_e,2)
endif
if n_elements(wrange) ne 2 then begin
  if n_elements(wrange) eq 0 then cc='WRANGE not defined' else $
  	cc='WRANGE must be a 2-element vector'
  stop,cc
endif
if n_elements(abund) lt 30 then begin
  if not keyword_set(hint) then hint='Grevesse et al.'
  abund=getabund(hint)
  print,'Using abundances from '+hint
endif
if not keyword_set(NH) then begin
  NH=0.
  print,'setting H-column density, NH='+strtrim(NH,2)
endif else begin
  if not keyword_set(fH2) then fH2=0.26
  if not keyword_set(He1) then He1=0.1*NH
  if not keyword_set(HeII) then HeII=0.01*NH
  if not keyword_set(Fano) then Fano=0
endelse
if n_elements(VEM) eq 0 then begin
  VEM=1d20 & print,'setting volume emission measure, VEM='+strtrim(VEM,2)
endif
if n_elements(dist) eq 0 then begin
  dist=0. & print,'setting distance to star, DIST='+strtrim(dist,2)+' pc'
endif
if n_elements(radius) eq 0 then begin
  radius=0. & print,'setting radius of star, RADIUS='+strtrim(radius,2)+' cm'
endif
if n_elements(elem) eq 0 then begin
  elem='ALL' & print,'Elements to consider: '+elem
endif
if n_elements(effar) eq 0 or n_elements(wvlar) eq 0 then begin
  print,'WARNING: effective area is not defined' & help,EFFAR,WVLAR
  print,''
endif
wmin=min(wrange,max=wmax)
lnlst=strtrim(elem,2)+':'+strtrim(wmin,2)+'-'+strtrim(wmax,2)+':'+lindir

;	have any of the crucial settings changed since last run?
if keyword_set(olindir) then begin
  if lindir ne olindir then linstr=1
endif
if keyword_set(ocondir) then begin
  if condir ne ocondir then constr=1
endif
if n_elements(owrange) eq 2 then begin
  if owrange(0) ne wrange(0) or owrange(1) ne wrange(1) then begin
    linstr=1 & constr=1
  endif
endif
if keyword_set(olnlst) then begin
  ok='ok' & nlst=n_elements(lnlst)
  if nlst ne n_elements(olnlst) then ok='not ok' else $
   for i=0,nlst-1 do if lnlst(i) ne olnlst(i) then ok='not ok'
  if ok ne 'ok' then linstr=1
endif
if keyword_set(oconroot) then begin
  if conroot ne oconroot then constr=1
endif
if keyword_set(on_e) then begin
  if n_e ne on_e then begin
    linstr=1 & constr=1
  endif
endif
if keyword_set(ochifil) then begin
  if chifil ne ochifil then linstr=1
endif
if keyword_set(ochidir) then begin
  if chidir ne ochidir then linstr=1
endif
;	save new settings
olindir=lindir & ocondir=condir & owrange=wrange & olnlst=lnlst
oconroot=conroot & on_e=n_e & ochifil=chifil & ochidir=chidir

;	read in emissivity databases if needed
if n_tags(linstr) eq 0 then begin
  linstr=rd_list(lnlst,sep=':',/incieq,chidir=chidir,chifil=chifil,n_e=n_e)
endif
if n_tags(constr) eq 0 then begin
  ;	CANNOT HANDLE CHANGES IN ELEM .. YET
  constr=1 & tmp=rd_cont(conroot,n_e=n_e,wrange=wrange,dbdir=condir,$
	abund=abund,fcstr=constr)
endif

;	temperatures at which to compute output
if n_elements(tlog) eq 0 then tlog=linstr.LOGT

;	all set?
help,'Parameters in use:',$
  n_e,elem,lindir,chidir,chifil,linstr,condir,conroot,constr,wrange,tlog,$
  VEM,abund,effar,wvlar,dist,radius,NH,fH2,He1,HeII

;	compute fluxes
flux=starflux(linstr,constr,tlog=tlog,DEM=VEM,abund=abund,wrange=wrange,$
	NH=NH,fH2=fH2,He1=He1,HeII=HeII,effar=effar,wvlar=wvlar,$
	dist=dist,lflux=lflux,cflux=cflux)

;	report
print,'outputs are in variables FLUX(TLOG), LFLUX(TLOG), and CFLUX(TLOG)'

end
