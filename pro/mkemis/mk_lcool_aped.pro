;+
;script	MK_LCOOL_APED.PRO
;	IDL program to set up the call to EXAPED and write out the
;	emissivities from the linelist.fits file
;
;usage
;	wmn=0. & wmx=3000.
;	atomdb='/data/fubar/SCAR/atomdb/atomdb_v2.0.2/'	;path to ATOMDB
;	apeddir='/data/fubar/SCAR/emissivity/apedD'	;output directory
;	.run mk_lcool_aped
;
;vinay kashyap (Jul02)
;-

;	initialize
if not keyword_set(wmn) then wmn=0.
if not keyword_set(wmx) then wmx=3000.
if not keyword_set(atomdb) then begin
  ivar=0 & defsysv,'!ATOMDB',exists=ivar
  if ivar eq 0 then ATOMDB='/data/fubar/SCAR/atomdb/atomdb_v2.0.2/' else $
    jnk=execute('ATOMDB=!ATOMDB')
endif
if not keyword_set(apeddir) then apeddir='/data/fubar/SCAR/emissivity/apedD'
if not keyword_set(verbose) then verbose=10
if n_elements(logT) eq 0 then logT=findgen(81)*0.05+4.

;	extract APED emissivities
exaped,'apec_v1.20',emis,Tlog,wvl,edens,Z,ion,desig,econf,jon,src,$
	/llist,atomdb=atomdb,logT=logT,verbose=verbose

;	make a structure out of it and write it out
nedens=n_elements(edens)
ow=where(abs(wvl) ge wmn and abs(wvl) le wmx,mow)
if mow eq 0 then message,'No lines selected'
for i=0,nedens-1 do begin
  linstr=mk_linstr(emis[*,ow,i],logT=logT,wvl=wvl[ow],Z=Z[ow],ion=ion[ow],$
	jon=jon[ow],desig=desig[*,ow],econf=econf[*,ow],src=src[ow],$
	verbose=verbose)
  wrt_ln_generic,linstr,apeddir,alog10(edens[i]),verbose=verbose
endfor

end
