;+
;script	MK_LCOOL_CHIANTI
;	IDL program to call WRT_LN_CHIANTI repeatedly for various pressures
;	and densities to generate a database of line intensities for all
;	available CHIANTI ions
;
;usage
;	set the variables WMN and WMX to fix wavelength range in [Ang]
;	set the variables NE_MIN, NE_MAX, D_NE [log10[cm^-3]] to set
;	range of densities, and P_MIN, P_MAX, D_P [log10[cm^3 K]] to
;	set range of pressures at which to compute the emissivities
;	set CHIDIR as the path to the CHIANTI installation
;	set OUTDIR as the path to where the outputs should be stored
;
;	if initale has been run, CHIDIR is set to !CHIDIR, and OUTDIR
;	is set to !TOPDIR/emissivity/chianti
;
;	invoke as
;	.run mk_lcool_chianti
;	
;
;vinay kashyap (Dec96)
;  updated to make it more robust to other people's setups (VK; Jul13)
;-

if not keyword_set(wmn) then wmn=0.
if not keyword_set(wmx) then wmx=3000.
if not keyword_set(chidir) then begin
  ivar=0 & defsysv,'!CHIDIR',exists=ivar
  if ivar eq 0 then chidir='/data/fubar/SCAR/CHIANTI/dbase/' else $
    jnk=execute('chidir=!CHIDIR')
endif
if not keyword_set(outdir) then begin
  ivar=0 & defsysv,'!TOPDIR',exists=ivar
  if ivar eq 0 then outdir='/data/fubar/SCAR/emissivity/chianti' else $
    jnk=execute("outdir=filepath('chianti',subdir='emissivity',root_dir="+!TOPDIR+")")
endif
if not keyword_set(ne_min) then ne_min=8.0
if not keyword_set(ne_max) then ne_max=15.0
if not keyword_set(d_ne) then d_ne=1.0
if not keyword_set(p_min) then p_min=13.0
if not keyword_set(p_max) then p_max=20.0
if not keyword_set(d_p) then d_p=1.0

;	this part computes and writes out the constant-density emissivities
;for logD=ne_min,ne_max,d_ne do begin
;  n_e = 10.^(logD)
;  print,n_e
;  wrt_ln_chianti,logP,outdir=outdir,chidir=chidir,wmn=wmn,wmx=wmx,n_e=n_e
;endfor

;	this part computes and writes out the constant-pressure emissivities
for logP=p_min,p_max,d_p do begin
  print,logP
  wrt_ln_chianti,logP,outdir=outdir,chidir=chidir,wmn=wmn,wmx=wmx
endfor

end
