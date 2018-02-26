pro wrt_ln_chianti,logP,logT,Z=Z,n_e=n_e,outdir=outdir, _extra=e
;+
;procedure	wrt_ln_chianti
;		writes out line intensities computed from a CHIANTI-type
;		database in a form that is easier to access.
;
;usage
;	wrt_ln_chianti,logP,logT,Z=Z,N_e=N_e,outdir=outdir,$
;	wmn=wmn,wmx=wmx,istate=istate,chidir=chidir
;
;parameters
;	logP	[INPUT; required] log10(Pressure [cm^-3 K])
;	logT	[INPUT] array of log10(Temperature [K]) at which to compute
;		line intensities.  default: logT=findgen(81)*0.05+4.
;		* best keep delta(logT) constant
;
;keywords
;	Z	[INPUT; default: ALL] element, ionic state (e.g.:'He','Fe XII')
;	n_e	[INPUT] if set, computes line emissivities at the
;		specified electron density [cm^-3]
;		* overrides LOGP
;		* appends a "D" to OUTDIR
;	outdir	[INPUT; default: emisschianti] directory tree under which to
;		store the output files
;		* WARNING: do not append a "/" (will crash if N_E is set)
;
;	_extra	[INPUT ONLY] allows passing defined keywords (WMN,WMX,N_E,
;		ISTATE,CHIDIR) to LINE_CHIANTI
;
;side effects
;	deposits (possibly large) files in OUTDIR
;	* ATOM_wvl	wavelengths of all the transitions
;	* ATOM_tem	temperatures
;	* ATOM_##.#	line intensities I(tem,wvl)@logP
;	* ATOM_ion	ionic states corresponding to wvl
;	* ATOM_jon	ionic states that matter to the ion balance
;	* ATOM_src	source of line information
;			2: CHIANTI
;	* ATOM_lvl	level designations of each transition
;	* ATOM_ecn	e-configurations of each level
;
;restrictions
;	requires CHIANTI-style database
;	requires LINE_CHIANTI and associated subroutines
;	input LOGP must be uniquely writeable in format f4.1 (e.g.,14.5)
;
;history
;	vinay kashyap (Nov96)
;	forced output filenames to not have spaces, changed output
;	  file names and format, esp. of e-configurations (VK; Dec96)
;	changed ATOM_src CHIANTI code from 0 to 2 (VK; Jan96)
;	added keyword N_E (VK; Jun97)
;	modified to account for dielectronic recombination lines in
;	 CHIANTIv3 (VK; OctMM)
;-

np=n_params()
if np eq 0 then begin
  print,'Usage: wrt_ln_chianti,logP,logT,Z=Z,N_e=N_e,outdir=outdir'
  print,'  writes out line intensities computed from a CHIANTI-type database'
  print,'  in a form that is easier to access -- see RD_LINE()'
  print,'  also accepts defined keywords WMN,WMX,ISTATE,CHIDIR for LINE_CHIANTI'
  return
endif

;	check inputs
;
if n_elements(logP) eq 0 then logP=15 & pres=10.^(logP)	;logP not defined
;
if n_elements(logT) eq 0 then logT=findgen(81)*0.05+4	;logT not defined
;
atom=['H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al','Si','P',$
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
if keyword_set(z) then atom=[strtrim(z,2)]		;element defined
nz=n_elements(atom)
;
if not keyword_set(outdir) then outdir='emisschianti'	;where to dump output
odir=strtrim(outdir,2) & if keyword_set(n_e) then odir=odir+'D'
;
fnum=logP(0) & if keyword_set(n_e) then fnum=alog10(n_e(0))
if fnum ge 10 or fnum le 0 then ffil=string(fnum,'(f4.1)') else $
	ffil='0'+string(fnum,'(f3.1)')
;
wfil='wvl' & tfil='tem' & ifil='ion' & jfil='jon' & sfil='src'
lfil='lvl' & efil='ecn'

;if OUTDIR doesn't exist, create it
;	clearly, this bit works only on UNIX
fil=findfile(odir,count=nfil)
if nfil eq 0 then spawn,'mkdir '+odir

for iz=0,nz-1 do begin				;for each element
  message,'	working on '+atom(iz),/info
  ;	get line cooling intensities [1e-30 erg cm^3/s]
  fx=line_chianti(atom(iz),pres,logT,ang,trans,ikey=ikey,jkey=jkey,$
	econf=econf,n_e=n_e,_extra=e)
  ;	if any lines were found...
  if ikey(0) ne 0 then begin
    fpre=strcompress(atom(iz),/remove_all)		;file prefixes
    message,'	writing to '+odir+'/'+fpre+'_*',/info
    ;	open files
    openw,uw,odir+'/'+fpre+'_'+wfil,/get_lun		;wavelengths
    openw,ut,odir+'/'+fpre+'_'+tfil,/get_lun		;temperatures
    openw,ui,odir+'/'+fpre+'_'+ifil,/get_lun		;ionic state
    openw,uj,odir+'/'+fpre+'_'+jfil,/get_lun		;ion balance state
    openw,uf,odir+'/'+fpre+'_'+ffil,/get_lun		;intensities
    openw,us,odir+'/'+fpre+'_'+sfil,/get_lun		;data source
    ;        odir+'/'+fpre+'_'+lfil			;level designations
    ;        odir+'/'+fpre+'_'+efil			;e configurations
    ;	write
    writeu,uw,n_elements(ang),ang
    writeu,ut,n_elements(logT),logT
    writeu,ui,n_elements(ikey),ikey
    writeu,uj,n_elements(jkey),jkey
    writeu,uf,n_elements(logT),n_elements(ang),fx
    writeu,us,n_elements(ikey),0*ikey+2
    tmp_trans=trans & save,file=odir+'/'+fpre+'_'+lfil,tmp_trans
    tmp_econf=econf & save,file=odir+'/'+fpre+'_'+efil,tmp_econf
    ;	close files
    close,uw & free_lun,uw & close,ut & free_lun,ut
    close,ui & free_lun,ui & close,uj & free_lun,uj
    close,us & free_lun,us & close,uf & free_lun,uf
  endif
endfor						;iz=0,nz-1

return
end
