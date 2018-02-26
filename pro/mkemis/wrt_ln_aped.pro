pro wrt_ln_aped,apedfil,logT,outdir=outdir,apeddir=apeddir,abund=abund,$
	_extra=e
;+
;procedure	wrt_ln_aped
;	writes out line intensities computed from an APED-type
;	database in a form that is easier to access.
;
;usage
;	wrt_ln_aped,apedfil,logT,outdir=outdir,apeddir=apeddir
;
;parameters
;	apedfil	[INPUT; required] FITS file, output of APEC
;		* passed w/o comment to RD_APED
;	logT	[INPUT] log_10(T[K]) at which to store the emissivities
;		* default is FINDGEN(81)*0.05+4
;
;keywords
;	outdir	[INPUT; default: aped] directory tree under which to
;		store the output files
;		* WARNING: do not append a "/" (will crash if N_E is set)
;	apeddir	[INPUT] directory containing APEDFIL
;		* passed w/o comment to RD_APED
;	abund	[INPUT] APED includes abundances, so "undo" by dividing
;		by the abundances
;		* default is to use Anders & Grevesse
;
;	_extra	[JUNK] here only to avoid crashing the program
;
;side effects
;	deposits (possibly large) files in OUTDIR
;	* ATOM_wvl	wavelengths of all the transitions
;	* ATOM_tem	temperatures
;	* ATOM_##.#	line intensities I(tem,wvl)@logP
;	* ATOM_ion	ionic states corresponding to wvl
;	* ATOM_src	source of line information
;			3: APED
;	* ATOM_lvl	level designations of each transition
;	* ATOM_ecn	e-configurations of each level
;
;restrictions
;	requires APED-style datafile
;	subroutines: RD_APED, KILROY, INICON
;
;history
;	vinay kashyap (AprMM)
;-

np=n_params()
if np eq 0 then begin
  print,'Usage: wrt_ln_aped,apedfil,logT,outdir=outdir,apeddir=apeddir'
  print,'  writes out line intensities computed from an APED-type database'
  print,'  in a form that is easier to access -- see RD_LINE()'
  return
endif

;	check LOGT
if n_elements(logT) eq 0 then logT=findgen(81)*0.05+4.

;	check OUTDIR
odir='apedD'	;where to dump output -- the default
if keyword_set(outdir) then odir=strtrim(outdir[0],2)
;	if OUTDIR doesn't exist, create it
;	clearly, this bit works only on UNIX
fil=findfile(odir,count=nfil)
if nfil eq 0 then spawn,'mkdir '+odir

;	check ABUND
if n_elements(abund) lt 30 then abnd=getabund('anders & grevesse') else $
	abnd=abund

;	read in the APED datafile
emiss=rd_aped(apedfil,wvls,temps,edens,Zs,ions,levs,$
	apeddir=apeddir,xtnames=xtnames)
;restore,'aped.sav',/verbose
if emiss[0] eq -1 then begin
  message,strtrim(apedfil[0],2)+': could not be read',/info
  return
endif

;	initialize
wfil='wvl' & tfil='tem' & ifil='ion' & sfil='src' & lfil='lvl' & efil='ecn'
nkTs=n_elements(temps) & neds=n_elements(edens) & nwvl=n_elements(wvls)
nT=n_elements(logT)
inicon,atom=atom
Tlogs=temps

;	convert [ph cm^3/s] -> [1e-23 ergs cm^3/s]
h=6.626176e-27 & c=2.9979e10	;[ergs s] & [cm/s]
nrg=h*c*1e8/abs(wvls)		;[ergs/ph]
ph2erg=nrg*1e23			;[1e-23 ergs/ph]

;	write out the files for each available atom
uZ=Zs[uniq(Zs,sort(Zs))] & nZ=n_elements(uZ)

for iz=0,nZ-1 do begin		;{for each element
  oz=where(Zs eq uz[iz],moz)
  if moz eq 0 then message,'BUG!'
  ww=wvls[oz] & ii=ions[oz] & ll=levs[*,oz] & ss=intarr(moz)+3

  ;	figure out the file name prefixes
  fpre=strtrim(atom[uz[iz]-1],2)
  ;	open files to store WVLS,TEMPS,IONS,SOURCE,LEVELS,ECONF
  openw,uw,odir+'/'+fpre+'_wvl',/get_lun,/swap_if_little_endian
  openw,ut,odir+'/'+fpre+'_tem',/get_lun,/swap_if_little_endian
  openw,ui,odir+'/'+fpre+'_ion',/get_lun,/swap_if_little_endian
  openw,us,odir+'/'+fpre+'_src',/get_lun,/swap_if_little_endian
  ;	openr,ul,odir+'/'+fpre+'_lvl',/get_lun,/swap_if_little_endian
  ;	openr,ue,odir+'/'+fpre+'_ecn',/get_lun,/swap_if_little_endian
  writeu,uw,moz,float(ww)
  writeu,ut,nT,float(logT)
  writeu,ui,moz,long(ii)
  writeu,us,moz,long(ss)
  tmp_trans=strtrim(ll,2) & save,file=odir+'/'+fpre+'_lvl',tmp_trans
  tmp_econf=strarr(2,moz) & save,file=odir+'/'+fpre+'_ecn',tmp_econf
  close,uw & free_lun,uw
  close,ut & free_lun,ut
  close,ui & free_lun,ui
  close,us & free_lun,us

  ;	write out the emissivities
  for iden=0,neds-1 do begin		;{for each density
    em=dblarr(nkTs,moz) & emis=dblarr(nT,nwvl)	;(from) -> (to)
    ;	get the (from)s (and convert to ergs while yer at it)
    for it=0,nkTs-1 do em[it,*]=emiss[oz,it,iden]*ph2erg
    ;	make the (to)s
    for iw=0L,moz-1L do begin		;{for each wavelength
      ee=reform(em[*,iw]) & eemax=max(ee)
      eee=(interpol(ee,Tlogs,logT) > 0) < eemax
      emis[*,iw]=eee[*]
    endfor				;IW=0,NWVL-1}

    ;	figure out file name suffixes
    fnum=alog10(edens[iden])
    if fnum ge 10 or fnum le 0 then ffil=string(fnum,'(f4.1)') else $
	ffil='0'+string(fnum,'(f3.1)')

    ;	and dump to files
    openw,uf,odir+'/'+fpre+'_'+ffil,/get_lun,/swap_if_little_endian
    writeu,uf,nT,moz,double(emis)
    close,uf & free_lun,uf

  endfor				;IDEN=0,NEDS-1}
endfor				;IZ=0,NZ-1}

return
end
