function rd_cont,root,pres=pres,logP=logP,n_e=n_e,wrange=wrange,dbdir=dbdir,$
	abund=abund,metal=metal,wvl=wvl,logT=logT,fcstr=fcstr,verbose=verbose,$
	chdwvl=chdwvl,twoph=twoph,noff=noff,nofb=nofb,cocofil=cocofil,$
	xlogT=xlogT, _extra=e
;+
;function	rd_cont
;	routine to read in the output of WRT_CONT_* and return continuum
;	emissivities [1e-23 erg cm^3/s/Ang] for specified criteria
;
;syntax
;	emis=rd_cont(root,pres=pres,logP=logP,n_e=n_e,wrange=wrange,$
;	dbdir=dbdir,abund=abund,metal=metal,wvl=wvl,logT=logT,fcstr=fcstr,$
;	chdwvl=chdwvl,/twoph,/noff,/nofb,cocofil=cocofil,/xlogT)
;
;parameters
;	root	[INPUT; required] prefix of files containing the emissivities
;		* if set to 'chianti', will carry out a real-time computation
;		  of the continua using CHIANTI routines FREEFREE, FREEBOUND,
;		  and TWO_PHOTON
;		* if set to 'apec', will read in the emissivities from the
;		  APEC continuum database table, which includes both regular
;		  continuum and the so-called pseudo-continuum
;keywords
;	pres	[INPUT; default: 1e15] pressure [cm^-3 K]
;	logP	[INPUT; default: 15] log(pressure [cm^-3 K])
;		* LOGP takes precedence over PRES
;		* if no match is found, uses nearest available dataset
;	n_e	[INPUT] electron density [cm^-3].  if set,
;		* overrides LOGP
;		* appends a "D" to DBDIR
;	wrange	[INPUT; default: all] wavelength range
;		* if 2-element array, wrange==[wmin,wmax]
;		* if 1-element array or scalar, wrange==[wr(0),wr(0)]
;	dbdir	[INPUT; default: "cont"] directory to look for files in
;		* $<whatever> : look for environment variable.  if missing,
;		  following special case is hardcoded --
;		  -- $CONT : see variable CONT_DB
;		* ignored when ROOT='chianti' or ROOT='apec'
;	abund	[INPUT; default: anders & grevesse] abundances
;	metal	[INPUT; default: 1] metallicity
;		* assumes ABUND to be of metallicity 1
;		* if set, scales abundances from ABUND appropriately
;	wvl	[OUTPUT] bin beginning wavelengths and final ending
;		wavelength [Angstroms]
;	logT	[OUTPUT] log10(Temperatures [K])
;	fcstr	[OUTPUT] if specified in call on input, will contain all
;		the output variables
;			{CONT_INT,LOGT,WVL,REH,TKEV,EKEV,MIDWVL}
;		in that order in a structure.
;	chdwvl	[INPUT] bin size to use for ROOT='chianti'
;		* default is 0.005 Angstrom
;       twoph	[INPUT] set this to _include_ two-photon continuum
;       noff    [INPUT] set this to _exclude_ free-free continuum
;       nofb    [INPUT] set this to _exclude_ free-bound continuum
;		* TWOPH, NOFF, NOFB are ignored unless ROOT='chianti'
;	cocofil	[INPUT] full path to the location of the APEC
;		continuum emissivities file
;		* default is to use '$APED', which points to
;		  /data/atomdb/apec_coco.fits
;		* used only if CEROOT points to APEC
;	xlogT	[INPUT] if set, returns the emissivity in whatever
;		native logT grid as is in the databases
;		* the default is to rebin it with REBINX() and put
;		  it on the grid logT=findgen(81)*0.05+4
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to RD_IONEQ(),
;		FREE_FREE_CH, FREE_BOUND_CH, and TWO_PHOTON_CH.
;
;restrictions
;	* WRT_CONT_* must have been run first
;	* cannot handle logT arrays that change from file to file
;	* requires subroutines
;	  INICON
;	  GETABUND
;	  SYMB2ZION [LAT2ARAB]
;	  SETSYSVAL [LEGALVAR]
;         RD_IONEQ
;         FREE_FREE_CH  [ITOH_CH,SUTHERLAND_CH]
;         FREE_BOUND_CH
;         TWO_PHOTON_CH [RD_CHIANTI]
;	  RD_COCO
;	  GETPOADEF
;	  from CHIANTI: READ_IONEQ, READ_GFFGU, BILINEAR, READ_IP,
;	    READ_KLGFB, CONCAT_DIR, GET_IEQ, FREEBOUND_ION, ZION2FILENAME,
;	    FILE_EXIST, DATATYPE, FILE_STAT, DATA_CHK, SINCE_VERSION,
;	    READ_FBLVL, IS_DIR, EXIST, GET_DELIM, OS_FAMILY, TAG_EXIST,
;	    IDL_RELEASE, VERNER_XS, KARZAS_XS
;	  
;
;history
;	vinay kashyap (Dec97; modified from RD_LINE)
;	allowed dbdir to be "$CONT" or something that may be read in
;	  via an environmental variable (VK; Dec98)
;	added missing file close statement (VK; 1999Feb)
;	bug was preventing reading nearest available database (VK; Mar99)
;	bug in last bug correction, wasn't reading element files (VK; Apr99)
;	allowed FCSTR output if set in command line (VK; AugMM)
;	changes suggested by Antonio Maggio: *_DB defaults (VK; NovMM)
;	added byte-order check; added calls to SWAP_ENDIAN() as suggested
;	  by Antonio Maggio; changed CONT_DB to CONT_CDB environment
;	  variable; added call to SETSYSVAL; now deciphers default
;	  emissivity directory location from where this file is located
;	  (VK; FebMMI)
;	added keyword VERBOSE (VK; Aug01)
;	corrected bug where first bin was being chopped off (VK; Jan04)
;	check to see whether input abundances are already ratios (VK; Apr04)
;	added option to use chianti routines to compute continuum, added
;	  TWO_PHOT keyword to toggle inclusion of two photon continuum when
;	  calculating via CHIANTI to save time (LL; Feb05) 
;	added NOFF and NOFB keywords to toggle inclusion of CHIANTI free-free
;	  and free-bound continuum, allowed system variable !XUVTOP to be set
;	  if not present (LL; Jun05) 
;	changed keyword TWO_PHOT to TWOPH, added keyword CHDWVL, and
;	  cleaned up interface and internals (VK; Jun05)
;	bug fix: $CONT/cont was missing the "/" on some systems (VK; May06)
;	added option to read in APEC_COCO files and added keyword XLOGT
;	  (VK; Apr09)
;	fixed complication with abund being input as ratio (VK; Nov13)
;	added call to GETPOADEF (VK; Aug15)
;-

;	usage
if n_elements(root) ne 1 then begin
  print,'Usage: fc=rd_cont(root,pres=pres,logP=logP,n_e=n_e,wrange=wrange,$'
  print,'       dbdir=dbdir,abund=abund,metal=metal,wvl=wvl,logT=logT,$'
  print,'       chdwvl=chdwvl,fcstr=fcstr,/twoph,/noff,/nofb,cocofil=cocofil,$
  print,'       /xlogT)'
  print,'  read in (or compute) continuum spectra'
  return,dblarr(1,1)-1
endif

;	initialize
elem=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']				;elements from 1-30
;	verbosity
vv=0
ivar=0 & defsysv,'!VERBOSE',exists=ivar	;if !VERBOSE exists
if ivar ne 0 then setsysval,'VERBOSE',vv,/getval
if n_elements(verbose) gt 0 then vv=long(verbose[0])
;	databases
;zTOPDIR=filepath('',root_dir='data',subdirectory=['fubar','SCAR'])
;zTOPDIR='/data/fubar/SCAR/'
zTOPDIR=getpoadef()
ivar=0 & defsysv,'!TOPDIR',exists=ivar	;if !TOPDIR exists
if ivar ne 0 then setsysval,'TOPDIR',zTOPDIR,/getval else begin
  tmp=routine_info('RD_CONT',/source,/functions)
  if n_tags(tmp) ne 0 then begin
    scardir=tmp.PATH
    jvar=strpos(scardir,'pro')
    if jvar ge 0 then begin
      if strpos(scardir,'rd_cont.pro',jvar) ge 0 then $
	zTOPDIR=strmid(scardir,0,jvar)
  		;we know where RD_CONT.PRO is
    endif
  endif
  if vv ge 10 then message,'Default emissivity directory is: '+zTOPDIR,/info
endelse
CONT_DB=getenv('CONT_CDB')
if CONT_DB eq '' then $
  CONT_DB=filepath('cont',root_dir=filepath('emissivity',root_dir=zTOPDIR))
	;CONT_DB=zTOPDIR+filepath('',root_dir='emissivity')+'cont'

;	check keywords
if keyword_set(pres) then pr=pres else pr=1e15		;pressure
if keyword_set(logP) then pr=10.^(logP)			;log10(pressure)
logP=alog10(pr)
case n_elements(wrange) of				;wavelength range
  1: wr=[wrange(0),wrange(0)]			;pick only one
  2: wr=wrange					;range
  else: wr=[0.,10000.]				;all possible
endcase
if not keyword_set(dbdir) then dbdir='cont'		;directory with files
;	remove trailing '/' and/or trailing 'D'
c=dbdir(0)
if strmid(c,strlen(c)-1,1) eq '/' then c=strmid(c,0,strlen(c)-1)
if strmid(c,strlen(c)-1,1) eq 'D' then c=strmid(c,0,strlen(c)-1)
if strmid(c,0,1) eq '$' then begin
  if c eq '$CONT' then c=CONT_DB else begin
    cdb=getenv(strmid(c,1,strlen(c)-1))	;environment variable
    if cdb eq '' then begin		;($WHATEVER was not set
      c0=strupcase(c) & c1=CONT_DB	;default
      if strpos(c0,'CONT',1) gt 0 then c1=CONT_DB
      message,c+': not defined; using '+c1,/info
      c=c1
    endif				;CDB='')
  endelse
endif
dbasedir=c
;	append 'D' if necessary
if keyword_set(n_e) then dbasedir=c+'D'
;if keyword_set(n_e) then dbasedir=dbdir+'D' else dbasedir=dbdir
if vv gt 0 then print,'reading from '+dbasedir
defabu=getabund('anders & grevesse')
if n_elements(abund) lt 30 then abnd=defabu else abnd=abund	;abundances
if n_elements(metal) ne 1 then met=1. else met=metal(0)	;metallicity
abnd(2:*)=abnd(2:*)*met				;include metallicity
;	check C, N, O, and Fe abundances to figure out whether
;	abundances already appear to be ratios
if abnd(6-1) gt 0.05 or $
  abnd(7-1) gt 0.01 or $
  abnd(8-1) gt 0.01 or $
  abnd(26-1) gt 0.005 then begin
  message,'input abundances appear to be ratios',/informational
  factor=abnd
  abnd=factor*defabu	;bug fix 2014-nov-13
endif else factor=abnd/defabu

case root of				;{special cases

  'apec': begin
    if not keyword_set(cocofil) then cocofil='$APED'
    jnk=rd_coco(cocofil,wrange=wrange,abund=abnd,$
    cemis=emis,wvl=wvl,wmid=wmid,logT=logT,verbose=vv, _extra=e)
    nw=n_elements(wvl)
    if not keyword_set(xlogT) then begin
      oldlogT=logT & logT=findgen(81)*0.05+4.
      emis=(rebinx(emis,oldlogT,logT)>0)
    endif
  end

  'chianti': begin			;(compute with CHIANTI

    xivar = 0 & defsysv,'!xuvtop',exists=xivar ;if !XUVTOP exists 
    if xivar eq 0 then begin 
      defsysv, '!chidir', exists=chivar 
      if chivar  eq 1 then defsysv,'!xuvtop',!chidir else $ 
        defsysv,'!xuvtop',filepath('',root_dir=zTOPDIR,subdirectory=['CHIANTI','dbase'])
    endif

    dwvl = 0.005 ; set the default wavelength binning
    if keyword_set(chdwvl) then dwvl=chdwvl[0]
    if dwvl le 0 then dwvl=0.005
    nw = (wr(1)-wr(0))/dwvl ; get number of wavelength bins
    wvl = findgen(nw+1)*dwvl + wr(0) ; create wavelength grid 
    wmid=0.5*(wvl(1:*)+wvl)

    ivar = 0 & defsysv,'!LOGT',exists=ivar 
    if ivar ne 0 then logt = !logt else logt=findgen(81)*0.05+4.0
    nt=n_elements(logt)

    ;	define outputs (note that actual output will be the
    ;	transpose of these -- i.e., fltarr(nt,nw))
    ffint=fltarr(nw,nt)
    fbint=fltarr(nw,nt)
    tpint=fltarr(nw,nt)

    ioneq = rd_ioneq(findgen(31)+1,logt,_extra=e) 
    ;   accnt for  ion and z arrays  exchanged btwn poa and chi in ioneq
    chioneq = transpose(ioneq,[0,2,1])
    if keyword_set(noff) and keyword_set(nofb) and not keyword_set(twoph) then begin 
      message, $ 
      'No continuum requested, please check keywords NOFF,NOFB, and TWOPH',/info 
      return, ffint+fbint+tpint
    endif 
    if not keyword_set(noff) then $    
      freefree_ch,10.^logt,wmid,ffint,ioneq=chioneq,ieq_logt=logt,$
                  abund=[abnd,0],/no_setup,_extra=e
    if not keyword_set(nofb) then $
      freebound_ch,10.^logt,wmid,fbint,ioneq=chioneq,ieq_logt=logt,$
                  abund=[abnd,0],/no_setup,_extra=e
    if keyword_set(twoph) then begin 
      if keyword_set(n_e) then begin 
        two_photon_ch,10.^logt,wmid,tpint,ioneq=chioneq,ieq_logt=logt,$ 
             abund=[abnd,0],edensity=n_e,/no_setup,_extra=e 
      endif else begin 
        message,'WARNING: N_E not set, TWO_PHOTON not included.',/info
      endelse
    endif
    emis = transpose(tpint+fbint+ffint)
    inicon, atom=atm
    atmn = findgen(30)+1
    ; unit conversions 1e-40 ergs/s/sr/ang -> 1e-23 ergs/s/ang
    emis = 4.*!pi*emis/1d17   ;4.*!pi*emis/1d17
  end					;CHIANTI)

  else: begin 				;(read from database

    ;	output arrays
    wvl=[-1.] & logT=wvl

    ;	find the right suffixes for pressure/density
    fnum=logp(0) & if keyword_set(n_e) then fnum=alog10(n_e(0))
    if fnum ge 0 and fnum lt 10 then fsfx='0'+string(fnum,'(f3.1)')
    if fnum lt 0 or fnum ge 10 then fsfx=string(fnum,'(f4.1)')

    ;	find out if there are any data files at all
    nfil=n_elements(fils) & fils=findfile(dbasedir+'/'+root+'*_'+fsfx,count=nfil)

    if nfil lt 4 then begin     ;(no files present?
    ;	find the nearest suffix
      filz=findfile(dbasedir+'/'+root+'_*.?',count=kfil)
      if kfil eq 0 then begin
        message,'	no spectra available for '+dbasedir+'/'+root+'_*',/info
        return,dblarr(1,1)-1.
      endif
      sufx=strarr(kfil)
      for i=0,kfil-1 do sufx(i)=strmid(filz(i),strlen(dbasedir+'/'+root+'_'),4)
      sufx=float(sufx) & tmp=min(abs(sufx-fnum),jj)
      message,'dataset with suffix '+fsfx+' not available.  using '+$
        string(sufx(jj),'(f4.1)')+' instead',/info
      fnum=sufx(jj)		;reset to nearest available
      if fnum ge 0 and fnum lt 10 then fsfx='0'+string(fnum,'(f3.1)')
      if fnum lt 0 or fnum ge 10 then fsfx=string(fnum,'(f4.1)')
      ;
      ;	ok, and files with this new suffix do exist..
      fils=findfile(dbasedir+'/'+root+'*_'+fsfx,count=nfil)
      if nfil lt 4 then message,'BUG!'
    endif                       ;NFIL<4)

    ;	need to figure out whether to swap byte order or not
    if n_elements(set_byte_swap) eq 0 then begin
      openr,uw,dbasedir+'/'+root+'_wvl',/get_lun
      nw=0L & readu,uw,nw
      if nw ge 16777216 or nw lt 0 then set_byte_swap=1 else set_byte_swap=0
      if float(strmid(!version.RELEASE,0,3)) lt 5.3 then $
        set_byte_swap=(-1)*set_byte_swap
      close,uw & free_lun,uw
    ;
    ;Antonio Maggio points out that keyword /SWAP_ENDIAN is not valid
    ;for earlier OPENRs, so the *function* SWAP_ENDIAN() must be used.
    ;however, SWAP_ENDIAN() is very slow for large double precision
    ;arrays, so we will do that only for older IDL versions, hence the
    ;multiplication by -1 above.
    ;
    endif

    ;	first read the wavelength array...
    wfil=dbasedir+'/'+root+'_wvl'
    if set_byte_swap gt 0 then openr,uw,wfil,/get_lun,/swap_endian else $
      openr,uw,wfil,/get_lun
    nw=0L & readu,uw,nw
    if set_byte_swap lt 0 then nw=swap_endian(nw)
    wvl=dblarr(nw) & readu,uw,wvl
    if set_byte_swap lt 0 then wvl=swap_endian(wvl)
    close,uw & free_lun,uw

    ;	... and the temperatures array ...
    tfil=dbasedir+'/'+root+'_tmp'
    if set_byte_swap gt 0 then openr,ut,tfil,/get_lun,/swap_endian else $
      openr,ut,tfil,/get_lun
    nt=0L & readu,ut,nt
    if set_byte_swap lt 0 then nt=swap_endian(nt)
    logT=fltarr(nt) & readu,ut,logT
    if set_byte_swap lt 0 then logT=swap_endian(logT)
    close,ut & free_lun,ut

    ;	... and the "parameter history" array ...
    sfil=dbasedir+'/'+root+'_src'
    if set_byte_swap gt 0 then openr,us,sfil,/get_lun,/swap_endian else $
      openr,us,sfil,/get_lun
    c1='' & src=['']
    while not eof(us) do begin
      readf,us,c1
      if set_byte_swap lt 0 then c1=swap_endian(c1) & src=[src,c1]
    endwhile & if n_elements(src) gt 1 then src=src(1:*)
    close,us & free_lun,us

    ;	... and the base emissivities.
    efil=dbasedir+'/'+root+'_'+fsfx
    if set_byte_swap gt 0 then openr,ue,efil,/get_lun,/swap_endian else $
      openr,ue,efil,/get_lun
    mt=nt & mw=nw & emis=dblarr(nt,nw-1)
    readu,ue,mt & readu,ue,mw & readu,ue,emis
    if set_byte_swap lt 0 then mt=swap_endian(mt)
    if set_byte_swap lt 0 then mw=swap_endian(mw)
    if set_byte_swap lt 0 then emis=swap_endian(emis)
    close,ue & free_lun,ue

    ;	and finally look through the element list
    for i=0,n_elements(elem)-1 do begin ;{the rootZ_##.# files
      atm=elem(i) & efil=dbasedir+'/'+root+atm+'_'+fsfx
      oo=where(fils eq efil,moo)
      if moo gt 0 then begin  ;(emissivities exist
        if set_byte_swap gt 0 then openr,ue,efil,/get_lun,/swap_endian else $
          openr,ue,efil,/get_lun
        emZ=dblarr(nt,nw-1)
        readu,ue,mt,mw,emZ
        if set_byte_swap lt 0 then mt=swap_endian(mt)
        if set_byte_swap lt 0 then mw=swap_endian(mw)
        if set_byte_swap lt 0 then emZ=swap_endian(emZ)
        close,ue & free_lun,ue
        emis=emis+emZ*factor(i) ;add 'em up!
        oi=where(wvl gt 20 and wvl lt 50)
      ;	or should it be		emis=emis+emZ*abnd(i)		??
                            ;print,efil,total(emZ),factor(i)
      endif                   ;MOO>0)
    endfor                      ;I=0,N_ELEMENTS(ELEM)-1}

  end				;read from database)

endcase				;ROOT}

;	select in wavelength range
;ww=wvl(0:nw-1-1)		;use bin-beginning values
ww=wvl(0:nw-1)		;use bin-beginning values
oo=where(ww ge wr(0) and ww le wr(1),moo)
if moo gt 0 then begin
  if oo[0] ne 0L then oo=[oo[0]-1L,oo]
  wvl=wvl([oo,oo(moo-1)+1]) & emis=emis(*,oo)
endif else begin
  wvl=[-1.] & emis=dblarr(nt,1)-1.
  message,'	No bins in this wavelength range',/info
  return,emis
endelse

;	get TkeV and EkeV
TkeV=10.^(logT)*1.380662e-16/1.6021892e-9		;[K]*[erg/K]/[erg/keV]
EkeV=12.3985/wvl					;[Ang]->[keV]

;	get mid-bin values
wmid=0.5*(wvl(1:*)+wvl)

if arg_present(fcstr) then begin
  fcstr=create_struct('CONT_INT',emis,'LOGT',logt,'WVL',wvl,$
	'TkeV',TkeV,'EkeV',EkeV,'midWVL',wmid)
endif

return,emis
end
