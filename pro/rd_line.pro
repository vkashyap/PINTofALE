function rd_line,atom,pres=pres,logP=logP,n_e=n_e,wrange=wrange,dbdir=dbdir,$
	wvl=wvl,logT=logT,z=z,ion=ion,jon=jon,src=src,desig=desig,econf=econf,$
	fstr=fstr,allah=allah,help=help,verbose=verbose, _extra=e
;+
;function	rd_line
;	supremely flexible routine to read in the output of WRT_LINE
;	and return line cooling emissivities [1e-23 erg cm^3/s] for
;	specified criteria in an array of size [NLOGT,NWVL]
;
;syntax
;	emis=rd_line(atom,pres=pres,logP=logP,n_e=n_e,wrange=wrange,$
;	dbdir=dbdir,wvl=wvl,logT=logT,z=z,ion=ion,jon=jon,src=src,$
;	desig=desig,econf=econf,fstr=fstr,/allah,/help,verbose=verbose)
;
;parameters
;	atom	[INPUT] element whose line intensities are to be read.
;		* can specify ionic state (e.g., 'FeXX')
;		* may be an array (e.g., ['He', 'c 5', 'N V', 's_8'])
;		* default -- ALL
;
;keywords
;	pres	[INPUT; default: 1e15] pressure [cm^-3 K]
;	logP	[INPUT; default: 15] log(pressure [cm^-3 K])
;		* LOGP takes precedence over PRES
;		* if no match is found, uses nearest available dataset
;	n_e	[INPUT] electron density [cm^-3].  if set,
;		* overrides LOGP
;		* appends a "D" to DBDIR
;		* NOTE: PRES, LOGP, and N_E are initalized in the
;		  system variables !GASPR, !LOGPR and !EDENS.
;		  HOWEVER, these system variables are *not* used as
;		  the defaults; these keywords must be explicitly set
;		  in order to take effect.
;	wrange	[INPUT; default: all] wavelength range
;		* if 2-element array, wrange==[wmin,wmax]
;		* if 1-element array or scalar, wrange==[wr[0],wr[0]]
;	dbdir	[INPUT; default: !LDBDIR; hard-coded default: "$CHIANTI"]
;		line database directory
;		* if an array, only the first element is used
;		* WARNING: if N_E is set, a "D" is appended at the end (but
;		  after stripping any pre-existing "D"s, mind!)
;		* $<whatever> : look for environment variable.
;		  some special cases are:
;		  -- $SPEX : from environment variable SPEX_LDB
;		  -- $CHIANTI : from environment variable CHIANTI_LDB
;		  -- $APED : from environment variable APED_LDB
;		  -- $SCAR : from environment variable SCAR_LDB
;		  if the environmental variable is undefined, then some
;		  special cases are hardcoded:
;		  -- SPEX : TOPDIR/emissivity/spex
;		  -- CHIANTI : TOPDIR/emissivity/chianti
;		  -- APED : TOPDIR/emissivity/aped
;		  -- SCAR : TOPDIR/emissivity/xinti
;	wvl	[OUTPUT] wavelengths [Angstroms]
;	logT	[OUTPUT] log10(Temperatures [K])
;	z	[OUTPUT] atomic numbers of elements generating given WVL
;	ion	[OUTPUT] ionic state of element generating given WVL
;	jon	[OUTPUT] the ionic state that matters, as far as ion-balance
;		is concerned
;	src	[OUTPUT] reference to source of line info
;		* 1: SPEX, 2: CHIANTI, 3: APED
;	desig	[I/O] level designations for lower & upper level
;		for each transition
;		* read in only if DESIG is present in call
;	econf	[I/O] e-configurations for lower & upper level
;		* read in only if ECONF is present in call
;	fstr	[OUTPUT] if present in the calling sequence, will contain all
;		the output variables
;			LINE_INT,LOGT,WVL,Z,ION,DESIG,ECONF,SRC,JON
;		in that order in a structure.
;	allah	[INPUT] if set, assumes all files are available in correct
;		forms, and does no checking.
;		* quicker, but UNFORGIVING
;		* faith is not an option
;	help	[INPUT] if set, prints out the calling sequence and quits
;	verbose	[INPUT] higher the number, greater the chatter
;		* default is !VERBOSE
;	_extra	[INPUT] junk -- here only to avoid crashing the program
;
;restrictions
;	* WRT_LINE must have been run first
;	* requires subroutines
;	  SYMB2ZION
;	  LAT2ARAB
;	  INICON
;	  SETSYSVAL
;	  IS_KEYWORD_SET
;	  GETPOADEF
;	* cannot (as yet) handle logT arrays that change from atom to atom
;
;usage examples
;	* f=rd_line(atom,pres=pres,logP=logP,wrange=wrange,dbdir=dbdir,$
;	    wvl=wvl,logT=logT,z=z,ion=ion,/allah)
;	* f=rd_line('He',logP=20,wr=[200.,500.],wvl=wvl,logT=logT,desig=desig,$
;	  econf=econf,fstr=fstr,/allah)
;	* f=rd_line(wr=16.776,wvl=wvl,logT=logT,z=z,ion=ion,desig=desig)
;	* econf=1 & f=rd_line(wvl=wvl,logt=logt,z=z,ion=ion,econf=econf,/allah)
;	  oo=sort(abs(wvl)) & surface,alog10(f[*,oo]+1),logt,abs(wvl[oo]),/ylog
;
;history
;	vinay kashyap (Nov96)
;	added keywords ECONF and FSTR, changed I/O (VK; Dec96)
;	corrected bug in reading in SRC (VK; Jan96)
;	added keyword N_E (VK; Jun97)
;	corrected bug that was returning nonsense labels when Z+ION was
;	  specified (VK; Jul97)
;	now handles trailing "/"s and/or "D"s correctly; accepts "0"
;	  as synonym for all atoms (VK; Jun98)
;	accepts '$SPEX', '$CHIANTI', and '$SCAR' as hardcoded default
;	  values for DBDIR (VK; Oct98)
;	accepts DBDIR defined via environmental variable; trapped DBDIR
;	  being input as vector; added keyword HELP (VK; Nov98)
;	"All" is accepted as an atom (VK; Dec98)
;	"_"s in DBDIR causing incorrect selections: corrected
;	  added keyword JON, added read to ATOM_jon (VK; 1999May)
;	added keyword VERBOSE, allowed input of the form "ALL ion" (VK; MarMM)
;	added /SWAP_IF_LITTLE_ENDIAN to OPENR calls (VK; AprMM)
;	added $APE to list of automatic *_DB (VK; AprMM)
;	changed to full compatibility with IDL5+ ("()"->"[]"); allowed
;	  read in of DESIG, ECONF and output of FSTR if present in call;
;	  bug correction if ATOM goes past Zn (VK; AugMM)
;	CHIANTI 3.0 has occassional NaNs in emissivity -- put in hack to
;	  avoid them (VK; SepMM)
;	changes suggested by A.Maggio: strmid params, *_DB defaults
;	  (VK; NovMM)
;	changed *_DB to *_LDB environment variables; if ALLAH was set,
;	  DBDIR was not being prepended (VK; DecMM)
;	added call to SETSYSVAL; now deciphers default emissivity directory
;	  location from where this file is located; added byte-order
;	  check; added calls to SWAP_ENDIAN() as suggested by Antonio Maggio
;	  (VK; FebMMI)
;	changed defaults for DBDIR, VERBOSE to read from !LDBDIR, !VERBOSE
;	  (VK; Jul01)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	added call to GETPOADEF (VK; 05Aug2015)
;	bug fix: compilation error (VK; 10Sep2015)
;-

np=n_params(0)
if np eq 0 or keyword_set(help) then begin
  print,'Usage: f=rd_line(atom,pres=pres,logP=logP,n_e=n_e,wrange=wrange,$'
  print,'  dbdir=dbdir,wvl=wvl,logT=logT,z=z,ion=ion,jon=jon,src=src,$'
  print,'  desig=desig,econf=econf,fstr=fstr,/allah,/help,verbose=verbose)'
  if np eq 0 and keyword_set(help) then return,create_struct('LINE_INT',[-1],$
	'LOGT',0.,'WVL',0.,'Z',0,'ION',0,$
	'DESIG','','CONFIG','','SRC',0,'JON',0)
endif

;	initialize

;	verbosity
v=0
ivar=0 & defsysv,'!VERBOSE',exists=ivar	;if !VERBOSE exists
if ivar ne 0 then setsysval,'VERBOSE',v,/getval
if n_elements(verbose) gt 0 then v=long(verbose[0])

;elements
;
inicon,atom=atomsymb
atomsymb=['X',atomsymb]
elm=atomsymb[0:30]
;elm = [	'X','H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
;	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
;	'Ni','Cu','Zn']
if keyword_set(allah) and not keyword_set(atom) then atom=['He','C',$
	'N','O','Ne','Mg','Al','Si','S','Ar','Ca','Fe','Ni']
if keyword_set(atom) then elem=[strtrim(atom,2)] else elem=elm[1:*]
;	check whether input ATOM is saying "all"
iall=strpos(strlowcase(elem[0]),'all',0)
if iall ge 0 then begin
  ;tmp=strmid(elem[0],iall+3) ;AM pointed out this is illegal in IDL5-
  tmp=strmid(elem[0],iall+3,(strlen(elem[0])-(iall+3)))
  elem=elm[1:*]+tmp
endif
nz=n_elements(elem)

;databases
;
zTOPDIR=getpoadef()
;zTOPDIR='/data/fubar/SCAR/'
ivar=0 & defsysv,'!TOPDIR',exists=ivar	;if !TOPDIR exists
if ivar ne 0 then setsysval,'TOPDIR',zTOPDIR,/getval else begin
  ;tmp=routine_info('RD_LINE',/source,/functions) & scardir=tmp.PATH
  ;jvar=strpos(scardir,'pro/rd_line.pro')	;we know where RD_LINE.PRO is
  ;zTOPDIR=strmid(scardir,0,jvar-1)
  if v ge 10 then message,'Default emissivity directory is: '+zTOPDIR,/info
endelse
SPEX_LDB=getenv('SPEX_LDB') & if SPEX_LDB eq '' then $
	SPEX_LDB=zTOPDIR+'/emissivity/spex'
CHIANTI_LDB=getenv('CHIANTI_LDB') & if CHIANTI_LDB eq '' then $
	CHIANTI_LDB=zTOPDIR+'/emissivity/chianti'
SCAR_LDB=getenv('SCAR_LDB') & if SCAR_LDB eq '' then $
	SCAR_LDB=zTOPDIR+'/emissivity/xinti'
APED_LDB=getenv('APED_LDB') & if APED_LDB eq '' then $
	APED_LDB=zTOPDIR+'/emissivity/aped'

;	check keywords
if keyword_set(pres) then pr=pres else pr=1e15		;pressure
if keyword_set(logP) then pr=10.^(logP) 		;log10(pressure)
logP=alog10(pr)
case n_elements(wrange) of				;wavelength range
  1: wrange=[wrange[0],wrange[0]]		;pick only one
  2: 						;range
  else: wrange=[-1.,1e10]			;all possible
endcase
;zLDBDIR='$CHIANTI'	;hard-coded default
;ivar=0 & defsysv,'!LDBDIR',exists=ivar	;if !LDBDIR exists
;if ivar ne 0 then setsysval,'LDBDIR',zLDBDIR,/getval
zLDBDIR=getpoadef('LDBDIR')	;hard-coded default
if not keyword_set(dbdir) then dbdir=zLDBDIR
;	remove trailing '/' and/or trailing 'D'
c=dbdir[0]
if strmid(c,strlen(c)-1,1) eq '/' then c=strmid(c,0,strlen(c)-1)
if strmid(c,strlen(c)-1,1) eq 'D' then c=strmid(c,0,strlen(c)-1)
;	translate $SPEX, $CHIANTI, $SCAR
if strmid(c,0,1) eq '$' then begin
  if c eq '$SPEX' then c=SPEX_LDB else $
   if c eq '$CHIANTI' then c=CHIANTI_LDB else $
    if c eq '$APE' or c eq '$APED' then c=APED_LDB else $
     if c eq '$SCAR' then c=SCAR_LDB else begin
       cdb=getenv(strmid(c,1,strlen(c)-1))	;environment variable
       if cdb eq '' then begin		;($WHATEVER was not set
	 c0=strupcase(c) & c1=SCAR_LDB	;default
	 if strpos(c0,'SPEX',1) gt 0 then c1=SPEX_LDB else $
	   if strpos(c0,'CHIANTI',1) gt 0 then c1=CHIANTI_LDB else $
	    if strpos(c0,'APE',1) gt 0 then c1=APED_LDB
	 message,c+': not defined; using '+c1,/info
	 c=c1
       endif				;CDB='')
     endelse
endif
dbasedir=c
;	append 'D' if necessary
if keyword_set(n_e) then dbasedir=c+'D'
;	waste time or not?
if keyword_set(desig) or arg_present(desig) then rdlvl=1 else rdlvl=0
if keyword_set(econf) or arg_present(econf) then rdecn=1 else rdecn=0

;	chatter
if v gt 1 then message,'Reading from database '+dbasedir,/info

;	output arrays
wvl=[-1.] & Z=[0L] & ion=[0L] & jon=[0L] & src=[0L]
desig=[''] & econf=desig & logT=wvl

for iz=0,nz-1 do begin					;{for each element

  ;	chatter
  if v gt 0 then begin
    print,'' & message,'	Working on '+elem[iz],/info
  endif

  ;	first make sure we have the right prefix
  zz=elem[iz] & symb2zion,zz,j,kion,verbose=v & zz=atomsymb[j]

  ;	and the right suffixes for logP/logNe
  fnum=logp[0] & if keyword_set(n_e) then fnum=alog10(n_e[0])
  if fnum ge 10 or fnum le 0 then fsfx=string(fnum,'(f4.1)') else $
	fsfx='0'+string(fnum,'(f3.1)')

  ;	find out if there are any data files at all
  fils=dbasedir+'/'+zz+'_'+['wvl','tem','ion','lvl','ecn','src',fsfx]
  nfil=n_elements(fils)
  if not keyword_set(allah) then fils=findfile(dbasedir+'/'+zz+'_*',count=nfil)

  if nfil ge 7 then begin				;{yep, files exist

    kw=0 & kt=1 & ki=2 & kl=3 & ke=4 & ks=5 & kj=-1 & kf=[7] & nkf=1
    if not keyword_set(allah) then begin		;{search in DBASEDIR
      for i=0,nfil-1 do begin	;{
        if strpos(fils[i],'_wvl') gt 0 then kw=i else begin		;{
         if strpos(fils[i],'_tem') gt 0 then kt=i else begin		;(
          if strpos(fils[i],'_ion') gt 0 then ki=i else begin		;{
           if strpos(fils[i],'_jon') gt 0 then kj=i else begin		;{
            if strpos(fils[i],'_lvl') gt 0 then kl=i else begin		;(
	     if strpos(fils[i],'_ecn') gt 0 then ke=i else begin	;{
	      if strpos(fils[i],'_src') gt 0 then ks=i else kf=[kf,i]
	     endelse							;}
	    endelse							;)
	   endelse							;}
	  endelse							;}
	 endelse							;)
        endelse								;}
      endfor			;}
      kf=kf[1:*] & nkf=n_elements(kf)
    endif		;allahu akabar}
    ;	and the input files are...
    wfl=fils[kw] & tfl=fils[kt] & ifl=fils[ki]
    lfl=fils[kl] & efl=fils[ke] & sfl=fils[ks] & ffl=fils[kf]
    if kj ge 0 then jfl=fils[kj] else jfl=ifl	;by default use ATOM_ion

    ;	which dataset is nearest to supplied value of pressure/density?
    pp=fltarr(nkf)
    tmp=str_sep(ffl[0],'_') & npp=n_elements(tmp)
    for ip=0,nkf-1 do pp[ip]=float((str_sep(ffl[ip],'_'))[npp-1L])
    ;	need to look only at the last bit, because what if there
    ;	are other underscores in the pathname?
    ip=(where(abs(pp-fnum) eq min(abs(pp-fnum))))[0] & ffl=ffl[ip]

    ;	chatter
    if v ge 5 then begin
      message,'Reading from files:',/info
      print,wfl,tfl,ifl,jfl,sfl,ffl
      if rdlvl eq 1 then print,lfl
      if rdecn eq 1 then print,efl
    endif

    ;	ok, read 'em (and weep?)
    kilroy; was here
    ;	this digression needed to figure out whether to swap byte order
    ;	or not
    if n_elements(set_byte_swap) eq 0 then begin
      openr,uw,wfl,/get_lun
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
    if set_byte_swap gt 0 then begin
      openr,uw,wfl,/get_lun,/swap_endian
      openr,ut,tfl,/get_lun,/swap_endian
      openr,ui,ifl,/get_lun,/swap_endian
      openr,us,sfl,/get_lun,/swap_endian
      openr,uf,ffl,/get_lun,/swap_endian
      if kj ge 0 then openr,uj,jfl,/get_lun,/swap_endian
    endif else begin
      openr,uw,wfl,/get_lun
      openr,ut,tfl,/get_lun
      openr,ui,ifl,/get_lun
      openr,us,sfl,/get_lun
      openr,uf,ffl,/get_lun
      if kj ge 0 then openr,uj,jfl,/get_lun
    endelse
    nw=0L & readu,uw,nw
    if set_byte_swap lt 0 then nw=swap_endian(nw)
    ang=fltarr(nw) & readu,uw,ang		;wvl
    if set_byte_swap lt 0 then ang=swap_endian(ang)
    nt=0L & readu,ut,nt
    if set_byte_swap lt 0 then nt=swap_endian(nt)
    temp=fltarr(nt) & readu,ut,temp		;tem
    if set_byte_swap lt 0 then temp=swap_endian(temp)
    fx=dblarr(nt,nw) & readu,uf,nt,nw,fx			;##.#
    if set_byte_swap lt 0 then nw=swap_endian(nw)
    if set_byte_swap lt 0 then nt=swap_endian(nt)
    if set_byte_swap lt 0 then fx=swap_endian(fx)
    ;	HACK TO AVOID NaNs in CHIANTI3
    ofin=where(finite(fx) eq 0,mofin) & if mofin gt 0 then fx(ofin)=0.
    ikey=lonarr(nw) & readu,ui,nw,ikey				;ion
    if set_byte_swap lt 0 then nw=swap_endian(nw)
    if set_byte_swap lt 0 then ikey=swap_endian(ikey)
    jkey=ikey & if kj ge 0 then readu,uj,nw,jkey		;jon
    if set_byte_swap lt 0 then nw=swap_endian(nw)
    if set_byte_swap lt 0 then jkey=swap_endian(jkey)
    srs=lonarr(nw) & readu,us,nw,srs				;src
    if set_byte_swap lt 0 then nw=swap_endian(nw)
    if set_byte_swap lt 0 then srs=swap_endian(srs)
    if rdlvl eq 1 then restore,file=lfl	;tmp_trans		;lvl
    if rdecn eq 1 then restore,file=efl	;tmp_econf		;ecn
    kilroy; was here
    close,uw & free_lun,uw
    close,ut & free_lun,ut
    close,ui & free_lun,ui
    close,us & free_lun,us
    close,uf & free_lun,uf
    if kj ge 0 then close,uj & if kj ge 0 then free_lun,uj

    ;	select in wavelength range
    oo=where(abs(ang) ge wrange[0] and abs(ang) le wrange[1])
    if oo[0] ne -1 then begin
      ang=ang[oo] & fx=fx[*,oo]
      ikey=ikey[oo] & jkey=jkey[oo] & srs=srs[oo]
      kilroy; was here
      if rdlvl eq 1 then tran1=tmp_trans[*,oo] else tran1=strarr(2,1)
      if rdecn eq 1 then tran2=tmp_econf[*,oo] else tran2=strarr(2,1)
      kilroy; was here
    endif else begin
      ang=[0.] & fx=dblarr(nt,1) & ikey=[-1L] & jkey=[-1L] & srs=[0L]
      tran1=strarr(2,1) & tran2=tran1
      print,'' & message,'	No lines in this wavelength range',/info
    endelse

    ;	select according to ION
    if kion gt 0 then begin
      oo=where(ikey eq kion)
      if oo[0] ne -1 then begin
        ang=ang[oo] & fx=fx[*,oo] & ikey=ikey[oo] & jkey=jkey[oo] & srs=srs[oo]
        if rdlvl eq 1 then tran1=tran1[*,oo] else tran1=strarr(2,1)
        if rdecn eq 1 then tran2=tran2[*,oo] else tran2=strarr(2,1)
      endif else begin
	ang=[0.] & fx=dblarr(nt,1) & ikey=[-1L] & jkey=[-1L] & srs=[0L]
        tran1=strarr(2,1) & tran2=tran1
	print,'' & message,'	No lines for this ionic state',/info
      endelse
    endif

   ;	add to output
   if ikey[0] gt 0 then begin			;any entries left?
     nw=n_elements(ang)
     wvl=[wvl,ang]
     Z=[Z,intarr(nw)+j] & ion=[ion,ikey] & jon=[jon,jkey] & src=[src,srs]
     if not is_keyword_set(f) then f=reform(fx[*]) else f=[f,fx[*]]
     kilroy; was here
     if rdlvl eq 1 then desig=[desig,tran1[*]]
     if rdecn eq 1 then econf=[econf,tran2[*]]
   endif
   kilroy; was here.

  endif else begin		;if (nwf+ntf+nif+nlf)=4 && (nff>0)}
    message,'	No line information available for '+zz,/info
  endelse

endfor	;iz=0,nz-1}

;	and reconfigure the arrays
if n_elements(wvl) gt 1 then begin
  print,''
  if v gt 0 then kilroy; was here.
  wvl=wvl[1:*] & Z=Z[1:*] & ion=ion[1:*] & jon=jon[1:*] & src=src[1:*]
  nw=n_elements(wvl) & logT=temp
  f=reform(f,nt,nw)
  if v gt 0 then kilroy; was here
  if rdlvl eq 1 then desig=reform(desig[1:*],2,nw)
  if rdecn eq 1 then econf=reform(econf[1:*],2,nw)
  if v gt 0 then kilroy; was here.
endif else f=fltarr(1,1)-1.

if arg_present(fstr) then begin
  if v gt 0 then kilroy; was here.
  fstr=create_struct('LINE_INT',f,'LOGT',logt,'WVL',wvl,'Z',z,'ION',ion,$
	'DESIG',desig,'CONFIG',econf,'SRC',src,'JON',jon)
  if v gt 0 then kilroy; was here.
endif

print,''

return,f
end
