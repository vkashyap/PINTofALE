pro wrt_ln_generic,linstr,outdir,tag,abund=abund,excieq=excieq,eps=eps,$
	verbose=verbose, _extra=e
;+
;procedure	wrt_ln_generic
;	writes out line intensities to a PoA style database from
;	a PoA style line emissivities structure.
;
;usage
;	wrt_ln_generic,linstr,outdir,tag,abund=abund,excieq=excieq,$
;	eps=eps,chidir=chidir,eqfile=eqfile
;
;parameters
;	linstr	[INPUT; required] line emissivities structure, of
;		the same form output by RD_LINE()
;		* NOTE: the emissivities must be in units of
;		  [1e-23 ergs cm^3/s]
;	outdir	[INPUT; required] directory to deposit the database files
;		* WARNING: the PINTofALE convention is to automatically
;		  append a "D" while reading from emissivities derived
;		  at constant density.  so if the emissivities being
;		  written out were derived at constant density, then
;		  please makes sure that OUTDIR has an extra "D" at the
;		  end.  this is not checked for, because a priori this
;		  program can have no idea of the nature of the emissivities
;		  being input.
;	tag	[INPUT; required] a value (such as pressure, or density)
;		defining the emissivities.  see description of output
;		files below for where it is used.
;		* MUST BE writeable in the format (f4.1)
;		* if > 99, the alog10 value is used
;		* if < 0, forced to be 00.0
;
;keywords
;	abund	[INPUT] if given, assumes that the input emissivities
;		have been multiplied by these abundance values, and
;		that you don't want this factor stored in the database
;		files, and divides them out.
;	excieq	[INPUT] if set, assumes that the input emissivities have
;		been multiplied by the ion balance defined by the
;		keywords CHIDIR and EQFILE, and that you don't want them
;		to be included in the stored database, and divides them
;		out over the range where the ion balance is > EPS
;	eps	[INPUT] a small number
;		* default is 1e-6
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords CHIDIR and EQFILE
;		to RD_IONEQ()
;
;side effects
;	deposits (possibly large) files in OUTDIR
;	* ATOM_wvl	wavelengths of all the transitions
;	* ATOM_tem	temperatures
;	* ATOM_##.#	line emissivities [tem,wvl]_TAG
;	* ATOM_ion	ionic states corresponding to wvl
;	* ATOM_jon	ionic states that matter to the ion balance
;	* ATOM_src	source of line information
;	* ATOM_lvl	level designations of each transition
;	* ATOM_ecn	e-configurations of each level
;
;subroutines:
;	ARRAYEQ
;	INICON
;	SETSYSVAL, LEGALVAR()
;	RD_IONEQ(), READ_IONEQ
;
;history
;	vinay kashyap (Jun02)
;	Antonio Maggio discovered that ion, jon, and src files were
;	  writing out float rather than long (VK; Aug02)
;-

;	usage
ok='ok' & np=n_params() & nl=n_elements(linstr) & ml=n_tags(linstr)
nd=n_elements(outdir) & nt=n_elements(tag)
if np lt 3 then ok='Insufficient parameters' else $
 if nl eq 0 then ok='LINSTR is undefined' else $
  if ml eq 0 then ok='LINSTR must be a structure -- see RD_LINE()' else $
   if nd eq 0 then ok='OUTDIR is undefined' else $
    if nd gt 1 then ok='OUTDIR must be a scalar' else $
     if nt eq 0 then ok='TAG is undefined' else $
      if nt gt 1 then ok='TAG must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: wrt_ln_generic,linstr,outdir,tag,abund=abund,excieq=excieq,$'
  print,'       eps=eps,verbose=verbose,chidir=chidir,eqfile=eqfile'
  print,'  writes line emissivities into PoA-type database [see RD_LINE()]'
  if np ne 0 then message,ok,/info
  return
endif

;	verify inputs
lnam=tag_names(linstr) & k=0
kli=0 & klt=0 & kw=0 & kz=0 & ki=0 & kj=0 & kd=0 & ke=0 & ks=0
for i=0,ml-1 do begin
  if lnam[i] eq 'LINE_INT' then kli=i+1
  if lnam[i] eq 'LOGT' then klt=i+1
  if lnam[i] eq 'WVL' then kw=i+1
  if lnam[i] eq 'Z' then kz=i+1
  if lnam[i] eq 'ION' then ki=i+1
  if lnam[i] eq 'JON' then kj=i+1
  if lnam[i] eq 'DESIG' then kd=i+1
  if lnam[i] eq 'CONFIG' then ke=i+1
  if lnam[i] eq 'SRC' then ks=i+1
endfor
ok='ok'
if not keyword_set(kli) then ok='linstr.LINE_INT is missing' else $
 if not keyword_set(klt) then ok='linstr.LOGT is missing' else $
  if not keyword_set(kw) then ok='linstr.WVL is missing' else $
   if not keyword_set(kz) then ok='linstr.Z is missing' else $
    if not keyword_set(ki) then ok='linstr.ION is missing' else $
     if not keyword_set(kj) then ok='linstr.JON is missing' else $
      if not keyword_set(kd) then ok='linstr.DESIG is missing' else $
       if not keyword_set(ke) then ok='linstr.CONFIG is missing' else $
	if not keyword_set(ks) then ok='linstr.SRC is missing'
if ok ne 'ok' then begin
  message,ok,/info
  message,'LINSTR is not in standard format.  returning.',/info
  return
endif
;
szd=size(outdir) & nszd=n_elements(szd)
if szd[nszd-2] ne 7 then begin
  message,'OUTDIR must be a string. cannot write to this. returning',/info
  return
endif else odir=outdir[0]
;
fcode='00.0'
szt=size(tag) & nszt=n_elements(szt)
if szt[nszt-2] eq 7 then fcode=strmid(strtrim(tag[0],2),0,4) else begin
  if tag[0] lt 0 then fcode='00.0' else $
   if tag[0] gt 99 then fcode=string(strtrim(alog10(tag[0]),2),'(f4.1)') else $
    fcode=string(strtrim(tag[0],2),'(f4.1)')
endelse
fcode=strtrim(fcode,2) & for i=1,4-strlen(fcode) do fcode='0'+fcode
;
vv=0 & ivar=0
defsysv,'!VERBOSE',exists=ivar	;if !VERBOSE exists
if ivar ne 0 then setsysval,'VERBOSE',vv,/getval
if n_elements(verbose) gt 0 then vv=long(verbose[0])>0	;override !VERBOSE
;
epsilon=1e-6 & if keyword_set(eps) then epsilon=eps[0]

;	initialize
slash='/'
case !version.OS_FAMILY of
  'unix': slash='/'
  'windows': slash='\'
  'macos': slash=':'
  'vms': slash='\'
  else: slash='/'	;unknown OS, assume UNIX
endcase
set_byte_swap=0

;	if OUTDIR doesn't exist, create it
;	clearly, this bit works only on UNIX
fil=findfile(odir,count=nfil)
if !version.OS_FAMILY eq 'unix' then begin
  if nfil eq 0 then spawn,'mkdir '+odir
endif else begin
  message,odir+' : does not exist.  please create it first.',/info
  return
endelse

;	initialize output file extensions
inicon,atom=atom
wfil='wvl' & tfil='tem' & ifil='ion' & jfil='jon'
ffil=fcode & sfil='src' & lfil='lvl' & efil='ecn'

;	extract the variables from LINSTR
emis=linstr.(kli-1) & logT=linstr.(klt-1) & wvl=linstr.(kw-1)
Z=linstr.(kz-1) & ion=linstr.(ki-1) & jon=linstr.(kj-1)
desig=linstr.(kd-1) & econf=linstr.(ke-1) & src=linstr.(ks-1)
sze=size(emis) & netem=sze[0] & newvl=sze[1]
ntem=n_elements(logT) & nwvl=n_elements(wvl)
szdesig=size(desig) & szeconf=size(econf)
nszdesig=n_elements(szdesig) & nszeconf=n_elements(szeconf)
idesig=0 & if szdesig[nszdesig-1] eq 2*nwvl then idesig=1
ieconf=0 & if szeconf[nszeconf-1] eq 2*nwvl then ieconf=1

;	write out the files for each atom
uZ=Z[uniq(Z,sort(Z))] & nuZ=n_elements(uZ)
for iz=0,nuZ-1 do begin			;{for each atom

  elm=atom[uZ[iz]-1]
  wfile=odir+slash+elm+'_'+wfil
  tfile=odir+slash+elm+'_'+tfil
  ifile=odir+slash+elm+'_'+ifil
  jfile=odir+slash+elm+'_'+jfil
  ffile=odir+slash+elm+'_'+ffil
  sfile=odir+slash+elm+'_'+sfil
  lfile=odir+slash+elm+'_'+lfil
  efile=odir+slash+elm+'_'+efil

  ;	filter
  oZ=where(Z eq uZ[iZ],moZ)
  if moZ eq 0 then message,'BUG!'
  ww=wvl[oz] & ii=ion[oz] & jj=jon[oz] & ff=emis[*,oz] & ss=src[oz]
  if idesig eq 1 then dd=desig[*,oz] else dd=strarr(2,moz)
  if ieconf eq 1 then ee=econf[*,oz] else ee=strarr(2,moz)

  ;	check if these files already exist
  ;	(only need to check WFILE and FFILE)
  ;btw, all this jumping through hoops with set_byte_swap only needed for
  ;early versions of OPENR which don't accept the /swap_endian keyword.
  wfl=findfile(wfile,count=iwfl) & ffl=findfile(ffile,count=iffl)

  if iwfl ne 0 then begin
    if vv gt 0 then message,wfile+' : already exists',/info
    openr,uw,wfile,/get_lun
      nw=0L & readu,uw,nw
      if keyword_set(set_byte_swap) then nw=swap_endian(nw) else begin
        if nw ge 16777216 or nw lt 0 then begin
	  set_byte_swap=1 & nw=swap_endian(nw)
	endif
      endelse
      ang=fltarr(nw) & readu,uw,ang
      if keyword_set(set_byte_swap) then ang=swap_endian(ang)
    close,uw & free_lun,uw
    wvleqang=arrayeq(wvl,ang)
    if wvleqang ne 1 then begin	;(the two arrays differ
      if vv gt 0 then message,$
	'input wavelengths differ from existing file values',/info
      c='Y'
      if vv ge 5 then begin
	c='overwrite? [Y/n] ' & print,c,form='($,a)' & c=get_kbrd(1)
	if strlowcase(c) ne 'n' then c='Y'
      endif
      if c ne 'Y' then begin
	message,'Exiting without overwriting '+wfile,/info
	return
      endif else begin
	if vv gt 0 then begin
	  message,'Will overwrite the following files:',/info
	  print,wfile & print,tfile & print,ifile & print,jfile
	  print,sfile & print,lfile & print,efile
	endif
      endelse
    endif			;WVLEQANG.NE.1)
  endif

  if iffl ne 0 then begin
    if vv gt 0 then message,ffile+' : already exists',/info
    openr,uf,ffile,/get_lun
      nw=0L & nt=0L & readu,uf,nt,nw
      if keyword_set(set_byte_swap) then begin
        nw=swap_endian(nw) & nt=swap_endian(nt)
      endif else begin
        if nw ge 16777216 or nw lt 0 then begin
	  set_byte_swap=1 & nw=swap_endian(nw) & nt=swap_endian(nt)
	endif
      endelse
      ;em=fltarr(nt,nw) & readu,uf,em
      ;if keyword_set(set_byte_swap) then em=swap_endian(em)
    close,uf & free_lun,uf
    if nt ne ntem then begin	;(the temperature grid differs
      if vv gt 0 then message,$
	'temperature grid differs for input emissivities',/info
      c='Y'
      if vv ge 5 then begin
	c='overwrite? [Y/n] ' & print,c,form='($,a)' & c=get_kbrd(1)
	if strlowcase(c) ne 'n' then c='Y'
      endif
      if c ne 'Y' then begin
	message,'Exiting without overwriting '+ffile,/info
	return
      endif else begin
	if vv gt 0 then message,'Will overwrite '+ffile,/info
      endelse
    endif			;NT.NE.NTEM)
  endif

    ;	correct for abundances if required
    if n_elements(abund) ge uZ[iz] then emis=emis/abund[uZ[iz]-1]

    ;	correct for ion balance if required
    if keyword_set(excieq) then begin
      ieq=rd_ioneq(uZ[iz],logT,verbose=vv, _extra=e)
      for j=0,uZ[iZ] do begin	;{for each ionic state
	oi=where(jj eq j,moi)
	if moi gt 0 then begin
	  ieqj=reform(ieq[*,j])
	  ok1=where(ieqj gt eps,mok1)
	  ok2=where(ieqj le eps,mok2)
	  if mok1 gt 0 then ff[ok1,*]=ff[ok1,*]/ieqj[ok1]
	  if mok2 gt 0 then ff[ok2,*]=0.
	endif
      endfor			;J=0,UZ[IZ]}
    endif

    ;	OK, now write them out
    openw,uw,wfile,/get_lun,/swap_if_little_endian
    openw,ut,tfile,/get_lun,/swap_if_little_endian
    openw,ui,ifile,/get_lun,/swap_if_little_endian
    openw,uj,jfile,/get_lun,/swap_if_little_endian
    openw,us,sfile,/get_lun,/swap_if_little_endian
    ;	openw,ul,lfile,/get_lun,/swap_if_little_endian
    ;	openw,ue,efile,/get_lun,/swap_if_little_endian
    openw,uf,ffile,/get_lun,/swap_if_little_endian
      writeu,uw,moz,float(ww)
      writeu,ut,ntem,float(logT)
      writeu,ui,moz,long(ii)
      writeu,uj,moz,long(jj)
      writeu,us,moz,long(ss)
      tmp_trans=dd & save,file=lfile,tmp_trans
      tmp_econf=ee & save,file=efile,tmp_econf
      writeu,uf,ntem,moz,double(ff)
    close,uw & free_lun,uw
    close,ut & free_lun,ut
    close,ui & free_lun,ui
    close,uj & free_lun,uj
    close,us & free_lun,us
    ;	close,ul & free_lun,ul
    ;	close,ue & free_lun,ue
    close,uf & free_lun,uf

    if vv gt 0 then message,'writing to '+odir,/info
    if vv gt 5 then begin
      print,wfile & print,tfile & print,ifile & print,jfile
      print,ffile & print,sfile & print,lfile & print,efile
    endif

endfor					;IZ=0,NUZ-1}

return
end
