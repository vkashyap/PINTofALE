function getpoadef,varname,force=force,poarc=poarc,verbose=verbose,_extra=e
;+
;function	getpoadef
;	returns the default value for a common PINTofALE variable
;	usually for filenames, which are often dependent on where
;	the files are installed.
;
;	The way this works is, if a system variable is already set,
;	the default is taken from that.  If not, after figuring out
;	what the default is (usually by checking where this program
;	is located), the appropriate system variable is set to it.
;
;syntax
;	default=getpoadef(varname,/force,poarc=poarc,verbose=verbose)
;
;parameters
;	var	[INPUT] name of variable as a character string
;		* may be an array of strings
;		* special variable names are the only ones recognized
;		* if 'help' or '?', lists all the variables for which
;		  defaults can be set
;		* if set to 'all' or 'ALL', calls itself and sets all
;		* in some cases default might be a number, but if the
;		  input is an array, then the output will always be a
;		  string array
;	poarc	[INPUT] the name of the file containing default system
;		variable definitions, e.g., !ARDB/initale.par
;		* first looks for filename as given (so if not full path,
;		  then in $cwd), then looks for the basename in $HOME,
;		  and finally in !ARDB/initale.par
;		* if not defined, does not do any of the above
;		* if given, then the corresponding variable is set if
;		  it is defined in the input definitions file, unless
;		  FORCE is set, in which case it resets to hardcoded
;		  default
;		* uses spawn and grep liberally, so only works with
;		  !version.os_family='unix'
;
;keywords
;	force	[INPUT] if set, ignores whether any default (except TOPDIR)
;		has already been set and resets it anyway to the closest
;		thing this routine has to a factory default
;		* to reset TOPDIR, just do it on the command line, as
;		  !TOPDIR='/path/to/PINTofALE'
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2015aug)
;	now if input is scalar, returns scalar (VK; 2015sep)
;	added keyword POARC, allowed option 'ALL' to VARNAME, made VERBOSE
;	  a bit more useful, now updates !ABUND if !ABREF has been changed
;	  (VK; 2016nov)
;	fixed creash when using FORCE and !VERBOSE was already defined (VK; 2017nov)
;	makes it easier to pick up !TOPDIR in GDL (VK; 2018jan)
;-

;	usage
ok='ok' & np=n_params() & defsysv,'!TOPDIR',exists=iPoA & nv=n_elements(varname) & szv=size(varname,/type)
if np ne 0 then begin
  if nv eq 0 then ok='VARNAME is not defined' else $
  if szv ne 7 then ok='VARNAME must be char string scalar or vector'
endif
if ok ne 'ok' then begin
  print,'Usage: default=getpoadef(varname,/force,poarc=poarc,verbose=verbose)'
  print,'  returns default values for common PINTofALE variables'
  print,getpoadef('help')
  if nv ne 0 then message,ok,/informational
  return,-1L
endif
if np eq 0 and iPoA ne 0 then return,!TOPDIR

;	keywords
defsysv,'!VERBOSE',exists=ivar
if ivar eq 1 then jnk=execute('vv=!VERBOSE') else vv=0L
if keyword_set(verbose) then vv=long(verbose[0])>1L
if not keyword_set(force) then begin
  if ivar eq 0 then defsysv,'!VERBOSE',vv
endif else begin
  if ivar eq 1 then jnk=execute('!VERBOSE=10') else defsysv,'!VERBOSE',10
endelse

;	is a POARC given?
if keyword_set(poarc) then begin
  tmp=file_search(strtrim(poarc[0],2),count=nfil)
  if nfil gt 0 then initfil=tmp[0] else begin
    if !version.os_family eq 'unix' then begin
      spawn,'basename '+strtrim(poarc[0],2),tmp2
      tmp=file_search(tmp2[0],count=nfil2)
      if nfil2 gt 0 then initfil=tmp[0] else initfil=filepath('initale.par',root_dir=getpoadef('ARDB'))
    endif else initfil=0
  endelse
endif else initfil=0
if keyword_set(initfil) and vv gt 1 then begin
  if keyword_set(force) then message,'ignoring POARC='+initfil,/informational else $
    message,'picking up defaults from '+initfil,/informational
endif

;	define output
def=strarr(nv>1)

;	basic default
if iPoA eq 1 then topdir=!TOPDIR else begin	;(iPoA=0

  topdir=filepath('PINTofALE',root=getenv('HOME'),subdir='MEGA')	;hardcoded default
  defsysv,'!TOPDIR',exists=ivar
  if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
    spawn,'grep \!TOPDIR '+initfil,tmp
    if tmp[0] ne '' then begin & defsysv,'!TOPDIR','' & jnk=execute(tmp[0]) & endif
    jnk=execute("topdir=!TOPDIR")
  endif else begin
    ;figure out where this program is and thence where TOPDIR should be
    ;	NOTE: this does not work in GDL, so that part is hardcoded, sorry
    defsysv,'!GDL',exists=ivar
    if ivar eq 0 then begin
      help,/source,output=srcpaths
      ok=where(strpos(srcpaths,'getpoadef.pro') gt 0,mok)
      if mok gt 0 then begin
        cc=strsplit(srcpaths[ok[0]],' ',/extract)
        ccc=cc[1] & itree=strpos(ccc,'/pro/getpoadef.pro')
        topdir=filepath('',root_dir=strmid(ccc,0,itree))
      endif else begin
        message,'something went wrong'
      endelse
    endif else begin
      print,'GDL does not support the OUTPUT keyword to HELP.'
      print,'We will define the !TOPDIR variable here, just'
      print,'reset it to whatever you want and rerun this function.
      print,'Or set it in your GDL_STARTUP.'
    endelse
    defsysv,'!TOPDIR',topdir
    jnk=execute("topdir=!TOPDIR")
  endelse
  if vv gt 0 then message,'!TOPDIR='+!TOPDIR,/informational
  if nv eq 0 then def=topdir

endelse						;iPoA=0)

;	step through each VARNAME and figure out the default
for i=0L,nv-1L do begin		;{I=0,NV-1

  cc=strtrim(strupcase(varname[i]),2)
  case cc of		;{variable name

    'HELP': begin
      print,'variable names can be:'
      print,'VERSION TOPDIR LDBDIR CDBDIR ARDB CEROOT CHIDIR IONQEF ATOMDB ABREF CALDB'
      print,'or ALL'
    end

    'ALL': begin
      jnk=getpoadef('TOPDIR',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('LDBDIR',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('CDBDIR',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('ARDB',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('CEROOT',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('CHIDIR',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('IONEQF',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('ATOMDB',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('ABREF',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('CALDB',force=force,poarc=poarc,verbose=vv)
      jnk=getpoadef('VERSION',force=force,poarc=poarc,verbose=vv)
    end

    'TOPDIR': def[i]=getpoadef() ;no need to do anything, this is basic default, see above

    'VERSION': begin	;(version
      defsysv,'!PoA',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!PoA else begin
        def[i]='2.97 (2016nov24)'
	defsysv,'!PoA',def[i]
	if vv gt 0 then message,'PINTofALE version '+def[i],/informational
      endelse		;VERSION)
    end

    'LDBDIR': begin	;(!TOPDIR/emissivity/chianti
      defsysv,'!LDBDIR',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!LDBDIR '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!LDBDIR','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!LDBDIR',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!LDBDIR else begin
	def[i]=filepath('chianti',root_dir=!TOPDIR,subdir='emissivity')
        defsysv,'!LDBDIR',def[i]
      endelse
      if vv gt 0 then message,'!LDBDIR='+!LDBDIR,/informational
    end			;LDBDIR)

    'CDBDIR': begin	;(!TOPDIR/emissivity/cont
      defsysv,'!CDBDIR',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!CDBDIR '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!CDBDIR','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!CDBDIR',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!CDBDIR else begin
	def[i]=filepath('cont',root_dir=!TOPDIR,subdir='emissivity')
        defsysv,'!CDBDIR',def[i]
      endelse
      if vv gt 0 then message,'!CDBDIR='+!CDBDIR,/informational
    end			;CDBDIR)

    'ARDB': begin	;(!TOPDIR/ardb
      defsysv,'!ARDB',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!ARDB '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!ARDB','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!ARDB',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!ARDB else begin
	def[i]=filepath('ardb',root_dir=!TOPDIR)
        defsysv,'!ARDB',def[i]
      endelse
      if vv gt 0 then message,'!ARDB='+!ARDB,/informational
    end			;ARDB)

    'CEROOT': begin	;(cie
      defsysv,'!CEROOT',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!CEROOT '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!CEROOT','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!CEROOT',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!CEROOT else begin
	def[i]='cie'
        defsysv,'!CEROOT',def[i]
      endelse
      if vv gt 0 then message,'!CEROOT='+!CEROOT,/informational
    end			;CEROOT)

    'CHIDIR': begin	;(!TOPDIR/../CHIANTI/dbase
      defsysv,'!CHIDIR',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!CHIDIR '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!CHIDIR','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!CHIDIR',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!CHIDIR else begin
	ccc=!TOPDIR & cccl=strmid(ccc,strlen(ccc)-1,1)
	while cccl eq '/' do begin & ccc=strmid(ccc,0,strlen(ccc)-1) & cccl=strmid(ccc,strlen(ccc)-1,1) & endwhile	;strip trailing "/"s
        itree=strpos(ccc,'/',/reverse_search)
	def[i]=filepath('dbase',root_dir=strmid(ccc,0,itree),subdir='CHIANTI')
        defsysv,'!CHIDIR',def[i]
      endelse
      if vv gt 0 then message,'!CHIDIR='+!CHIDIR,/informational
    end			;CHIDIR)

    'IONEQF': begin	;(ioneq/chianti.ioneq
      defsysv,'!IONEQF',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!IONEQF '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!IONEQF','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!IONEQF',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!IONEQF else begin
	def[i]=filepath('chianti.ioneq',root_dir='ioneq')
        defsysv,'!IONEQF',def[i]
      endelse
      if vv gt 0 then message,'!IONEQF='+!IONEQF,/informational
    end			;IONEQF)

    'ATOMDB': begin	;(!TOPDIR/../atomdb
      defsysv,'!ATOMDB',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!ATOMDB '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!ATOMDB','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!ATOMDB',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!ATOMDB else begin
	ccc=!TOPDIR & cccl=strmid(ccc,strlen(ccc)-1,1)
	while cccl eq '/' do begin & ccc=strmid(ccc,0,strlen(ccc)-1) & cccl=strmid(ccc,strlen(ccc)-1,1) & endwhile	;strip trailing "/"s
        itree=strpos(ccc,'/',/reverse_search)
	def[i]=filepath('',root_dir=strmid(ccc,0,itree),subdir='atomdb')
        defsysv,'!ATOMDB',def[i]
      endelse
      if vv gt 0 then message,'!ATOMDB='+!ATOMDB,/informational
    end			;ATOMDB)

    'ABREF': begin	;(grevesse et al.
      defsysv,'!ABREF',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!ABREF '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!ABREF','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!ABREF',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!ABREF else begin
	def[i]='Grevesse et al.'
        defsysv,'!ABREF',def[i]
      endelse
      ;
      ;overwrite !ABUND for this !ABREF
      jnk=execute("tmp=!ABREF")
      defsysv,'!ABUND',getabund(tmp)
      if vv gt 0 then message,'!ABREF='+!ABREF,/informational
      if vv gt 0 then message,'!ABUND has been updated',/informational
    end			;ABREF)

    'CALDB': begin	;(!TOPDIR/caldb
      defsysv,'!CALDB',exists=ivar
      if ivar eq 0 and !version.os_family eq 'unix' and keyword_set(initfil) then begin
        spawn,'grep \!CALDB '+initfil,tmp
	if tmp[0] ne '' then begin & defsysv,'!CALDB','' & jnk=execute(tmp[0]) & endif
      endif
      defsysv,'!CALDB',exists=ivar
      if keyword_set(force) then ivar=0	;force override
      if ivar eq 1 then def[i]=!CALDB else begin
	def[i]=filepath('',root_dir=!TOPDIR,subdir='caldb')
        defsysv,'!CALDB',def[i]
      endelse
      if vv gt 0 then message,'!CALDB='+!CALDB,/informational
    end			;CALDB)

    else: begin	;(
      if vv gt 0 then message,cc+': variable name not understood',/informational
    end		;else)

  endcase		;case CC}

endfor				;I=0,NV-1}

if nv eq 1 then def=def[0]

return,def
end
