;+
;script	initale
;	extremely flexible initialization script for the Package for
;	Interactive Analysis of Line Emission (PINTofALE); defines some
;	useful variables as system variables.  Examples and scripts in
;	the standard distribution use these variables.
;
;	NOTE: Strictly speaking, none of these system variables are
;	_necessary_ in order to use basic PINTofALE functionality.
;	However, they set up the default values for some commonly
;	used keywords and usually must be explicitly included in
;	calls to PINTofALE routines.  Most scripts and documentation
;	examples do use these directly, and some routines set their
;	default values at run-time if the system variables are defined.
;	An additional use of this script is that it sets up the !PATH
;	to include PINTofALE, CHIANTI, and IDL-Astro if they are not
;	already included.
;
;syntax
;	.run initale
;	(or)
;	.run /full/path/initale
;
;system variables and hardcoded default values
;	!PoA	 current version number
;	!TOPDIR	 the top-level PINTofALE directory
;		 * first checks for environment variables
;			PINTofALE, PoA, SCARDIR, and SCAR
;		   in that order.  if none of them exist, then uses
;		   value deduced from the location of _this_ script.
;	!LDBDIR	 [!TOPDIR+'emissivity/chianti'] directory of line emissivities
;	!CDBDIR	 [!TOPDIR+'emissivity/cont'] directory of continuum emissivities
;	!CHIDIR	 [!TOPDIR+'CHIANTI'] path to CHIANTI installation
;		 NOTE: !TOPDIR is added to initial values of above 3 ONLY IF:
;		 -- !TOPDIR is not already contained in them
;		 -- the first character is not a "$"
;		 -- the first character is not a "/"
;		 -- for CHIANTIv4+, path name must include the "dbase" part
;	!ATOMDB  [!TOPDIR+'atomdb/'] directory of local installation of ATOMDB
;	!APECDIR where the IDL tools of APED are
;		 * first checks for environment variable APEC_DIR, and if
;		   that is missing, looks in succession in
;		   -- !ATOMDB+'apec_v11_idl/'
;		   -- !TOPDIR+'apec_v11_idl/'
;		   -- !TOPDIR+'atomdb/apec_v11_idl/'
;	!ARDB	 [!TOPDIR+'ardb/'] the Analysis Reference DataBase directory,
;		 contains sundry useful documents
;	!CEROOT	 ['cie'] root prefix for continuum emissivity files
;	!IONEQF	 ['ioneq/bryans_etal_09.ioneq'] path name relative to
;		   !CHIDIR of the ion-balance file to be used
;	!ABREF	 ['grevesse et al.'] abundance reference
;	!CALDB	 ['/data/caldb/'] directory containing instrument
;		    calibration data products
;	!METALS	 [0.0] [Fe/H] relative to !ABREF
;	!ABUND	 abundances of elements H..Zn corresponding to !ABREF and
;		 !METALS
;	!GASPR	 [1e15 cm^-3 K] gas pressure
;	!LOGPR	 [15] log10(gas pressure) -- if set, overrides !GASPR
;	!EDENS	 [1e10 cm^-3] electron density -- if set, overrides use of
;		   !GASPR and !LOGPR
;	!LOGT	 [findgen(81)*0.05+4] the default temperature grid
;	!DEM	 [dblarr(81)+1d12] a sample DEM(!LOGT)
;	!NH	 [1e18 cm^-2] H column density
;	!FH2	 [0.26] fraction of molecular H2 relative to HI
;	!HE1	 [1e17 cm^-2] He I column
;	!HEII	 [1e16 cm^-2] He II column
;	!WMIN	 [1.239854 Ang] minimum wavelength
;	!WMAX	 [3000.0 Ang] maximum wavelength
;	!VERBOSE [5] verbosity -- controls chatter
;	!ATOM	 atomic symbols
;	!AMU	 atomic weights
;	!FUNDAE	 a whole bunch of fundamental constants
;	!FIP	 firt ionization potentials [eV]
;	!ROMAN	 roman numerals, upto max(!ATOM)+1
;	!AA	 ['!3'+string(byte(197))+'!X'] the symbol for Angstrom
;	!ANGSTROM	An array of Angstrom symbols, should work for
;			every conceivable font choice.
;
;control variables
;	factory	 determines whether to reset all the system variables to
;		 the "factory default" or not.  set this variable to 0 or 1 
;		 prior to running the script.
;		 * if explicitly set to 0, then only initializes those
;		   variables which have not been defined yet
;		 * if PARFIL is defined, automatically gets unset
;		 * if set, or is not defined, then overwrites all variables
;		   with hardcoded defaults
;		 -- one may want to use FACTORY=0 when some of the variables
;		    have been set (e.g, !TOPDIR) during the session without
;		    an initial call to this script, or when a limited number
;		    of the variables are being reset via a parameter file.
;		 -- upon execution, if FACTORY is not set, it will be set to 0
;	parfil	 if set to a named parameter file that contains new definitions
;		 of the defined system variables.  They must be placed one on
;		 each line, and must be fully legal IDL statements.
;		* if not set, will automatically look in
;		  -- $HOME/.pintofale if UNIX, pintofale.par if Windows
;		     used to also look in !ARDB/initale.par, but with the
;		     advent of getpoadef(), no longer does so.  The file
;		     !ARDB/initale.par is still present, as a template for
;		     user customization
;		* DO NOT try to define the following in PARFIL:
;		   !ABUND, !ATOM, !AMU, !FUNDAE, !FIP, !ROMAN, !PoA
;	yCHIANTI if set then checks to see whether the CHIANTI distribution is
;		 available for use and includes it if it is not.
;	yATOMDB  if set then checks to see whether the ATOMDB IDL routines
;		 are accesible and adds them if they are
;
;calls subroutines
;	GETABUND
;	INICON
;	PEASECOLR
;
;restrictions
;	subroutine WHICH works only in UNIX, but it is used only if none
;	  of the environment variables PINTofALE, PoA, SCARDIR, or SCAR
;	  are set.
;	Unlike normal IDL variables, system variables, once defined,
;	  cannot change their size, though their values may be changed.
;	  Thus, care must be taken to match the expected sizes of the
;	  following variables to the hardcoded defaults if they are to
;	  be redefined within PARFIL:
;	  !ABUND, !LOGT, !DEM, !ATOM, !AMU, !FUNDAE, !FIP, !ROMAN
;
;history
;	Vinay Kashyap(OctMM)
;	set FACTORY to 0 if undefined; changed !ABREF default; added !ARDB;
;	  !TOPDIR now automagic (VK; DecMM)
;	added !PoA, !CALDB, !LOGT, !DEM, and !AA (VK; JanMMI)
;	fixed bug with !TOPDIR becoming an array if none of the env vars
;	  were defined; change prompt for large !VERBOSE (VK FebMMI/Dan Dewey)
;	changed call to WHICH to WHEREIS (VK; Dec2001)
;	commented out call to WHEREIS by using the output of HELP, and now
;	  adds !TOPDIR/pro to !PATH if it doesn't exist; also check for
;	  other packages such as IDL-ASTRO and CHIANTI; more robust to
;	  OS filename dependencies for v5.3+ due to extensive use of
;	  FILEPATH function; added !ATOMDB, !APECDIR, and yATOMDB
;	  (VK; Jun2002)
;	added call to PEASECOLR (VK; Jul02)
;	some people have IDL-ASTRO installed under 'astrolib' (VK; Aug02)
;	added !ANGSTROM; made !AA robust to font choices (VK; Dec'02)
;	suppose environment variable is defined, but incorrectly?  i.e.,
;	  what if $PINTofALE, $SCARDIR, etc. do not exist? (VK; May'03)
;	if IDL_DIR is not defined, assumes that IDL library is in anything
;	  that has '...idl/lib' (VK; Jul'03)
;	corrected case where if env variable defining TOPDIR had a trailing
;	  "/", would miss matching it in !PATH and not reset the order of
;	  directories (VK; MarMMIV)
;	made a few of the system variables read-only (VK; May04)
;	made windows compatible, by and large (VK; Apr05)
;	forced IDL-Astro to move ahead of Chianti's astron directory
;	  in !PATH (VK; FebMMVII)
;	removed early call to strsplit() and forced restoration of
;	  correctly compiled versions of STRSPLIT, UNIQ, and STR_2_ARR
;	  (VK; JunMMVII)
;	bugfix: wasn't setting TOPDIR correctly when it had to be determined
;	  from the location of this file (VK; OctMMIX)
;	changed default for !IONEQF from Mazzotta et al. to Bryans et al. (2009)
;	  (VK; MayMMXIII)
;	removed pre-5.3 code; changed default for !IONEQF from Bryans et al.
;	  to chianti.ioneq; now first looks for inits in ~/.pintofale (*nix) or
;	  $CWD\pintofale.par (windows)
;	  (VK; MayMMXIV)
;	bug fix: some variables were not being set if not defined in initale.par
;	cleaned up a bit, fixed bugs due to calls to getpoadef, and now no longer
;	  loads !ARDB/initale.par by default (VK; SepMMXV)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;	initialize hard-coded defaults
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
forward_function filepath,getabund,getpoadef
zslash='/'
case !version.OS_FAMILY of
  'unix': zslash='/'
  'Windows': zslash='\'
  ;'macos': zslash=':'
  'vms': zslash='\'
  else: message,!version.OS_FAMILY+': operating system not understood;'+$
	' assuming UNIX',/info,/noname
endcase

;	first things first -- get TOPDIR and make sure it is in !path and if not, add it
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;{	old code, left here only for reference
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;	figure out installation path
;tmp=file_search(ztopdir,count=nn) & if nn eq 0 then ztopdir=''	;reset the default and try to work it out from location of initale.pro
;if (findfile(ztopdir))[0] eq '' then ztopdir=''
;if not keyword_set(ztopdir) then begin
;  ;	surely we now know where *this* file resides
;  help,/source_files,output=scardir
;  zpath=filepath('initale.pro',root_dir='pro',subdirectory='scrypt')
;  for i=0,n_elements(scardir)-1 do begin
;    ivar=strpos(scardir[i],zpath,0)
;    if ivar ge 0 then begin
;      ;jvar=strpos(scardir[i],zslash,0)-1 > 0
;      ;ztopdir=strmid(scardir[i],jvar,ivar-jvar-1L)
;      ;	the following block uses strsplit,
;      ;	which conflicts with SSW and CHIANTI,
;      ;	so we will try something different:{
;      ;		cc=strsplit(scardir[i],/extract)
;      ;		ztopdir=strmid(cc[1],0,strpos(cc[1],zpath)-1)
;      ;	:}{
;      jvar=strpos(scardir[i],' ',/reverse_search)
;      ztopdir=strmid(scardir[i],jvar+1,strlen(scardir[i])-jvar-(strlen(zpath)+1))
;      ;	}
;    endif
;    ;if ivar ge 0 then begin
;    ;  case !version.OS_FAMILY of
;    ;	'unix': jvar=strpos(scardir[i],'/',0)
;    ;	'Windows': jvar=strpos(scardir[i],':',0)-1 > 0
;    ;	'macos': jvar=strpos(scardir[i],':',0)-1 > 0
;    ;	else: begin
;    ;	  message,!version.OS_FAMILY+': operating system not understood.'+$
;    ;		'  assuming UNIX',/info,/noname
;    ;	  jvar=strpos(scardir[i],'/',0)
;    ;	end
;    ;  endcase
;    ;  ztopdir=strmid(scardir[i],jvar,ivar-jvar-1L)
;    ;endif
;  endfor
;endif
;ztopdir=strtrim(ztopdir,2)
;;if not keyword_set(ztopdir) then begin
;;  scardir=whereis('initale.pro',/dironly)
;;  ivar=strpos(scardir[0],'pro/scrypt')  ;we know where INITALE.PRO resides
;;  ztopdir=strmid(scardir[0],0,ivar-1)
;;endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	old code above}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nn=0 & if keyword_set(ztopdir) then tmp=file_search(ztopdir,count=nn)
if nn eq 0 then ztopdir=''	;reset the default and try to work it out from location of initale.pro
if not keyword_set(ztopdir) then begin
  ;	surely we now know where *this* file resides
  help,/source,output=scardir
  zpath=filepath('initale.pro',root_dir='pro',subdirectory='scrypt')
  for i=0,n_elements(scardir)-1 do begin
    ivar=strpos(scardir[i],zpath,0)
    if ivar ge 0 then begin
      jvar=strpos(scardir[i],' ',/reverse_search)
      ztopdir=strmid(scardir[i],jvar+1,strlen(scardir[i])-jvar-(strlen(zpath)+1))
    endif
  endfor
endif
ztopdir=strtrim(ztopdir,2)
;	add TOPDIR/pro to !path
if strpos(!path,ztopdir) lt 0 then !path = expand_path('+'+filepath('',root=ztopdir,subdir='pro'))+':'+!path

;	now define everything that can be set using environmental variables or GETPOADEF()
defsysv,'!PoA',getpoadef('VERSION'),1
if not keyword_set(ztopdir) then ztopdir=getenv('PINTofALE')
if not keyword_set(ztopdir) then ztopdir=getenv('PoA')
if not keyword_set(ztopdir) then ztopdir=getenv('SCARDIR')
if not keyword_set(ztopdir) then ztopdir=getenv('SCAR')
if not keyword_set(ztopdir) then ztopdir=getpoadef('TOPDIR')
if not keyword_set(zldbdir) then zldbdir=getenv('LDBDIR')
if not keyword_set(zldbdir) then zldbdir=getpoadef('LDBDIR')
if not keyword_set(zcdbdir) then zcdbdir=getenv('CDBDIR')
if not keyword_set(zcdbdir) then zcdbdir=getpoadef('CDBDIR')
if not keyword_set(zardb) then zardb=getenv('ARDB')
if not keyword_set(zardb) then zardb=getpoadef('ARDB')
if not keyword_set(zceroot) then zceroot=getenv('CEROOT')
if not keyword_set(zceroot) then zceroot=getpoadef('CEROOT')
if not keyword_set(zchidir) then zchidir=getenv('CHIDIR')
if not keyword_set(zchidir) then zchidir=getpoadef('CHIDIR')
if not keyword_set(zioneqf) then zioneqf=getenv('IONEQF')
if not keyword_set(zioneqf) then zioneqf=getpoadef('IONEQF')
if not keyword_set(zatomdb) then zatomdb=getenv('ATOMDB')
if not keyword_set(zatomdb) then zatomdb=getpoadef('ATOMDB')
if keyword_set(zatomdb) then zapecdir=filepath('',root_dir=zatomdb,subdir='atomdb_idl-2.00')
tmp=file_search(zapecdir,count=nn) & if nn eq 0 then zapecdir=filepath('',root_dir=zatomdb,subdir='apec_v11_idl')
if not keyword_set(zabref) then zabref=getenv('ABREF')
if not keyword_set(zabref) then zabref=getpoadef('ABREF')
if not keyword_set(zcaldb) then zcaldb=getenv('CALDB')
if not keyword_set(zcaldb) then zcaldb=getpoadef('CALDB')

;	now define other variables of interest
zmetals=0.0
zgaspr=1e15 & zlogpr=15. & zedens=1e10		;[cm^-3]
zlogT=findgen(81)*0.05+4.			;[log[degK]]
zDEM=dblarr(n_elements(zlogT))+1d12		;[cm^-5, typically]
znh=1e18 & zfh2=0.26 & zhe1=1e17 & zheII=1e16	;[cm^-2]
zwmin=1.239854 & zwmax=3000.0			;[Ang]
zverbose=5
zaa=[	string(byte(197)),$
	'!3'+string(byte(197))+'!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!4!sA!r!u!9 %!n!X',$
	'!5!sA!r!u!9 %!n!X',$
	'!6!sA!r!u!9 %!n!X',$
	'!7!sA!r!u!9 %!n!X',$
	'!8!sA!r!u!9  %!n!X',$
	'!9!s3!r!u!9  %!n!X',$
	'!10!sA!r!u!9 %!n!X',$
	'!11!sA!r!u!9  %!n!X',$
	'!12!sA!r!u!9  %!n!X',$
	'!13!sA!r!u!9   %!n!X',$
	'!14!sA!r!u!9  %!n!X',$
	'!15!sA!r!u!9  %!n!X',$
	'!16!sA!r!u!9 %!n!X',$
	'!17!sA!r!u!9 %!n!X',$
	'!18!sA!r!u!9  %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!20!s@!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X',$
	'!3!sA!r!u!9 %!n!X' ]

;	check if the variables are already defined
;	if not defined, then define.
;	if defined, then redefine iff FACTORY=1
zfactory=1
if n_elements(factory) gt 0 and not keyword_set(factory) then zfactory=0
if keyword_set(parfil) then zfactory=0
;
ivar=0 & defsysv,'!PoA',exists=ivar
if (ivar eq 0) then defsysv,'!PoA',getpoadef('VERSION'),1
ivar=0 & defsysv,'!TOPDIR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!TOPDIR',zTOPDIR
ivar=0 & defsysv,'!LDBDIR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!LDBDIR',zLDBDIR
ivar=0 & defsysv,'!CDBDIR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!CDBDIR',zCDBDIR
ivar=0 & defsysv,'!CEROOT',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!CEROOT',zCEROOT
ivar=0 & defsysv,'!CHIDIR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!CHIDIR',zCHIDIR
ivar=0 & defsysv,'!ATOMDB',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!ATOMDB',zATOMDB
ivar=0 & defsysv,'!APECDIR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!APECDIR',zAPECDIR
ivar=0 & defsysv,'!ARDB',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!ARDB',zARDB
ivar=0 & defsysv,'!IONEQF',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!IONEQF',zIONEQF
ivar=0 & defsysv,'!ABREF',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!ABREF',zABREF
ivar=0 & defsysv,'!CALDB',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!CALDB',zCALDB
ivar=0 & defsysv,'!METALS',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!METALS',zMETALS
ivar=0 & defsysv,'!GASPR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!GASPR',zGASPR
ivar=0 & defsysv,'!LOGPR',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!LOGPR',zLOGPR
ivar=0 & defsysv,'!EDENS',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!EDENS',zEDENS
ivar=0 & defsysv,'!LOGT',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!LOGT',zLOGT
ivar=0 & defsysv,'!DEM',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!DEM',zDEM
ivar=0 & defsysv,'!NH',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!NH',zNH
ivar=0 & defsysv,'!FH2',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!FH2',zFH2
ivar=0 & defsysv,'!HE1',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!HE1',zHE1
ivar=0 & defsysv,'!HEII',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!HEII',zHEII
ivar=0 & defsysv,'!WMIN',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!WMIN',zWMIN
ivar=0 & defsysv,'!WMAX',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!WMAX',zWMAX
ivar=0 & defsysv,'!VERBOSE',exists=ivar
if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!VERBOSE',zVERBOSE
ivar=0 & defsysv,'!AA',exists=ivar
if ivar eq 0 then defsysv,'!AA','!3'+string(byte(197))+'!X'
ivar=0 & defsysv,'!ANGSTROM',exists=ivar
if ivar eq 0 then defsysv,'!ANGSTROM',zaa

;	override any with PARFIL?
zparfil=0
;	force PARFIL to be read from
;		$HOME/.pintofale if UNIX, pintofale.par if Windows
;		!TOPDIR/ardb/initale.par
;	need to jump through EXECUTE to avoid calling !TOPDIR before it
;	may have been defined, and FILEPATH to make it OS independent
if not keyword_set(parfil) and zfactory eq 1 then begin
  if !version.OS_FAMILY eq 'unix' then homefil=getenv('HOME')+'/.pintofale' else homefil='pintofale.par'
  tmp=file_search(homefil,count=nrc)
  if nrc gt 0 then zparfil=tmp[0]
  ;	{used to fall back on !ARDB/initale.par
  ;	not necessary any more
  ;	else zparfil=filepath('initale.par',root_dir=!TOPDIR,subdirectory='ardb')
  ;	}
endif
;	read from user-defined PARFIL
if keyword_set(parfil) then zparfil=parfil
if keyword_set(zparfil) then begin
  fil=findfile(zparfil,count=nfil)
  if nfil eq 0 then begin
    if !VERBOSE gt 0 then message,'PARFIL is illegible; ignoring',/info,/noname
  endif else begin
    if !VERBOSE gt 1 then message,'reading from: '+fil[0],/info,/noname
    openr,upar,fil[0],/get_lun		;{open PARFIL for input
      while not eof(upar) do begin
        zline='' & readf,upar,zline
        ivar=execute(zline)	;beware that this will override ANY defaults set by getpoadef()
      endwhile
    close,upar & free_lun,upar		;close PARFIL}
  endelse
endif

;	check for IDL-ASTRO
if strpos(!PATH,'astron',0) lt 0 and strpos(!PATH,'astrolib',0) lt 0 then begin
  message,'WARNING: IDL-ASTRO library is not in !PATH',/info,/noname,/noprefix
endif else begin
  ivar=0 & defsysv,'!TEXTOUT',exists=ivar
  if ivar eq 0 then astrolib
endelse

;	check for and add CHIANTI if asked for
ivar=strpos(!CHIDIR,'dbase')
;if ivar lt 0 then ivar=strlen(!CHIDIR)-1L
if ivar lt 0 then ivar=strlen(!CHIDIR)
zchidir=strmid(!CHIDIR,0,ivar)
if strpos(!PATH,zCHIDIR,0) lt 0 then begin
  if keyword_set(yCHIANTI) then begin
    message,'Appending CHIANTI distribution to !PATH',/info,/noname
    if !version.os_family eq 'Windows' then junk=';' else $
     if !version.os_family eq 'vms' then junk=',' else junk=':'
    !path=expand_path(!path+junk+'+'+zchidir)
    junk=findfile(filepath('use_chianti.pro',root_dir=zchidir,subdirectory=['idl','SETUP']),count=ivar)
    if ivar ne 0 then use_chianti,!chidir	;for CHIANTI 4 and above
  endif
endif

;	check for and add ATOMDB if asked for
if strpos(!PATH,!APECDIR,0) lt 0 then begin
  if keyword_set(yATOMDB) then begin
    fil1=filepath('init_atomdb_idl.pro',root_dir=!ATOMDB,subdirectory='atomdb_idl-2.00')
    fil2=filepath('init_atomdb_idl.pro',root_dir=!APECDIR)
    fil3=filepath('init_atomdb_idl.pro',root_dir=!TOPDIR,subdirectory=['atomdb','atomdb_idl-2.00'])
    ;
    fil1o=filepath('start.pro',root_dir=!ATOMDB,subdirectory='apec_v11_idl')
    fil2o=filepath('start.pro',root_dir=!APECDIR)
    fil3o=filepath('start.pro',root_dir=!TOPDIR,subdirectory=['atomdb','apec_v11_idl'])
    ;
    fdir1=filepath('',root_dir=!ATOMDB,subdirectory='atomdb_idl-2.00')
    fdir2=filepath('',root_dir=!APECDIR)
    fdir3=filepath('',root_dir=!TOPDIR,subdirectory=['atomdb','atomdb_idl-2.00'])
    ;
    fdir1o=filepath('',root_dir=!ATOMDB,subdirectory='apec_v11_idl')
    fdir2o=filepath('',root_dir=!APECDIR)
    fdir3o=filepath('',root_dir=!TOPDIR,subdirectory=['atomdb','apec_v11_idl'])
    ;
    if not keyword_set(ivar1) then junk=findfile(fil1,count=ivar1)
    if not keyword_set(ivar2) then junk=findfile(fil2,count=ivar2)
    if not keyword_set(ivar3) then junk=findfile(fil3,count=ivar3)
    ;
    if not keyword_set(ivar1o) then junk=findfile(fil1o,count=ivar1o)
    if not keyword_set(ivar2o) then junk=findfile(fil2o,count=ivar2o)
    if not keyword_set(ivar3o) then junk=findfile(fil3o,count=ivar3o)
    ;
    if ivar1 ne 0 then zapecdir=fdir1 else begin
      if ivar2 ne 0 then zapecdir=fdir2 else $
       if ivar3 ne 0 then zapecdir=fdir3 else $
        if ivar1o ne 0 then zapecdir=fdir1o else $
         if ivar2o ne 0 then zapecdir=fdir2o else $
          if ivar3o ne 0 then zapecdir=fdir3o
      message,'WARNING: Resetting !APECDIR to '+zapecdir,/info,/noname
      !APECDIR=zapecdir
    endelse
    if ivar1 ne 0 or ivar2 ne 0 or ivar3 ne 0 or ivar1o ne 0 or ivar2o ne 0 or ivar3o ne 0 then begin
      message,'Appending APED distribution to !PATH',/info,/noname
      if !version.os_family eq 'Windows' then junk=';' else $
       if !version.os_family eq 'vms' then junk=',' else junk=':'
      !path=expand_path(!path+junk+'+'+!APECDIR)
    endif else begin
      message,!APECDIR+': does not exist -- not including ATOMDB',/info
    endelse
  endif
endif

;	if PINTofALE is not in !PATH, then add it
if strpos(!PATH,!TOPDIR+zslash+'pro',0) lt 0 and $
   strpos(!PATH,!TOPDIR+'pro',0) lt 0 then begin
  if !version.os_family eq 'Windows' then junk=';' else $
   if !version.os_family eq 'vms' then junk=',' else junk=':'
  !path=expand_path('+'+filepath('pro',root_dir=!TOPDIR)+junk+!path)
endif

;	rearrange the !path in the following order:
;		$IDL_DIR/lib, IDL-Astro,
;		PoA,
;		CHIANTI's astron, rest of CHIANTI,
;		everything else
junk=expand_path(!PATH,/array) & ipath=lindgen(n_elements(junk))
zlib=getenv('IDL_DIR') & if not keyword_set(zlib) then zlib='idl'
zlib=zlib+zslash+'lib'
zo1=where(strpos(junk,zlib,0) ge 0,izo1)
;	NOTE: this will force IDL to put the main IDL-Astro library
;	ahead of Chianti's bastardized version in the !PATH
;	BEWARE that this may cause Chianti routines to fail in
;	this environment.  If you need to use those Chianti routines,
;	then run IDL without initializing PINTofALE with this script.
zo1a=where(strpos(junk,'astron',0) ge 0 or strpos(junk,'astrolib',0) ge 0,izo1a)
zo2=where(strpos(junk,!TOPDIR+zslash+'pro',0) ge 0 or strpos(junk,!TOPDIR+'pro',0) ge 0,izo2)
zo3=where(strpos(junk,zchidir,0) ge 0 and strpos(junk,'astron',0) ge 0,izo3)
zo4=where(strpos(junk,zchidir,0) ge 0 and strpos(junk,'astron',0) lt 0,izo4)
if izo1 gt 0 then ipath[zo1]=-1
if izo1a gt 0 then ipath[zo1a]=-1
if izo2 gt 0 then ipath[zo2]=-2
if izo3 gt 0 then ipath[zo3]=-3
if izo4 gt 0 then ipath[zo4]=-4
zo5=where(ipath ge 0,izo5) & zpath=['']
if izo5 gt 0 then zpath=[junk[zo5],zpath]
if izo4 gt 0 then zpath=[junk[zo4],zpath]
if izo3 gt 0 then zpath=[junk[zo3],zpath]
if izo2 gt 0 then zpath=[junk[zo2],zpath]
if izo1a gt 0 then zpath=[junk[zo1a],zpath]
if izo1 gt 0 then zpath=[junk[zo1],zpath]
ivar=n_elements(zpath) & zpath=zpath[0L:ivar-2L]
if !version.os_family eq 'Windows' then junk=';' else $
 if !version.os_family eq 'vms' then junk=',' else junk=':'
!path=zpath[0] & for i=1L,ivar-2L do !path=!path+junk+zpath[i]

;skip this if running GDL
iGDL=0 & defsysv,'!GDL',exists=iGDL
if not keyword_set(iGDL) then begin
;	there are some conflicting routine names, for example in SSW
;	if these routines have been compiled, catch them and ruthlessly
;	recompile them
;	examples: STRSPLIT, STR_2_ARR, UNIQ
;resolve_routine,filepath('str_2_arr.pro',root_dir=!TOPDIR,subdirectory=['pro','misc'])
message,'recompiling some offending subroutines:',/informational,/noname
if strpos(!ARDB,!TOPDIR) eq 0 then $
  restore,filepath('strsplit_uniq.sav',root_dir=!ARDB),/verbose else $
  restore,filepath('strsplit_uniq.sav',root_dir=!TOPDIR,subdirectory=!ARDB),/verbose
	;override existing compilations of STRSPLIT and UNIQ
endif else message,'GDL cannot restore compiled programs,'+$
	' so SSW-recompiled routines like STRSPLIT and UNIQ cannot be overridden',/information

;	define abundances
zABUND=getabund(!ABREF) & zabund[2:*]=zabund[2:*]*10.^(!METALS)
defsysv,'!ABUND',zABUND
;ivar=0 & defsysv,'!ABUND',exists=ivar
;if (ivar eq 0) or (ivar eq 1 and zfactory eq 1) then defsysv,'!ABUND',zABUND

inicon,atom=zatom,amu=zamu,roman=zroman,fundae=zfundae,fip=zfip
ivar=0 & defsysv,'!ATOM',exists=ivar,1
if ivar eq 0 then defsysv,'!ATOM',zATOM
ivar=0 & defsysv,'!AMU',exists=ivar,1
if ivar eq 0 then defsysv,'!AMU',zAMU
ivar=0 & defsysv,'!ROMAN',exists=ivar,1
if ivar eq 0 then defsysv,'!ROMAN',zROMAN
ivar=0 & defsysv,'!FUNDAE',exists=ivar,1
if ivar eq 0 then defsysv,'!FUNDAE',zFUNDAE
ivar=0 & defsysv,'!FIP',exists=ivar,1
if ivar eq 0 then defsysv,'!FIP',zFIP

;	prepend !TOPDIR as needed
junk='ok'
;if strpos(!LDBDIR,'$') eq 0 then junk=!LDBDIR+' is in environment'
if strpos(!LDBDIR,'$') eq 0 then !LDBDIR=getenv(strmid(!LDBDIR,1))
if strpos(!LDBDIR,!TOPDIR) eq 0 then junk=!LDBDIR+' already contains '+!TOPDIR
if strpos(!LDBDIR,zslash) eq 0 then junk=!LDBDIR+' is full pathname'
if junk ne 'ok' then message,junk,/info,/noname else $
 !LDBDIR=filepath(!LDBDIR,root_dir=!TOPDIR)
junk='ok'
;if strpos(!CDBDIR,'$') eq 0 then junk=!CDBDIR+' is in environment'
if strpos(!CDBDIR,'$') eq 0 then !CDBDIR=getenv(strmid(!CDBDIR,1))
if strpos(!CDBDIR,!TOPDIR) eq 0 then junk=!CDBDIR+' already contains '+!TOPDIR
if strpos(!CDBDIR,zslash) eq 0 then junk=!CDBDIR+' is full pathname'
if junk ne 'ok' then message,junk,/info,/noname else $
 !CDBDIR=filepath(!CDBDIR,root_dir=!TOPDIR)
junk='ok'
;if strpos(!CHIDIR,'$') eq 0 then junk=!CHIDIR+' is in environment'
if strpos(!CHIDIR,'$') eq 0 then !CHIDIR=getenv(strmid(!CHIDIR,1))
if strpos(!CHIDIR,!TOPDIR) eq 0 then junk=!CHIDIR+' already contains '+!TOPDIR
if strpos(!CHIDIR,zslash) eq 0 then junk=!CHIDIR+' is full pathname'
if junk ne 'ok' then message,junk,/info,/noname else $
 !CHIDIR=filepath(!CHIDIR,root_dir=!TOPDIR)
junk='ok'
;if strpos(!ATOMDB,'$') eq 0 then junk=!ATOMDB+' is in environment'
if strpos(!ATOMDB,'$') eq 0 then !ATOMDB=getenv(strmid(!ATOMDB,1))
if strpos(!ATOMDB,!TOPDIR) eq 0 then junk=!ATOMDB+' already contains '+!TOPDIR
if strpos(!ATOMDB,zslash) eq 0 then junk=!ATOMDB+' is full pathname'
if junk ne 'ok' then message,junk,/info,/noname else $
 !ATOMDB=filepath(!ATOMDB,root_dir=!TOPDIR)
junk='ok'
;if strpos(!APECDIR,'$') eq 0 then junk=!APECDIR+' is in environment'
if strpos(!APECDIR,'$') eq 0 then !APECDIR=getenv(strmid(!APECDIR,1))
if strpos(!APECDIR,!TOPDIR) eq 0 then junk=!APECDIR+' already contains '+!TOPDIR
if strpos(!APECDIR,zslash) eq 0 then junk=!APECDIR+' is full pathname'
if junk ne 'ok' then message,junk,/info,/noname else $
 !APECDIR=filepath(!APECDIR,root_dir=!TOPDIR)
junk='ok'
;if strpos(!ARDB,'$') eq 0 then junk=!ARDB+' is in environment'
if strpos(!ARDB,'$') eq 0 then !ARDB=getenv(strmid(!ARDB,1))
if strpos(!ARDB,!TOPDIR) eq 0 then junk=!ARDB+' already contains '+!TOPDIR
if strpos(!ARDB,zslash) eq 0 then junk=!ARDB+' is full pathname'
if junk ne 'ok' then message,junk,/info,/noname else $
 !ARDB=filepath(!ARDB,root_dir=!TOPDIR)

;	set factory if previously unset
if n_elements(factory) eq 0 then factory=0

;	make sure things will work with 24-bit color displays
if !D.N_COLORS gt 256 then device,decomposed=0
;	and load color table with some useful colors
loadct,3	;standard "heat" at high color numbers
peasecolr	;specific colors at small color numbers

if !VERBOSE ge 10 then begin
  print,'' & message,'Defined system variables:',/informational & print,''
  print,"!TOPDIR='"+!TOPDIR+"'"
  print,"!LDBDIR='"+!LDBDIR+"'"
  print,"!CDBDIR='"+!CDBDIR+"'"
  print,"!CHIDIR='"+!CHIDIR+"'"
  print,"!ATOMDB='"+!ATOMDB+"'"
  print,"!APECDIR='"+!APECDIR+"'"
  print,"!ARDB='"+!ARDB+"'"
  print,"!CEROOT='"+!CEROOT+"'"
  print,"!IONEQF='"+!IONEQF+"'"
  print,"!ABREF='"+!ABREF+"'"
  print,"!CALDB='"+!CALDB+"'"
  print,"!METALS="+strtrim(!METALS,2)
  print,"!GASPR="+strtrim(!GASPR,2)
  print,"!LOGPR="+strtrim(!LOGPR,2)
  print,"!EDENS="+strtrim(!EDENS,2)
  print,"!NH="+strtrim(!NH,2)
  print,"!WMIN="+strtrim(!WMIN,2)
  print,"!WMAX="+strtrim(!WMAX,2)
  print,"!VERBOSE="+strtrim(!VERBOSE,2)
  print,"..."
  print,"!ABUND,!LOGT,!DEM,!ATOM,!AMU,!FUNDAE,!FIP,!ROMAN,!AA,!ANGSTROM,!FH2,!HE1,!HEII"
endif

;	declare
print,'' & print,'	PINTofALE v'+!PoA & print,''
if !VERBOSE gt 0 then begin
  print,'	Package for Interactive Analysis of Line Emission'
  print,'	Kashyap, V., \& Drake, J.J.\ 2000, Bull.Astr.Soc.India, 28, 475'
  if !VERBOSE gt 1 then $
    print,'	[ pintofale@head.cfa.harvard.edu ]'
endif

;	reset prompt for the nonwary
if !VERBOSE ge 5 then !PROMPT='PoA> '

;	clean up intermediate variables
;	(this trick of using temporary variables thanks to Pavel Romashkin)
;	(the alternative is to use
;		delvar,scardir,ztopdir,zldbdir,zcdbdir,zceroot,zchidir,$
;		zioneqf,zabref,zmetals,zabund,zgaspr,zlogpr,zedens,znh,$
;		zfh2,zhe1,zheII,zwmin,zwmax,zverbose,zatom,zamu,zroman,$
;		zfundae,zfip,zfactory,ivar,zparfil,fil,nfil,upar,zline,$
;		zardb,zcaldb,zlogT,zDEM, etc.
;	which however screws up the IDL environment if .RESET_SESSION is
;	used at any time afterwards.)
if keyword_set(scardir) then junk=temporary(scardir)
junk=temporary(zslash)
junk=temporary(ztopdir) & junk=temporary(zldbdir)
junk=temporary(zcdbdir) & junk=temporary(zceroot) & junk=temporary(zchidir)
junk=temporary(zapecdir) & junk=temporary(zatomdb)
junk=temporary(zioneqf) & junk=temporary(zabref) & junk=temporary(zmetals)
junk=temporary(zabund) & junk=temporary(zgaspr) & junk=temporary(zlogpr)
junk=temporary(zedens) & junk=temporary(znh) & junk=temporary(zfh2)
junk=temporary(zhe1) & junk=temporary(zheII) & junk=temporary(zwmin)
junk=temporary(zwmax) & junk=temporary(zverbose) & junk=temporary(zatom)
junk=temporary(zamu) & junk=temporary(zroman) & junk=temporary(zfundae)
junk=temporary(zfip) & junk=temporary(zfactory) & junk=temporary(ivar)
if keyword_set(jvar) then junk=temporary(jvar)
fil=0 & fil1=0 & fil2=0 & fil3=0 & fil4=0 & i=0 & ivar1=0 & ivar2=0 & ivar3=0 & ivar4=0
junk=temporary(fil) & junk=temporary(i)
junk=temporary(fil1) & junk=temporary(fil2) & junk=temporary(fil3) & junk=temporary(fil4)
junk=temporary(ivar1) & junk=temporary(ivar2) & junk=temporary(ivar3) & junk=temporary(ivar4)
junk=temporary(zparfil)
if keyword_set(nfil) then junk=temporary(nfil)
if keyword_set(upar) then junk=temporary(upar)
if keyword_set(zline) then junk=temporary(zline)
junk=temporary(zardb) & junk=temporary(zcaldb)
junk=temporary(zdem) & junk=temporary(zlogt)
junk=temporary(zpath) & junk=temporary(zlib)
junk=temporary(ipath) & junk=temporary(izo1)
junk=temporary(izo2) & junk=temporary(izo3) & junk=temporary(izo4)
junk=temporary(izo5) & junk=temporary(zo1) & junk=temporary(zo2)
junk=temporary(zo3) & junk=temporary(zo4) & junk=temporary(zo5)
junk=temporary(zaa)
junk=0

end
