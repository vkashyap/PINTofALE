;+
;script	PoA_constructor
;	this is a "front-end" to the initialization script INITALE of
;	PINTofALE, needed for older versions of IDL, for which system
;	variables need to be defined before running ANY procedure which
;	uses them.
;	
;	this script defines the system variables with default values
;	and then calls INITALE.PRO to complete the initialization.
;
;syntax
;	@PoA_constructor
;
;variables
;	for a complete list of the PINTofALE system variables, see
;	the description for INITALE.PRO .
;
;calls subroutines
;	INICON
;	INITALE
;
;restrictions
;	it is assumed that this is called only once per session, and
;	that it is called _first_
;
;history
;	Vinay Kashyap (FebMMI; based on suggestion by A.Maggio)
;	added !ATOMDB and !APECDIR (VK; Jun02)
;-

ivar=0 & defsysv,'!PoA',exists=ivar
if (ivar eq 0) then defsysv,'!PoA',2.61
ivar=0 & defsysv,'!TOPDIR',exists=ivar
if (ivar eq 0) then defsysv,'!TOPDIR','/foo/bar/PINTofALE/'
ivar=0 & defsysv,'!LDBDIR',exists=ivar
if (ivar eq 0) then defsysv,'!LDBDIR','emissivity/chianti'
ivar=0 & defsysv,'!CDBDIR',exists=ivar
if (ivar eq 0) then defsysv,'!CDBDIR','emissivity/cont'
ivar=0 & defsysv,'!CEROOT',exists=ivar
if (ivar eq 0) then defsysv,'!CEROOT','cie'
ivar=0 & defsysv,'!CHIDIR',exists=ivar
if (ivar eq 0) then defsysv,'!CHIDIR','CHIANTI/dbase/'
ivar=0 & defsysv,'!ATOMDB',exists=ivar
if (ivar eq 0) then defsysv,'!ATOMDB','/data/atomdb/'
ivar=0 & defsysv,'!APECDIR',exists=ivar
if (ivar eq 0) then defsysv,'!APECDIR','apec_v11_idl'
ivar=0 & defsysv,'!ARDB',exists=ivar
if (ivar eq 0) then defsysv,'!ARDB','ardb'
ivar=0 & defsysv,'!IONEQF',exists=ivar
if (ivar eq 0) then defsysv,'!IONEQF','ioneq/mazzotta_etal.ioneq'
ivar=0 & defsysv,'!ABREF',exists=ivar
if (ivar eq 0) then defsysv,'!ABREF','grevesse et al.'
ivar=0 & defsysv,'!CALDB',exists=ivar
if (ivar eq 0) then defsysv,'!CALDB','/data/caldb/'
ivar=0 & defsysv,'!METALS',exists=ivar
if (ivar eq 0) then defsysv,'!METALS',0.0
ivar=0 & defsysv,'!GASPR',exists=ivar
if (ivar eq 0) then defsysv,'!GASPR',1e15
ivar=0 & defsysv,'!LOGPR',exists=ivar
if (ivar eq 0) then defsysv,'!LOGPR',15.
ivar=0 & defsysv,'!EDENS',exists=ivar
if (ivar eq 0) then defsysv,'!EDENS',1e10
ivar=0 & defsysv,'!LOGT',exists=ivar
if (ivar eq 0) then defsysv,'!LOGT',findgen(81)*0.05+4.
ivar=0 & defsysv,'!DEM',exists=ivar
if (ivar eq 0) then defsysv,'!DEM',dblarr(81)+1d12
ivar=0 & defsysv,'!NH',exists=ivar
if (ivar eq 0) then defsysv,'!NH',1e18
ivar=0 & defsysv,'!FH2',exists=ivar
if (ivar eq 0) then defsysv,'!FH2',0.26
ivar=0 & defsysv,'!HE1',exists=ivar
if (ivar eq 0) then defsysv,'!HE1',1e17
ivar=0 & defsysv,'!HEII',exists=ivar
if (ivar eq 0) then defsysv,'!HEII',1e16
ivar=0 & defsysv,'!WMIN',exists=ivar
if (ivar eq 0) then defsysv,'!WMIN',1.239854
ivar=0 & defsysv,'!WMAX',exists=ivar
if (ivar eq 0) then defsysv,'!WMAX',3000.0
ivar=0 & defsysv,'!VERBOSE',exists=ivar
if (ivar eq 0) then defsysv,'!VERBOSE',5
ivar=0 & defsysv,'!ABUND',exists=ivar
if (ivar eq 0) then defsysv,'!ABUND',getabund('grevesse et al.')
;
inicon,atom=zatom,amu=zamu,roman=zroman,fundae=zfundae,fip=zfip
;
ivar=0 & defsysv,'!ATOM',exists=ivar
if ivar eq 0 then defsysv,'!ATOM',zATOM
ivar=0 & defsysv,'!AMU',exists=ivar
if ivar eq 0 then defsysv,'!AMU',zAMU
ivar=0 & defsysv,'!ROMAN',exists=ivar
if ivar eq 0 then defsysv,'!ROMAN',zROMAN
ivar=0 & defsysv,'!FUNDAE',exists=ivar
if ivar eq 0 then defsysv,'!FUNDAE',zFUNDAE
ivar=0 & defsysv,'!FIP',exists=ivar
if ivar eq 0 then defsysv,'!FIP',zFIP
ivar=0 & defsysv,'!AA',exists=ivar
if ivar eq 0 then defsysv,'!AA',string(byte(197))

.run initale
