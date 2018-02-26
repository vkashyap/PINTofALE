pro initstuff,atom,rom,funcon
;+
;procedure	initstuff
;	return arrays of atomic symbols, appropriate number of roman numerals,
;	fundamental and other fun constants
;
;syntax
;	initstuff,atom,rom,funcon
;
;parameters
;	atom	[OUTPUT] atomic symbols, as many as there are elements in
;		abundance array
;	rom	[OUTPUT] roman numerals, upto max(ATOM)+1
;	funcon	[OUTPUT] structure containing some useful physical and
;		astronomical constants
;		* just do help,funcon,/str -- should be obvious
;
;keywords	NONE
;
;requires
;	CREATE_STRUCT
;
;history
;	vinay kashyap (Aug98)
;	removed call to GETABUND and replaced ABUND and HINT by FUNCON
;	  (VK; Apr99)
;-

message,'OBSOLETE; use INICON instead',/informational

;	usage
np=n_params()
if np eq 0 then begin
  print,'Usage: initstuff,atom,rom,funcon'
  print,'  strictly utility routine returns atomic symbols, roman numerals'
  print,'  physical and astronomical constants'
  return
endif

;	elements from 1-99
Z00=[	'H', 'He','Li','Be','B', 'C', 'N', 'O', 'F' ]		;001-009
Z01=[	'Ne','Na','Mg','Al','Si','P', 'S', 'Cl','Ar','K' ]	;010-019
Z02=[	'Ca','Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu']	;020-029
Z03=[	'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ]	;030-039
Z04=[	'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In']	;040-049
Z05=[	'Sn','Sb','Te','I', 'Xe','Cs','Ba','La','Ce','Pr']	;050-059
Z06=[	'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm']	;060-069
Z07=[	'Yb','Lu','Hf','Ta','W', 'Re','Os','Ir','Pt','Au']	;070-079
Z08=[	'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac']	;080-089
Z09=[	'Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es']	;090-099
atom=[	Z00, Z01, Z02, Z03, Z04, Z05, Z06, Z07, Z08, Z09]
nZ=n_elements(atom)

;	roman numerals from 1-100
rom=[	'I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[	rom,'X'+rom,'XX'+rom,'XXX'+rom, 'XXXXI','XXXXII','XXXXIII','XXXXIV',$
	'XXXXV','XXXXVI','XXXXVII','XXXXVIII']
rom=[rom,'IL','L','L'+rom,'IC','C']
nrom=n_elements(rom)

;	physical and astronomical constants

c=2.99792458d10  	& c_h='Speed of light [cm/s]'
s=create_struct('c',c) & SH=create_struct('c',c_h)

h=6.626176d-27  	& h_h="Planck's constant [erg/s]"
s=create_struct(s,'h',h) & SH=create_struct(SH,'h',h_h)

G=6.672e-8      	& G_h='Gravitational constant [cm^3/gm/s^2]'
s=create_struct(s,'G',G) & SH=create_struct(SH,'G',G_h)

Me=9.109558d-28  	& Me_h='Electron mass [gm]'
Mp=1.67248d-24  	& Mp_h='Proton mass [gm]'
s=create_struct(s,'Me',Me,'Mp',Mp) & SH=create_struct(SH,'Me',Me_h,'Mp',Mp_h)

e=1.6021892d-19  	& e_h='electron charge [C]'
esu=4.803d-10      	& esu_h='electron charge [ESU]'
s=create_struct(s,'e',e,'esu',esu) & SH=create_struct(SH,'e',e_h,'esu',esu_h)

kB=1.380662d-16  	& kB_h="Boltzmann's constant [erg/K]"
s=create_struct(s,'kB',kB,'k',kB) & SH=create_struct(SH,'kB',kB_h,'k',kB_h)

loga=alog10(8./15.)+5.*alog10(!DPI)+4.*alog10(kB)-3.*alog10(c)-3.*alog10(h)
a=10.^(loga)       	& a_h='Radiation Pressure constant [erg/cm^3/deg^4]'
sigma=a*c/4.    	& sigma_h='Stefan-Boltzmann constant [erg/cm^2/deg^4/s]'
wein=0.289780D  	& wein_h='Wein Displacement Law constant [cm K]'
s=create_struct(s,'a',a,'sigma',sigma,'wein',wein)
SH=create_struct(SH,'a',a_h,'sigma',sigma_h,'wein',wein_h)

NA=6.0249d23    	& NA_h='Avogadro number [/mole]'
s=create_struct(s,'NA',NA) & SH=create_struct(SH,'NA',NA_h)

atm=1013250.D   	& atm_h='1 Atmosphere [dynes/cm^2]'
s=create_struct(s,'atm',atm) & SH=create_struct(SH,'atm',atm_h)

eVwav=12379.7d-8	& eVwav_h='1 eV in wave numbers [/cm]'
degeV=8.6173468d-05	& degeV_h='1 deg K in eV [eV]'
radian=3600.*180./!DPI 	& radian_h='1 radian [arcsec]'
s=create_struct(s,'eVwav',eVwav,'degeV',degeV,'radian',radian)
SH=create_struct(SH,'eVwav',eVwav_h,'degeV',degeV_h,'radian',radian_h)

Ryd=109677.58*eVWav   	& Ryd_h='Rydberg Constant for H [eV]'
s=create_struct(s,'Ryd',Ryd) & SH=create_struct(SH,'Ryd',Ryd_h)

RBohr=2.*alog10(h)-alog10(4.*!DPI^2)-alog10(Me)-2.*alog10(esu)
RBohr=10.^(RBohr)	& RBohr_h='Bohr Radius [cm]'
s=create_struct(s,'RBohr',RBohr) & SH=create_struct(SH,'RBohr',RBohr_h)

day=24.*3600.+3.*60.+56.555   	& day_h='mean solar day [sec]'
year=365.24219878D*day 	& year_h='Equinoctial Year [sec]'
yr=365.25636556D*day	& yr_h='Sidereal Year [sec]'
s=create_struct(s,'day',day,'year',year,'yr',yr)
SH=create_struct(SH,'day',day_h,'year',year_h,'yr',yr_h)

pc=3.26*c*yr    	& pc_h='1 parsec [cm]'
ly=c*yr         	& ly_h='1 Light Year [cm]'
AU=1.495985e13  	& AU_h='Astronomical Unit [cm]'
s=create_struct(s,'pc',pc,'ly',ly,'AU',AU)
SH=create_struct(SH,'pc',pc_h,'ly',ly_h,'AU',AU_h)

Msun=1.989d33   	& Msun_h='Mass of Sun [gm]'
Lsun=3.826d33       	& Lsun_h='Luminosity of Sun [erg/s]'
Rsun=6.969e10      	& Rsun_h='Radius of Sun [cm]'
Tsun=5770.0      	& Tsun_h='Effective Temperature of Sun [K]'
s=create_struct(s,'Msun',Msun,'Lsun',Lsun,'Rsun',Rsun,'Tsun',Tsun)
SH=create_struct(SH,'Msun',Msun_h,'Lsun',Lsun_h,'Rsun',Rsun_h,'Tsun',Tsun_h)

Mgeo=5.977d27   	& Mgeo_h='Mass of Earth [gm]'
Rgeo=6371.23    	& Rgeo_h='mean radius of Earth [km]'
s=create_struct(s,'Mgeo',Mgeo,'Rgeo',Rgeo)
SH=create_struct(SH,'Mgeo',Mgeo_h,'Rgeo',Rgeo_h)

;RgeoEQ=6378.17  	& RgeoEQ_h='Equatorial radius of Earth [km]'
;RgeoPO=6356.90  	& RgeoPO_h='Polar radius of Earth [km]'
;s=create_struct(s,'RgeoEQ',RgeoEQ,'RgeoPO',RgeoPO)
;SH=create_struct(SH,'RgeoEQ',RgeoEQ_h,'RgeoPO',RgeoPO_h)

;Mmoon=Mgeo/81.33	& Mmoon_h='Mass of Moon [gm]'
;Rmoon=1737.9    	& Rmoon_h='Radius of Moon [km]'
;Dmoon=384404.0  	& Dmoon_h='mean distance of Moon from Earth [km]'
;s=create_struct(s,'Mmoon',Mmoon,'Rmoon',Rmoon,'Dmoon',Dmoon)
;SH=create_struct(SH,'Mmoon',Mmoon_h,'Rmoon',Rmoon_h,'Dmoon',Dmoon_h)

;DBode=0.4+[1.,0.3*2.^(findgen(9))] & DBode_h="Bode's Law distances [AU]'
;Dplanet=[0.39,0.72,1,1.52,2.9,5.2,9.55,19.2,30.1,39.5]
;Dplanet_h='distance of planets from Sun [AU]'
;s=create_struct(s,'DBode',DBode,'Dplanet',Dplanet)
;SH=create_struct(SH,'DBode',DBode_h,'Dplanet',Dplanet_h)

funcon=create_struct('HELP',SH,s)

return
end
