pro inicon,atom=atom,aname=aname,amu=amu,fip=fip,fundae=fundae,roman=roman,$
	help=help
;+
;procedure	inicon
;	return atomic symbols, names, ionization potentials, fundamental
;	and other fun constants, atomic weights, roman numerals
;
;syntax
;	inicon,atom=atom,aname=aname,amu=amu,fip=fip,fundae=fundae,$
;	roman=roman,/help
;
;parameters	NONE
;
;keywords
;	atom	[OUTPUT] atomic symbols
;	aname	[OUTPUT] names of elements
;	amu	[OUTPUT] structure containing atomic weights in AMU (1H=1)
;		* AMU.(Z) is an array, with first element the atomic weight
;		  for terrestrial composition, then the AMUs of the isotopes
;		  in order of abundance (i.e., most abundant is 2nd element,
;		  second most abundant is 3rd element, etc. -- but not
;		  completely implemented)
;	fip	[OUTPUT] first ionization potentials [eV]
;		* set to -1 if not known
;	fundae	[OUTPUT] structure containing some useful physical and
;		astronomical constants
;		* just do help,fundae,/str -- should be obvious
;	roman	[OUTPUT] roman numerals, upto max(ATOM)+1
;		* beware that CHIANTI has a function named ROMAN()
;		  which returns numerals that are offset by 1 compared
;		  to this variable.
;
;		NOTE: in IDL versions prior to 5, the keywords must be
;		pre-defined if output is desired
;
;requires
;	CREATE_STRUCT
;
;history
;	modified from INITSTUFF (VK; 99May)
;	bug correction -- ROMAN wasn't working (VK; 99Jul)
;	added Jupiter (VK; MIM.XI)
;	added keV-to-Ang (VK; MMII.I)
;	added planetary data (VK; MMII.IX)
;	changed AMU tag names from symbol to full names (VK; MMII.X)
;	corrected AMU value (JD; IMVIM.I)
;	added JouleV, ergeV, and arcsr (VK; IMVIM.V)
;	corrected Bode's Law and extended to 12 "planets" (VK; MMV.VII)
;	changed Planets field to "Classical Planets", moved Pluto to
;	  separate Plutoid field (VK; MMVIII.VII)
;	added Thompson cross section (per JJD req), Compton wvl
;	  (VK; MMXI.XI)
;	new IAU definition for A.U. (VK; MMXII.IX)
;	added new Pluto and Charon measurements from New Horizons (VK; MMXV.VII)
;	updated h to match new SI definition (Cho 2018, Science 362, 625; VK MMXVII.XI)
;-

;	usage
idlvers=float(!version.release)
if idlvers ge 5 then begin
  if arg_present(atom)   then get_atm=1 else get_atm=0
  if arg_present(aname)  then get_anm=1 else get_anm=0
  if arg_present(amu)    then get_amu=1 else get_amu=0
  if arg_present(fip)    then get_fip=1 else get_fip=0
  if arg_present(fundae) then get_fun=1 else get_fun=0
  if arg_present(roman)  then get_rom=1 else get_rom=0
endif else begin
  if keyword_set(atom)   then get_atm=1 else get_atm=0
  if keyword_set(aname)  then get_anm=1 else get_anm=0
  if keyword_set(amu)    then get_amu=1 else get_amu=0
  if keyword_set(fip)    then get_fip=1 else get_fip=0
  if keyword_set(fundae) then get_fun=1 else get_fun=0
  if keyword_set(roman)  then get_rom=1 else get_rom=0
endelse
if get_atm+get_anm+get_amu+get_fip+get_fun+get_rom eq 0 then help=1
if keyword_set(help) then begin
  print,'Usage: inicon,atom=atom,aname=aname,amu=amu,fip=fip,fundae=fundae,$'
  print,'       roman=roman,/help'
  print,'  strictly utility routine returns atomic symbols, names, weights,'
  print,'  ionization potentials, some useful constants, roman numerals'
endif

;	elements from 1-109
if get_atm gt 0 then begin				;(ATOM
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
  Z10=[	'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt']	;100-109
  atom=[ Z00, Z01, Z02, Z03, Z04, Z05, Z06, Z07, Z08, Z09, Z10]
  nZ=n_elements(atom)
endif							;ATOM)

;	elements from 1-109
if get_anm gt 0 then begin				;(ANAME
  Z00=[	'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',$
	'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine']
  Z01=[	'Neon', 'Sodium', 'Magnesium', 'Aluminum', 'Silicon',$
	'Phosphorus', 'Sulfur', 'Chlorine', 'Argon', 'Potassium']
  Z02=[	'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium',$
	'Manganese', 'Iron', 'Cobalt', 'Nickel', 'Copper']
  Z03=[	'Zinc', 'Gallium', 'Germanium', 'Arsenic', 'Selenium',$
	'Bromine', 'Krypton', 'Rubidium', 'Strontium', 'Yttrium']
  Z04=[	'Zirconium', 'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium',$
	'Rhodium', 'Palladium', 'Silver', 'Cadmium', 'Indium']
  Z05=[	'Tin', 'Antimony', 'Tellurium', 'Iodine', 'Xenon',$
	'Cesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium']
  Z06=[	'Neodymium', 'Promethium', 'Samarium', 'Europium', 'Gadolinium',$
	'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium']
  Z07=[	'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten',$
	'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold']
  Z08=[	'Mercury', 'Thallium', 'Lead', 'Bismuth', 'Polonium',$
	'Astatine', 'Radon', 'Francium', 'Radium', 'Actinium']
  Z09=[	'Thorium', 'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',$
	'Americium', 'Curium', 'Berkelium', 'Californium', 'Einsteinium']
  Z10=[	'Fermium', 'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium',$
	'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium']
  aname=[Z00, Z01, Z02, Z03, Z04, Z05, Z06, Z07, Z08, Z09, Z10]
endif							;ANAME)

;	roman numerals from 1-100
if get_rom gt 0 then begin				;(ROMAN
  rom0=['I','II','III','IV','V','VI','VII','VIII','IX','X']
  rom1=[rom0,'X'+rom0,'XX'+rom0,'XXX'+rom0,'XXXXI','XXXXII','XXXXIII',$
	'XXXXIV','XXXXV','XXXXVI','XXXXVII','XXXXVIII']
  roman=[rom1,'IL','L','L'+rom1,'IC','C','C'+rom0]
  nrom=n_elements(roman)
endif							;ROMAN)

;	atomic weights
if get_amu gt 0 then begin				;(AMU
  ;	from http://www.webelements.com/webelements/elements/text/isot/
  if n_elements(aname) eq 0 then begin
    if idlvers lt 5 then aname=1 & inicon,aname=aname
  endif
  k=1 & m=[1.0079, 1,2] & amu=create_struct(aname(k-1L),m)
  k=2 & m=[4.0026, 4,3] & amu=create_struct(amu,aname(k-1L),m)
  k=3 & m=[6.941, 7,6] & amu=create_struct(amu,aname(k-1L),m)
  k=4 & m=[9.01218, 9] & amu=create_struct(amu,aname(k-1L),m)
  k=5 & m=[10.81, 11,10] & amu=create_struct(amu,aname(k-1L),m)
  k=6 & m=[12.011, 12,13] & amu=create_struct(amu,aname(k-1L),m)
  k=7 & m=[14.0067, 14,15] & amu=create_struct(amu,aname(k-1L),m)
  k=8 & m=[15.9994, 16,18,17] & amu=create_struct(amu,aname(k-1L),m)
  k=9 & m=[18.998403, 19] & amu=create_struct(amu,aname(k-1L),m)
  k=10 & m=[20.1797, 20,22,21] & amu=create_struct(amu,aname(k-1L),m)
  k=11 & m=[22.98977, 23] & amu=create_struct(amu,aname(k-1L),m)
  k=12 & m=[24.305, 24,26,25] & amu=create_struct(amu,aname(k-1L),m)
  k=13 & m=[26.98154, 27] & amu=create_struct(amu,aname(k-1L),m)
  k=14 & m=[28.0855, 28,29,30] & amu=create_struct(amu,aname(k-1L),m)
  k=15 & m=[30.97376, 31] & amu=create_struct(amu,aname(k-1L),m)
  k=16 & m=[32.066, 32,34,33,36] & amu=create_struct(amu,aname(k-1L),m)
  k=17 & m=[35.453, 35,37] & amu=create_struct(amu,aname(k-1L),m)
  k=18 & m=[39.948, 40,36,38] & amu=create_struct(amu,aname(k-1L),m)
  k=19 & m=[39.0983, 39,41,40] & amu=create_struct(amu,aname(k-1L),m)
  k=20 & m=[40.078, 40,44,42,48,43,46] & amu=create_struct(amu,aname(k-1L),m)
  k=21 & m=[44.95591, 45] & amu=create_struct(amu,aname(k-1L),m)
  k=22 & m=[47.867, 48,46,47,49,50] & amu=create_struct(amu,aname(k-1L),m)
  k=23 & m=[50.9415, 51,50] & amu=create_struct(amu,aname(k-1L),m)
  k=24 & m=[51.9961, 52,53,50,54] & amu=create_struct(amu,aname(k-1L),m)
  k=25 & m=[54.93805, 55] & amu=create_struct(amu,aname(k-1L),m)
  k=26 & m=[55.845, 56,54,57,58] & amu=create_struct(amu,aname(k-1L),m)
  k=27 & m=[58.9332, 59] & amu=create_struct(amu,aname(k-1L),m)
  k=28 & m=[58.6934, 58,60,62,61,64] & amu=create_struct(amu,aname(k-1L),m)
  k=29 & m=[63.546, 63,65] & amu=create_struct(amu,aname(k-1L),m)
  k=30 & m=[65.39, 64,66,68,67,70] & amu=create_struct(amu,aname(k-1L),m)
  k=31 & m=[69.723, 69,71] & amu=create_struct(amu,aname(k-1L),m)
  k=32 & m=[72.61, 74,72,70,73,76] & amu=create_struct(amu,aname(k-1L),m)
  k=33 & m=[74.9216, 75] & amu=create_struct(amu,aname(k-1L),m)
  k=34 & m=[78.96, 80,78,76,82,77,74] & amu=create_struct(amu,aname(k-1L),m)
  k=35 & m=[79.904, 79,81] & amu=create_struct(amu,aname(k-1L),m)
  k=36 & m=[83.80, 84,86,82,83,80,78] & amu=create_struct(amu,aname(k-1L),m)
  k=37 & m=[85.4678, 85,87] & amu=create_struct(amu,aname(k-1L),m)
  k=38 & m=[87.62, 88,86,87,84] & amu=create_struct(amu,aname(k-1L),m)
  k=39 & m=[88.9058, 89] & amu=create_struct(amu,aname(k-1L),m)
  k=40 & m=[91.22, 90,94,92,91,96] & amu=create_struct(amu,aname(k-1L),m)
  k=41 & m=[92.9064, 93] & amu=create_struct(amu,aname(k-1L),m)
  k=42 & m=[95.94, 98,96,95,92,100,97,94] & amu=create_struct(amu,aname(k-1L),m)
  k=43 & m=[98., 98] & amu=create_struct(amu,aname(k-1L),m)
  k=44 & m=[101.07, 102,104,101,99,100,96,98] & amu=create_struct(amu,aname(k-1L),m)
  k=45 & m=[102.9055, 103] & amu=create_struct(amu,aname(k-1L),m)
  k=46 & m=[106.42, 106,108,105,110,104,102] & amu=create_struct(amu,aname(k-1L),m)
  k=47 & m=[107.868, 107,109] & amu=create_struct(amu,aname(k-1L),m)
  k=48 & m=[112.41, 114,112,111,110,113,116,106,108] & amu=create_struct(amu,aname(k-1L),m)
  k=49 & m=[114.818, 115,113] & amu=create_struct(amu,aname(k-1L),m)
  k=50 & m=[118.71, 120,118,116,119,117,124,122,112,114,115] & amu=create_struct(amu,aname(k-1L),m)
  k=51 & m=[121.76, 121,123] & amu=create_struct(amu,aname(k-1L),m)
  k=52 & m=[127.60, 130,128,126,125,124,122,123,120] & amu=create_struct(amu,aname(k-1L),m)
  k=53 & m=[126.9045, 127] & amu=create_struct(amu,aname(k-1L),m)
  k=54 & m=[131.29, 132,129,131,134,146,130,128,124,126] & amu=create_struct(amu,aname(k-1L),m)
  k=55 & m=[132.90545, 133] & amu=create_struct(amu,aname(k-1L),m)
  k=56 & m=[137.33, 138,137,136,135,134,130,132] & amu=create_struct(amu,aname(k-1L),m)
  k=57 & m=[138.9055, 139,138] & amu=create_struct(amu,aname(k-1L),m)
  k=58 & m=[140.116, 140,142,138,136] & amu=create_struct(amu,aname(k-1L),m)
  k=59 & m=[140.90765, 141] & amu=create_struct(amu,aname(k-1L),m)
  k=60 & m=[144.24, 142,144,146,143,145,148,150] & amu=create_struct(amu,aname(k-1L),m)
  k=61 & m=[145., 145] & amu=create_struct(amu,aname(k-1L),m)
  k=62 & m=[150.36, 152,154,147,149,148,150,144] & amu=create_struct(amu,aname(k-1L),m)
  k=63 & m=[151.964, 153,151] & amu=create_struct(amu,aname(k-1L),m)
  k=64 & m=[157.25, 158,160,156,157,155,154,152] & amu=create_struct(amu,aname(k-1L),m)
  k=65 & m=[158.92534, 159] & amu=create_struct(amu,aname(k-1L),m)
  k=66 & m=[162.50, 164,162,163,161,160,158,156] & amu=create_struct(amu,aname(k-1L),m)
  k=67 & m=[164.9303, 165] & amu=create_struct(amu,aname(k-1L),m)
  k=68 & m=[167.26, 166,168,167,170,164,162] & amu=create_struct(amu,aname(k-1L),m)
  k=69 & m=[168.9342, 169] & amu=create_struct(amu,aname(k-1L),m)
  k=70 & m=[173.04, 174,172,173,171,176,170,168] & amu=create_struct(amu,aname(k-1L),m)
  k=71 & m=[174.967, 175,176] & amu=create_struct(amu,aname(k-1L),m)
  k=72 & m=[178.49, 180,178,177,179,176,174] & amu=create_struct(amu,aname(k-1L),m)
  k=73 & m=[180.9479, 181,180] & amu=create_struct(amu,aname(k-1L),m)
  k=74 & m=[183.84, 184,186,182,183,180] & amu=create_struct(amu,aname(k-1L),m)
  k=75 & m=[186.207, 187,185] & amu=create_struct(amu,aname(k-1L),m)
  k=76 & m=[190.23, 192,190,189,188,187,186,184] & amu=create_struct(amu,aname(k-1L),m)
  k=77 & m=[192.217, 193,191] & amu=create_struct(amu,aname(k-1L),m)
  k=78 & m=[195.078, 195,194,196,198,192,190] & amu=create_struct(amu,aname(k-1L),m)
  k=79 & m=[196.9665, 197] & amu=create_struct(amu,aname(k-1L),m)
  k=80 & m=[200.59, 202,200,199,201,198,204,196] & amu=create_struct(amu,aname(k-1L),m)
  k=81 & m=[204.3833, 205,203] & amu=create_struct(amu,aname(k-1L),m)
  k=82 & m=[207.2, 208,206,207,204] & amu=create_struct(amu,aname(k-1L),m)
  k=83 & m=[208.98038, 209] & amu=create_struct(amu,aname(k-1L),m)
  k=84 & m=[210., 210] & amu=create_struct(amu,aname(k-1L),m)
  k=85 & m=[210., 210] & amu=create_struct(amu,aname(k-1L),m)
  k=86 & m=[222., 222] & amu=create_struct(amu,aname(k-1L),m)
  k=87 & m=[223., 223] & amu=create_struct(amu,aname(k-1L),m)
  k=88 & m=[226., 226] & amu=create_struct(amu,aname(k-1L),m)
  k=89 & m=[227., 227] & amu=create_struct(amu,aname(k-1L),m)
  k=90 & m=[232.038, 232] & amu=create_struct(amu,aname(k-1L),m)
  k=91 & m=[231.03588, 231] & amu=create_struct(amu,aname(k-1L),m)
  k=92 & m=[238.0289, 238,235,234] & amu=create_struct(amu,aname(k-1L),m)
  k=93 & m=[237., 237] & amu=create_struct(amu,aname(k-1L),m)
  k=94 & m=[244., 244] & amu=create_struct(amu,aname(k-1L),m)
  k=95 & m=[243., 243] & amu=create_struct(amu,aname(k-1L),m)
  k=96 & m=[247., 247] & amu=create_struct(amu,aname(k-1L),m)
  k=97 & m=[247., 247] & amu=create_struct(amu,aname(k-1L),m)
  k=98 & m=[251., 251] & amu=create_struct(amu,aname(k-1L),m)
  k=99 & m=[252., 252] & amu=create_struct(amu,aname(k-1L),m)
  k=100 & m=[257., 257] & amu=create_struct(amu,aname(k-1L),m)
  k=101 & m=[258., 258] & amu=create_struct(amu,aname(k-1L),m)
  k=102 & m=[259., 259] & amu=create_struct(amu,aname(k-1L),m)
  k=103 & m=[262., 262] & amu=create_struct(amu,aname(k-1L),m)
  k=104 & m=[261., 261] & amu=create_struct(amu,aname(k-1L),m)
  k=105 & m=[262., 262] & amu=create_struct(amu,aname(k-1L),m)
  k=106 & m=[266., 266] & amu=create_struct(amu,aname(k-1L),m)
  k=107 & m=[264., 264] & amu=create_struct(amu,aname(k-1L),m)
  k=108 & m=[269., 269] & amu=create_struct(amu,aname(k-1L),m)
  k=109 & m=[268., 268] & amu=create_struct(amu,aname(k-1L),m)
  for i=k+1,n_elements(aname)-1 do amu=create_struct(amu,$
	aname(i-1),[2.7*i,fix(2.4*i+0.5)])
endif							;AMU)

;	first ionization potentials [eV]
if get_fip gt 0 then begin				;(FIP
  f00=[	13.5984, 24.5874, 5.3917, 9.3227, 8.2980, 11.2603,$
	14.5341, 13.6181, 17.4228]				;001-009
  f01=[	21.5646, 5.1391, 7.6462, 5.9858, 8.1517, 10.4867,$
	10.3600, 12.9676, 15.7596, 4.3407]			;010-019
  f02=[	6.1132, 6.5615, 6.8281, 6.7462, 6.7665, 7.4340,$
	7.9024, 7.8810, 7.6398, 7.7264]				;020-029
  f03=[	9.3942, 5.9993, 7.8994, 9.7886, 9.7524, 11.8138,$
	13.9996, 4.1771, 5.6949, 6.2171]			;030-039
  f04=[	6.6339, 6.7589, 7.0924, 7.28, 7.3605, 7.4589,$
	8.3369, 7.5762, 8.9938, 5.7864]				;040-049
  f05=[	7.3439, 8.6084, 9.0096, 10.4513, 12.1298, 3.8939,$
	5.2117, 5.5769, 5.5387, 5.473]				;050-059
  f06=[	5.5250, 5.582, 5.6436, 5.6704, 6.1501, 5.8638,$
	5.9389, 6.0215, 6.1077, 6.1843]				;060-069
  f07=[	6.2542, 5.4259, 6.8251, 7.5496, 7.8640, 7.8335,$
	8.4382, 8.9670, 8.9587, 9.2255]				;070-079
  f08=[	10.4375, 6.1082, 7.4167, 7.2856, 8.417, 9.5,$			;FIP for At from Sobel'man
	10.7485, 4.0727, 5.2784, 5.17]				;080-089
  f09=[	6.3067, 5.89, 6.1941, 6.2657, 6.0262, 5.9738,$
	5.9915, 6.1979, 6.2817, 6.42]				;090-099
  f10=[	6.50, 6.58, 6.65, 4.9, 6.0, -1,$
	-1, -1, -1, -1]						;100-109
  fip=[f00, f01, f02, f03, f04, f05, f06, f07, f08, f09, f10]
endif							;FIP)

;	physical and astronomical constants
if get_fun gt 0 then begin				;(FUNDAE

  c=2.99792458d10  	& c_h='Speed of light (defines SI meter) [cm/s]'
  s=create_struct('c',c) & SH=create_struct('c',c_h)

  h=6.626176d-27  	& h_h="Planck's constant [erg s]"
  h=6.626207015d-27	& h_h="Planck's constant (defines SI kg) [erg s]"
  s=create_struct(s,'h',h) & SH=create_struct(SH,'h',h_h)

  G=6.672e-8      	& G_h='Gravitational constant [cm^3/gm/s^2]'
  s=create_struct(s,'G',G) & SH=create_struct(SH,'G',G_h)

  Me=9.109558d-28  	& Me_h='Electron mass [gm]'
  Mp=1.67248d-24  	& Mp_h='Proton mass [gm]'
  Mn=1.674929d-24	& Mn_h='Neutron mass [gm]'
  h_1=1.66054d-24	& h_1_h='Atomic Mass Unit [gm]'
  s=create_struct(s,'Me',Me,'Mp',Mp,'Mn',Mn,'AMU',h_1)
  SH=create_struct(SH,'Me',Me_h,'Mp',Mp_h,'Mn',Mn_h,'AMU',h_1_h)

  e=1.6021892d-19  	& e_h='electron charge [C]'
  esu=4.803d-10      	& esu_h='electron charge [ESU]'
  eps0=8.854187d-12	& eps0_h='permittivity of vacuum [F/m]'
  s=create_struct(s,'e',e,'esu',esu,'eps0',eps0)
  SH=create_struct(SH,'e',e_h,'esu',esu_h,'eps0',eps0_h)

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

  ;kevAng=12.39852066		& kevAng_h='keV*Ang (1e8*h*c/(e*1e10))'
  kevAng=1d8*h*c/(e*1d10)	& kevAng_h='keV*Ang (1e8*h*c/(e*1e10))'
  eVwav=12379.7d-8		& eVwav_h='1 eV in wave numbers [/cm]'
  degeV=8.6173468d-05		& degeV_h='1 deg K in eV [eV] (degK*degev/1e3->keV, or degev*1e3/keV->degK)'
  JouleV=e			& JouleV_h='1 eV in Joule = Coulomb*meter [J]'
  ergeV=JouleV*1e7		& ergeV_h='1 eV in ergs [erg]'
  s=create_struct(s,'kevAng',kevAng,'eVwav',eVwav,'degeV',degeV,'JouleV',JouleV,'ergeV',ergeV)
  SH=create_struct(SH,'kevAng',kevAng_h,'eVwav',eVwav_h,'degeV',degeV_h,'JouleV',JouleV_h,'ergeV',ergeV_h)

  sigthom=6.6524586e-25	& sigthom_h='Thompson cross-section for electron [cm^2]'
  compton=h*1e8/(Me*c)	& compton_h='Compton wavelength for electron, h/mc [Angstrom]'
  	;note: CODATA value is 0.024263102175±33 Ang
  s=create_struct(s,'Thompson',sigthom,'Compton',compton) & SH=create_struct(SH,'Thompson',sigthom_h,'Compton',compton_h)

  radian=3600.*180./!DPI	& radian_h='1 radian [arcsec]'
  arcsr=(1./radian)^2		& arcsr_h='1 arcsec^2 in steradians [sr]'
  s=create_struct(s,'radian',radian,'arcsr',arcsr)
  SH=create_struct(SH,'radian',radian_h,'arcsr',arcsr_h)

  pi='3.14159265358979323846'	& pi_h='PI to 20 decimals'
  ee=2.718281828459D		& ee_h='e'
  gg=0.5772156649D		& gg_h="Euler's constant, gamma"
  s=create_struct(s,'pi',pi,'euler',ee,'gamma',gg)
  SH=create_struct(SH,'pi',pi_h,'euler',ee_h,'gamma',gg_h)

  Ryd=109677.58*eVWav   & Ryd_h='Rydberg Constant for H [eV]'
  s=create_struct(s,'Ryd',Ryd) & SH=create_struct(SH,'Ryd',Ryd_h)

  RBohr=2.*alog10(h)-alog10(4.*!DPI^2)-alog10(Me)-2.*alog10(esu)
  RBohr=10.^(RBohr)	& RBohr_h='Bohr Radius [cm]'
  s=create_struct(s,'RBohr',RBohr) & SH=create_struct(SH,'RBohr',RBohr_h)

  day=24.*3600.+3.*60.+56.555  	& day_h='mean solar day [sec]'
  year=365.24219878D*day 	& year_h='Equinoctial Year [sec]'
  yr=365.25636556D*day		& yr_h='Sidereal Year [sec]'
  s=create_struct(s,'day',day,'year',year,'yr',yr)
  SH=create_struct(SH,'day',day_h,'year',year_h,'yr',yr_h)

  pc=3.26*c*yr    	& pc_h='1 parsec [cm]'
  ly=c*yr         	& ly_h='1 Light Year [cm]'
  ;AU=1.495985e13  	&
  AU=14959787070000.D	;new IAU definition
  AU_h='Astronomical Unit [cm]'
  s=create_struct(s,'pc',pc,'ly',ly,'AU',AU)
  SH=create_struct(SH,'pc',pc_h,'ly',ly_h,'AU',AU_h)

  Msun=1.989d33   	& Msun_h='Mass of Sun [gm]'
  Lsun=3.826d33       	& Lsun_h='Luminosity of Sun [erg/s]'
  Rsun=6.969e10      	& Rsun_h='Radius of Sun [cm]'
  Tsun=5770.0      	& Tsun_h='Effective Temperature of Sun [K]'
  Psun=25.38		& Psun_h='Rotational period of Sun [day]'
  s=create_struct(s,'Msun',Msun,'Lsun',Lsun,'Rsun',Rsun,'Tsun',Tsun,'Psun',Psun)
  SH=create_struct(SH,'Msun',Msun_h,'Lsun',Lsun_h,'Rsun',Rsun_h,'Tsun',Tsun_h,'Psun',Psun_h)

  Mjup=1.8986d30	& Mjup_h='Mass of Jupiter [gm]'
  Rjup=6.9911e9 	& Rjup_h='volumetric radius of Jupiter at 1 bar [cm]'
  s=create_struct(s,'Mjup',Mjup,'Rjup',Rjup)
  SH=create_struct(SH,'Mjup',Mjup_h,'Rjup',Rjup_h)

  Mgeo=5.977d27   	& Mgeo_h='Mass of Earth [gm]'
  Mgeo=5.972d27	;http://www.cnn.com/2000/TECH/space/05/01/earth.weight.ap/
  Rgeo=6371.23    	& Rgeo_h='mean radius of Earth [km]'
  s=create_struct(s,'Mgeo',Mgeo,'Rgeo',Rgeo)
  SH=create_struct(SH,'Mgeo',Mgeo_h,'Rgeo',Rgeo_h)

  RgeoEQ=6378.17  	& RgeoEQ_h='Equatorial radius of Earth [km]'
  RgeoPO=6356.90  	& RgeoPO_h='Polar radius of Earth [km]'
  s=create_struct(s,'RgeoEQ',RgeoEQ,'RgeoPO',RgeoPO)
  SH=create_struct(SH,'RgeoEQ',RgeoEQ_h,'RgeoPO',RgeoPO_h)

  Mmoon=Mgeo/81.33	& Mmoon_h='Mass of Moon [gm]'
  Rmoon=1737.9    	& Rmoon_h='Radius of Moon [km]'
  Dmoon=384404.0  	& Dmoon_h='mean distance of Moon from Earth [km]'
  s=create_struct(s,'Mmoon',Mmoon,'Rmoon',Rmoon,'Dmoon',Dmoon)
  SH=create_struct(SH,'Mmoon',Mmoon_h,'Rmoon',Rmoon_h,'Dmoon',Dmoon_h)

  DBode=0.4+[0.,0.3*2.^(findgen(12))] & DBode_h="Bode's Law distances [AU]'
  ;Dplanet=[0.39,0.72,1,1.52,2.9,5.2,9.55,19.2,30.1,39.5]
  ;Dplanet_h='distance of planets from Sun [AU]'
  s=create_struct(s,'DBode',DBode)
  SH=create_struct(SH,'DBode',DBode_h)

  ;	http://nssdc.gsfc.nasa.gov/planetary/planetfact.html
  Planets=['Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune']
  Planets_h='The Classical Planets'
  MPlanet=1d27*[0.3302,4.8685,5.9736,0.64185,1898.6,568.46,86.832,102.43]
  MPlanet_h='Masses of the planets [gm]'
  RPlanet=[2439.7,6051.8,6371.0,3390.,69911.,58232.,25362.,24624.]*1d5
  RPlanet_h='Volumetric mean radius [cm]'
  DPlanet=[57.91,108.21,149.60,227.92,778.57,1433.53,2872.46,4495.06]*1d11
  DPlanet_h='semi-major axis of orbit [cm]'
  PPlanet=[87.969,224.701,365.256,686.980,4332.589,10759.22,30685.4,60189.]
  PPlanet_h='sidereal orbit period [days]'
  RotPlanet=[1407.6,-5832.5,23.9345,24.6229,9.9250,10.656,-17.24,16.11]
  RotPlanet_h='sidereal rotation period [hr]'
  DayPlanet=[4222.6,2802.0,24.,24.6597,9.9259,10.656,17.24,16.11]
  DayPlanet_h='length of day [hrs]'
  OPlanet=[0.01,177.36,23.45,25.19,3.13,26.73,97.77,28.32]
  OPlanet_h='obliquity to orbit [deg]'
  ePlanet=[0.2056, 0.0067, 0.0167, 0.0935, 0.0489, 0.0565, 0.0457, 0.0113]
  ePlanet_h='orbit eccentricity'
  s=create_struct(s,'Planets',Planets,'MPlanet',MPlanet,'RPlanet',RPlanet,$
	'DPlanet',DPlanet,'PPlanet',PPlanet,'RotPlanet',RotPlanet,'DayPlanet',$
	DayPlanet,'OPlanet',OPlanet,'ePlanet',ePlanet)
  SH=create_struct(SH,'Planets',Planets_h,'MPlanet',MPlanet_h,'RPlanet',RPlanet_h,$
	'DPlanet',DPlanet_h,'PPlanet',PPlanet_h,'RotPlanet',RotPlanet_h,$
	'DayPlanet',DayPlanet_h,'OPlanet',OPlanet_h,'ePlanet_h',ePlanet_h)

  Plutoids=['Pluto','Charon','Eris','Makemake']
  Plutoids_h='The Plutoids'
  MPlutoid=1d25*[1.305, 0.152, 1.67, 0.4]	;Makemake goes from 0.24-0.54
  MPlutoid_h='Masses of the Plutoids [gm]'
  ;RPlutoid=1e6*[1.195, 1.3, 0.75]
  RPlutoid=1e8*[1.185, 0.6, 1.16, 0.715]	;Pluto & Charon from New Horizons
  RPlutoid_h='estimated radius [cm]'
  DPlutoid=1e12*[5.906376272d, 5.906376272d, 14.60d, 685.03d]
  DPlutoid_h='semi-major axis of orbit [cm]'
  PPlutoid=[90613.3055, 203600., 113183.]
  PPlutoid_h='orbital period [days]'
  RotPlutoid=[-153.2928, -153.2928, 8., !values.F_NAN]	;Pluto & Charon tidally locked
  RotPlutoid_h='sidereal rotation period [hr]'
  DayPlutoid=[153.2820, 153.2820, 8., !values.F_NAN]	;Pluto & Charon tidally locked
  DayPlutoid_h='length of day [hrs]'
  OPlutoid=[122.53, 122.53, !values.F_NAN, !values.F_NAN]
  OPlutoid_h='obliquity to orbit [deg]'
  ePlutoid=[0.24880766, 0.24880766, 0.44177, 0.159]
  ePlutoid_h='orbit eccentricity'
  s=create_struct(s,'Plutoids',Plutoids,'MPlutoid',MPlutoid,'RPlutoid',RPlutoid,$
	'DPlutoid',DPlutoid,'PPlutoid',PPlutoid,'RotPlutoid',RotPlutoid,'DayPlutoid',$
	DayPlutoid,'OPlutoid',OPlutoid,'ePlutoid',ePlutoid)
  SH=create_struct(SH,'Plutoids',Plutoids_h,'MPlutoid',MPlutoid_h,'RPlutoid',RPlutoid_h,$
	'DPlutoid',DPlutoid_h,'PPlutoid',PPlutoid_h,'RotPlutoid',RotPlutoid_h,$
	'DayPlutoid',DayPlutoid_h,'OPlutoid',OPlutoid_h,'ePlutoid_h',ePlutoid_h)

  fundae=create_struct('HELP',SH,s)
endif							;FUNDAE)

return
end
