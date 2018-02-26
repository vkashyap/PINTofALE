function popsol,temp,pres,chianti,n_e=n_e,delec=delec,z=z,jon=jon,$
	radtemp=radtemp,dilute=dilute,pe_ratio=pe_ratio,$
	n_levels=n_levels,data_str=data_str
;+
;function	popsol
;	returns the level populations at given temperature and
;	pressure for a set of ions of a given element as an array
;	of size [N(TEMP),N(PRES or DENSITY),N(N_LEVELS)]
;
;usage
;	pop=popsol(temp,pres,chianti,n_e=n_e,/delec,$
;	radtemp=radtemp,dilute=dilute,pe_ratio=pe_ratio,$
;	n_levels=n_levels,data_str=data_str)
;
;parameters
;	temp	[INPUT; required] Temperature [K]
;	pres	[INPUT; required] pressure [cm^-3 K]
;	chianti	[INPUT; required] structure that contains all the relevant
;		information gleaned from the CHIANTI database (type
;			DBHELP=RD_CHIANTI() & HELP,DBHELP,/STRUCTURE
;		for a description of what this structure contains)
;
;keywords
;	n_e	[INPUT] electron density in cm^-3
;		* if set, overrides PRES/TEMP
;	delec	[INPUT] if set, assumes these are DR lines
;	z	[INPUT] atomic number of element
;	jon	[INPUT] ionic state of species
;		* Z and JON are used only if DELEC is set
;	radtemp	[INPUT] The blackbody radiation field temperature
;		* default is 6000 K
;	dilute	[INPUT] a dilution factor for photoionization, based
;		on the distance above a star's surface -- see CHIANTI
;		routine R2W()
;		* if not set, assumed to be 0
;	pe_ratio [INPUT] proton-to-electron ratio
;		* if not given, assumed to be 1.0
;	n_levels [INPUT] number of levels in the model
;		* (from CHIANTI/idl/low_level/pop_solver.pro)
;		This allows the number of levels in the model to be reduced.
;		E.g., if the full model contains 100 levels, one could set
;		N_LEVELS=50. This can be useful if one is interested in
;		looking at the effects of cascading from higher levels.
;	data_str [OUTPUT] some potentially useful data
;		* (from CHIANTI/idl/low_level/pop_solver.pro)
;		If POP_SOLVER is called for just 1 temperature and density,
;		then the individual data arrays for each of the physical
;		processes can be output through DATA_STR.  This allows the
;		user to check for the dominant processes affecting the
;		population of a given level. DATA_STR is a structure with
;		the following tags:
;			.aa	A-values
;			.aax	Photoexcitation/stimulated emission
;			.cc	Electron rate coefficients
;			.ccp	Proton rate coefficients
;		Each tag is a 2D array, and is such that, e.g.,
;		aa[0,20] corresponds to an excitation, while aa[20,0]
;		is a de-excitation.
;
;A note about the default value of PE_RATIO:
;	The original version of this program,
;		CHIANTI/dbase/idl/low_level/pop_solver.pro,
;	set the default to be 0.83, which is the limit for a high-T
;	plasma with cosmic abundances.  Later versions include a
;	calculated value passed in via a common block.  Here we
;	exclude this factor because we want to store the emissivities
;	without reference to their ion balance and abundance, on which
;	this factor depends.  The intensity of a line transition from
;	upper level j to lower level i is in general
;		dI(lam_ij) = (1/4pi)(hc/lam_ij) N_j A_ji dV [ergs/s/sr]
;	where lam_ij is the wavelength of the transition, A_ji is the
;	Einstein coefficient (the spontaneous transition probability),
;	dV is the emitting volume, and N_j is the number density of the
;	ions in level j.  We can write N_j as the product
;		N_j = (N_j/N_XIon) (N_XIon/N_X) (N_X/N_H) (N_H/N_e) N_e
;	where N_XIon is the number density of the ions,
;	N_j/N_XIon is the population fraction of the ions at level j,
;	N_XIon/N_X is the ion fraction,
;	N_X/N_H is the relative abundance of element X wrt H, and
;	N_H/N_e is the H density relative to the free-electron density.
;
;	In PINTofALE, we calculate and store the emissivities as
;		e_ji = (hc/lam_ij) A_ji (1/N_e) (N_j/N_XIon)
;	and write the line intensity as
;		f(lam_ij) = e_ji (N_XIon/N_X) (N_X/N_H) EM (N_H/N_e)
;	where EM = N_e^2 V is the emission measure, and the factor
;	(1/4pi) is always ignored and the factor (hc/lam_ij) is usually
;	taken out again.  This calculation is carried out in LINEFLX(),
;	which provides for the inclusion of the ion balance, abundance,
;	and the proton-to-electron ratio.  Hence the use of a default of
;	1.0 for PE_RATIO here.
;
;subroutines
;	CHIANTI routine DESCALE_ALL
;	NENH
;
;history
;	POP_SOLVER.PRO written by Peter Young
;{..
;..	Ver 1, PRY 29-Mar-99
;..	Ver 2, PRY 30-Apr-99, added call to get_prot_rates
;..	Ver 3, PRY 15-Dec-99, added deu to upsilon common block in order 
;..		to be consistent with the main Chianti routines.
;..	Ver 4, PRY 9-May-00, corrected problem with threshold when dealing 
;..		with sparse matrices. Basically values less than 1.e-30 in 
;..		the c-matrix were being set to zero and giving rise to 
;..		NaN's in certain circumstances.
;..       Ver.5, PRY 14-Jul-00, changed elvl common block to the elvlc common 
;..               block which is now the Chianti standard. Also, when 
;..               descaling upsilons, the routine now uses the Delta-E from 
;..               the .splups file.
;..       Ver.6, PRY 9-Aug-00, changed routine to deal better with the 
;..               dielectronic recombination files
;..       Ver.7, PRY 17-Aug-00, routine does not call LINBCG now if radtemp 
;..               is non-zero.
;..       Ver.8, PRY 29-Aug-00, the sparse matrix section has been disabled.
;..       Ver.9, PRY 12-Nov-01, calls routine proton_dens() to calculate the 
;..               proton to electron ratio.
;..       Ver.10, PRY, 6-Dec-01, corrected bug when there are more levels 
;..               in .splups file than in .elvlc file (ZnXXV).
;..       Ver.11, PRY, 11-Jul-02, removed ION keyword
;..       Ver.12, PRY, 9-Aug-02, within the equation solving section, I've set 
;..               the population of the ground level (rather the n_level level) 
;..               to 1, and this seems to stop negative populations appearing 
;..               in extreme conditions.
;..       Ver.12, PRY, 21-Aug-02, changed exp(-1/1/a) to exp(-a) in electron
;..               excitation section which caused a hang-up in some 
;..               circumstances. Also, the routine now uses vector ECMC 
;..               (combined experimental and theoretical energies) in 
;..               determining if a level lies above or below another level. 
;..               Previously only used the observed energy vector. Also, the 
;..               exponential in the electron excitation section now uses the 
;..               (accurate) .elvlc energy separation rather than the .splups 
;..               energy separation, which can cause significant (~20-30%) 
;..               differences in level populations of high-lying levels at 
;..               low temperatures.
;..       Ver.13, PRY, 10-Sep-02, corrected bug for proton rates. The excitation 
;..               and de-excitation rates were being swapped.
;..
;..       V. 14  4-Oct-2003  Giulio Del Zanna (GDZ).
;..               -removed all COMMON blocks (note that only proton_dens.pro has
;..                one: COMMON elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref)
;..               -only the essential information input is passed to the routine
;..                via a new input structure.
;..               -fixed a bug, that affected all the satellite lines, and was 
;..                introduced in v.12,  included in  CHIANTI v.4.0.
;..                basically the ionization potential was not subtracted when
;..                calculating the Delta E in the exponential.
;..
;.. VERSION     : 14, 4-Oct-2003 
;..}
;	modified to eliminate common blocks and the call to SETUP_ION,
;	  and altered the I/O interface to mesh with PINTofALE, changed
;	  call to PROTON_DENS() to call NENH() and changed default PE_RATIO
;	  to 1.0 (Vinay Kashyap, Apr03)
;	modified to account for v14, for CHIANTI v4.2 (VK; Apr04)
;-

;	usage
ok='ok' & np=n_params(0)
nt=n_elements(temp) & npr=n_elements(pres) & nch=n_tags(chianti)
if np lt 3 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='Temperature: not defined' else $
  if npr eq 0 and not keyword_set(n_e) then ok='Pressure: not defined' else $
   if nch eq 0 then ok='CHIANTI structure is undefined'
if ok ne 'ok' then begin
  print,'Usage: pop=popsol(temperature,pressure,chianti_structure,$'
  print,'       n_e=edensity,radtemp=radiation_temperature,$'
  print,'       dilute=dilution_factor,pe_ratio=proton_to_electron_ratio,$'
  print,'       n_levels=n_levels,data_str=data_str)'
  print,'  returns level populations at given temperature and pressure'
  if np ne 0 then message,ok,/info
  return,-1L
endif


;------------------------------------------------------
;{VK	the following bit is how POP_SOLVER is set up
;------------------------------------------------------
;;these are required:
;;-------------------
;
; gname = input.gname
; jj= input.jj
; ecm= input.ecm
; ecmth= input.ecmth
; wvl= input.wvl
; a_value= input.a_value
; splstr= input.splstr
;
;;optional ones:
;;-------------------
;
;IF tag_exist(input, 'radtemp') THEN  radtemp= input.radtemp
;IF tag_exist(input, 'dilute') THEN  dilute= input.dilute
;IF tag_exist(input, 'prot_struc') THEN  prot_struc= input.prot_struc
;IF tag_exist(input, 'pe_ratio') THEN  pe_ratio= input.pe_ratio
;
;
;;COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;;COMMON wgfa, wvl,gf,a_value
;;COMMON upsilon, splstr
;;COMMON radiative, radtemp,dilute
;;COMMON proton, prot_struc, pe_ratio
;
;IF n_params(0) LT 4 THEN BEGIN
;   print,' use>  pop_solver,input, temperature,density,populations, $'
;   print,'                   [n_levels=n_levels, data_str=data_str] '
;   return
;ENDIF
;
;convertname,gname,iz,ion
;ion2spectroscopic,gname,snote, dielectronic=dielectronic
;
;if dielectronic then begin
;read_ip,concat_dir(concat_dir(!xuvtop, 'ip'), 'chianti.ip'),ionpot,ipref
;ip=ionpot(iz-1,ion-1)
;endif else ip=0
;
;
;xne = DOUBLE(xne)
;t = DOUBLE(t)
;;
;; need the following to turn t into an array if it only has 1 element
;;
;IF n_elements(t) EQ 1 THEN BEGIN
;  t0=t
;  t=dblarr(1)
;  t[0]=t0
;ENDIF
;------------------------------------------------------
;VK:	thus ends POP_SOLVER setup}{below is PoA translation
;------------------------------------------------------

;	inputs
t=[double(temp)]
if keyword_set(n_e) then xne=double(n_e) else xne=double(pres/t)

;	decode CHIANTI variables
l1=chianti.lev1 & l2=chianti.lev2
w1=chianti.wvl & g1=chianti.gf & a1=chianti.a
nlev=max([l1,l2]) & wvl=fltarr(nlev,nlev) & gf=wvl & a_value=wvl
;
o1=where(w1 eq 0,mo1) & o2=where(w1 ne 0,mo2)
wvl(l1-1,l2-1)=abs(w1) & gf(l1-1,l2-1)=g1
if mo1 gt 0 then begin
  a_value(l1(o1)-1,l2(o1)-1) = a1(o1)
  if mo2 gt 0 then a_value(l1(o2)-1,l2(o2)-1)=a_value(l1(o2)-1,l2(o2)-1)+a1(o2)
endif else begin
  if mo2 gt 0 then a_value(l1(o2)-1,l2(o2)-1)=a1(o2)
endelse
;
ecm=chianti.ecmobs & ecmth=chianti.ecmthr & jj=chianti.j & mult=2*jj+1
splstr=chianti.splstr
prot_struc=chianti.psplstr
;
deu=chianti.deryd
t_type=chianti.trans & c_ups=chianti.c_ups & splups=chianti.spups
ionpot=chianti.ip
if not keyword_set(delec) then ip=0 else begin
  iz=1 & if keyword_set(Z) then iz=fix(Z[0])
  ion=1 & if keyword_set(jon) then ion=fix(jon[0])
  ip=ionpot(iz-1,ion-1)
endelse

;	special keywords
;	RADTEMP
if not keyword_set(radtemp) then radtemp=6000.D
;	DILUTE
if not keyword_set(dilute) then dilute=0.
;	PE_RATIO
if not keyword_set(pe_ratio) then pe_ratio=dblarr(nt)+1.0
;------------------------------------------------------
;VK:	thus ends PoA setup}
;------------------------------------------------------

;------------------------------------------------------
;{VK: from here on down is uninterrupted POP_SOLVER.PRO
;	{except for the bit about calling PROTON_DENS()
;------------------------------------------------------

ecmc=ecm
ind=where(ecm EQ 0.)
IF ind[0] NE -1 THEN ecmc[ind]=ecmth[ind]

mult=2.*jj+1.
;
hck=1.98648d-16/1.38062d-16
ryd2cm=109737.31534d
;
n_elvl=n_elements(ecm)
wsize=size(wvl)
n_wgfa1=wsize(1)
n_wgfa2=wsize(2)
usize=max([splstr.lvl1,splstr.lvl2])
;
IF N_ELEMENTS(n_levels) EQ 0 THEN n_levels=min([n_elvl,n_wgfa1,n_wgfa2,usize])
;
;
c=dblarr(n_levels,n_levels)
d=dblarr(n_levels,n_levels)
b=dblarr(n_levels)

diag=dblarr(n_levels)

nt=N_ELEMENTS(t)       ; no. of temperatures
nxne=N_ELEMENTS(xne)   ; no. of densities


pop=dblarr(nt,nxne,n_levels)

;;------------------------------[]
; The arrays
;
; e.g., aa(0,19) will be zero (no 0 -> 19 A value)
;       aa(19,0) will be non-zero
;
;       qq(0,19) electron excitation
;       qq(19,0) electron de-excitation
;;------------------------------[]

ident=make_array(n_levels,val=1.)  ; use for making de arrays
ecmn=ecmc(0:n_levels-1)        ; n_levels
den=ecmn#ident & dem=ident#ecmn

aat=a_value(0:n_levels-1,0:n_levels-1)        ; transpose of aa

aa=DOUBLE(TRANSPOSE(aat))

aax=c

IF N_ELEMENTS(dilute) EQ 0 THEN dilute=0.

;;------------------------------------------------------------------[-]
; The following loads up the photoexcitation (pexc) and stimulated 
; emission (stem) arrays)
;
IF dilute NE 0. THEN BEGIN
  stem=c & pexc=c & ede=c
 ;
  multn=mult[0:n_levels-1]      ; in case mult and ecm are bigger than 
 ;
  mm=TRANSPOSE(multn#(1/multn))      ; needed for photoexcitation
 ;
  dd=ABS(den-dem)*hck/radtemp
  ede=exp(dd) - 1.
 ;
  ind=WHERE( (aat NE 0.) AND (ede NE 0.) )
  pexc[ind]=aat[ind]*dilute*mm[ind]/ede[ind]
 ;
  ind=where( (aa NE 0.) AND (ede NE 0.) )
  stem[ind]=aa[ind]*dilute/ede[ind]
 ;
  aax=pexc+stem
ENDIF
;;------------------------------------------------------------------[-]

;________________
; create a ppr array for the proton rates
;
ppr=MAKE_ARRAY(n_levels,n_levels,nt,/double)
IF n_tags(prot_struc) NE 0 THEN BEGIN
 ;
 ;------------------------------------------------------
 ;{VK: replace call to PROTON_DENS() by NENH()
 ;------------------------------------------------------
 ;
  IF (n_elements(pe_ratio) NE nt) THEN BEGIN
    print,'%POP_SOLVER: WARNING, pe_ratio size does not match temp'
    print,n_elements(pe_ratio),nt
 ;   pe_ratio=proton_dens(alog10(t))
    pe_ratio=1./(nenh(alog10(t),/nproton)>1.)
  ENDIF
 ;
 ;------------------------------------------------------
 ;VK: replaced call to PROTON_DENS() by NENH()}
 ;	and from here on out, back to uninterrupted POP_SOLVER}
 ;------------------------------------------------------
 ;
 ;
  FOR i=0,n_elements(prot_struc)-1 DO BEGIN
    l1=prot_struc[i].lvl1-1
    l2=prot_struc[i].lvl2-1
    de=ABS(prot_struc[i].de)
    descale_all,t,prot_struc,i,prate
    IF ecmc(l1) LT ecmc(l2) THEN BEGIN
      ppr[l1,l2,*]=prate*pe_ratio
      ppr[l2,l1,*]=prate*pe_ratio*mult[l1]/mult[l2]* $
           exp(de*13.61/8.617/10.^(-5)/t)
    ENDIF ELSE BEGIN
      ppr[l2,l1,*]=prate*pe_ratio
      ppr[l1,l2,*]=prate*pe_ratio*mult[l2]/mult[l1]* $
           exp(de*13.61/8.617/10.^(-5)/t)
    ENDELSE
  ENDFOR
ENDIF
;________________




;______________
; Create a qq array for electron rates
;
qq=MAKE_ARRAY(n_levels,n_levels,nt,/double)
;
l1=splstr.lvl1-1
l2=splstr.lvl2-1

;GDZ- added ip
kte=(hck*abs(ecmc[l1]-(ecmc[l2]-ip))) # (1d0/t)
;******************************************


; kte=(hck*abs(ecmc[l1]-ecmc[l2])) # (1d0/t)


ind_pos=where(ecmc[l2] GT ecmc[l1])
ind_neg=where(ecmc[l2] LT ecmc[l1])
xx=dblarr(n_elements(splstr),nt)
yy=xx
;
; xx and yy contain all factors in the expression for the rate coefficient, 
; except for the upsilon. They can be generated using array operations - the 
; upsilons need a for loop.
;
IF ind_neg[0] NE -1 THEN BEGIN
  xx[ind_neg,*]=(8.63d-6/(mult[l1[ind_neg]])#(1./sqrt(t)))
  yy[ind_neg,*]=8.63d-6* exp(-kte[ind_neg,*]) * $
       1./( mult[l2[ind_neg]] # sqrt(t) )
ENDIF
IF ind_pos[0] NE -1 THEN BEGIN
  yy[ind_pos,*]=(8.63e-6/mult[l2[ind_pos]]) # (1./sqrt(t))
  xx[ind_pos,*]=8.63e-6* exp(-kte[ind_pos,*]) * $
       1./(mult[l1[ind_pos]] # sqrt(t))
ENDIF
;
; this is the for loop for the upsilons
;
FOR i=0,n_elements(splstr)-1 DO BEGIN
  IF (l1[i] LE n_levels-1) AND (l2[i] LE n_levels-1) THEN BEGIN
    descale_all,t,splstr,i,ups
    qq[l1[i],l2[i],*]=xx[i,*]*ups
    qq[l2[i],l1[i],*]=yy[i,*]*ups
  ENDIF
ENDFOR
;______________



FOR j=0,nt-1 DO BEGIN             ;; loop over temperatures

;
;  We are solving the equation Ax=0 where A is the matrix c, and x is a 
;  vector containing the level populations. There are an infinity of 
;  solutions to Ax=0 of the form kx', 
;  where k is real. This is the case here, and so we throw away one of 
;  the linear equations, and replace it with an equation setting the 
;  ground level population to 1. Once the linear equations have been 
;  solved, the populations are renormalised such that the sum of the level 
;  populations is 1.
;  The equation we throw away is the first row in A.
;  Note with the normalisation, we are solving Ax=b, where all elements 
;  of b are zero, except the first which is one.
;
;  Most of the c arrays we are dealing with in CHIANTI are sparse, with 
;  many zero elements. IDL has special routines for dealing with sparse 
;  matrices, and we use this below. We define frac_not_zero to denote the 
;  number of non-zero elements in c. If this is less than a specified number, 
;  then we switch to using the sparse matrix routine SPRSIN.
;

FOR i=0,nxne-1 DO BEGIN            ;; loop over densities
 ;
  c=aa + xne(i)*qq[*,*,j] + xne(i)*ppr[*,*,j] + aax    ; creates the c matrix
 ;
  ind=where(c NE 0.)
  frac_not_zero=float(n_elements(ind))/float(n_elements(c)) 
 ;
  diag = -total(c,2)
  c[findgen(n_levels),findgen(n_levels)] = diag
 ;
  c(*,0)=0d0                    ; set this row to zeros...
  c(0,0)=1d0                    ; ...except for ground level
 ;
  b(0)=1d0                      ; b is zero except for last element
 ;
 ;
 ; Currently (v4 of CHIANTI) the sparse matrix solver (SPRSIN) is not used 
 ; in CHIANTI so the IF statement below is set to always use the INVERT 
 ; routine.
 ;
  IF frac_NOT_zero GT 0. THEN BEGIN
    pp1=invert(transpose(temporary(c)),status)#b
  ENDIF ELSE BEGIN
    pp1=LINBCG(SPRSIN(temporary(c),thresh=1d-40),b,b)
  ENDELSE
 ;
  pp=pp1/total(pp1)     ; force total population to be 1
  pop[j,i,*]=pp         ; temp first, then dens

ENDFOR

ENDFOR

IF nt EQ 1 AND nxne EQ 1 THEN $
     data_str={aa: aa, aax: aax, cc: xne[0]*qq[*,*,0], ccp: xne[0]*ppr[*,*,0]}

;-------------------------------------------------------------
;VK: above was uninterrupted POP_SOLVER.PRO}
;-------------------------------------------------------------

return,pop
end
