function popsol3,temp,pres,chianti,n_e=n_e,delec=delec
;+
;function	popsol3
;	returns the level populations at given temperature and
;	pressure for a set of ions of a given element
;
;usage
;	pop=popsol3(temp,pres,chianti,n_e=n_e)
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
;	delec	[INPUT] a flag to denote dielectronic recombination
;		transitions
;
;requires
;	DESCALE_UPS_VK, modified version of CHIANTI routine DESCALE_UPS
;
;history
;..	POP_SOLVER.PRO is a faster version of Ken Dere's original POPULATE.PRO
;..	PROGRAMMING NOTES:
;..	The following methods have been used to speed up the routine.
;..	(i) the radiative data is now manipulated with array operations 
;..	rather than for-loops
;..	(ii) the for-loops for the collisional data are now performed over 
;..	half the indices used previously
;..	(iii) for sparse c-arrays (defined to be those with >50% zeros), 
;..	the LINBCG idl routine is used to solve the linear equations. This 
;..	is quicker than the INVERT routine used previously.
;..	(iv) Vectors of temperatures and densities can be input instead 
;..	of just single values. The electron excitation array only has to 
;..	be loaded for different temperatures, and so considerable time 
;..	savings can be made if a number of densities are input at once.
;..	Other changes include:
;..	* the routine does not solve for level populations, but for 
;..	max(A_ji)/q_1j/N_e * N_j which is ~1 for low densities.
;..	Ver 1, PRY (Peter Young, CfA, pyoung@cfa) 29-Mar-99
;..	Ver 2, PRY 30-Apr-99, added call to get_prot_rates
;..	Ver 3, PRY 15-Dec-99, added deu to upsilon common block in order 
;..		to be consistent with the main Chianti routines.
;..	Ver 4, PRY 9-May-00, corrected problem with threshold when dealing 
;..		with sparse matrices. Basically values less than 1.e-30 in 
;..		the c-matrix were being set to zero and giving rise to 
;..		NaN's in certain circumstances.
;..	Ver.5, PRY 14-Jul-00, changed elvl common block to the elvlc common 
;..		block which is now the Chianti standard. Also, when 
;..		descaling upsilons, the routine now uses the Delta-E from 
;..		the .splups file.
;
;	modified to eliminate common blocks and calls to GET_PROT_RATES and
;	  SETUP_ION, and mesh with SCAR (Vinay Kashyap, SepMM)
;	added keyword DELEC (VK; OctMM)
;	changed name to POPSOL3 because now it doesn't work in v4 (VK; Jun02)
;-

message,'OBSOLETE!',/informational

;	usage
np=n_params(0)
if np lt 3 then begin
  print,'Usage: pop=popsol3(temperature,pressure,chianti_structure,$'
  print,'       n_e=edensity,delec=dielectric_recombinations)'
  print,'  returns level populations at given temperature and pressure'
  return,-1L
endif

;	inputs
t=temp & if keyword_set(n_e) then xne=n_e else xne=pres/t

;	decode CHIANTI variables
l1=chianti.lev1 & l2=chianti.lev2
w1=chianti.wvl & g1=chianti.gf & a1=chianti.a
nlev=max([l1,l2]) & wvl=fltarr(nlev,nlev) & gf=wvl & a_value=wvl
wvl(l1-1,l2-1)=w1 & gf(l1-1,l2-1)=g1 & nl1=n_elements(l1)
for j=0,nl1-1 do a_value(l1(j)-1,l2(j)-1)=a_value(l1(j)-1,l2(j)-1)+a1(j)
;
ecm=chianti.ecmobs & jj=chianti.j & mult=2*jj+1
;
deu=chianti.deryd
t_type=chianti.trans & c_ups=chianti.c_ups & splups=chianti.spups

;COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;COMMON wgfa, wvl,gf,a_value
;COMMON upsilon,t_type,deu,c_ups,splups
;COMMON radiative, radtemp,dilute
;COMMON proton, prot_struc

;-------------------------------------------------------------
;{VK: from here on down is uninterrupted POP_SOLVER.PRO, with 3 exceptions
;-------------------------------------------------------------

mult=2.*jj+1.
;
hck=1.98648e-16/1.38062e-16
ryd2cm=109737.31534d
;
n_elvl=n_elements(ecm)
wsize=size(wvl)
n_wgfa1=wsize(1)
n_wgfa2=wsize(2)
usize=size(splups)
n_ups1=usize(1)
n_ups2=usize(2)
;
IF N_ELEMENTS(n_levels) EQ 0 THEN n_levels=min([n_elvl,n_wgfa1,n_wgfa2,n_ups1,n_ups2])
;
;
c=dblarr(n_levels,n_levels)
d=dblarr(n_levels,n_levels)
b=dblarr(n_levels)

diag=dblarr(n_levels)

nt=N_ELEMENTS(t)       ; no. of temperatures
nxne=N_ELEMENTS(xne)   ; no. of densities

pp=make_array(n_levels,value=1.,/double)   ; inital guess for pp. Used by 
                                           ; linbcg.pro
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
ecmn=ecm(0:n_levels-1)        ; n_levels
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


qq=c

FOR j=0,nt-1 DO BEGIN             ;; loop over temperatures

;______________
; Create a qq array
;
FOR m=0,n_levels-1 DO BEGIN
  FOR n=m+1,n_levels-1 DO BEGIN
    decm=ecm(n)-ecm(m)
    deryd = ABS(deu[n,m])
   ;
    IF decm LT 0. THEN BEGIN
      IF t_type(n,m) GT 0 THEN BEGIN
         temp=t(j)
         xt=temp/(hck*ryd2cm*deryd)
	 ;VK #1: replace DESCALE_UPS by DESCALE_UPS_VK
         ;descale_ups,n,m,xt,ups
         descale_ups_vk,n,m,xt,ups,t_type,c_ups,splups
         qq(m,n)=8.63e-6*ups/(mult(m)*sqrt(temp))                 ; de-exc
         qq(n,m)=8.63e-6*ups*exp(-1./xt)/(mult(n)*sqrt(temp))     ; exc
      ENDIF
    ENDIF 
    IF decm GT 0. THEN BEGIN
      IF t_type(n,m) GT 0 THEN BEGIN
         temp=t(j)
         xt=temp/(hck*ryd2cm*deryd)
	 ;VK #2: replace DESCALE_UPS by DESCALE_UPS_VK
         ;descale_ups,n,m,xt,ups
         descale_ups_vk,n,m,xt,ups,t_type,c_ups,splups
         qq(n,m)=8.63e-6*ups/(mult(n)*sqrt(temp))                 ; de-exc
         qq(m,n)=8.63e-6*ups*exp(-1./xt)/(mult(m)*sqrt(temp))     ; exc
      ENDIF
    ENDIF
  ENDFOR
ENDFOR

;________________
; create a ppr array for the proton rates
;
ppr=MAKE_ARRAY(n_levels,n_levels,/double)
;VK #3: COMMENT OUT THIS BLOCK
;	IF KEYWORD_SET(proton) THEN BEGIN
;	  IF (prot_struc(0).lvl1 NE -1) THEN BEGIN
;	    ppr=get_prot_rates(t(j),prot_struc,n_levels)
;	  ENDIF
;	ENDIF
;________________


mqqj=DBLARR(n_levels)
maaj=mqqj
fact=mqqj


FOR i=0,nxne-1 DO BEGIN            ;; loop over densities

  c=aa + xne(i)*qq + xne(i)*ppr + aax    ; creates the c matrix

 ;
 ; Solving the level balance equations. For v.low densities, there is a huge 
 ; difference between the ground level population (1) and excited level 
 ; populations (<1.e-15). I am actually going to solve for 
 ; scaled level populations.
 ; The form of the scaling varies depending on the conditions under 
 ; consideration. The `max' commands below pick out the dominant 
 ; contributor to the level population, and puts them in mqqj and maaj. 
 ; These are then used to form the scaling vector `fact'. The reasons for 
 ; the names mqqj and maaj stem from the fact that if radiative decay and 
 ; electron excitation are the dominant processes, then the scaling 
 ; factor is of the form  A_j0/q_0j/N_e.
 ;
  FOR k=0,n_levels-1 DO BEGIN
    IF qq(0,k) EQ 0. THEN qq_gr=1.e-14 ELSE qq_gr=qq(0,k)
    mqqj(k)=MAX([qq_gr,aax(0,k)/xne(i)])
    IF qq(k,0) EQ 0. THEN qq_gr=1.e-14 ELSE qq_gr=qq(k,0)
    maaj(k)=MAX([REFORM(aa(k,*)),qq_gr*xne(i)])
  ENDFOR

  fact=maaj/mqqj/xne(i)      ; create fact
  ind=where(fact eq 0.)      ; This should only be for ground level
  IF ind(0) NE -1 THEN fact(ind)=1.

  ; Create the diagonal elements of c
  ;
  FOR k=0,n_levels-1 DO BEGIN
    diag(k) = - TOTAL( c(k,*) )
    c(k,k)=diag(k)
  ENDFOR

  FOR k=0,n_levels-1 DO c(*,k)=c(*,k)/fact

 ;  For a matrix A, there an infinity of solutions to Ax=0 of the form kx', 
 ;  where k is real. This is the case here, and so we throw away one of 
 ;  the linear equations, and replace it with an equation setting the 
 ;  ground level population to 1. Once the linear equations have been 
 ;  solved, the populations are renormalised such the sum of the level 
 ;  populations is 1.
 ;  Note with the normalisation, we are solving Ax=b, where all elements 
 ;  of b are zero, except the last which is one.
 ;
  c(*,n_levels-1)=0.          ; set this row to zeros...
  c(0,n_levels-1)=1.          ; ...except for ground level
 ;
  b(n_levels-1)=1.            ; b is zero except for last element
 ;
 ;
 ; frac_zero is the fraction of zeros in the c array. If this is greater 
 ; than 0.5, then I will use LINBCG to solve the linear equations
 ;
  ind=where(c eq 0.)
  frac_zero=float(n_elements(ind))/float(n_elements(c))

  IF (frac_zero GT 0.5) AND (MIN(ABS(c)) GE 1.e-30) THEN BEGIN
    IF pp(n_levels-1) NE 1. THEN pp=pp1
    pp1=LINBCG(SPRSIN(c,thresh=1.e-30),b,pp)
  ENDIF ELSE BEGIN
    pp1=invert(transpose(c),status)#b   ; very marginally faster this way
  ENDELSE

  pp=pp1/fact         ; re-scale populations
  pp=pp/total(pp)     ; force total population to be 1
  pop(j,i,*)=pp       ; temp first, then dens

ENDFOR

ENDFOR

;-------------------------------------------------------------
;VK: above was (almost) uninterrupted POP_SOLVER.PRO}
;-------------------------------------------------------------

return,pop
end
