function popsol,temp,pres,chianti,n_e=n_e,$
	radtemp=radtemp,dilute=dilute,pe_ratio=pe_ratio
;+
;function	popsol
;	returns the level populations at given temperature and
;	pressure for a set of ions of a given element
;
;usage
;	pop=popsol(temp,pres,chianti,n_e=n_e,
;	radtemp=radtemp,dilute=dilute,pe_ratio=pe_ratio)
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
;	radtemp	[INPUT] The blackbody radiation field temperature
;		* default is 6000 K
;	dilute	[INPUT] a dilution factor for photoionization, based
;		on the distance above a star's surface -- see CHIANTI
;		routine R2W()
;		* if not set, assumed to be 0
;	pe_ratio [INPUT] proton-to-electron ratio
;		* if not given, assumed to be 1.0
;
;A note about the default value of PE_RATIO:
;	The original version of this program, POP_SOLVER.PRO, sets
;	the default to be 0.83, which is the limit for a high-T
;	plasma with cosmic abundances.  Here we exclude this factor
;	because we want to store the emissivities without reference
;	to their ion balance and abundance, on which this factor
;	depends.  The intensity of a line transition from upper level j
;	to lower level i is in general
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
;
;history
;	POP_SOLVER.PRO written by Peter Young
;{..
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
;..	(iv) to make LINBCG work well, I need to scale the c-matrix and 
;..	level populations in order to make the scaled populations close 
;..	to 1
;..	(v) Vectors of temperatures and densities can be input instead 
;..	of just single values. The electron excitation array only has to 
;..	be loaded for different temperatures, and so considerable time 
;..	savings can be made if a number of densities are input at once.
;..	LINBCG routine
;..	--------------
;..	  This is IDL's routine for solving linear equations with a sparse 
;..	matrix. A difference with INVERT is that it requires an initial 
;..	guess for the solution, and studying the performance of the 
;..	routine when dealing with the dielectronic recombination (DR) 
;..	files has revealed that this initial guess can lead to problems.
;..	  For some of the DR files (fe_22d is an example), solving the 
;..	level balance equations leads to many levels having zero 
;..	population. If my initial guess at the (scaled) level population 
;..	was all 1's, then using LINBCG gave final populations which were 
;..	not 0's, but very small numbers, some of which were negative. 
;..	Another option was to set the initial guess to 1 for the ground 
;..	level, and zero for all other levels. This solved the zero 
;..	population problem, but the iteration for one ion (fe_24d) became 
;..	very long (20secs for just one inversion), which was corrected by 
;..	going back to all 1's for the initial guess.
;..	  My solution for this was to check the t_type array. By summing 
;..	along columns, I arrive at a vector whose elements will be non-zero 
;..	if there are transitions which excite or de-excite this level, and 
;..	zero otherwise. I use this to set the initial guess for LINBCG: 1's 
;..	if the t_type vector are non-zero; 0's if they're zero. This 
;..	solves both of the problem cases above.
;..	PROTON RATES
;..	------------
;..	  To include the proton rates, it is necessary to have the 
;..	proton-to-electron ratio. This needs to be calculated before the 
;..	call to pop_solver, and the resulting ratio(s) passed through 
;..	'pe_ratio' in the common block 'proton'. If the ratio is not passed 
;..	then a default value of 0.83 is assumed.
;..	  Note that there is no keyword to switch off proton rates (i.e., 
;..	no /NOPROT keyword). To switch off proton rates, it is necessary 
;..	to set pstr=-1. This should be done by the calling routine.
;..	Ver 1, PRY 29-Mar-99
;..	Ver 2, PRY 30-Apr-99, added call to get_prot_rates
;..	Ver 3, PRY 15-Dec-99, added deu to upsilon common block in order 
;..		to be consistent with the main Chianti routines.
;..	Ver 4, PRY 9-May-00, corrected problem with threshold when dealing 
;..		with sparse matrices. Basically values less than 1.e-30 in 
;..		the c-matrix were being set to zero and giving rise to 
;..		NaN's in certain circumstances.
;..	Ver.5, PRY 14-Jul-00, changed elvl common block to the elvlc common 
;..	        block which is now the Chianti standard. Also, when 
;..	        descaling upsilons, the routine now uses the Delta-E from 
;..	        the .splups file.
;..	Ver.6, PRY 9-Aug-00, changed routine to deal better with the 
;..	        dielectronic recombination files
;..	Ver.7, PRY 17-Aug-00, routine does not call LINBCG now if radtemp 
;..	        is non-zero.
;..	Ver.8, PRY 29-Aug-00, the sparse matrix section has been disabled.
;..	Ver.9, PRY 12-Nov-01, calls routine proton_dens() to calculate the 
;..	        proton to electron ratio.
;..	Ver.10, PRY, 6-Dec-01, corrected bug when there are more levels 
;..	        in .splups file than in .elvlc file (ZnXXV).
;..	CONTACT:
;..	Peter Young, CfA, pyoung@cfa.harvard.edu
;..}
;
;	modified to eliminate common blocks and the call to SETUP_ION,
;	  and altered the I/O interface to mesh with PINTofALE, and
;	  changed default PE_RATIO to 1.0 (Vinay Kashyap, Jun02)
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
  print,'       dilute=dilution_factor,pe_ratio=proton_to_electron_ratio)'
  print,'  returns level populations at given temperature and pressure'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
t=temp & if keyword_set(n_e) then xne=double(n_e) else xne=double(pres/t)

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
ecm=chianti.ecmobs & jj=chianti.j & mult=2*jj+1
splstr=chianti.splstr
prot_struc=chianti.psplstr
;
deu=chianti.deryd
t_type=chianti.trans & c_ups=chianti.c_ups & splups=chianti.spups

;	special keywords
;	RADTEMP
if not keyword_set(radtemp) then radtemp=6000.D
;	DILUTE
if not keyword_set(dilute) then dilute=0.
;	PE_RATIO
if not keyword_set(pe_ratio) then pe_ratio=dblarr(nt)+1.0

;COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;COMMON wgfa, wvl,gf,a_value
;COMMON upsilon, splstr
;COMMON radiative, radtemp,dilute
;COMMON proton, prot_struc, pe_ratio

;------------------------------------------------------
;{VK: from here on down is uninterrupted POP_SOLVER.PRO
;------------------------------------------------------

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


;;; pp is no longer needed, so I've commented out the following

;pp_set = TOTAL(splups,3)
;pp_set = TOTAL(pp_set,2)
;pp=MAKE_ARRAY(n_levels,value=0.,/double)   ; inital guess for pp. Used by 
                                           ; linbcg.pro
;ind = WHERE(pp_set NE 0)
;pp[ind] = 1.

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

;________________
; create a ppr array for the proton rates
;
ppr=MAKE_ARRAY(n_levels,n_levels,nt,/double)
IF n_tags(prot_struc) NE 0 THEN BEGIN
 ;
  IF (n_elements(pe_ratio) NE nt) THEN BEGIN
    print,'%POP_SOLVER: WARNING, pe_ratio size does not match temp'
    print,n_elements(pe_ratio),nt
    pe_ratio=dblarr(nt)+0.83        ;--- assume default value
  ENDIF
 ;
  FOR i=0,n_elements(prot_struc)-1 DO BEGIN
    l1=prot_struc[i].lvl1-1
    l2=prot_struc[i].lvl2-1
    de=ABS(prot_struc[i].de)
    descale_all,t,prot_struc,i,prate
    IF ecm(l1) GT ecm(l2) THEN BEGIN
      ppr[l1,l2,*]=prate*pe_ratio
      ppr[l2,l1,*]=prate*pe_ratio*mult[l2]/mult[l1]* $
           exp(de*13.61/8.617/10.^(-5)/t)
    ENDIF ELSE BEGIN
      ppr[l2,l1,*]=prate*pe_ratio
      ppr[l1,l2,*]=prate*pe_ratio*mult[l2]/mult[l1]* $
           exp(de*13.61/8.617/10.^(-5)/t)
    ENDELSE
  ENDFOR
ENDIF
;________________


qq=MAKE_ARRAY(n_levels,n_levels,nt,/double)

;______________
; Create a qq array
;
l1=splstr.lvl1-1
l2=splstr.lvl2-1
de=ABS(splstr.de)
kte=(1./(hck*ryd2cm*de)) # t
ind_pos=where(ecm[l2] GT ecm[l1])
ind_neg=where(ecm[l2] LT ecm[l1])
xx=dblarr(n_elements(splstr),nt)
yy=xx
IF ind_neg[0] NE -1 THEN BEGIN
  xx[ind_neg,*]=(8.63d-6/(mult[l1[ind_neg]])#(1./sqrt(t)))
  yy[ind_neg,*]=8.63d-6* exp(-1./kte[ind_neg,*]) * $
       1./( mult[l2[ind_neg]] # sqrt(t) )
ENDIF
IF ind_pos[0] NE -1 THEN BEGIN
  yy[ind_pos,*]=(8.63e-6/mult[l2[ind_pos]]) # (1./sqrt(t))
  xx[ind_pos,*]=8.63e-6* exp(-1./kte[ind_pos,*]) * $
       1./(mult[l1[ind_pos]] # sqrt(t))
ENDIF
FOR i=0,n_elements(splstr)-1 DO BEGIN
  IF (l1[i] LE n_levels-1) AND (l2[i] LE n_levels-1) THEN BEGIN
    descale_all,t,splstr,i,ups
    qq[l1[i],l2[i],*]=xx[i,*]*ups
    qq[l2[i],l1[i],*]=yy[i,*]*ups
  ENDIF
ENDFOR



FOR j=0,nt-1 DO BEGIN             ;; loop over temperatures



; FOR m=0,n_levels-1 DO BEGIN
;   FOR n=m+1,n_levels-1 DO BEGIN
;     decm=ecm(n)-ecm(m)
;     deryd = ABS(deu[n,m])
;    ;
;     IF decm LT 0. THEN BEGIN
;       IF t_type(n,m) GT 0 THEN BEGIN
;          temp=t(j)
;          xt=temp/(hck*ryd2cm*deryd)
;          descale_ups,n,m,xt,ups
;          qq(m,n)=8.63e-6*ups/(mult(m)*sqrt(temp))                 ; de-exc
;          qq(n,m)=8.63e-6*ups*exp(-1./xt)/(mult(n)*sqrt(temp))     ; exc
;       ENDIF
;     ENDIF 
;     IF decm GT 0. THEN BEGIN
;       IF t_type(n,m) GT 0 THEN BEGIN
;          temp=t(j)
;          xt=temp/(hck*ryd2cm*deryd)
;          descale_ups,n,m,xt,ups
;          qq(n,m)=8.63e-6*ups/(mult(n)*sqrt(temp))                 ; de-exc
;          qq(m,n)=8.63e-6*ups*exp(-1./xt)/(mult(m)*sqrt(temp))     ; exc
;       ENDIF
;     ENDIF
;   ENDFOR
; ENDFOR


mqqj=DBLARR(n_levels)
maaj=mqqj
fact=mqqj+1.d
leave = mqqj
enter = mqqj

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
;  For a matrix A, there an infinity of solutions to Ax=0 of the form kx', 
;  where k is real. This is the case here, and so we throw away one of 
;  the linear equations, and replace it with an equation setting the 
;  ground level population to 1. Once the linear equations have been 
;  solved, the populations are renormalised such the sum of the level 
;  populations is 1.
;  Note with the normalisation, we are solving Ax=b, where all elements 
;  of b are zero, except the last which is one.
;

FOR i=0,nxne-1 DO BEGIN            ;; loop over densities
 ;
  c=aa + xne(i)*qq[*,*,j] + xne(i)*ppr[*,*,j] + aax    ; creates the c matrix
 ;
  ind=where(c NE 0.)
  frac_not_zero=float(n_elements(ind))/float(n_elements(c))
 ;
 ;
 ; The following checks to see whether c will be treated as a sparse matrix, 
 ; or just a normal one.
 ;
  IF (frac_NOT_zero LE 0.) AND (MIN(ABS(c[ind])) GE 1.d-30) $
     AND (dilute EQ 0.d) THEN BEGIN
   ;
   ; The following gets the factors which will be used for the scaling.
   ;
    FOR k=0,n_levels-1 DO BEGIN
       IF qq(0,k) EQ 0.d THEN qq_gr=1.d-14 ELSE qq_gr=qq(0,k)
       mqqj(k)=MAX([qq_gr,aax(0,k)/xne(i)])
       IF qq(k,0) EQ 0.d THEN qq_gr=1.d-14 ELSE qq_gr=qq(k,0)
       maaj(k)=MAX([REFORM(aa(k,*)),qq_gr*xne(i)])
    ENDFOR
   ;
   ; fact is a vector containing the scaling factors for each level
   ;
    indfact = where(mqqj NE 0.d)
    fact[indfact]=maaj[indfact]/mqqj[indfact]/xne[i]
   ;
    diag = -total(c,2)
    c[findgen(n_levels),findgen(n_levels)] = diag
   ;
    FOR k=0,n_levels-1 DO c[*,k]=c[*,k]/fact
   ;
    c(*,n_levels-1)=0.d         ; set this row to zeros...
    c(0,n_levels-1)=1.d         ; ...except for ground level
   ;
    b(n_levels-1)=1d0           ; b is zero except for last element
   ;
   ;
   ; If pop_solver has been given several temps or densities, then it's 
   ; better to use the populations from the previous temp or density as 
   ; the initial guess.
   ;
;    IF n_elements(pp1) NE 0 THEN pp=pp1
    pp1=LINBCG(SPRSIN(c,thresh=1d-30),b,pp)
   ;
    pp1=pp1/fact         ; re-scale populations
   ;
  ENDIF ELSE BEGIN
   ;
    diag = -total(c,2)
    c[findgen(n_levels),findgen(n_levels)] = diag
   ;
    c(*,n_levels-1)=0d0         ; set this row to zeros...
    c(0,n_levels-1)=1d0         ; ...except for ground level
   ;
    b(n_levels-1)=1d0           ; b is zero except for last element
   ;
    pp1=invert(transpose(c),status)#b
   ;
  ENDELSE

  pp=pp1/total(pp1)     ; force total population to be 1
  pop[j,i,*]=pp       ; temp first, then dens

ENDFOR

ENDFOR

;-------------------------------------------------------------
;VK: above was uninterrupted POP_SOLVER.PRO}
;-------------------------------------------------------------

return,pop
end
