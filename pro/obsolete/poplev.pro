function poplev,temp,pres,chianti,n_e=n_e
;+
;function	poplev
;		returns the level populations at given temperature and
;		pressure for a set of ions of a given element
;
;usage
;	pop=poplev(temp,pres,chianti,n_e=n_e)
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
;
;requires
;	DESCALE2_UPS, a CHIANTI routine
;
;history
;	originally written as POPULATE.PRO, a CHIANTI routine,
;		by Ken Dere (v2: Mar96)
;	modified to eliminate common blocks, mesh with SCAR, and
;	fltarr->dblarr (to remove the -ve population levels)
;		by Vinay Kashyap (Nov96)
;-

message,'OBSOLETE: use POPSOL() instead',/info

np=n_params(0)
if np lt 3 then begin
  print,'Usage: pop=poplev(temperature,pressure,chianti_structure,n_e=edensity)'
  print,'  returns level populations at given temperature and pressure'
  return,-1L
endif

;inputs
t=temp & if keyword_set(n_e) then xne=n_e else xne=pres/t

;decode CHIANTI parameters
;
l1=chianti.lev1 & l2=chianti.lev2
w1=chianti.wvl & g1=chianti.gf & a1=chianti.a
nlev=max([l1,l2]) & wvl=fltarr(nlev,nlev) & gf=wvl & a_value=wvl
wvl(l1-1,l2-1)=w1 & gf(l1-1,l2-1)=g1 & nl1=n_elements(l1)
for j=0,nl1-1 do a_value(l1(j)-1,l2(j)-1)=a_value(l1(j)-1,l2(j)-1)+a1(j)
;
ecm=chianti.ecmobs & jj=chianti.j & mult=2*jj+1
;
t_type=chianti.trans & c_ups=chianti.c_ups & splups=chianti.spups

;-------------------------------------------------------------
;VK: from here on down is POPULATE.PRO, with 3 exceptions
;-------------------------------------------------------------
hck=1.98648e-16/1.38062e-16
;
n_elvl=n_elements(ecm)
wsize=size(wvl)
n_wgfa1=wsize(1)
n_wgfa2=wsize(2)
usize=size(splups)
n_ups1=usize(1)
n_ups2=usize(2)
;
n_levels=min([n_elvl,n_wgfa1,n_wgfa2,n_ups1,n_ups2])
;
;
;----------------------------------------------------
;VK: changed fltarr to dblarr
;----------------------------------------------------
c=dblarr(n_levels,n_levels)
d=dblarr(n_levels,n_levels)
b=dblarr(n_levels)
pop=dblarr(n_levels)
;----------------------------------------------------
;
;  population rates for A values
;
for m=0,n_levels-1 do begin
  for n=0,n_levels-1 do begin
      decm=ecm(n)-ecm(m)
    if(decm gt 0.) then begin
;        n level higher than m
;        contribution of level n to population of level m
;
      c(n,m)=c(n,m)+a_value(m,n)
;
   endif else if(decm lt 0.) then begin
;       n level lower than m
      c(m,m)=c(m,m)-a_value(n,m)
   endif
endfor  ; n
endfor  ;m
;
;
;  population rates for collision values
;
for m=0,n_levels-1 do begin
  for n=0,n_levels-1 do begin
      decm=ecm(n)-ecm(m)
    if(decm lt 0.) then begin
;        n level lower than m
;        contribution of level n to population of level m
;
      if(t_type(n,m) gt 0) then begin
         xt=t/(hck*abs(decm))
	 ;----------------------------------------------------
	 ;VK: call DESCALE2_UPS to avoid common blocks
	 ;----------------------------------------------------
         ;descale_ups,n,m,xt,ups
	 cc=c_ups(n,m) & tt=t_type(n,m) & ss=reform(splups(n,m,*))
         descale2_ups,xt,cc,tt,ss,ups
	 ;----------------------------------------------------
         c(m,m)=c(m,m)-xne*8.63e-6*ups/(mult(m)*sqrt(t))
         c(n,m)=c(n,m)+xne*8.63e-6*ups*exp(-1./xt)/(mult(n)*sqrt(t))
      endif
   endif else if(decm gt 0.) then begin
;        n level higher than m
;
      if(t_type(m,n) gt 0) then begin
        xt=t/(hck*abs(decm))
	;----------------------------------------------------
	;VK: call DESCALE2_UPS to avoid common blocks
	;----------------------------------------------------
        ;descale_ups,m,n,xt,ups
	cc=c_ups(m,n) & tt=t_type(m,n) & ss=reform(splups(m,n,*))
        descale2_ups,xt,cc,tt,ss,ups
	;----------------------------------------------------
        c(m,m)=c(m,m)-xne*8.63e-6*ups*exp(-1./xt)/(mult(m)*sqrt(t))
        c(n,m)=c(n,m)+xne*8.63e-6*ups/(mult(n)*sqrt(t))
      endif
    endif
  endfor  ; n
endfor  ;m
;
;  to set total population = 1.
;
m=n_levels-1
  for n=0,n_levels-1 do begin
    c(n,m)=1.  
endfor
b(n_levels-1)=1.
;
;
pop=invert(transpose(c))#b     ; very marginally faster this way

return,pop
end
