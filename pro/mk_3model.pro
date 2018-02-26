function mk_3model,x,p,w,h,pder,group=group,delp=delp,type=type,asis=asis,$
	allcomp=allcomp, _extra=e
;+
;function	mk_3model
;	generate and return a compound model using the supplied
;	parameters of components made up of 3-parameter functions
;
;syntax
;	model=mk_3model(x,p,w,h,pder,group=group,delp=delp,type=type,/asis,$
;	/allcomp,/fwhm,/normflx,missing=missing,betap=betap,vrot=vrot)
;
;parameters
;	x	[INPUT; required] points where model is computed
;	p	[INPUT; required] positions of components
;	w	[INPUT; required] widths of components
;	h	[INPUT; required] heights of components
;	pder	[OUTPUT] contains 2D array (NX,3*NPAR) of partial derivatives
;		of each component wrt parameters
;
;keywords
;	group	[INPUT] grouping index of gaussians
;	delp	[INPUT] offsets from first position value in group
;		NOTE: GROUP and DELP are ignored if ASIS is set
;		NOTE: if GROUP and DELP are not given, ASIS is set
;	type	[INPUT; default='G*aussian'] what type of curve to generate
;		for the model.  other possibilities are
;			L*orentzian
;			B*eta-profile (set "Beta=<value>")
;			SL*ant (set "slant=<angle>,<beta>")
;			R*otation-convolved Gaussian (set "rot=<velocity>")
;			P*ower-law (set "pw=<break>" or "pe=<break>")
;			I*dentity
;			U*ser_defined (not yet implemented), etc.
;			SI*nusoid (set "sin=<phase>")
;		* if not defined, uses "Gauss"
;		* if insufficiently defined, uses last defined as default
;		  for the rest
;	asis	[INPUT] the default action is to verify that the p are
;		consistent with GROUP and DELP and to correct as necessary.
;		if ASIS is set, this is not done.
;	allcomp	[INPUT] if set, returns the sum of components and also each
;		component in a separate column of a 2D array of size
;		(N_ELEMENTS(P)+1,N_ELEMENTS(X))
;	_extra	[INPUT] pass defined variables to subroutine:-
;		MK_GAUSS: MISSING, FWHM, NORMFLX
;		MK_LORENTZ: BETAP, MISSING, NORMFLX
;		MK_SLANT: ANGLE, BETAP, MISSING, NORMFLX
;		MK_ROGAUSS: VROT, MISSING, FWHM, NORMFLX
;		MK_IDENT: GIGO
;		MK_SINUSOID: PHASE, MISSING
;		MK_MIXDECAY: YBKG
;
;usage summary
;	* call as a function
;	* returns composite model by adding up components
;	* input component parameters as sequence of arrays of
;	  positions, widths, and heights
;	* input grouping info by 2 additional arrays
;	* specify type of model with keyword TYPE
;	* partial derivatives are computed only if ALL 5 parameters are
;	  supplied in call
;	* example TYPEs:
;	  'Gauss', 'gaussian', 'gau', 'g'
;	  'Lorentz', 'lorentz', 'lorentzian', 'lor', 'l'
;	  'beta=1.8', 'BETA=2.4', 'B=2', 'b 1.8', 'b:1.8', 'b	1.8'
;	  'slant=45,1.8','SLANT=89,2.4', 's=10','s=,2', 's:3,2'
;	  'rotational velocity=1e-4', 'rot=30', 'r:30e5', 'R 1e-4'
;	  'pe=1', 'pw=12.3985', 'pw@12.4'
;	  'sin', 'sin=90'
;	  'mixdecay', 'm'
;	  'identity', 'i'
;
;subroutines
;	MK_GAUSS
;	MK_LORENTZ
;	MK_SLANT
;	MK_ROGAUSS
;	MK_POWLAM
;	MK_SINUSOID
;	MK_MIXDECAY
;	KILROY
;
;history
;	vinay kashyap (Oct96)
;	added parameter to return partial derivatives (VK; Nov96)
;	added _EXTRA, added call to MK_LORENTZ (VK; Oct98)
;	allowed setting BETAP for MK_LORENTZ via TYPE definition (VK; Dec98)
;	added MK_SLANT calls (VK; MarMM)
;	added MK_ROGAUSS, MK_POWLAM calls (VK; JunMM)
;	changed name from MK_MODEL to MK_3MODEL (VK; Aug01)
;	converted array notation to IDL 5; added MK_IDENT (VK; Apr02)
;	changed keyword NORM to NORMFLX (VK; Oct02)
;	added call to MK_SINUSOID (VK; Sep04)
;	added call to MK_MIXDECAY (VK; Apr17)
;-

np=n_params(0)
if np lt 4 then begin
  print, 'Usage: model=mk_3model(x,p,w,h,pder,group=group,delp=delp,type=type,$'
  print, '       /asis,/allcomp,/fwhm,/normflx,missing=missing,betap=betap,$'
  print, '       vrot=vrot,phase=phase)'
  print, '  generate compound model of multiple 3-parameter curves'
  if np eq 0 then return,-1L else return,0*x
endif

;how many models, how many bins?
z=p & nmod=n_elements(p) & nbin=n_elements(x)

;figure out which curve
mtyp=strarr(nmod)+'g'		;default is Gaussian
mtyp2=strarr(nmod)
modtyp=strarr(nmod)+'gauss'
ntype=n_elements(type)
if ntype eq 0 then type='gauss'
for i=0,nmod-1 do mtyp[i]=$
	strmid(strtrim(strlowcase(type[[i]]),2),0,1)
for i=0,nmod-1 do mtyp2[i]=$
	strmid(strtrim(strlowcase(type[[i]]),2),1,1)
for i=0,nmod-1 do modtyp[i]=type[[i]]

;bother with groups and all?
ngp=n_elements(group)
if ngp gt 0 then grp=group else grp=indgen(nmod)+1
ndp=n_elements(delp)
if ndp ne ngp then delp=grp*0.

ig=grp[uniq(grp,sort(grp))] & ng=n_elements(ig)
if keyword_set(allcomp) then m=fltarr(ng+1,nbin) else m=fltarr(nbin)

;do consistency check
if not keyword_set(asis) then begin
  for i=0,ng-1 do begin
    oo=where(ig[i] eq grp,moo)
    if moo gt 1 then begin
      pp=z[oo] & dp=delp[oo] & zz=pp[0]+dp
      gotcha=total(abs(zz-pp))
      c1='Input positions inconsistent with grouping...correcting'
      ;if gotcha gt 0 and keyword_set(type) then message,c1,/info
      if gotcha gt 0 then kilroy,dot='*'	;make a mark
      z[oo]=zz[*]
    endif
  endfor
endif

if np ge 5 then pder=fltarr(nbin,3*nmod)

for i=0,nmod-1 do begin
  pos=z[i] & wdt=w[i] & hgt=h[i] & jg=grp[i]
  typ=mtyp[i]
  case typ of
    'g': begin				;{GAUSSIAN
      if np lt 5 then mm=mk_gauss(x,pos,wdt,hgt, _extra=e) else $
        mm=mk_gauss(x,pos,wdt,hgt,pd, _extra=e)
    end					;GAUSSIAN}
    'l': begin				;{modified Lorentzian
      if np lt 5 then mm=mk_lorentz(x,pos,wdt,hgt, _extra=e) else $
	mm=mk_lorentz(x,pos,wdt,hgt,pd, _extra=e)
    end					;LORENTZIAN}
    'b': begin				;{MODIFIED BETA-PROFILE
      tip=modtyp[i]
      c=str_sep(tip,'=')                   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,' ')   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,'	') & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,':')   & nc=n_elements(c)
      if nc le 1 then betap=1. else betap=float(c[1])
      if np lt 5 then mm=mk_lorentz(x,pos,wdt,hgt,betap=betap,_extra=e) else $
	mm=mk_lorentz(x,pos,wdt,hgt,pd,betap=betap,_extra=e)
    end					;BETA-PROFILE}
    's': begin				;{ASYMMETRICAL BETA-PROFILE/SINUSOID
      tip=modtyp[i]
      if mtyp2[i] eq 'l' then begin	;(SLANT
        angle=0. & betap=1.	;the defaults
        c=str_sep(tip,'=')                   & nc=n_elements(c)
        if nc eq 1 then c=str_sep(tip,' ')   & nc=n_elements(c)
        if nc eq 1 then c=str_sep(tip,'	') & nc=n_elements(c)
        if nc eq 1 then c=str_sep(tip,':')   & nc=n_elements(c)
        if nc gt 1 then begin
	  cc=str_sep(c[1],',')                       & ncc=n_elements(cc)
	  if ncc eq 1 then cc=str_sep(c[1],' ')      & ncc=n_elements(cc)
	  if ncc eq 1 then cc=str_sep(c[1],'	') & ncc=n_elements(cc)
	  if strtrim(cc[0],2) ne '' then angle=float(cc[0])
	  if ncc gt 1 then if strtrim(cc[1],2) ne '' then betap=float(cc[1])
        endif
        if np lt 5 then $
	  mm=mk_slant(x,pos,wdt,hgt,angle=angle,betap=betap,_extra=e) else $
	  mm=mk_slant(x,pos,wdt,hgt,pd,angle=angle,betap=betap,_extra=e)
      endif				;SLANT)
      if mtyp2[i] eq 'i' then begin	;(SINUSOID
        c=str_sep(tip,'=')                   & nc=n_elements(c)
        if nc eq 1 then c=str_sep(tip,' ')   & nc=n_elements(c)
        if nc eq 1 then c=str_sep(tip,'	')   & nc=n_elements(c)
        if nc eq 1 then c=str_sep(tip,':')   & nc=n_elements(c)
	if nc le 1 then phase=0. else phase=float(c[1])
        if np lt 5 then mm=mk_sinusoid(x,pos,wdt,hgt,phase=phase,_extra=e) $
	  else mm=mk_sinusoid(x,pos,wdt,hgt,pd,phase=phase,_extra=e)
      endif				;SINUSOID)
    end					;SLant/SInusoid}
    'r': begin				;{ROTATIONALLY CONVOLVED GAUSSIAN
      tip=modtyp[i]
      c=str_sep(tip,'=')                   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,' ')   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,'	') & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,':')   & nc=n_elements(c)
      if nc le 1 then vrot=0. else vrot=float(c[1])
      if np lt 5 then mm=mk_rogauss(x,pos,wdt,hgt,vrot=vrot,_extra=e) else $
	mm=mk_rogauss(x,pos,wdt,hgt,pd,vrot=vrot,_extra=e)
    end					;ROGAUSS}
    'p': begin				;{BROKEN POWER-LAW
      pa1=pos & pa2=wdt & pn=hgt
      tip=modtyp[i]
      c=str_sep(tip,'=')                   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,' ')   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,'	') & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,':')   & nc=n_elements(c)
      if nc eq 1 then c=str_sep(tip,'@')   & nc=n_elements(c)
      if nc le 1 then pbreak=1. else pbreak=float(c[1])
      tip2=strmid(strtrim(c[0],2),0,2)
      ;if tip2 ne 'pw' and tip2 ne 'pe' then begin
      if tip2 ne 'pw' then begin
	message,'Power-law form not understood.  using PL[w]',/info
	tip2='pw'
      endif
      if tip2 eq 'pw' then begin	;(PL=PL(wavelength)
	if np lt 5 then mm=mk_powlam(x,pa1,pa2,pn,pbreak=pbreak,_extra=e) $
	  else mm=mk_powlam(x,pa1,pa2,pn,pd,break=pbreak,_extra=e)
      endif else begin			;)(PL=PL(energy)
	message,'PL(E) not implemented'
	;if np lt 5 then mm=mk_pownrg(x,pa1,pa2,pn,pbreak=pbreak,_extra=e) $
	;  else mm=mk_pownrg(x,pa1,pa2,pn,pd,break=pbreak,_extra=e)
      endelse				;PL)
    end					;BROKEN POWER-LAW}
    'm': begin				;{MIXDECAY
      x0=pos & tau0=wdt & y0=hgt
      if np lt 5 then mm=mk_mixdecay(x,x0,tau0,y0,_extra=e) else $
        mm=mk_mixdecay(x,x0,tau0,y0,pd,_extra=e)
    end					;MIXDECAY}
    'i': begin				;{IDENTITY GIGO
      if np lt 5 then mm=mk_ident(x,hgt,pos,wdt,_extra=e) else $
	mm=mk_ident(x,hgt,pos,wdt,pd,_extra=e)
    end					;IDENTITY GIGO}
    else : begin
      message,'sorry, '+strtrim(type,2)+' is not implemented yet',/info
      stop
      return,0*x
    end
  endcase
  if keyword_set(allcomp) then begin
    m[0,*]=m[0,*]+mm[*]
    m[jg,*]=m[jg,*]+mm[*]
  endif else m=m+mm
  if np ge 5 then begin
    pder[*,3*i]=pd[*,0] & pder[*,3*i+1]=pd[*,1] & pder[*,3*i+2]=pd[*,2]
  endif
endfor

return,m
end
