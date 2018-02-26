pro libmodel,x,a,f,pder,anew,type=type,verbose=verbose,renorm=renorm,$
	betap=betap,angle=angle,vrot=vrot,x0=x0,pbreak=pbreak,$
	phase=phase,$
	angrise=angrise,angfall=angfall,ybase=ybase,yoff=yoff,$
	_extra=e
;+
;procedure	libmodel
;	computes F(X;A) and dF/dA for a library of model functions.
;	procedure written to be compatible with FIT_LEVMAR, and also
;	with IDL's CURVEFIT, GHRS-IDL's WFIT, etc.
;
;syntax
;	libmodel,x,a,f,pder,anew,type=type,verbose=verbose,$
;	missing=missing,/norm,/fwhm,betap=betap,angle=angle,vrot=vrot,$
;	Xo=Xo,pbreak=pbreak,NH=NH,wvlar=wvlar,effar=effar,exptime=exptime,$
;	dellam=dellam,X0=X0,angrise=angrise,angfall=angfall,ybase=ybase,$
;	yoff=yoff
;
;parameters
;	x	[INPUT; required] points at which to generate models
;	a	[INPUT; required] array of parameter values
;		* should match the number expected from TYPE.
;		* no checks are made to verify that the user is not shooting
;		  hirself in the foot.
;	f	[OUTPUT; required] output f=f(x;a)
;	pder	[OUTPUT; optional] partial derivatives for each parameter
;	anew	[OUTPUT; optional] return a new set of A in which the
;		normalization parameters have all been adjusted to account
;		for the supplied RENORM
;
;keywords
;	type	[INPUT; required] string array denoting the model functions
;		to be called.
;		* this should rightly be a parameter, because it is
;		  a very necessary part of the program.  however, most
;		  fitting functions do not allow that degree of flexibility
;		  and so we are forced to use this hack.
;		* some models may be multiplicative (e.g., 'absorb'), and
;		  their location in the TYPE array matters.  a multiplicative
;		  model multiplies the model that is defined after it.  so
;		  to apply the same absorption model to a powerlaw and a
;		  gaussian, set
;		    TYPE=['absorb','power','absorb','gauss']
;		  and then use TIES=['a6=a1'],FREEZE=[6] to make sure the same 
;		  absorption is applied to both components
;	betap	[I/O] keyword passed w/o checking to MK_LORENTZ and MK_SLANT
;	angle	[I/O] keyword passed w/o checking to MK_SLANT
;	vrot	[I/O] keyword passed w/o checking to MK_ROGAUSS
;	x0	[I/O] keyword passed w/o checking to MK_POWLAM,MK_BKNPOWER
;	pbreak	[I/O] keyword passed w/o checking to MK_POWLAM,MK_BKNPOWER
;	phase	[I/O] keyword passed w/o checking to MK_SINUSOID
;	angrise	[I/O] keyword passed w/o checking to MK_LCECLIPSE
;	angfall	[I/O] keyword passed w/o checking to MK_LCECLIPSE
;	ybase	[I/O] keyword passed w/o checking to MK_LCECLIPSE
;	yoff	[I/O] keyword passed w/o checking to MK_LCECLIPSE
;	verbose	[INPUT; default=0] controls chatter
;	renorm	[INPUT; default=1] renormalize the output by multiplying
;		by this factor
;		* the adjusted normalizations for each component that
;		  correspond to this renormalization are returned in ANEW
;
;	_extra	[INPUT] pass defined keywords to subroutines
;		MK_GAUSS: MISSING, FWHM, NORM
;		MK_LORENTZ: BETAP, MISSING, NORM
;		MK_SLANT: ANGLE, BETAP, MISSING, NORM
;		MK_ROGAUSS: VROT, FWHM, NORM, MISSING
;		MK_POWLAM: XO,PBREAK,NH,WVLAR,EFFAR,EXPTIME,DELLAM
;		MK_BKNPOWER: XO,NOBREAK,VERBOSE
;		MK_POLY1D: X0
;		MK_SINUSOID: PHASE, MISSING
;		MK_SINE: MISSING
;		MK_LCECLIPSE: ANGRISE, ANGFALL, YBASE, YOFF
;		MK_ABSORB: FH2,HE1,HEII,FANO,IKEV,WAM,BAM,MAM,NOHEH,ABUND
;		MK_EFOLD: VERBOSE
;		MK_BBANG: VERBOSE
;		MK_BBKEV: VERBOSE
;		MK_RECIPROCALPOW: VERBOSE
;		MK_SPLINE1D: XLOC,USEINIT,TENSION,VERBOSe
;
;subroutines
;	MK_GAUSS
;	MK_LORENTZ
;	MK_SLANT
;	MK_ROGAUSS
;	MK_POWLAM
;	MK_BKNPOWER
;	MK_POLY1D
;	MK_LCECLIPSE
;	MK_ABSORB
;	MK_EFOLD
;	MK_BBANG
;	MK_BBKEV
;	MK_RECIPROCALPOW
;	MK_SINUSOID
;	MK_SINE
;
;history
;	vinay kashyap (Aug01; based on X3MODEL and MK_3MODEL)
;	added call to MK_SINUSOID (VK; Sep04)
;	added call to MK_LCECLIPSE, allowed return of extra parameters
;	  via keywords (VK; May07)
;	added call to MK_BKNPOWER via "pow" and "pe"; added MK_ABSORB
;	  via "absorb"; allowed for multiplicative as well as additive
;	  models (VK; Jul08)
;	added call to MK_EFOLD via "ef", "efold", and "exp" and to
;	  MK_BBANG and MK_BBKEV via "bbang" and "bbkev" (VK; Aug08)
;	added call to MK_RECIPROCALPOW via "reci/pow" and MK_SPLINE1D
;	  via "spline=KNOTS"; added parameter ANEW and keyword RENORM
;	  (VK; Sep08)
;	added call to MK_SINE (VK; Jan13)
;-

;	usage
ok='ok' & np=n_params(0) & nx=n_elements(x) & na=n_elements(a)
nt=n_elements(type)
if np lt 3 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is undefined' else $
  if na eq 0 then ok='A is undefined' else $
   if nt eq 0 then ok='TYPE is required'
if ok ne 'ok' then begin
  print,'Usage: libmodel,x,a,f,pder,anew,type=type,verbose=verbose,$'
  print,'       missing=missing,/norm,/fwhm,betap=betap,angle=angle,$'
  print,'       vrot=vrot,Xo=Xo,pbreak=pbreak,NH=NH,wvlar=wvlar,$'
  print,'       effar=effar,exptime=exptime,dellam=dellam,X0=X0,phase=phase'
  print,'       angrise=angrise,angfall=angfall,ybase=ybase,yoff=yoff'
  print,'  evaluate F and partial(F)/partial(A) at X for parameters A'
  print,'  TYPE may be any of:'
  print,'  	+Gaussian (g; gauss; Gauss; etc -- 3 parameters: mean, sig, peak)'
  print,'  	+Lorentzian (l; lorentz; Lorentzian; etc -- 3 parameters: mean, width, peak)'
  print,'  	+Beta-profile (b; beta=value; Beta:value; etc -- 4 parameters: mean, width, peak)'
  print,'  	+asymmetric beta (sl; slant=angle,beta; slant:angle,beta; etc -- 5 parameters: mean, width, peak, angle, beta)'
  print,'  	+sinusoidal (sinu; sinu=phase, sinusoid:phase; etc. -- 4 parameters: ampl, freq, dcoffset, phase)'
  print,'  	+sine (sine; -- 3 parameters: ampl, freq, phase)'
  print,'  	+spline (spline=N; spline:N, sp@N; etc. -- N parameters: spline knots at c0..cN-1)'
  print,'  	+rotation broadened Gaussian (r; rot=value; etc -- 4 parameters: mean, sig, peak, vrot)'
  print,'  	+power-law in wavelength (pw; pw=break; pw@break; etc -- 4 parameters: a1, a2, norm, break)'
  print,'  	+power-law in general (pow; power; pe; pe=break; pe@break; etc -- 4 parameters: norm, a1, a2, break)'
  print,'  	+reciprocal power-law (reciproc; recipow; rp=join; ep@join; etc -- 4 parameters: norm, a1, a2, xjoin)'
  print,'   	+polynomial (pol=N; poly:N; polynomial=N; etc -- N parameters: coefficients c0..cN-1)'
  print,'  	+eclipse (ec; eclipse; etc -- 7 parameters: midpt,hwidth,depth,angfall,angrise,ybase,yoff)'
  print,'	+efold (ef; efold; exp; exponential; etc., 3 parameters: norm, scale, offset)'
  print,'	+bbkev (bbk; bbkev; 2 parameters: norm, temperature)'
  print,'	+bbang (bba; bbang; 2 parameters: norm, temperature)'
  print,'  	*absorb (a; ab; absorb; etc -- 1 parameter: NH)'
  if np ne 0 then message,ok,/info
  return
endif

;	verbosity
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1

;	output
f=fltarr(nx) & pder=fltarr(nx,na)

;	call each model function in sequence
k0=-1L & k1=0L	;parameter count
;
multfactor=fltarr(nx)+1.0	;multiplicative factor
;
anew=a & allnorm=1.0 & if keyword_set(renorm) then allnorm=abs(renorm[0])
;
for i=0L,nt-1L do begin			;{stepping through the TYPEs
  cc=strtrim(type[i],2) & c1=strlowcase(strmid(cc,0,1))
  addmodel=1
  case c1 of				;{usually enough to tell

    'a': begin			;{ISMabsorption
      addmodel=0
      npar=1	;number of parameters (NH)
      k0=k0+1L & k1=k0+npar-1L
      if na-1 lt k1 then begin
	message,'insufficient parameters specified for MK_ABSORB()',/info
	return
      endif
      NH=a[k0]
      if np lt 3 then mm=mk_absorb(x,NH, _extra=e) else $
       mm=mk_absorb(x,NH,pd, _extra=e)
    end				;ISMabsorption}

    'g': begin			;{Gaussian
      npar=3	;number of parameters (pos,wdt,hgt)
      k0=k0+1L & k1=k0+npar-1L
      if na-1L lt k1 then begin
	message,'insufficient parameters specified for MK_GAUSS()',/info
	return
      endif
      pos=a[k0] & wdt=a[k0+1] & hgt=a[k0+2]*allnorm & anew[k0+2]=hgt
      if np lt 4 then mm=mk_gauss(x,pos,wdt,hgt, _extra=e) else $
	mm=mk_gauss(x,pos,wdt,hgt,pd, _extra=e)
    end				;Gaussian}

    'l': begin			;{Lorentzian
      npar=3	;number of parameters (pos,wdt,hgt)
      k0=k0+1L & k1=k0+npar-1L
      if na-1L lt k1 then begin
	message,'insufficient parameters specified for MK_LORENTZ()',/info
	return
      endif
      pos=a[k0] & wdt=a[k0+1] & hgt=a[k0+2]*allnorm & anew[k0+2]=hgt
      if np lt 4 then mm=mk_lorentz(x,pos,wdt,hgt, _extra=e) else $
	mm=mk_lorentz(x,pos,wdt,hgt,pd, _extra=e)
    end				;Lorentzian}

    'b': begin			;{modified beta-profile or blackbody
      mtyp='unknown'
      if strpos(cc,'be') ge 0 then mtype='beta'
      if strpos(cc,'bb') ge 0 then mtype='blackbody'
      if strpos(cc,'blackbody') ge 0 then mtype='blackbody'
      case mtype of
	'beta': begin		;{modified beta-profile
          npar=4	;number of parameters (pos,wdt,hgt,betap)
          c=str_sep(cc,'=')                     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')     & nc=n_elements(c)
          if nc gt 1 then begin
	    betap=float(c[1]) & npar=3	;because BETAP is fixed
          endif else betap=1.	;the default
          k0=k0+1L & k1=k0+npar-1L
          if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_LORENTZ()',/info
	    return
          endif
          pos=a[k0] & wdt=a[k0+1] & hgt=a[k0+2]*allnorm & if npar eq 4 then betap=a[k0+3]
	  anew[k0+2]=hgt
          if np lt 4 then mm=mk_lorentz(x,pos,wdt,hgt,betap=betap,_extra=e) else $
	    mm=mk_lorentz(x,pos,wdt,hgt,pd,betap=betap, _extra=e)
	end			;MK_LORENTZ(BETAP=BETAP)}
	'blackbody': begin		;{different types of blackbody functions
	  if strpos(cc,'kev') ge 0 then begin	;(blackbody in keV
	    npar=2
	    k0=k0+1L & k1=k0+npar-1L
	    if na-1L lt k1 then begin
	      message,'insufficient parameters specified for MK_BBKEV()',/info
	      return
	    endif
	    norm=a[k0]*allnorm & temperature=a[k0+1]
	    anew[k0]=norm
	    if np lt 4 then mm=mk_bbkev(x,norm,temperature, verbose=vv) else $
	      mm=mk_bbkev(x,norm,temperature,pd, verbose=vv)
	  endif					;BBKEV)
	  if strpos(cc,'ang') ge 0 then begin	;(blackbody in Ang
	    npar=2
	    k0=k0+1L & k1=k0+npar-1L
	    if na-1L lt k1 then begin
	      message,'insufficient parameters specified for MK_BBANG()',/info
	      return
	    endif
	    norm=a[k0]*allnorm & temperature=a[k0+1]
	    anew[k0]=norm
	    if np lt 4 then mm=mk_bbang(x,norm,temperature, verbose=vv) else $
	      mm=mk_bbang(x,norm,temperature,pd, verbose=vv)
	  endif					;BBANG)
	end				;blackbody functions}
	else: message,cc+' : model not understood',/informational
      endcase

    end				;modified beta-profile or blackbody}

    's': begin
      mtyp='unknown'
      if strpos(cc,'sl') ge 0 then mtyp='slant'
      if strpos(cc,'sinu') ge 0 then mtyp='sinusoid'
      if strpos(cc,'sine') ge 0 then mtyp='sine'
      if strpos(cc,'sp') ge 0 then mtyp='spline'
      case mtyp of			;{slant, sine, sinusoid, or spline

	'slant': begin			;{asymmetrical beta-profile
          npar=5	;number of parameters (pos,wdt,hgt,betap,angle)
          c=str_sep(cc,'=')                     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')     & nc=n_elements(c)
          if nc gt 1 then begin
	    cb=str_sep(c[1],',') & ncb=n_elements(cb)
	    if ncb eq 1 then cb=str_sep(c[1],' ') & ncb=n_elements(cb)
	    if ncb eq 1 then cb=str_sep(c[1],'	') & ncb=n_elements(cb)
	    if strtrim(cb[0],2) ne '' then begin
	      angle=float(cb[0]) & npar=4	;because ANGLE is fixed
	    endif else angle=0.
	    if ncb gt 1 then begin
	      if strtrim(cb[1],2) ne '' then betap=float(cb[1])
	      npar=npar-1	;because BETAP is also fixed
	    endif ;else betap=1.
          endif
          k0=k0+1L & k1=k0+npar-1L
          if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_SLANT()',/info
	    return
          endif
          pos=a[k0] & wdt=a[k0+1] & hgt=a[k0+2]*allnorm
	  anew[k0+2]=hgt
          if npar eq 4 then angle=a[k0+3]
          if npar eq 5 then betap=a[k0+4]
          if np lt 4 then $
	    mm=mk_slant(x,pos,wdt,hgt,angle=angle,betap=betap, _extra=e) else $
	    mm=mk_slant(x,pos,wdt,hgt,pd,angle=angle,betap=betap, _extra=e)
        end				;asymmetrical beta-profile}

	'sinusoid': begin		;{sinusoid
	  npar=4	;number of parameters (ampl,freq,offset,phase)
          c=str_sep(cc,'=')                     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')     & nc=n_elements(c)
          if nc gt 1 then begin
	    phase=float(c[1]) & npar=3	;because PHASE is fixed
          endif else phase=0.	;the default
          k0=k0+1L & k1=k0+npar-1L
          if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_SINUSOID()',/info
	    return
          endif
          ampl=a[k0]*allnorm & freq=a[k0+1] & offset=a[k0+2]
	  anew[k0]=ampl
	  if npar eq 4 then phase=a[k0+3]
          if np lt 4 then $
	    mm=mk_sinusoid(x,ampl,freq,offset,phase=phase, _extra=e) else $
	    mm=mk_sinusoid(x,ampl,freq,offset,pd,phase=phase, _extra=e)
	end				;sinusoid}

	'sine': begin			;{sine
	  npar=3	;number of parameters (ampl,freq,offset,phase)
          k0=k0+1L & k1=k0+npar-1L
          if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_SINUSOID()',/info
	    return
          endif
          ampl=a[k0]*allnorm & freq=a[k0+1] & phase=a[k0+2]
	  ;anew[k0]=ampl
	  if np lt 3 then $
	    mm=mk_sine(x,ampl,freq,phase, _extra=e) else $
	    mm=mk_sine(x,ampl,freq,phase,pd, _extra=e)
	end				;sine}

	'spline': begin			;{1D spline
          c=str_sep(cc,'=')                         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'@')         & nc=n_elements(c)
          if nc gt 1 then begin
	    npar=long(c[1])
	    k0=k0+1L & k1=k0+npar-1L
	    if na-1L lt k1 then begin
	      message,'insufficient parameters specified for '+$
		'MK_SPLINE1D()',/info
	    endif
	    aa=a[k0:k1]*allnorm & anew[k0:k1]=aa
	    if np lt 3 then mm=mk_spline1d(x,aa, _extra=e) else $
		mm=mk_spline1d(x,aa,pd, _extra=e)
	  endif else begin
	    ;no other way to know!
	    message,cc+': unknown number of coefficients for MK_SPLINE1D()',/info
	  endelse

	end				;spline1d}

	else: message,cc+': model not understood',/informational

      end				;slant, sinusoid, spline}
    end

    'r': begin			;{rogauss or reciprocalpow
      mtyp='unknown'
      if strpos(cc,'ro') ge 0 then mtyp='rogauss'
      if strpos(cc,'re') ge 0 or strpos(cc,'rp') ge 0 or $
        (strpos(cc,'rec') ge 0 and strpos(cc,'pow') ge 0) then mtyp='recipow'
      case mtyp of
	'rogauss': begin	;{Rotationally convolved Gaussian
          npar=4	;number of parameters (pos,wdt,hgt,vrot)
          c=str_sep(cc,'=')                     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')     & nc=n_elements(c)
          if nc gt 1 then begin
	    vrot=float(c[1]) & npar=3	;because VROT is fixed
          endif else vrot=0.	;the default
          k0=k0+1L & k1=k0+npar-1L
          if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_ROGAUSS()',/info
	    return
          endif
          pos=a[k0] & wdt=a[k0+1] & hgt=a[k0+2]*allnorm & if npar eq 4 then vrot=a[k0+3]
	  anew[k0+2]=hgt
          if np lt 4 then mm=mk_rogauss(x,pos,wdt,hgt,vrot=vrot,_extra=e) else $
	    mm=mk_rogauss(x,pos,wdt,hgt,pd,vrot=vrot, _extra=e)
	end			;Rotationally convolved Gaussian}

	'recipow': begin	;{reciprocal power law
	  npar=4	;number of parameters (norm,g1,g2,xjoin)
          c=str_sep(cc,'=')                     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')     & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'@')     & nc=n_elements(c)
	  if nc gt 1 then begin
	    xjoin=float(c[1]) & npar=3	;because XJOIN is fixed
	  endif
	  k0=k0+1L & k1=k0+npar-1L
	  if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_RECIPROCALPOW()',/info
	    return
	  endif
	  norm=a[k0]/allnorm & g1=a[k0+1] & g2=a[k0+2] & if npar eq 4 then xjoin=a[k0+3]
	  anew[k0]=norm
	  if np lt 4 then mm=mk_reciprocalpow(x,norm,g1,g2,xjoin,_extra=e) else $
	    mm=mk_reciprocalpow(x,norm,g1,g2,xjoin,pd,_extra=e)
	end			;reciprocal power law}

        else: message,cc+': model type not understood',/info

      endcase

    end				;rogauss or reciprocalpow}

    'p': begin			;{polynomial or power-law
      mtyp='unknown'
      if strpos(cc,'pl') ge 0 or strpos(cc,'power') ge 0 or $
	 strpos(cc,'pw') ge 0 or strpos(cc,'pe') ge 0 or $
	 strpos(cc,'pow') ge 0 then mtyp='power-law'
      if strpos(cc,'pol') ge 0 or strpos(cc,'poly') ge 0 then mtyp='polynomial'

      case mtyp of

	'power-law': begin		;{power-law
	  ;message,cc+': sorry, not implemented yet',/info
	  npar=4	;number of parameters (gamma1,gamma2,norm,pbreak)
          c=str_sep(cc,'=')                         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'@')         & nc=n_elements(c)
	  if nc gt 1 then begin
	    pbreak=float(c[1]) & npar=3	;because PBREAK is fixed
	  endif else pbreak=1.	;the default
	  k0=k0+1L & k1=k0+npar-1L
	  if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_POWLAM()/MK_BKNPOWER()',/info
	    return
	  endif
	  c2=strmid(cc,0,2) & ptyp='pw'
	  if c2 ne 'pw' and c2 ne 'pe' and c2 ne 'po' then begin
	    message,'Power-law form not understood; using MK_BKNPOWER()',/info
	    ptyp='po'
	  endif else ptyp=c2
	  if ptyp eq 'pw' then begin
	    pa1=a[k0] & pa2=a[k0+1] & pn=a[k0+2]*allnorm & anew[k0+2]=pn
	    if npar eq 4 then pbreak=a[k0+3]
	    if np lt 4 then $
		mm=mk_powlam(x,pa1,pa2,pn,pbreak=pbreak,_extra=e) else $
		mm=mk_powlam(x,pa1,pa2,pn,pd,pbreak=pbreak,_extra=e)
	  endif else begin
	    pn=a[k0]*allnorm & pa1=a[k0+1] & pa2=a[k0+2] & anew[k0]=pn
	    if npar eq 4 then pbreak=a[k0+3]
	    if np lt 4 then $
		mm=mk_bknpower(x,pn,pa1,pa2,pbreak,_extra=e) else $
		mm=mk_bknpower(x,pn,pa1,pa2,pbreak,pder,_extra=e)
	  endelse
	end				;power-law}

	'polynomial': begin
          c=str_sep(cc,'=')                         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,' ')         & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,'	')  & nc=n_elements(c)
          if nc eq 1 then c=str_sep(cc,':')         & nc=n_elements(c)
          if nc gt 1 then begin
	    npar=long(c[1])
	    k0=k0+1L & k1=k0+npar-1L
	    if na-1L lt k1 then begin
	      message,'insufficient parameters specified for '+$
		'MK_POLY1D()',/info
	    endif
	    aa=a[k0:k1] & aa[0]=aa[0]*allnorm & anew[k0]=aa[0]
	    if np lt 3 then mm=mk_poly1d(x,aa, _extra=e) else $
		mm=mk_poly1d(x,aa,pd, _extra=e)
	  endif else begin
	    ;no other way to know!
	    message,cc+': unknown number of coefficients for MK_POLY1D()',/info
	  endelse
	end

        else: message,cc+': model type not understood',/info

      endcase
    end				;polynomial or power-law}

    'e': begin			;{eclipse || efold
      mtyp='unknown'
      if strpos(cc,'ec') ge 0 or strpos(cc,'eclipse') ge 0 then mtyp='eclipse'
      if strpos(cc,'ef') ge 0 or strpos(cc,'efold') ge 0 or $
	 strpos(cc,'exp') ge 0 or strpos(cc,'exponential') ge 0 then mtyp='efold'
      case mtyp of

	'eclipse': begin	;{eclipse light curve
          npar=7
          k0=k0+1L & k1=k0+npar-1L
          if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_LCECLIPSE()',/info
	    return
          endif
          if na ge 1 then midpt=a[k0]
          if na ge 2 then hwidth=a[k0+1]
          if na ge 3 then depth=a[k0+2]
          if na ge 4 then angfall=a[k0+3]
          if na ge 5 then angrise=a[k0+4]
          if na ge 6 then ybase=a[k0+5]
          if na ge 7 then yoff=a[k0+6]
	  ;	note: RENORM has no effect on this model
          if np lt 4 then mm=mk_lceclipse(x,midpt,hwidth,depth,$
	      angfall=angfall,angrise=angrise,ybase=ybase,yoff=yoff,_extra=e) else $
	    mm=mk_lceclipse(x,midpt,hwidth,depth,pd,$
	      angfall=angfall,angrise=angrise,ybase=ybase,yoff=yoff,_extra=e)
	end			;eclipse light curve}

	'efold': begin		;{exponential decay
	  npar=3
	  k0=k0+1L & k1=k0+npar-1L
	  if na-1L lt k1 then begin
	    message,'insufficient parameters specified for MK_EFOLD()',/info
	    return
	  endif
	  if na ge 1 then norm=a[k0]*allnorm
	  if na ge 2 then scale=a[k0+1]
	  if na ge 3 then yoff=a[k0+2]
	  anew[k0]=norm
	  if np lt 4 then mm=mk_efold(x,norm,scale,yoff, verbose=vv) else $
	    mm=mk_efold(x,norm,scale,yoff,pd, verbose=vv)
	end			;exponential decay}

	else: message,cc+': model type not understood',/informational

      endcase
    end				;eclipse || efold}
    else: message,cc+': model type not understood',/info

  endcase				;case C1}

  if addmodel eq 1 then begin
    f=f+multfactor*mm
    multfactor=0.*f+1.0	;reset multiplicative factor
  endif else begin
    multfactor=mm	;multiply the next model by this factor
  endelse
  ;if addmodel eq 1 then f=f+mm else f=f*mm
  if np eq 4 then for j=0,npar-1 do pder[*,k0+j]=pd[*,j]
  k0=k1 & k1=k1+1L
endfor					;I=0,NT-1}

return
end
