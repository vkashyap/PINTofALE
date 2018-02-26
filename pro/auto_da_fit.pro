function auto_da_fit, x,y,ysig=ysig,pos=pos,wdt=wdt,flx=flx,thaw=thaw,effar=effar,effwgrid=effwgrid,$
perrp=perrp,perrm=perrm,perrc=perrc,werrp=werrp,werrm=werrm,werrc=werrc,epithet=epithet,$
ferrp=ferrp,ferrm=ferrm,ferrc=ferrc,crng=crng,type=type,algo_type=algo_type,conlev=conlev,$
comment=comment,ties=ties,_extra=e

;+
;function    auto_da_fit 
;       command line-based function to fit Gaussians and Lorentzians (or
;       combinations, or other 3-parameter functions) to lines in a spectrum 
;       and determine line fluxes. This is done automatically at the command
;       line, not interactivley with a GUI as in FITLINES() 
;
;       Given a set of lines, AUTO_DA_FIT() will identify which lines are
;       blended (hardcoded as 5*wdt from pos). Individual lines and blended 
;       groups are then fit successively and errors are computed for each.
;       Parameter freezing and thawing can be controlled with keyword THAW.
;       Otherwise, all parameters are considered free. A second blended group 
;       identificatin is performed after the preliminary fit to ensure that 
;       there are no components mistakingly left out because of an initial
;       wdt value that is too small. 
;       
;       It is HIGHLY RECOMMENDED that  routines LINEID(), ID_TO_FITPAR() are
;       used to prepare inputs. It is also advisable to visually inspect fits 
;       using FITLINES() after automatic fits are made. See USAGE EXAMPLES
;
;       returns a structure of the form, identical to FITLINES():  
;          {POS,PERRP,PERRM,PERRC,$
;	   FLX,FERRP,FERRM,FERRC,$
;	   WDT,WERRP,WERRM,WERRC,$
;	   THAW,TYPE,TIES,EPITHET,CONLEVX,CONLEV,CONSIGX,CONSIG,COMMENT}
;	All these are also available as output individually as keywords.
;	See keyword descriptions for what the variables mean.
;
;syntax  
;     fitstr = auto_da_fit(x,y,ysig=ysig,pos=pos,wdt=wdt,flx=flx,thaw=thaw,$
;     effar=effar,effwgrid=effwgrid,perrp=perrp,perrm=perrm,perrc=perrc,$
;     werrp=werrp,werrm=werrm,werrc=werrc,ties=ties,epithet=epithet,$
;     ferrp=ferrp,ferrm=ferrm,ferrc=ferrc,crng=crng, algo_type=algo_type,$
;     conlev=conlev)
;
;parameters
;	x	[INPUT; required] absissa, e.g., wavelength or energy
;	y	[INPUT; required] count spectrum Y(X)
;		* size of Y must match that of X
;keywords
;	ysig	[INPUT] errors on Y
;		* default is sqrt(abs(Y)+0.75)+1.
;	pos	[I/O] best-fit line positions
;		* length of vector shows number of components
;	perrp	[I/O] 1-sided error on POS (POS+PERRP is upper bound)
;	perrm	[I/O] 1-sided error on POS (POS-PERRM is lower bound)
;	perrc	[I/O] DCHI threshold used in computing PERR?
;	flx	[I/O] best-fit fluxes in the lines
;	ferrp	[I/O] 1-sided error on FLX (FLX+FERRP is upper bound)
;	ferrm	[I/O] 1-sided error on FLX (FLX-FERRM is lower bound)
;	ferrc	[I/O] DCHI threshold used in computing FERR?
;	wdt	[I/O] best-fit widths (sigma, core-radius, etc.) of the lines
;	werrp	[I/O] 1-sided error on WDT (WDT+WERRP is upper bound)
;	werrm	[I/O] 1-sided error on WDT (WDT-WERRM is lower bound)
;	werrc	[I/O] DCHI threshold used in computing WERR?
;		* on input, FLX, WDT, and ?ERR? are forced to match the length
;		  of POS: excess elements are thrown away, and insufficient
;		  lengths are made up by padding with first elements, 1s, etc.
;		* on output ?ERRM are identical to ?ERRP >unless< the
;		  projected errors have been calculated using ERORS, in
;		  which case the values may be different
;		* on output, places where ?ERRC contain 0's are those where
;		  the computed errors are 1-sigma formal curvature errors
;	thaw	[I/O] integer array signaling frozen (0) or thawed parameter (1)
;		* refers to sequences of {POS,WDT,FLX} for each component in
;		  a 3-parameter model (cf. X3MODEL) -- whatever goes for the
;		  appropriate user-defined function.
;		* length must match 3*length(POS)
;		* default is to freeze defined input, thaw all new additions
;	type	[I/O] type of model ('gauss', 'lorentz', etc.)
;		* default is 'gauss'
;	epithet	[I/O] label for each component
;		* labels, obtained, for example, with IDLABEL
;		* Merriam-Webster> a characterizing word or phrase accompanying
;		  or occurring in place of the name of a person or thing
;	ties	[I/O] any ties between parameters?
;	conlev	[I/O] the continuum that was taken out of the spectrum
;		* CONLEV must match the size of X and Y else ignored
;	consig	[I/O] error on CONLEV
;		* default for CONSIG is sqrt(abs(CONLEV)+0.75)+1.
;		* NOTE: CONLEV and CONSIG are compressed using SPLAC in
;		  the output that gets returned via the structure.  That
;		  structure therefore also has appropriate abscissae
;		  CONLEVX and CONSIGX.
;	comment	[OUTPUT] descriptive string
;	_extra	[INPUT ONLY] use this to pass defined keywords to subroutines
;		PICKRANGE: XSIZE, YSIZE, WID, DYNRNG
;		LINEREM: POSve, NEGve
;		FIT_LEVMAR: ITMAX, CHITHR, DUMB
;		LEVMARQ: JUMPUP, JUMPDN, SVDTHR
;		MK_3MODEL: MISSING
;		MK_SLANT: ANGLE, BETAP
;		ERORS: VERBOSE
;		STAMPLE: NUTHIN, NODAY, NOTIME, NOUSER, NOPACK, STACOL,
;			 STASIZ, STATHK
;		MK_LORENTZ: BETAP, MISSING, NORM
;               MCERROR
;		IS_KEYWORD_SET
;
;usage example
;       
;       PoA> ;      restore V711 Tau spectrum 
;       PoA> restore, !ARDB+'demo_v711Tau.save',/v ; restore V711 Tau spectrum 
;       PoA> 
;       PoA> ;      restore PoA line ID structure 
;       PoA> restore, !ARDB+'example_2.save',   /v ; restore PoA line ID structure    
;       PoA> 
;       PoA> ;      prepare inputs with ID_TO_FITPAR()
;       PoA> id_to_fitpar, lineid_idstr, pars, ties, thaw, pos=pos,wdt=wdt,$
;       PoA>  epithet=epithet ,shackle='flux,position',flx=flx,ties=ties, leeway = 'position,width',$ 
;       PoA>  perr = 0.6,/numer, lsfwdt = '0.01',werr=0.003,nozero=1
;       PoA>
;       PoA> ;      call AUTO_DA_FIT() 
;       PoA> fstr = auto_da_fit( LAM_M1P_V711TAU,SPC_M1P_V711TAU, pos=pos, wdt=wdt ,flx=flx                 
;       PoA>           effar=EFFAR_M1P,effwgrid=WVLAR_M1P ,crng = [2.4,3.4], verbose = 0,$             
;       PoA>           funcs = 'x3model',function_name='x3model',conlev=conlev ,/normflx,$            
;       PoA>           type = replicate('beta=2.5',n_elements(pos)),ties=ties ,thaw = thaw ,mciter = 10.0) 
;       PoA> 
;       PoA> ;      visually inspect with FITLINES() 
;       PoA> newstr = fitlines(LAM_M1P_V711TAU,SPC_M1P_V711TAU,oldstr=fstr,/dumb) 
;       
;history          
;        liwei lin (Jul05) some code lifted directly from FITLINES()
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	changed name from AUTO_DA_FET to AUTO_DA_FIT (VK; Apr2006)
;	some bug fixes (VK; Apr2009)
;-

;            check parameters 
forward_function mpfit,mpfitfun

;     	ALGO
;  algo_type='LevMarq+SVD'			;default
;  algval=['LevMarq','Curvefit','MPFIT','tifsrhg']  
;  if keyword_set(algo) then fitalgo=strtrim(algval(algo(0)),2) else fitalgo=algo_type     
;  i0=-1L
;  i00=strpos(strlowcase(fitalgo),'lev',0) & i0 = i0 > i00
;  i01=strpos(strlowcase(fitalgo),'mar',0) & i0 = i0 > i01
;  i02=strpos(strlowcase(fitalgo),'lm',0) & i0 = i0 > i02
;  i03=strpos(strlowcase(fitalgo),'svd',0) & i0 = i0 > i03
;  i10=strpos(strlowcase(fitalgo),'idl',0) & i0 = i0 > i10
;  i11=strpos(strlowcase(fitalgo),'curve',0) & i0 = i0 > i11
;  i12=strpos(strlowcase(fitalgo),'mpfit',0) & i0 = i0 > i12 
;  if i00 ge 0 or i01 ge 0 or i02 ge 0 or i03 ge 0 then algo_type='LevMarq+SVD'
;  if i10 ge 0 or i11 ge 0 then algo_type='IDL-Curvefit'
;  if i12 ge 0 then algo_type ='MPFIT'
;  if algo_type eq 'IDL-Curvefit' and nf ne 0 then begin
;    message,algo_type+' does not understand parameter freezing',/info
;  endif
;  if i0 lt 0 then message,$
;  	fitalgo+': not understood.  reverting to '+algo_type,/info
;  if v ge 1 then message,'using fitting algorithm: '+algo_type,/info


np = n_elements(pos) 
if is_keyword_set(thaw) then begin 
    nth = n_elements(thaw)  
    if nth ne np*3 then begin 
       message, 'Length of keyword THAW  must match 3*length(pos)',/info 
       return, 0L 
    endif  
    if total(thaw) le 0 then thaw = fltarr(nth)+1 
endif else begin 
    thaw = fltarr(np*3.0)+1 
endelse

;      estimate cotinuum level 
areff = interpol(effar, effwgrid, x) 
oo = where(x lt crng(1) and x gt crng(0)) 
conlev = total(y(oo))*areff/total(areff(oo)) 
if not is_keyword_set(consig) then consig = sqrt(abs(CONLEV)+0.75)+1.
dwvl = mean(x[1:*]-x)


;       prep for fitting
nw=n_elements(pos)  & aa=fltarr(nw*3)  
 tmp=indgen(nw*3) & tmp2=indgen(nw) & tmp3 = fltarr(nw*3)+1.0
aa(3.0*tmp2+0L)=pos  & aa(3.0*tmp2+1L)=wdt   & aa(3.0*tmp2+2L)=flx 
aaf = fltarr(nw,3)  & errfu = fltarr(nw,3) & errfl = fltarr(nw,3)
lst = fltarr(nw)    & xmn = min(x)     & xmx = max(x) 
i = 0  & iter = 0
oy = y-conlev >0
dwvl = mean(x[1:*]-x)
xmin = min(x) & xmax = max(x) 
if not keyword_set(algo_type) then algo_type = 'LevMarq+SVD'

;       loop through each line and get the best-fit
while total(lst) ne nw do begin 
    ;       establish subset for fit (there must be a better way)
    coda: 
    xmn=pos(i)-5*wdt(i) > xmin
    xmx=pos(i)+5*wdt(i) < xmax

    oo = where(pos gt xmn and pos lt xmx, noo) 
    roe= where(pos le xmn or pos ge xmx, rnoe) 
    if total(noo) gt 1 then begin 
       oe = oo & noe = noo & noo = noo-1
       While noe gt noo  do begin 
          noo = noe 
          mini = min(pos(oo),ii) & ii = oo(ii)
          maxi = max(pos(oo),xi) & xi = oo(xi) 
          xmn=mini-5*wdt(ii) > xmin
          xmx=maxi+5*wdt(xi) < xmax             
          oe = where(pos gt xmn and pos lt xmx, noe) 
          roe= where(pos le xmn or pos ge xmx, rnoe) 
       endwhile
       oo = oe            
    endif  
    ;        establish xrange for fit 
    fo = where(x gt xmn and x lt xmx)
    xx = x(fo) & yy = oy(fo)
    ;        establish free parameters 
    tmp3 = fltarr(nw*3)+1.0
    tmpfrz = ([tmp(roe*3.+0L),tmp(roe*3.+1L),tmp(roe*3.+2L)])
    tmp3([tmp(roe*3.+0L),tmp(roe*3.+1L),tmp(roe*3.+2L)]) = 0 
    tf  = where(thaw eq 0,ntf) 
    if ntf gt 0 then  tmp3(tf) = tmp3(tf)-1.0 
    freeze = where(tmp3 le 0) 
    iter   = iter+1 
    ;        if  subset is indeed different after coda then re-fit  
    if iter eq 2 then if noe eq old_noe then goto,coda2
    ;old_noe = noe 
    if keyword_set(noe) then old_noe = noe 
    ;        fit by algorithm of choice (default is fit_levmar) 
    if algo_type eq 'LevMarq+SVD' then fit_levmar,xx,yy,$
        aa,yfunc,freeze=freeze,erra=erra,chisq=chisq,sig=sigy,/dumb,/normflx,$
          funcs='x3model',function_name='x3model',ties=ties, _extra=e  
    if algo_type eq 'IDL-Curvefit' then yfunc=curve_fit(xx,yy,ywt,aa,erra,chisq=chisq,_extra=e)
    if algo_type eq 'MPFIT' then begin
         if keyword_set(mpfunc) then userfuncs=mpfunc else begin 
         message, 'mpfunc keyword required for mcerror calculation with MPFIT',/info
         return, ostr
         if v ge 5 then quiet = 0 else quiet = 1
         a=MPFITFUN(userfuncs,xx,yy,sigy,aa,PERROR=erra,yfit=yfit,$
	    	parinfo=parinfo,functarg=functarg,quiet=quiet, _extra=e)
         endelse
    endif
    ;        compute errors with algorithm of choice (default is mcerror)         
      mcerror,xx,yy,aa,erru,errl,ysig=zsig,freeze=freeze,dchi=dchisq,$
  	    algo=algo,yfunc=ymod,erra=erra,$
  	    ties=ties,type=type,mpfunc=userfuncs,x_seg=x_seg, rmfstr=rmfstr,$
              parinfo=parinfo,functarg=functarg,mcsims=mcsims,/dumb,$
              funcs='x3model',function_name='x3model', _extra=e,/normflx
    ;        re-asses subset considering newly fit parameters
    stop
    if iter lt 2 then goto,coda 
    coda2: 
    ;          update list 
    lst(oo) = 1  & zi = where(lst eq 0) & i = zi(0) & iter = 0
endwhile 

;      prepare  outputs
  pos = aa[3.0*tmp2+0L] & wdt = aa[3.0*tmp2+1L] & flx = aa[3.0*tmp2+2L]/dwvl 
  perrp = erru[3.0*tmp2+0L] & perrm = errl[3.0*tmp2+0L] & perrc = 1.0
  werrp = erru[3.0*tmp2+1L] & werrm = errl[3.0*tmp2+1L] & errc = 1.0
  ferrp = erru[3.0*tmp2+2L]/dwvl - aa[3.0*tmp2+2L]/dwvl
  ferrm = aa[3.0*tmp2+2L]/dwvl - errl[3.0*tmp2+2L]/dwvl
  ferrc = 1.0 & werrc = 1.0 & conlev = conlev & conlevx = x

if not keyword_set(comment) then comment = 'Preliminary fit using AUTO_DA_FIT()' 
if not keyword_set(epithet) then epithet = '['+strtrim(pos,2)+']' 
ostr=$
	{POS:pos, PERRP:perrp, PERRM:perrm, PERRC:perrc,$
	 FLX:flx, FERRP:ferrp, FERRM:ferrm, FERRC:ferrc,$
	 WDT:wdt, WERRP:werrp, WERRM:werrm, WERRC:werrc,$
	 THAW:thaw, TYPE:type, TIES:ties, EPITHET:epithet,$
	 CONLEVX:conlevx, CONLEV:conlev,CONSIGX:x,CONSIG:consig, COMMENT:comment}

return, ostr
end 
