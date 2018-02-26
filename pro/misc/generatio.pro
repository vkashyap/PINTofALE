pro generatio,fx,rcode,rx,fxerr=fxerr,rxerr=rxerr,verbose=verbose,$
	dfx_mul=dfx_mul,dfx_add=dfx_add, _extra=e
;+
;procedure	generatio
;	a generalized mechanism to compute flux ratios and errors
;	from given input fluxes and a well-defined definition of
;	how the ratios are constructed
;
;	all of the input fluxes are assumed to be independent
;	of each other and uncorrelated
;
;syntax
;	generatio,fx,rcode,rx,fxerr=fxerr,rxerr=rxerr,verbose=v,$
;	dfx_mul=dfx_mul,dfx_add=dfx_add
;
;parameters
;	fx	[INPUT; required] input fluxes
;	rcode	[INPUT] string array describing which of the
;		input fluxes must be considered only as ratios
;		* basic format is: "sP#[,sP#[,...]]" where
;		  -- "#" is an integer flag describing the ratio being
;		     constructed
;		  -- "P" is a positional descriptor and can take on
;		     values N (for numerator) or D (for denominator)
;		  -- "s" stands for the sign with which the flux
;		     goes into the ratio "+" or "-"
;		* e.g., to construct a simple ratio FX[2]/FX[1],
;		  RCODE=['','+D1','+N1']
;		* e.g., to construct two hardness ratios
;			FX[3]/FX[1] and
;			(FX[3]-FX[1])/(FX[3]+FX[1]),
;		  RCODE=['','+D1,-N2,+D2','','+N1,+N2,+D2']
;		* if size is incompatible with FX, then no action
;		  is taken.
;	rx	[OUTPUT] array of ratios constructed using FX and RCODE
;
;keywords
;	fxerr	[INPUT] 1-sigma errors on FX
;		* size must match that of FX.  otherwise,
;		  -- if single-element, assumed to represent
;		     -- a fractional error if >0 and <1
;		     -- a percentage error if >1 and <100
;		     -- an abs(constant) error if >100 or <0
;		  -- ignored otherwise
;	rxerr	[OUTPUT] 1-sigma errors on RX, computed only if FXERR
;		is given and is legal
;	verbose	[INPUT] controls chatter
;	dfx_mul	[INPUT] multiplicative factor by which to jiggle FX
;		while computing partial derivatives to propagate
;		errors (default is 0.05)
;		* if scalar, taken to be same factor for all FX
;		* if vector, must be of same size as FX
;	dfx_add	[INPUT] additive factor by which to jiggle FX
;		while computing partial derivatives to propagate
;		errors (default is 0.05*FXERR)
;		* if scalar, taken to be same offset for all FX
;		* if vector, must be of same size as FX
;	_extra	[JUNK] here only to avoid crashing the program
;
;example
;	;this makes the ratios FX[0]/FX[1] and (FX[1]-FX[0])/(FX[1]+FX[0])
;	fx=[10.,5.] & fxerr=0.1 & rcode=['+N1,-N2,+D2','+D1,+N2,+D2']
;	generatio,fx,rcode,rx,fxerr=fxerr,rxerr=rxerr
;	for i=0,n_elements(rx)-1 do print,rx[i],' +- ',rxerr[i]
;
;	Warning: The numeral suffix in the ratio specification is
;	only as a placeholder to determine uniqueness and does _not_
;	translate to the index in the output.  For example, try:
;	fx=[1.,2.,1.,3.,1.,4.,1.,5.,1.,6.]
;	rcode=['+D1','+N1','+D2','+N2','+N3','+D3','+N4','+D4','+D5','+N5']
;	generatio,fx,rcode,rx,verbose=100
;	print,rcode & print,rx
;	generatio,fx,reverse(rcode),rx_reverse,verbose=100
;	print,reverse(rcode) & print,rx_reverse
;
;subroutines
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Nov'02)
;	cosmetic changes (VK; Dec'02)
;	won't spit out error messages if RCODE is 'X', though the FX
;	  corresponding to that is ignored (VK; Apr'03)
;	added warning about ratio sequence -- it's a feature, not a bug;
;	  corrected bug in case of sequence number exceeding 9 (VK; Dec'03)
;	now lets specifying just numerator or just denominator (VK,LL; Jun'04)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	changed DFX_ADD and DFX_MUL to be arrays of same size as FX
;-

;	usage
ok='ok' & np=n_params() & nfx=n_elements(fx) & nrc=n_elements(rcode)
szrc=size(rcode,/tname)
if np le 1 then ok='Insufficient parameters' else $
 if nfx eq 0 then ok='FX undefined' else $
  if szrc ne 'STRING' then ok='RCODE must be a string array'
if ok ne 'ok' then begin
  print,'Usage: generatio,fx,rcode,rx,fxerr=fxerr,rxerr=rxerr,verbose=v,$'
  print,'       dfx_mul=dfx_mul,dfx_add=dfx_add'
  print,'  compute ratios from input fluxes FX'
  print,'  RCODE format is: sP#[,sP#[,...]], e.g., "[-N1,+D1]"'
  if np ne 0 then message,ok,/informational
  return
endif

;	check inputs
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
;
ok='ok'
if nrc eq 0 then ok='RCODE undefined; returning' else $
 if nrc ne nfx then ok='RCODE scrambled; returning'
if ok ne 'ok' then begin
  if vv gt 1 then message,ok,/informational & return
endif

;  parse RCODE
idxrat=0
;	first find out how many ratios we are dealing with
for i=0L,nrc-1L do begin	;{for each flux value
  cc=strtrim(rcode[i],2)
  ss=str_sep(cc,',') & nss=n_elements(ss)
  for j=0L,nss-1L do begin	;{for each use of flux in some ratio
    c1=strmid(ss[j],0,1)
    c2=strlowcase(strmid(ss[j],1,1))
    c3=strmid(ss[j],2,10000)
    ok='ok'
    if c1 eq '' then ok='' else $
     if strlowcase(c1) eq 'x' then ok='X' else $
      if c1 ne '+' and c1 ne '-' then $
      ok='sP#: format must start with sign' else $
       if c2 eq '' or c2 eq ' ' then ok='' else $
        if c2 ne 'n' and c2 ne 'd' then $
        ok='sP#: format must specify Numerator or Denominator' else $
	 if c3 eq '' or c3 eq ' ' then ok='' else $
          if c3 eq '' then $
          ok='sP#: format must end with unique numerical index'
    if ok eq 'ok' then begin	;(correct format
      rr=c1+'F'+strtrim(i,2)
      if is_keyword_set(idxrat) then begin	;(not the first ratio
	k=long(c3) & oo=where(idxrat eq k,moo)
	if moo eq 0 then begin		;(new ratio
 	  idxrat=[idxrat,k] & nidx=n_elements(idxrat)
	  if c2 eq 'n' then begin
	    ratnum=[ratnum,rr] & ratden=[ratden,'']
	  endif else begin
	    ratden=[ratden,rr] & ratnum=[ratnum,'']
	  endelse
	endif else begin		;new)(existing ratio
	  if c2 eq 'n' then ratnum[oo[0]]=ratnum[oo[0]]+rr else $
	  	ratden[oo[0]]=ratden[oo[0]]+rr
	endelse				;existing ratio)
      endif else begin			;not first)(first ratio
        k=long(c3) & idxrat=[k]
	if c2 eq 'n' then begin
	  ratnum=[rr] & ratden=['']
	endif else begin
	  ratden=[rr] & ratnum=['']
	endelse
      endelse					;whichth ratio)
      if vv gt 50 then print,rr,j,i
    endif else begin		;format OK)(not OK
      if ok ne '' and ok ne 'X' then message,ok,/informational
    endelse			;format)
  endfor			;J=0,NSS-1}
endfor				;I=0,NRC-1}

;	create new variables F0..(NFX-1)
for i=0L,nfx-1L do jnk=execute('F'+strtrim(i,2)+'=fx['+strtrim(i,2)+']')

;	and make the ratios
nrat=n_elements(ratnum)
if nrat eq 0 then begin
  if vv gt 1 then message,'No ratios to compute; returning.',/informational
  return
endif
rx=fltarr(nrat) & rxerr=rx
for i=0L,nrat-1L do begin
  cc='rx['+strtrim(i,2)+']='
  if strtrim(ratnum[i],2) ne '' then cc=cc+'('+ratnum[i]+')' else cc=cc+'1.'
  if strtrim(ratden[i],2) ne '' then cc=cc+'/('+ratden[i]+')'
  if vv gt 50 then print,cc & jnk=execute(cc)
  if vv gt 200 then stop,'Halting.  type .CON to continue'
endfor

if arg_present(fxerr) then begin		;(propagate errors?
  nfxe=n_elements(fxerr)
  if nfxe eq nfx then fxe=fxerr else begin	;(just what is FXERR?
    if nfxe eq 1 then begin		;(scale the errors
      if fxerr[0] lt 0 then fxe=0.*fx+abs(fxerr[0]) else $
       if fxerr[0] lt 1 then fxe=fx*fxerr[0] else $
	if fxerr[0] lt 100 then fxe=fx*fxerr[0]/100. else $
	 fxe=0.*fx+fxerr[0]
    endif else begin			;NFXE=1)(NFXE.ne.1
      if vv gt 1 then message,'FXERR incompatible with FX; ignoring',$
	/informational
      return
    endelse				;NFXE ? 1)
  endelse					;FXERR)
  ;
  ;	compute partial derivatives
  ;	by brute force because we don't know the exact form of the
  ;	ratio until run time.  so compute the value of the ratio
  ;	at two points upstream and downstream of each Fx, and
  ;	compute delta(ratio)/delta(flux) as the partials
  dmul=0.05+fltarr(nfx) & dadd=0.05*fxerr
  if n_elements(dfx_mul) eq 1 then dmul[*]=dfx_mul[0] else $
   if n_elements(dfx_mul) eq nfx then dmul=dfx_mul
  if n_elements(dfx_add) eq 1 then dadd[*]=dfx_add[0] else $
   if n_elements(dfx_add) eq nfx then dadd=dfx_add
  drx_dfx=fltarr(nrat,nfx) & fxp=fltarr(nfx) & fxm=fltarr(nfx)
  for i=0L,nfx-1L do jnk=execute('FXP['+strtrim(i,2)+']=(1.+dmul[i])*fx['+$
	strtrim(i,2)+']+dadd[i]')
  for i=0L,nfx-1L do jnk=execute('FXM['+strtrim(i,2)+']=(1.-dmul[i])*fx['+$
	strtrim(i,2)+']-dadd[i]')
  for i=0L,nrat-1L do begin	;{for each ratio in question
    for j=0L,nfx-1L do begin	;{for each flux value in the list
      ;this seems like overkill.. for k=0L,nfx-1L do jnk=execute('F'+strtrim(k,2)+'=fx['+strtrim(k,2)+']')
      ;	get ratio at upstream point
      jnk=execute('F'+strtrim(j,2)+'=FXP['+strtrim(j,2)+']')
      cc='RXP=' & if strtrim(ratnum[i],2) ne '' then cc=cc+'('+ratnum[i]+')' else cc=cc+'1.'
      if strtrim(ratden[i],2) ne '' then cc=cc+'/('+ratden[i]+')'
      jnk=execute(cc)
      jnk=execute('F'+strtrim(j,2)+'=FXM['+strtrim(j,2)+']')
      ;	get ratio at downstream point
      cc='RXM=' & if strtrim(ratnum[i],2) ne '' then cc=cc+'('+ratnum[i]+')' else cc=cc+'1.'
      if strtrim(ratden[i],2) ne '' then cc=cc+'/('+ratden[i]+')'
      jnk=execute(cc)
      ;jnk=execute('RXM=('+ratnum[i]+')/('+ratden[i]+')')
      ;	partial derivative
      drx_dfx[i,j]=(rxp-rxm)/(fxp[j]-fxm[j])
      ;	reset Fj
      jnk=execute('F'+strtrim(j,2)+'=fx['+strtrim(j,2)+']')
    endfor			;J=0,NFX-1}
    rxerr[i]=sqrt(total(drx_dfx[i,*]^2*fxerr^2))
  endfor			;I=0,NRAT-1}
endif						;propagate errors)

return
end
