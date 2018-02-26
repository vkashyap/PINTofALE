pro id_to_fitpar,idstr,pars,ties,thaw,pos=pos,wdt=wdt,flx=flx,lsfwdt=lsfwdt,$
	shackle=shackle,leeway=leeway,perr=perr,werr=werr,ferr=ferr,epithet=epithet,$
	numer=numer,nozero=nozero,_extra=e
;+
;procedure	id_to_fitpar
;	given the ID structure resulting from LINEID, figures out
;	initial fitting parameters for use in, e.g., FITLINES, FIT_LEVMAR,
;	etc.
;
;syntax
;	id_to_fitpar,idstr,pars,ties,thaw,pos=pos,wdt=wdt,flx=flx,$
;	lsfwdt=lsfwdt,shackle=shackle,leeway=leeway,perr=perr,werr=werr,$
;	ferr=ferr
;
;parameters
;	idstr	[INPUT; required] ID structure containing the observed
;		wavelengths and appropriate matches
;	pars	[OUTPUT] initial values of all the parameters set up in
;		a single array for consumption by, e.g., FIT_LEVMAR or
;		X3MODEL.
;	ties	[OUTPUT] constraints on the parameters
;	thaw	[OUTPUT] integer array of same size as PARS, signaling
;		frozen (0) or thawed (1) parameter
;		* all parameters that are SHACKLEd are frozen
;
;keywords
;	pos	[OUTPUT] initial-guess positions of all lines
;	wdt	[OUTPUT] initial-guess widths of all lines
;	flx	[OUTPUT] initial-guess fluxes of all lines
;       epithet [OUTPUT] species label (e.g. 'OVII' or 'FeXIX') string array
;	lsfwdt	[INPUT] width of the instrument's line spread function, and
;		the minimum value of WDT
;		* if not defined, assumed to be same as PERR
;		  (on the theory that your position determination is uncertain
;		  on the same order as the LSF width)
;	shackle	[INPUT] specify how to handle multiple ID components:-
;		"position" => line centers are fixed relative to location
;		  of line of maximum strength
;		"width" => freeze all WDT at maximum determined value OR
;		  value of WERR if given
;		"flux" => relative strengths of fluxes are held fixed
;		* if array, only the first element is considered
;		* if integer (as in /shackle), all 3 of above are set
;		* any two (or all 3) may be specified simply as "flux & pos",
;		  "width and flux", "position,width,flux", etc.
;		* all affected parameters are frozen
;	leeway	[INPUT] specify how to handle multiple ID components:-
;		"position" => relative locations of line centers are held
;		  to within a narrow range of position of line of maximum
;		  strength, defined by PERR
;		"width" => all WDT are constrained to be > minimum determined
;		  width, and < maximum determined width, and this range may
;		  be further expanded by WERR
;		"flux" => relative strengths of fluxes are held to within
;		  a narrow range defined by FERR
;		* obviously, SHACKLE takes precedence if both are set
;		* if array, only the first element is used
;		* if integer (e.g., /leeway), all of above are set
;		* unlike SHACKLE, all affected parameters are kept THAWed.
;	perr	[INPUT; default=1.0] slippage allowed in POS
;		* setting this is crucial in case the IDs include lines
;		  from higher orders or unknown IDs, because this will
;		  tell the program that deviations larger than PERR from
;		  the location of the observed line are suspect.
;		* if LSFWDT is not set, PERR constitutes a >lower limit< to WDT
;	werr	[INPUT] allowed excess variation in range of WDT
;		* default is 0.5*PERR[0]
;		* subtracted from min and added to max if LEEWAY requires it
;	ferr	[INPUT; default=0.1] allowed fractional error in FLX
;		* has an effect only if LEEWAY is set
;		* NOTE: for PERR,WERR,FERR, if size does not match that of
;		  POS,WDT,FLX respectively OR the number of components in
;		  IDSTR, only the first element is used
;       numer   [INPUT; default = 0] if set, ties specified via leeway, perr,
;                 werr, and ferr will be specified numerically rather than 
;                 symbolically
;       nozero  [INPUT; default = 0] if set, then fluxes are constained to 
;                 be positive. if set negative, then fluxes are constrained 
;                 to be negative 
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Nov98)
;	added keyword LSFWDT (VK; Mar99)
;	changed call to INITSTUFF to INICON (VK; 99May)
;       added EPITHET keyword (LL; Feb04) 
;       added NUMER keyword (LL; Jul05) 
;       BUGFIX: for multiplet component leeway width specification (LL; Jul 05) 
;       had used lower bound for position instead of width 
;       added NOZERO keyword (LL; Jul05) 
;	BUGFIX: for multiplet component leeway position specification
;	  was written with the wrong sign (VK; May06)
;-

;	usage
nid=n_tags(idstr)
if nid eq 0 then begin
  print,'Usage: id_to_fitpar,idstr,pars,ties,thaw,pos=pos,wdt=wdt,flx=flx,$'
  print,'       lsfwdt=lsfwdt,shackle=shackle,leeway=leeway,perr=perr,$'
  print,'       werr=werr,ferr=ferr'
  print,'  generate initial guesses for fit parameters for lines'
  return
endif

;	initialize
pos=[0.] & wdt=[0.] & flx=[0.] & ties='' & pars=[pos,wdt,flx]
atom=1 & rom=1 & inicon,atom=atom,roman=rom
atom=['X',atom] & rom=['',rom]

;	check input
idnam=tag_names(idstr) & ok='ok'
if idnam(0) ne 'WVL' then ok='Input structure in unknown form' else $
 if idnam(1) eq 'WVL_COMMENT' then ok='No data in structure'
if ok ne 'ok' then begin
  message,ok,/info & return
endif

;	disentangle input
owvl=abs(idstr.WVL) & nwvl=n_elements(owvl)	;Nwvl is the total # of lines
mwvl=0L & for i=0,nwvl-1 do mwvl=mwvl+n_elements(idstr.(i+1).WVL)
					;MWVL is the total # of IDs
;	output goes into..
pos=fltarr(mwvl) & wdt=pos & flx=pos & labl=strarr(nwvl)
pars=fltarr(3*mwvl)		;all the parameters
thaw=intarr(3*mwvl)+1		;all parameters thawed unless otherwise set
idx=lonarr(mwvl)		;matching the IDs to IDSTR components
k=0L
;	(the following figures out how many IDs in this component
for i=0L,nwvl-1L do begin
  idi=idstr.(i+1) & mw=n_elements(idi.WVL)
  for j=0L,mw-1L do idx(k+j)=i & k=k+mw
endfor		;IDX keeps track of which IDs belong to which line)

;	establish the allowed slippage
npe=n_elements(perr) & nwe=n_elements(werr) & nfe=n_elements(ferr)
pe=0*pos+0.1 & we=0.5*pe & fe=0*flx+0.1		;hard-coded defaults
if npe ge 1 then pe(*)=perr(0)			;override default
we=0.5*pe					;the "real" default
if nwe ge 1 then we(*)=werr(0)			;override default
if nfe ge 1 then fe(*)=ferr(0)			;override default
if npe eq mwvl then pe=perr			;one for each match
if nwe eq mwvl then we=werr			;one for each match
if nfe eq mwvl then fe=ferr			;one for each match
if npe eq nwvl or nwe eq nwvl or nfe eq nwvl then begin	;(one per component
  for i=0L,nwvl-1L do begin
    oo=where(idx eq i,moo) & if moo eq 0 then message,'BUG!'
    if npe eq nwvl then pe(oo)=perr(i)
    if nwe eq nwvl then we(oo)=werr(i)
    if nfe eq nwvl then fe(oo)=ferr(i)
  endfor
endif							;one per component)

;	instrumental width
wdtlsf=pe & if keyword_set(lsfwdt) then wdtlsf(*)=lsfwdt(0)


;	guess the parameter values
k=0L
for i=0,nwvl-1 do begin			;{for each observed wvl
  idi=idstr.(i+1) & ww=idi.wvl & mw=n_elements(ww)
  zz=idi.Z & jon=idi.ION & ff=idi.FLUX
  ;
  delP=owvl(i)-abs(ww)		;need this to compute possible WDT
  oi=where(idx eq i,moi) & if moi ne mw then message,'bug!'
  for j=0,mw-1 do begin			;{for each ID, correct for high orders
    if abs(delP(j)) gt 3*pe(oi(j)) then begin	;(here's a problem..
      deltaP=abs(delP(j)) & wwj=abs(ww(j)) & l=1L
      while deltaP gt 3*pe(oi(j)) do begin
	wwj=wwj*l & deltaP=owvl(i)-wwj & l=l+1L
      endwhile
      delP(j)=abs(deltaP) < 3*pe(oi(j))		;in case deltaP overshoots
    endif					;DELP too large)
  endfor				;J=0,MW-1}
  ww=owvl(i)-delP		;reset WW
  ;
  sig1=total(delP^2)/mw & sig1=sqrt(sig1)		;(2 ways of..
  sig2=total(delP^2)					;..getting..
  for j=1,mw-1 do sig2=sig2+total((ww(j-1)-ww)^2)	;..estimates of..
  for j=0,mw-1 do sig2=sig2/(mw-j) & sig2=sqrt(sig2)	;..expected..
  sig=sig1 & if sig lt sig2 then sig=sig2		;..widths)
  ;
  for j=0,mw-1 do begin
    k3=3L*k
    if sig lt wdtlsf(k) then sig=wdtlsf(k)
    pos(k)=ww(j) & wdt(k)=sig & flx(k)=ff(j)
    pars(k3:k3+3L-1L)=[pos(k),wdt(k),flx(k)]
    k=k+1L
  endfor
endfor					;I=0,NWVL-1}

;	and now to figure out the ties

;	see if we are allowed any leeway
pleek=0 & wleek=0 & fleek=0		;by default, no leeway
if keyword_set(leeway) then begin			;{LEEWAY
  szl=size(leeway) & nszl=n_elements(szl)
  if szl(nszl-2) eq 7 then begin		;(LEEWAY is a string
    i=strpos(strlowcase(leeway(0)),'po',0) & if i ge 0 then pleek=1
    i=strpos(strlowcase(leeway(0)),'wi',0) & if i ge 0 then wleek=1
    i=strpos(strlowcase(leeway(0)),'fl',0) & if i ge 0 then fleek=1
  endif else begin				;)(not a string
    ;	right now, for want of anything better to do, set all
    pleek=1 & wleek=1 & fleek=1
  endelse					;leeway type)
  if keyword_set(numer) then begin 
    lpk = strcompress(string(pos),/remove_all)   
    lwk = strcompress(string(wdt),/remove_all)
    lfk = strcompress(string(flx),/remove_all) 
  endif 
endif							;LEEWAY}

;	or the parameters must be shackled
pshaq=0 & wshaq=0 & fshaq=0		;by default, no shackles
if keyword_set(shackle) then begin			;{SHACKLE
  szl=size(shackle) & nszl=n_elements(szl)
  if szl(nszl-2) eq 7 then begin		;(SHACKLE is a string
    i=strpos(strlowcase(shackle(0)),'po',0) & if i ge 0 then pshaq=1
    i=strpos(strlowcase(shackle(0)),'wi',0) & if i ge 0 then wshaq=1
    i=strpos(strlowcase(shackle(0)),'fl',0) & if i ge 0 then fshaq=1
  endif else begin				;)(not a string
    ;	right now, for want of anything better to do, set all
    pshaq=1 & wshaq=1 & fshaq=1
  endelse					;shackle type)
endif							;SHACKLE}

;	ok, NOW step through the parameters and figure out the ties
k=0L
wmin=min(wdt,max=wmax)
for i=0,nwvl-1 do begin			;{for each component in IDSTR
  k3=3L*k & ow=where(idx eq i,mow)
  pp=pos(ow) & ww=wdt(ow) & ff=flx(ow)
  tt='' & dw=owvl(i)-ww & fmx=max(ff,imx)
  ;	set the leeway
  kk=k3+3L*imx
  if keyword_set(numer) then begin 
    pkk=lpk(k3/3L) & wkk=lwk(k3/3L) & fkk=lfk(k3/3L) 
  endif 
  pk='A'+strtrim(kk,2) & wk='A'+strtrim(kk+1,2) & fk='A'+strtrim(kk+2,2)

  if pleek ne 0 then begin
    if not keyword_set(numer) then begin 
      tt=pk+'= (('+pk+' < ('+pk+'+'+strtrim(abs(pe(k+imx)),2)+')) > ('+pk+$
      '-'+strtrim(abs(pe(k+imx)),2)+'))'
    endif else begin 
      tt=pk+'= (('+pk+' < ('+pkk+'+'+strtrim(abs(pe(k+imx)),2)+')) > ('+pkk+$
      '-'+strtrim(abs(pe(k+imx)),2)+'))'
    endelse
  endif & if tt ne '' then ties=[ties,tt] & tt=''

  if wleek ne 0 then begin
    if not keyword_set(numer) then begin
      tt=wk+'= (('+wk+' < ('+wk+'+'+strtrim(wmax+abs(we(k+imx)),2)+$
	 ')) > ('+wk+'-'+strtrim(wmin+abs(we(k+imx)),2)+'))'
      tt=wk+'= (('+wk+' < ('+strtrim(wmax+abs(we(k+imx)),2)+$
	 ')) > ('+strtrim(wmin-abs(we(k+imx)),2)+'))'
    endif else begin 
      tt=wk+'= (('+wk+' < ('+wkk+'+'+strtrim(wmax+abs(we(k+imx)),2)+$
	 ')) > ('+wkk+'-'+strtrim(wmin+abs(we(k+imx)),2)+'))'
      tt=wk+'= (('+wk+' < ('+strtrim(wmax+abs(we(k+imx)),2)+$
	 ')) > ('+strtrim(wmin-abs(we(k+imx)),2)+'))'
    endelse
  endif & if tt ne '' then ties=[ties,tt] & tt=''
  if fleek ne 0 then begin
    if not keyword_set(numer) then begin     
      tt=fk+'= (('+fk+' < ('+fk+'+'+strtrim(abs(fe(k+imx)),2)+')) > ('+fk+$
	'-'+strtrim(abs(fe(k+imx)),2)+'))'
    endif else begin 
      tt=fk+'= (('+fkk+' < ('+fkk+'+'+strtrim(abs(fe(k+imx)),2)+')) > ('+fkk+$
	'-'+strtrim(abs(fe(k+imx)),2)+'))'
    endelse
  endif & if tt ne '' then ties=[ties,tt] & tt=''
  if mow gt 0 then begin		;(there are multiple IDs
    for j=0,mow-1 do begin
      kj=k3+3L*j
      if keyword_set(numer) then begin 
         pkkj=lpk(k3/3L) & wkkj=lwk(k3/3L) & fkkj=lfk(k3/3L)
      endif 
        pkj='A'+strtrim(kj,2) & wkj='A'+strtrim(kj+1,2) & fkj='A'+strtrim(kj+2,2)
      if j ne imx then begin		;(haven't already dealt with above
	;	(leeway
	if pleek ne 0 and pshaq eq 0 then begin
	  dup=pp(imx)-pp(j)-sqrt(pe(k)^2+pe(imx)^2)
	  ddn=pp(imx)-pp(j)+sqrt(pe(k)^2+pe(imx)^2)
	  if keyword_set(numer) then begin 
            tt=pkj+'= (('+pkj+' < ('+pkk+'-('+strtrim(dup,2)+'))) > ('+$
	  	pkk+'-('+strtrim(ddn,2)+')))'
          endif else begin 
            tt=pkj+'= (('+pkj+' < ('+pk+'-('+strtrim(dup,2)+'))) > ('+$
                 pk+'-('+strtrim(ddn,2)+')))'
          endelse
	endif & if tt ne '' then ties=[ties,tt] & tt=''
	if wleek ne 0 and wshaq eq 0 then begin
	  wup=wmax+abs(we(j)) & wdn=wmin-abs(we(j))

          if keyword_set(numer) then begin 
            tt=wkj+'= (('+wkj+' < ('+strtrim(wup,2)+')) > ('+strtrim(wdn,2)+'))'
          endif else begin 
            tt=wkj+'= (('+wkj+' < ('+strtrim(wup,2)+')) > ('+strtrim(wdn,2)+'))'
          endelse
	endif & if tt ne '' then ties=[ties,tt] & tt=''
	if fleek ne 0 and fshaq eq 0 then begin
	  frup=abs(ff(imx)-fe(imx))
	  if frup gt 0 then frup=abs(ff(j)+fe(j))/frup else frup=1.+fe(j)
	  frdn=abs(ff(imx)+fe(imx))
	  if frdn gt 0 then frdn=abs(ff(j)-fe(j))/frdn else frdn=1.-fe(j)
          if keyword_set(numer) then begin 
              tt=fkj+'= (('+fkkj+' < ('+fkk+'*('+strtrim(frup,2)+'))) > ('+$
	      fkk+'*('+strtrim(frdn,2)+')))'
          endif else begin 
              tt=fkj+'= (('+fkj+' < ('+fk+'*('+strtrim(frup,2)+'))) > ('+$
	      fk+'*('+strtrim(frdn,2)+')))'
          endelse         
	endif & if tt ne '' then ties=[ties,tt] & tt=''
	;	leeway)(shackle
	if pshaq ne 0 then begin
	  tmp=pp(imx)-pp(j)
	  thaw(kj)=0
	  tt=pkj+' = '+pk+'- ('+strtrim(tmp,2)+')'
	endif & if tt ne '' then ties=[ties,tt] & tt=''
	if wshaq ne 0 then begin
	  tmp=ww(j) & if tmp le 0 then tmp=wmax
	  thaw(kj+1L)=0
	  tt=wkj+' = '+strtrim(tmp,2)
	endif & if tt ne '' then ties=[ties,tt] & tt=''
	if fshaq ne 0 then begin
	  tmp=ff(imx) & if tmp ne 0 then tmp=ff(j)/tmp else tmp=0
	  thaw(kj+2L)=0
	  tt=fkj+' = '+fk+'* ('+strtrim(tmp,2)+')'
	endif & if tt ne '' then ties=[ties,tt] & tt=''
	;	shackle)
      endif				;J.NE.IMX)
    endfor
  endif 				;MOW>0)
  k=k+mow
endfor					;I=0,NWVL-1}

if keyword_set(nozero) then begin 
       kk = 'A'+strtrim(lindgen(n_elements(idx))*3L+2L,2)
   if nozero lt 0  then begin 
       tt = kk+'='+kk+' < '+'0' 
   endif 
   if nozero gt 0  then begin 
       tt = kk+' = '+kk+' > '+'0' 
   endif
       ties = [ties,tt]    
endif 

ntie=n_elements(ties)
if ntie gt 1 then begin
  otie=lindgen(ntie-1L)+1L & ties=ties(otie)
endif


;	prepare epithet keyword
if arg_present(epithet) then begin 
   epithet = 'jnk'  
   for j = 0, nwvl-1 do begin 
      idz = idstr.(j+1) & nw = n_elements(idz.wvl)
      for qq = 0, nw -1 do begin 
         epithet=[epithet,atom(idz.z[qq])+rom(idz.ion[qq])]
      endfor
   endfor
   epithet = epithet[1:*]
endif

return
end
