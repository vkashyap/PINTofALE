function lsd,wrange,fluxes,wvls,elem=elem,edens=edens,DEM=DEM,tlog=tlog,$
	dbdir=dbdir,ceiling=ceiling,flor=flor,ratmax=ratmax,flxmax=flxmax,$
	outZ=outZ,outIon=outIon,outIDX=outIDX, _extra=e
;+
;function	lsd
;	returns list of density sensitive lines in given wavelength range
;
;syntax
;	ld=lsd(wrange,fluxes,wvls,elem=elem,edens=edens,DEM=DEM,tlog=tlog,$
;	dbdir=dbdir,ceiling=ceiling,flor=flor,ratmax=ratmax,flxmax=flxmax,$
;	outZ=outZ,outIon=outIon,outIDX=outIDX, tol=tol,chifil=chifil,$
;	chidir=chidir,eqfile=eqfile,/temp,abund=abund,/noph,effar=effar,$
;	wvlar=wvlar,/ikev)
;
;parameters
;	wrange	[INPUT; required] wavelength range in which to search for
;		density sensitive lines
;		* if only one element is given, set the range to be 10% of
;		  it on either side
;		* if the single element is a string, apply RD_LIST rules
;		  to decipher it ("WVL +- dW", "WMIN,WMAX", "WMIN-WMAX"
;		  are all allowed)
;	fluxes	[OUTPUT] 2D array of predicted fluxes for lines that are
;		found to be density sensitive, FLUXES(EDENS,WVLS)
;	wvls	[OUTPUT] array of wavelengths of lines that are found to
;		be density sensitive
;
;keywords
;	elem	[INPUT] if set, confine attention to just those elements and
;		ionic species specified.
;		* default is to look at all available lines
;	edens	[I/O] array of electron densities [cm^-3] at which to
;		check for density sensitivity.
;		* default is [1e8,1e14]
;		* if not specified, insufficiently specified, or otherwise
;		  meaningless, uses default
;	DEM	[INPUT] input array of differential emission measures [cm^-5]
;		* default is to use 1e12 [cm^-5] at specified TLOG
;	tlog	[INPUT] log(temperatures) at which DEM is defined
;		* if not specified, set to 6.0, unless DEM has more steps,
;		  in which case interpolated to the 4..8 range
;	dbdir	[INPUT] line database directory to use to search for the
;		density sensitive lines
;		* default is $CHIANTI
;	ceiling	[INPUT] maximum relative variation the fluxes of a given
;		line must have before it is accepted as density sensitive
;		* default is 2
;		* if 0 < CEILING < 1, it is assumed that the *reciprocal*
;		  of CEILING is given (e.g., CEILING=0.5 is translated to 2)
;		* if CEILING < 0, assumed to be percentages (e.g., CEILING=-1
;		  is translated to 1/0.01 == x100)
;	flor	[INPUT] if set, ignores any line which has a peak flux
;		a factor FLOR below the flux of the strongest of the
;		lines filtered through CEILING.
;		* default is 0
;		* if FLOR > 1, assumed that the reciprocal of FLOR is given,
;		  i.e., FLOR=100 translates to 0.01
;		* if FLOR < 0, taken to be the absolute flux threshold, in
;		  whatever units it is being calculated in
;		* NOTE on difference between CEILING and FLOR:-
;		  CEILING works on ratios, and FLOR works on the intensities
;	ratmax	[OUTPUT] the maximum deviations found for each of the
;		density-sensitive lines found
;	flxmax	[OUTPUT] the peak fluxes for each of the density sensitive
;		line found
;	outZ	[OUTPUT] atomic numbers
;	outIon	[OUTPUT] atomic numbers
;	outIDX	[OUTPUT] indices which passed the test
;		* be careful in how you use this.  a common mistake is to
;		  assume that it matches with whatever line emissivity
;		  database that has been read in.  that might very well be
;		  true, but is definitely not guaranteed.
;		
;	_extra	[INPUT] use this to pass defined keywords to subroutines:-
;		FOLD_IONEQ: CHIFIL
;		RD_IONEQ: CHIDIR, EQFILE
;		LINEFLX: TEMP, ABUND, NOPH, EFFAR, WVLAR, IKEV
;		ARRAYEQ: TOL
;
;example usage
;	ls=lsd('10-120',dbdir='$CHIANTI',edens=[1e9,1e12],ceiling=2,ratmax=rmx)
;	for i=0,n_elements(ls)-1 do print,ls(i),1./rmx(i),i
;	;
;	eden=10.^(findgen(6)+8)
;	ls=lsd('113.35?0.1',flx,elem='Fe20',dbdir='$CHIANTI',edens=eden)
;	plot,eden,flx
;
;subroutines
;	MK_DEM
;	RD_LINE
;	FOLD_IONEQ
;	    RD_IONEQ
;	        READ_IONEQ (CHIANTI subroutine)
;	LINEFLX
;	ARRAYEQ
;	KABOOM
;	INICON
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Nov98)
;	changed default EDENS to [1e8,1e14]; added keywords DEM,TLOG,OUTZ,
;	  OUTION; added call to MK_DEM; changed keyword FLOOR to FLOR; changed
;	  behavior of FLOR<0 (VK; Dec98)
;	modified ion balance calcs (VK; 99May)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	added keyword outIDX (VK; JanMMI)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;-

;	usage
szw=size(wrange) & nszw=n_elements(szw) & nw=szw(nszw-1L)
if szw(nszw-2) eq 0 or nw gt 2 then begin
  print,'Usage: ld=lsd(wrange,fluxes,wvls,elem=elem,edens=edens,DEM=DEM,$'
  print,'       tlog=tlog,dbdir=dbdir,ceiling=ceiling,flor=flor,ratmax=rmx,$'
  print,'       flxmax=fmx,outZ=outZ,outIon=outIon,outIDX=outIDX, tol=tol,$'
  print,'       chifil=chifil,chidir=chidir,eqfile=eqfile,/temp,abund=abund,$'
  print,'       /noph,effar=effar,wvlar=wvlar,/ikev)'
  print,'  return list of density sensitive lines'
  return,''
endif

;	check input
ok='ok'
if nw gt 2 then ok='cannot understand WRANGE -- too many elements?'
if szw(nszw-2) eq 7 then begin			;(string input
  if nw gt 1 then message,'using only first element of WRANGE',/info
  c=strtrim(wrange(0),2)
  j=strpos(c,'[',0) & if j ge 0 then c=strmid(c,j+1,strlen(c)-j-1)
  j=strpos(c,']',0) & if j ge 0 then c=strmid(c,0,j)
  j=strpos(c,'(',0) & if j ge 0 then c=strmid(c,j+1,strlen(c)-j-1)
  j=strpos(c,')',0) & if j ge 0 then c=strmid(c,0,j)
  j=strpos(c,'{',0) & if j ge 0 then c=strmid(c,j+1,strlen(c)-j-1)
  j=strpos(c,'}',0) & if j ge 0 then c=strmid(c,0,j)
  j=0						; W (exact match)
  if strpos(c,'+-',0) ge 0 then j=1		; W +- dW
  if strpos(c,'-',0) ge 0 and j ne 1 then j=2	; Wmin - Wmax
  if strpos(c,',',0) ge 0 then j=3		; Wmin , Wmax
  if strpos(c,'?',0) ge 0 then j=4		; W ? dW (same as W +- dW)
  case j of					;{for each case
    1: begin
      cc=str_sep(c,'+-') & ncc=n_elements(cc)
      w=float(cc(0)) & if ncc eq 1 then dw=0.1*w else dw=float(cc(1))
      w1=w-dw & w2=w+dw
    end
    2: begin
      cc=str_sep(strtrim(c,2),'-') & ncc=n_elements(cc)
      w1=float(cc(0)) & if ncc eq 1 then w2=w1+1. else w2=float(cc(1))
    end
    3: begin
      cc=str_sep(strtrim(c,2),',') & ncc=n_elements(cc)
      w1=float(cc(0)) & if ncc eq 1 then w2=w1+1. else w2=float(cc(1))
    end
    4: begin
      cc=str_sep(strtrim(c,2),'?') & ncc=n_elements(cc)
      w=float(cc(0)) & if ncc eq 1 then dw=0.1*w else dw=float(cc(1))
      w1=w-dw & w2=w+dw
    end
    else: begin
      w=float(c) & w1=0.9*w & w2=1.1*w
    endelse
  endcase					;J}
endif						;WRANGE is a string)
if szw(nszw-2) le 5 and szw(nszw-2) ge 2 then begin	;(numeric input
  if nw eq 1 then begin
    w=0.+wrange(0) & w1=0.9*w & w2=1.1*w
  endif else begin
    w1=0.+wrange(0) & w2=0.+wrange(1)
  endelse
endif							;WRANGE=[Wmin,Wmax])
if w1 gt w2 then begin
  w=w1 & w1=w2 & w2=w			;switch if necessary
endif
if szw(nszw-2) gt 7 or szw(nszw-2) lt 2 then ok='WRANGE: Unknown data type'
if ok ne 'ok' then begin
  message,ok,/info & return,''
endif

;	check keywords
if not keyword_set(elem) then elem=0
;
nden=n_elements(edens)
if nden lt 2 then begin
  edens=[1e8,1e14] & nden=2L
endif
;
if not keyword_set(dbdir) then dbdir='$CHIANTI'
;
cthr=0.5
if keyword_set(ceiling) then begin
  if ceiling(0) gt 1 then cthr=1./ceiling(0)
  if ceiling(0) gt 0 and ceiling(0) le 1 then cthr=ceiling(0)
  if ceiling(0) lt 0 then cthr=abs(ceiling(0))/100.
endif
if cthr gt 1 then cthr=1./cthr
;
fthr=1e-6
if keyword_set(flor) then begin
  if flor(0) gt 0 then fthr=flor(0)
  if flor(0) lt 0 then fthr=1e-6
endif
;if fthr gt 1 then fthr=1./fthr

;	need some initialization
atom=1 & rom=1 & inicon,atom=atom,roman=rom
nlin=0L

;	extract the emissivities, fold in ion balance, get line fluxes
for i=0L,nden-1L do begin		;{for each density
  ff=rd_line(elem,n_e=edens(i),wrange=[w1,w2],dbdir=dbdir,$
	wvl=wvl,logT=logT,Z=Z,ion=ion,jon=jon, _extra=e)
  if ff(0) le -1 then goto,skip		;{no lines read..
  wvl=abs(wvl)
  ff=fold_ioneq(ff,Z,jon,logT=logT, _extra=e)
  if keyword_set(DEM) then D_E_M=DEM
  if n_elements(D_E_M) ne n_elements(logT) then begin
    message,'interpolating DEM to fit the temperature grid',/info
    D_E_M=mk_dem('interpolate',logT=logT,indem=DEM,pardem=tlog)
  endif
  fx=lineflx(ff,logT,wvl,Z,DEM=D_E_M, _extra=e)

  if not keyword_set(ii) then begin	;(check consistency at EDENS(i)
    ii=1			;first pass
    nlin=n_elements(fx)
    owvl=wvl & outZ=Z & oution=ion		;save them..
    fmax=fx
  endif else begin
    ;	subsequent passes
    nln=n_elements(fx) & aeq=arrayeq(wvl,owvl,_extra=e)
    if nln ne nlin or aeq ne 1 then begin		;(problem..
      message,'holes exist in the line database',/info
      kaboom,/flash,/beep
      ffx=0*fx & k=0L
      if nln gt nlin then begin		;(more lines than before..
	for j=0L,nlin-1L do begin
	  if abs(owvl(j)-wvl(k)) lt 1e-3 and outZ(j) eq Z(k) and oution(j) eq ion(k) then begin
	    ffx(j)=fx(k) & k=k+1L
	  endif
	endfor
      endif else begin			;)(fewer lines than before..
	for j=0L,nln-1L do begin
	  if abs(owvl(k)-wvl(j)) lt 1e-3 and outZ(k) eq Z(j) and oution(k) eq ion(j) then begin
	    ffx(k)=fx(j) & k=k+1L
	  endif
	endfor
      endelse				;NLN v/s NLIN)
      fx=ffx
    endif						;..problem)
    fmax = fmax > fx	;locate the maxima
  endelse				;EDENS consistency check)

  ;	form flux array..
  if not is_keyword_set(flux) then flux=fx else flux=[flux,fx]
  if not is_keyword_set(DenE) then DenE=edens(i) else DenE=[DenE,edens(i)]

  skip:	;get here directly if no lines were read..}
endfor					;I=0,NDEN-1}

;	simple error checks
if nlin eq 0 then return,''				;nothing was found
nfx=n_elements(flux) & if nfx lt nlin then return,''	;nothing was found
nd=n_elements(DenE) & if nd lt 2 then return,''		;nothing was found
ok=where(fmax gt 0,mok) & if mok eq 0 then return,''	;nothing varied
nd2=nfx/nlin & if nd ne nd2 then message,'bug'

;	now find the density sensitive lines
frat=fltarr(nlin)+1.
for i=0L,nd-1L do begin		;{check all lines at each density
  k1=i*nlin & k2=k1+nlin-1L
  frat1=flux(k1:k2) & frat(ok)=frat1(ok)/fmax(ok) < (frat(ok))
endfor				;I=0,ND-1}
;
ot=where(frat le cthr,mot)
if mot eq 0 then begin
  message,'no density sensitive lines were found',/info
  return,''
endif
if keyword_set(flor) then begin
  fmaxmax=max(fmax(ot))
  if fthr lt 1 then ot=where(frat lt cthr and fmax ge fthr*fmaxmax,mot) else $
    ot=where(frat lt cthr and fmax ge abs(flor(0)),mot)
  if mot eq 0 then message,'no lines were found for this flux threshold',/info
endif

;	the output
sep='	'
ls=strarr(mot)
for i=0L,mot-1L do ls(i)=atom(outZ(ot(i))-1)+' '+rom(oution(ot(i))-1)+sep+$
	strtrim(string(owvl(ot(i)),'(g12.5)'),2)+sep+dbdir
ratmax=frat(ot)
flxmax=fmax(ot)
outZ=outZ(ot)
outIon=outIon(ot)
outIDX=ot

np=n_params()
if np gt 1 then begin
  fluxes=fltarr(nd,mot)
  for i=0L,nd-1L do begin
    k1=i*nlin & k2=k1+nlin-1L
    fx=flux(k1:k2) & fluxes(i,*)=fx(ot)
  endfor
endif
if np gt 2 then wvls=owvl(ot)

return,ls
end
