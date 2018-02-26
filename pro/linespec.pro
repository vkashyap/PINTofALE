function linespec,atom,wmin=wmin,wmax=wmax,nbin=nbin,wlog=wlog,ww=ww,$
	ikev=ikev,fstr=fstr, _extra=e
;+
;function	linespec
;	returns a spectrum of atomic lines.  output will be in units
;	of [erg,ph]/s/cm^2/[keV,Ang]
;	
;	if IKEV is set, WMIN,WMAX,WW must be input in [keV] and the
;	output will have units .../keV
;
;	if NOPH is set, output will have units of erg/...
;
;	if EFFAR and WVLAR are not set, output will have units of .../cm^2
;
;syntax
;	spec=linespec(atom,wmin=wmin,wmax=wmax,nbin=nbin,/wlog,ww=ww,/ikev,$
;	fstr=fstr,pres=pres,logP=logP,n_e=n_e,desig=desig,econf=econf,$
;	/allah,dbdir=dbdir,dem=dem,abund=abund,effar=effar,wvlar=wvlar,$
;	/noph,chifil=chifil,chidir=chidir,eqfile=eqfile,verbose=verbose)
;
;parameters
;	atom	[INPUT; default: all] element(s) whose line intensities are
;		to be read.
;		* can specify ionic state (e.g., 'FeXX')
;		* may be an array (e.g., ['He', 'c 5', 'N V', 's_8'])
;
;keywords
;	wmin	[I/O] minimum value of wavelength range
;		* default is 1.23985A (10 keV)
;		* overwritten if WW is set
;	wmax	[I/O] maximum value of wavelength range
;		* default is 900A
;		* overwritten if WW is set
;	nbin	[I/O] number of bins in spectrum
;		* default is 100
;		* overwritten if WW is set
;	wlog	[I/O] if set, binning will be logarithmic
;	ww	[I/O] if set as an array on input, assumed to be the
;		bin boundaries --
;		* NBIN will be reset to n_elements(WW)-1 on output
;		* WMIN will be reset to min([WW(0),WW(nbin)]) on output
;		* WMAX will be reset to max([WW(0),WW(nbin)]) on output
;	ikev	[INPUT] if set, WMIN,WMAX,WW are all assumed to be in keV
;		* it is also caught at this stage and is not passed down to
;		  LINEFLX
;	fstr	[I/O] if defined on input, will return all the spectral
;		lines info LINE_INT,LOGT,WVL,Z,ION,DESIG,ECONF in that order
;		in a structure.
;		* ion balance is included in LINE_INT, but abundances
;		  are not.
;		* if defined as a structure on input, RD_LINE and FOLD_IONEQ
;		  will not be called, so use with caution!
;	_extra	[INPUT] allows specifying defined keywords to subroutines
;		called by this program
;		* RD_LINE: PRES,LOGP,N_E,DESIG,ECONF,ALLAH,DBDIR,VERBOSE
;		* FOLD_IONEQ: CHIFIL,VERBOSE
;		* RD_IONEQ: CHIDIR,EQFILE
;		* LINEFLX: DEM,ABUND,EFFAR,WVLAR,NOPH
;
;subroutines
;	RD_LINE [KILROY, SYMB2ZION [LAT2ARAB]]
;	FOLD_IONEQ [WHEE, GETLOCMAX, RD_IONEQ [READ_IONEQ (CHIANTI)]]
;	LINEFLX [WHEE]
;	HASTOGRAM [KILROY]
;
;history
;	vinay kashyap (Jan97)
;	included ion balance into FSTR (VK;Feb97)
;	bug -- crashes when FSTR not specified (VK;Mar97)
;	speeded up spectrum accumulation by adding call to HASTOGRAM (VK;Jan98)
;	modified ion balance calcs (VK; 99May)
;	changed keyword KEV to IKEV; DESIG and ECONF must be explicitly set
;	  to 0 to avoid reading them in in RD_LINE; streamlined IKEV behavior
;	  (VK; JanMMI)
;-

;	usage
if n_params(0) eq 0 then begin
  print,'Usage: spec=linespec(atom,wmin=wmin,wmax=wmax,nbin=nbin,/wlog,$'
  print,'       ww=ww,/ikev,fstr=fstr)'
  print,'  returns a spectrum generated from atomic lines.  also accepts'
  print,'  keywords PRES,LOGP,N_E,DESIG,ECONF,ALLAH,DBDIR,VERBOSE (RD_LINE);'
  print,'  DEM,ABUND,EFFAR,WVLAR,NOPH (LINEFLX); CHIDIR,CHIFIL (FOLD_IONEQ);'
  print,'  EQFILE (RD_IONEQ)'
endif

;	check keywords
nww=n_elements(ww)
if nww gt 1 then begin
  nbin=nww-1 & wmin=min([ww(0),ww(nbin)],max=wmax)	;reset
endif else begin
  if not keyword_set(wmin) then wmin=1.23985		;if wmin not set
  if not keyword_set(wmax) then wmax=900.		;if wmax not set
  if not keyword_set(nbin) then nbin=100L		;if nbin not set
  if wmin le 0 then wmin=1.23985			;something wrong?
  if wmax le 0 then wmax=900.				;something wrong?
  if wmin gt wmax then begin
    ;	exchange WMIN & WMAX
    wmin=wmin-wmax & wmax=wmin+wmax & wmin=wmax-wmin
    ;	not doing X=WMIN & WMIN=WMAX & WMIN=X because WMAX will not be
    ;	so much larger than WMIN that I don't expect floating point errors
    ;	will be significant.
  endif
  if not keyword_set(wlog) then begin
    dw=float((wmax-wmin))/float(nbin)
    ww=findgen(nbin+1)*dw+wmin				;linear binning
  endif else begin
    w0=alog10(wmin) & w1=alog10(wmax) & dw=(w1-w0)/nbin
    ww=findgen(nbin+1)*dw+w0 & ww=10.^(ww)		;log binning
  endelse
endelse
w0=wmin & w1=wmax
if keyword_set(ikev) then begin
  w0=12.3985/wmax & w1=12.3985/wmin
  ww=12.3985/ww & ow=sort(ww) & ww=ww(ow)
endif

;	initialize
spec=fltarr(nbin)

;	read in the lines
szf=size(fstr)
if szf(n_elements(szf)-2) ne 8 then begin	;FSTR not structure
  ;	read in database
  ff=rd_line(atom,wrange=[w0,w1],wvl=wvl,logT=logT,Z=Z,ion=ion,jon=jon,$
	fstr=fstr,/desig,/econf, _extra=e)
  ;	include ion balance
  ff=fold_ioneq(ff,z,jon,logt=logt, _extra=e)
  if n_tags(fstr) gt 0 then fstr.(0)=ff
  wvl=abs(wvl)		;CHIANTI theoretical wavelengths are -ve
endif else begin				;FSTR be structure
  ff=fstr.LINE_INT
  logt=fstr.LOGT
  wvl=abs(fstr.WVL)
  z=fstr.Z
  ion=fstr.ION
endelse
if z(0) eq 0 then return,spec			;no lines found?

;	compute line fluxes (comes out as [ph/s/cm^2] or [ph/s])
f=lineflx(ff,logt,wvl,z, _extra=e)

;	add to spectrum
if keyword_set(ikev) then begin
  wave=12.3985/wvl
  ;also make sure outputs go back in correct order
  ww=12.3985/ww & ow=sort(ww) & ww=ww(ow)
  wmin=12.3985/wmin & wmax=12.3985/wmax
endif else wave=wvl
dw=ww(1:*)-ww
spec=hastogram(wave,ww,wts=f)/dw

;;for iw=0,nbin-1 do begin
;;  whee,spoke=spoke,/moveit			;patience...
;;  w0=ww(iw) & w1=ww(iw+1) & dw=abs(w1-w0)
;;  ;in the general case, where the binning may not even be logarithmic, 
;;  ;the easy solution (IW=(WAVE-WMIN)/NBIN) won't work
;;  oo=where(abs(wave) ge w0 and abs(wave) lt w1,moo)
;;  if moo gt 0 then spec(iw)=total(f(oo))/dw
;;  ;;for i=0,moo-1 do begin
;;  ;;  ;doing it this way instead of SPEC(IW)=TOTAL(F(OO))/DW allows future
;;  ;;  ;generalization to complicated (i.e., non-delta function!) line profiles
;;  ;;  spec(iw)=spec(iw)+f(oo(i))/dw
;;  ;;endfor
;;endfor

return,spec
end
