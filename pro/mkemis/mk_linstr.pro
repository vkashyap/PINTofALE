function mk_linstr,emis,logT=logT,wvl=wvl,Z=Z,ion=ion,jon=jon,$
	desig=desig,econf=econf,src=src,verbose=verbose, _extra=e
;+
;function	mk_linstr
;	create a standard emissivity structure by putting together
;	the various inputs.  this emissivity structure can then be
;	used in PINTofALE exactly as though it had been read in via
;	RD_LINE().
;
;syntax
;	linstr=mk_linstr(emis,logT=logT,wvl=wvl,Z=Z,ion=ion,jon=jon,$
;	desig=desig,econf=econf,src=src,verbose=verbose)
;
;parameters
;	emis	[INPUT; required] emissivities
;		* if 2D, assumed to be an array of size (NLOGT,NWVL)
;		* if 1D, makes a semi-intelligent guess as to whether
;		  it is an array corresponding to LOGT or WVL
;		  -- first checks with LOGT and if size matches, that's it
;		  -- else checks with WVL and
;		  -- if that doesn't match then quits with an error
;
;keywords
;	logT	[INPUT] log10(Temperatures [K])
;		* default is a regular grid going from 4..8
;		* if any of the values are > 100, then it is assumed
;		  that the input is in [K] and the log10 is taken
;		* if size does not match first dimension of 2D EMIS,
;		  then ignored and default is used
;	wvl	[INPUT] wavelengths [Angstroms]
;		* if size does not match 2nd dimension of EMIS, then
;		  quits with an error
;	Z	[INPUT] atomic numbers of elements generating given WVL
;		* missing values will be set to 1
;	ion	[INPUT] ionic state of element generating given WVL
;		* missing values are set to 1
;	jon	[INPUT] the ionic state that matters as far as ion-balance
;		is concerned
;		* missing values are set to ION
;	src	[INPUT] reference to source of line info
;		* e.g., 1=SPEX, 2=CHIANTI, 3=APED
;		* if not set, assumed to be -99
;		* if given but in insufficient numbers, the value of the
;		  first element is extrapolated out
;	desig	[INPUT] level designations for lower & upper level
;	econf	[INPUT] e-configurations for lower & upper level
;		* DESIG and ECONF are expected to be array of size (2,NWVL)
;		  and are filled out with blank strings if parts are missing
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] ignore -- here only to prevent crashing the program
;
;history
;	vinay kashyap (Jun02)
;-

;	usage
ok='ok' & np=n_params()
linstr=create_struct('LINE_INT',[-1],'LOGT',0.,'WVL',0.,'Z',0,'ION',0,$
	'DESIG','','CONFIG','','SRC',0,'JON',0)
nem=n_elements(emis) & szem=size(emis) & nszem=n_elements(szem)
if np eq 0 then ok='Insufficient parameters' else $
 if nem eq 0 then ok='EMIS is undefined' else $
  if szem[0] gt 2 then ok='EMIS can at most be 2D'
if ok ne 'ok' then begin
  print,'Usage: linstr=mk_linstr(emis,logT=logT,wvl=wvl,Z=Z,ion=ion,jon=jon,$'
  print,'       desig=desig,econf=econf,src=src,verbose=verbose)'
  if np ne 0 then message,ok,/info
  return,linstr
endif

;	check inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
nt=n_elements(logT) & nw=n_elements(wvl) & nz=n_elements(Z)
ni=n_elements(ion) & nj=n_elements(jon) & ns=n_elements(src)
szd=size(desig) & nszd=n_elements(szd) & ndes=n_elements(desig)
sze=size(econf) & nsze=n_elements(sze) & necf=n_elements(econf)

;	EMIS, LOGT, WVL
if szem[0] eq 2 then begin	;(EMIS is 2D
  nemt=szem[1] & nemw=szem[2]
  tlog=findgen(nemt)*(4./(nemt-1L))+4.
  if nt ne nemt then begin	;(LOGT size does not match EMIS
    if nt gt 0 and vv gt 0 then message,$
	'LOGT incompatible with EMIS; using the default',/info
  endif else tlog=[logT*1.0]	;NT v NEMT)
  if nw ne nemw then begin
    if nw gt 0 then begin
      if vv gt 0 then message,'WVL incompatible with EMIS; exiting',/info
      return,linstr
    endif
  endif else wave=[wvl*1.0]
  line_int=emis
endif else begin		;2D)(EMIS is 1D
  if nem eq nt then begin	;(number of elements matches size of LOGT
    line_int=reform(emis,nT,1)
    tlog=[logT]
    if nw eq 1 then wave=wvl[0] else begin
      if nw eq 0 then begin
	if vv gt 0 then message,'WVL is missing; assuming 0.0',/info
	wave=[0.0]
      endif else begin
	if vv gt 0 then message,$
		'WVL array incompatible with EMIS; using only first element',$
		/info
	wave=[wvl[0]*1.0]
      endelse
    endelse
  endif else begin		;LOGT)(number doesn't match LOGT)
    if nem eq nw then begin	;(#elements matches size of WVL
      line_int=reform(emis,1,nw)
      if nt eq 1 then tlog=logT[0] else begin
	if nt eq 0 then begin
	  if nt eq 0 then begin
	    if vv gt 0 then message,'LOGT is missing; setting to 4.0',/info
	    tlog=[4.0]
	  endif
	endif else begin
	  if vv gt 0 then message,$
		'LOGT array incompatible with EMIS; using only first element',$
		/info
	  tlog=[logT[0]*1.0]
	endelse
      endelse
    endif else begin		;WVL)(# doesn't match WVL
      message,'neither LOGT nor WVL matches EMIS; exiting',/info
      return,linstr
    endelse			;not WVL)
  endelse			;not LOGT)
endelse				;1D)
nt=n_elements(tlog) & nw=n_elements(wave)

;	Z
zz=intarr(nw)+1 & if nz ge nw then zz=Z[0L:nw-1L] else zz[0L:nz-1L]=Z[*]
;	ION
ii=intarr(nw)+1 & if ni ge nw then ii=ion[0L:nw-1L] else ii[0L:ni-1L]=ion[*]
;	JON
jj=ii & if nj ge nw then jj=jon[0L:nw-1L] else jj[0L:nj-1L]=jon[*]
;	SRC
ss=intarr(nw)-99 & if ns ge nw then ss=src[0L:nw-1L] else ss[0L:ns-1L]=src[*]
;	DESIG
dd=strarr(2,nw)
if ndes gt 0 then begin
  if ndes eq nw then dd[0,*]=desig[*] else $
   if ndes eq 2*nw then dd=reform(desig,2,nw) else $
    for i=0L,szd[2]-1L do dd[*,i]=desig[*,i]
endif
;	ECONF
ee=strarr(2,nw)
if necf gt 0 then begin
  if necf eq nw then ee[0,*]=econf[*] else $
   if necf eq 2*nw then ee=reform(econf,2,nw) else $
    for i=0L,sze[2]-1L do ee[*,i]=econf[*,i]
endif

;	make the output
linstr=create_struct('LINE_INT',line_int,'LOGT',tlog,'WVL',wave,$
	'Z',zz,'ION',ii,'JON',jj,$
	'DESIG',desig,'CONFIG',econf,'SRC',src)

return,linstr
end
