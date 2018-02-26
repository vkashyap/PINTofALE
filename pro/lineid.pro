function lineid,lamda,spec,elem,wrng=wrng,wshift=wshift,locmax=locmax,$
	best=best,topX=topX,batch=batch,stor=stor,oldid=oldid,omatch=omatch,$
	markw=markw,feature=feature,wrange=wrange,temp=temp,verbose=verbose,$
	_extra=e
;+
;function	lineid
;	allows identification of spectral features and returns ID
;	information in a structure (see structure variable IDHLP
;	for complete description of fields)
;
;syntax
;	id=lineid(lamda,spec,elem,wrng=wrng,wshift=wshift,locmax=locmax,$
;	best=best,topX=topX,batch=batch,stor=stor,oldid=oldid,omatch=omatch,$
;	markw=markw,feature=feature,verbose=verbose,$
;	n_e=n_e,logP=logP,pres=pres,dbdir=dbdir,/allah,chifil=chifil,$
;	eqfile=eqfile,chidir=chidir,dem=dem,abund=abund,effar=effar,$
;	wvlar=wvlar,sigmay=sigmay,xsize=xsize,ysize=ysize,wid=wid,$
;	dynrng=dynrng,markp=markp,marks=marks,marko=marko,$
;	/nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$
;	stathk=stathk)
;
;parameters
;	lamda	[INPUT; required] wavelength(s) at which spectrum is defined
;		* interactive if LAMDA is a multi-element array, batch mode
;		  otherwise (also see keyword BATCH)
;	spec	[INPUT; optional] the spectrum
;		* default is 10+random
;		* if SPEC is a character string or array, taken to be ELEM
;		* if LAMDA and SPEC do not match in size, taken to be atomic
;		  numbers
;	elem	[INPUT; optional] confine matches to this particular element
;		* overrides anything inferred "via" SPEC
;
;keywords
;	wrng	[INPUT] interval in which to search for matches
;		* default is WRNG=0.01
;		* if scalar and between 0 and 1, look for matches in the range
;		  	[W-W*WRNG,W+W*WRNG]
;		* if scalar and > 1, or scalar and < 0, abs(WRNG) is assumed
;		  to be the actual step size and is used as is, i.e.,
;		  	W+abs(WRNG)*[-1,1]
;		* if 1-element vector, then range is W+abs(WRNG(0))*[-1,1]
;		* if multi-element vector, then [W-WRNG[0],W+WRNG[1]]
;		* NOTE: in interactive mode the range may be altered, but
;		  the >behavior< of the range is set by the initial values.
;		* also note that WRNG is not changed on output.
;	wshift	[INPUT] amount by which to shift the wavelength scale
;		* WRNG will be reset to WRNG+WSHIFT
;	locmax	[INPUT] if set, finds matches at all local maxima
;		* search window width (2*abs(LOCMAX)+1)>3
;	best	[INPUT] if set, returns "best" match -- the strongest
;		possible nearest line
;		* obviously, there is weighting involved
;		* BEST.GE.0: pick largest of Flx[w]/(w-w0)^BEST
;		* otherwise: return the nearest match(es)
;	topX	[I/O; default: 20] returns at most TOPX matches
;		* MUST be scalar
;		* in interactive mode, presents only TOPX for matching choice
;		* if 0 < TOPX < 1, uses this as a >threshold<, selecting only
;		  those lines whose theoretical predicted fluxes are this
;		  fraction or higher than the maximum flux in the set.
;		* if TOPX < 0, the integer value is used as the upper limit
;		  to number of matches, and the fractional value is used as
;		  the threshold.
;		* may get altered within LINEID_MENU
;		* examples:
;		  TOPX = 0 brings up 20 matches
;		  TOPX = 10 brings up 10 matches
;		  TOPX = 0.6 brings up all lines with fluxes > 0.6 of max
;		  TOPX = -40.01 brings up 40 lines with fluxes > 1% of max
;	batch	[INPUT] if set, operates in batch mode (i.e., no widget-based
;		selections)
;		* batch mode is also forced if:-
;		  o N_ELEMENTS(LAMDA)=1
;		  o LOCMAX is set
;		  o TOPX is set to 1
;		  o !D.NAME is not 'X'
;	stor	[I/O] if present in call, will store the output of RD_LINE and
;		FOLD_IONEQ as a structure.  if a structure, then MUST be of
;		the form {LINE_INT,LOGT,WVL,Z,ION,DESIG,CONFIG}, because
;		then RD_LINE and FOLD_IONEQ will be bypassed!
;		* use with caution!!
;	oldid	[INPUT] if set and is in the same format as the output of
;		LINEID, starts off from this point.  i.e., things can be
;		added or modified from here on.
;	omatch	[INPUT] if OLDID is given *and* OMATCH is set, forces each
;		wavelength selection to match the nearest matched wavelength
;		in OMATCH.
;		* if OMATCH < 0, won't verify, but will force a match if
;		  new location is within abs(OMATCH), and will force a new
;		  match otherwise
;		* if OMATCH > 0, will ask
;	markw	[INPUT] if set, places marks on the plot at these wavelengths
;	feature	[INPUT] labels for MARKW
;		* passed with minimal modification to PICKRANGE
;	verbose	[INPUT] controls chatter
;	wrange	[IGNORED] keyword for RD_LINE, will be trapped and ignored.
;	temp	[IGNORED] keyword for LINEFLX, will be trapped and ignored.
;	_extra	[INPUT ONLY] pass predefined input keywords to subroutines:
;		pres	(RD_LINE) pressure [cm^-3 K]
;		logP	(RD_LINE) log10(pressure [cm^-3 K])
;		n_e	(RD_LINE) electron density [cm^-3]
;		dbdir	(RD_LINE) line emissivity database
;		allah	(RD_LINE) which elements to choose
;		chifil	(FOLD_IONEQ) set to read CHIANTI ion balance files
;		eqfile	(FOLD_IONEQ) name of ion balance file
;		chidir	(RD_IONEQ) location of CHIANTI installation
;		DEM	(LINEFLX) Differential Emission Measure distribution
;		abund	(LINEFLX) abundances
;		effar	(LINEFLX) effective area [cm^2]
;		wvlar	(LINEFLX) wavelength grid for EFFAR
;		sigmaY	(PICKRANGE) error-bars on SPEC
;		xsize	(PICKRANGE) window size
;		ysize	(PICKRANGE) window size
;		wid	(PICKRANGE) window ID
;		dynrng	(PICKRANGE) dynamic range in log scale
;		markp	(PICKRANGE) PSYM for (FEATURE,MARKW)
;		marks	(PICKRANGE) SYMSIZE for (FEATURE,MARKW)
;		marko	(PICKRANGE) ORIENTATION for (FEATURE,MARKW)
;		oxr	(PICKRANGE) final output XRANGE
;		oyr	(PICKRANGE) final output YRANGE
;		nuthin	(STAMPLE) don't stamp the plots
;		noday	(STAMPLE) don't stamp the plots with day
;		notime	(STAMPLE) don't stamp the plots with time
;		nouser	(STAMPLE) don't stamp the plots with username
;		nopack	(STAMPLE) don't stamp the plots with package name
;		stacol	(STAMPLE) stamp color
;		stasiz	(STAMPLE) stamp size
;		stathk	(STAMPLE) stamp thickness
;
;description
;	read in database of line cooling intensities and for each observed
;	wavelength which requires a match, find the set of matches.  if
;	interactive, allow selection of multiple matches.  once done,
;	reorganize into structure of
;	{WVL: observed wavelength;
;	  {corresponding
;	    WVL: matched wavelengths,
;	    Z: atomic numbers,
;	    ION: ionic states,
;	    LABL: e-configurations or whatever is available,
;	    FLUX: theoretical fluxes for all new or modified matches
;	    	(fluxes in OLDID are not touched unless ID is changed),
;	    FLUXERR: errors on FLUX (set to zero),
;	    LOGT: temperatures at which emissivities are defined
;	    EMIS: line emissivities for all IDs (ion-balance included)
;	    NOTES: scribbles about this ID, for human eyes only
;	  }
;	}
;
;restrictions
;	requires IDL v5 or greater
;	interactive selections require X-windows
;	requires subroutines
;	    RD_LINE [KILROY, SYMB2ZION [LAT2ARAB]]
;	    FOLD_IONEQ [WHEE, GETLOCMAX, RD_IONEQ [READ_IONEQ]]
;	    LINEFLX [WHEE]
;	    GETLOCMAX
;	    PICKRANGE
;	    LINEID_MENU [CAT_ID]
;	    CREATE_STRUCT
;	    KABOOM
;	    INICON
;	    STR_2_ARR
;
;usage examples
;	* lid=lineid(12.3985/7.83,['Fe','Ni'],wr=0.1*[1,1],topX=50)
;	* stor=1 & lid=lineid(x,y,wr=[1],/allah,logP=20,stor=stor)
;	  lid=lineid(x,y,wr=[1],/allah,logP=20,stor=stor)
;	* lid=lineid(x,y,'Fe',logP=20,econf=0)
;
;history
;	vinay kashyap (Dec96)
;	bug corrections: ion balance in STOR, format sizes in menu input;
;	  added ID history info to PICKRANGE and menu; added keyword OLDID
;	  (VK; Feb97)
;	added keyword OMATCH; allowed expeditious return when no matches
;	  found (VK; Jun98)
;	modified behavior of STOR to eliminate excess variables (VK; Oct98)
;	added LOGT and EMIS (EMIS includes ion balance) to the output; bug
;	  correction re STOR not being specified; corrected bug with
;	  BEST actually returning the worst; added keywords BATCH, MARKW,
;	  FEATURE; added calls to KABOOM, INITSTUFF, STR2ARR; TOPX now acts
;	  as threshold too; added support for multiple grating orders
;	  (VK; Nov98)
;	changed behavior of WRANGE; added renormalization of fluxes; bug fix
;	  older matches not echoing in menu window; added field FLUXERR to
;	  ID structure (VK; Dec98)
;	added field NOTES to ID structure (VK; Mar99)
;	modified ion balance calcs (VK; 99May)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	changed to full IDL5+ compatibility ("()"->"[]"); stor does not
;	  have to be defined on input, just be present in the call; added
;	  dummy keywords WRANGE and TEMP (VK; AugMM)
;	added keyword VERBOSE, and now doesn't default to random if NX eq NY+1
;	  (i.e., bin boundaries are given rather than mid-bin values: VK; Nov01)
;	changed call to STR_2_ARR from STR2ARR (VK; AprMMV)
;	asks for user confirmation in batch mode only if VERBOSE>10 (VK; JunMMXVI)
;-

;	initialize
iwav=0 & jwav=0 & desig=1 & econf=1 & code=0
atom=1 & rom=1 & inicon,atom=atom,roman=rom & atom=['X',atom] & rom=['',rom]

;	on-line help
w1=create_struct(name='W_'+strtrim(1,2),$
	'WVL','matching wavelengths [A] (fltarr(N_1))',$
	'Z','corresponding atomic numbers (intarr(N_1))',$
	'ION','corresponding ion states (intarr(N_1))',$
	'LABL','electron configurations & designations (strarr(2,N_1))',$
	'FLUX','observed or computed ([ph/s{/cm^2}]) fluxes (fltarr(N_1))',$
	'FLUXERR','Errors on FLUX, set to zero if theoretical (fltarr(N_1))',$
	'LOGT','temperatures at which EMIS are defined',$
	'EMIS','line emissivities [ergs cm^3/s]',$
	'NOTES','scribble while you work..')
wi=create_struct(name='W_i',$
	'WVL','...','Z','...','ION','...','LABL','...','FLUX','...',$
	'FLUXERR','...','LOGT','...','EMIS','...','NOTES','(optional)')
wN=create_struct(name='W_N',$
	'WVL','matching wavelengths [A] (fltarr(N_N))',$
	'Z','corresponding atomic numbers (intarr(N_N))',$
	'ION','corresponding ion states (intarr(N_N))',$
	'LABL','electron configurations & designations (strarr(2,N_N))',$
	'FLUX','observed or computed ([ph/s{/cm^2}]) fluxes (fltarr(N_N))',$
	'FLUXERR','Errors on FLUX, set to zero if theoretical (fltarr(N_1))',$
	'LOGT','temperatures at which EMIS are defined',$
	'EMIS','line emissivities [ergs cm^3/s]')
idhlp=create_struct(name='LINEID_HELP','WVL',fltarr(3),$
	'WVL_COMMENT','observed wavelengths [A] (fltarr(N))',$
	'ID'+strtrim(1,2),w1, 'IDi',wi, 'IDN',wN)
;
if n_elements(lamda) eq 0 then begin
  print,'Usage: idstr=lineid(lamda,spec,elem,wrng=wrng,wshift=wshift,locmax=locmax,$'
  print,'       best=best,topX=topX,batch=batch,stor=stor,oldid=oldid,omatch=omatch,$'
  print,'       markw=markw,feature=feature,verbose=verbose,$'
  print,'       n_e=n_e,logP=logP,pres=pres,dbdir=dbdir,/allah,chifil=chifil,$'
  print,'       eqfile=eqfile,chidir=chidir,dem=dem,abund=abund,effar=effar,$'
  print,'       wvlar=wvlar,sigmay=sigmay,xsize=xsize,ysize=ysize,wid=wid,$'
  print,'       dynrng=dynrng,markp=markp,marks=marks,marko=marko,$'
  print,'       /nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$'
  print,'       stathk=stathk)'
  print,'  matches wavelengths with lines and returns info in structure'
  return,idhlp
endif

;	check inputs
x=[lamda]						;wavelength
nx=n_elements(x) & ny=n_elements(spec)
sz=size(spec) & ns=n_elements(sz) & type=sz(ns-2)
if n_elements(spec) ne 0 then y=[spec] else y=10.+randomu(seed,nx)>5.
if n_elements(elem) eq 0 and type eq 7 then begin
  y=10.+randomu(seed,nx)>5.				;constant, rough top
  elem=spec						;ELEM, not SPEC
endif
if nx ne ny and nx ne ny+1L then begin
  y=10.+randomu(seed,nx)>5.				;constant, rough top
  ;	assume atomic #s
  if n_elements(elem) eq 0 and ny ne 0 then begin
    print,'' & print,''
    message,'WARNING: LAMBDA and SPEC do not match.',/info
    message,'Assuming that you mean SPEC->ELEM and SPEC==RANDOM',/info
    print,'' & print,''
    elem=atom[long(spec)]
  endif
endif

;	check keywords
vv=0 & ivar=0 & defsysv,'!VERBOSE',exists=ivar	;if !VERBOSE exists
if ivar ne 0 then setsysval,'VERBOSE',vv,/getval
if n_elements(verbose) gt 0 then vv=long(verbose(0)) > 1
;
if not keyword_set(wrng) then wrng=0.01			;search area
if n_elements(best) eq 0 then best=-1			;find nearest
if n_elements(topX) eq 0 then topX=20			;at most XX
thresh=0.		;include all lines with flux > THRESH
if topX lt 1 then begin
  tmp=abs(topX) & thresh=tmp-fix(tmp)		;topX acts as a threshold
endif
if not keyword_set(wshift) then ws=0. else ws=wshift	;shift wvl scale
;
if keyword_set(oldid) then begin			;begin where ye left off
  if n_tags(oldid) gt 1 then begin			;structure, allright
    tname=tag_names(oldid)
    ;	output of LINEID, and not IDHLP either
    if tname[0] eq 'WVL' and tname[1] eq 'ID1' then begin
      outwvl=oldid.(0)				;wavelengths
      idwvl=create_struct(tname[1],oldid.(1))	;strip the IDs
      xw=[oldid.(1).wvl]
      for i=2,n_tags(oldid)-1 do begin
	idwvl=create_struct(idwvl,tname[i],oldid.(i))
        xw=[xw,oldid.(i).wvl]
      endfor
      jwav=n_tags(idwvl)
    endif
  endif
endif

;	check whether interactive or batch mode
if n_elements(batch) eq 0 then batch=0			;default is interactive
if nx eq 1 then batch=1					;only one point
if topX eq 1 then batch=1				;1 for all, all for 1
if keyword_set(locmax) then batch=1			;at all local maxima
if !d.name ne 'X' then begin				;not X-windows?  feh.
  message,'WARNING: current device is: '+!D.NAME,/info
  message,'forcing BATCH mode',/info
  batch=1
endif

;	figure out search areas
wr=float(wrng) & szr=size(wr) & frang=0
if szr[0] eq 0 then frang=1				;fractional ranges
if frang eq 1 then begin
  if wr[0] eq 0 then wr[*]=0.1				;default
  if wr[0] lt 0 or wr[0] gt 1 then begin
    frang=0 & wr=abs(wr[0])*[1,1]			;actual width
  endif
endif
if n_elements(wr) eq 1 then wr=abs(wr[0])*[1,1] else wr=[abs(wr[0]),abs(wr[1])]
if wr[1] lt wr[0] then wr=[wr[1],wr[0]]

;	interval of interest
w0=min(x,max=w1)
if frang eq 1 then begin
  w0=w0-w0*abs(wr[0]) & w1=w1+w1*abs(wr[1])
endif
if frang eq 0 then begin
  w0=w0-abs(wr[0]) & w1=w1+abs(wr[1])
endif

;	compute fluxes for all lines in interval of interest
nstor=n_elements(stor) & mstor=n_tags(stor)
if nstor eq 0 then stor=0
if mstor gt 0 then begin		;(verify that input is correct
  stornam=tag_names(stor) & ok='ok'
  if n_elements(stornam) lt 5 then ok='Structure incomplete?' else $
   if stornam[0] ne 'LINE_INT' then ok='missing emissivities' else $
    if stornam[1] ne 'LOGT' then ok='missing temperatures' else $
     if stornam[2] ne 'WVL' then ok='missing wavelengths' else $
      if stornam[3] ne 'Z' then ok='missing atomic numbers' else $
       if stornam[4] ne 'ION' then ok='missing ionic states'
  if ok ne 'ok' then mstor=0
endif
if mstor eq 0 then begin		;(first run
  if arg_present(stor) then begin	;(time-saver
    ff=rd_line(elem,wrange=[w0,w1],wvl=wvl,logT=logT,Z=Z,ion=ion,jon=jon,$
      desig=desig,econf=econf,fstr=stor,verbose=vv, _extra=e)
    ff=fold_ioneq(ff,z,jon,logt=logt,tmax=tmax,trng=trng,verbose=vv,_extra=e)
    stor.line_int=ff
  endif else begin			;)(space-saver
    ff=rd_line(elem,wrange=[w0,w1],wvl=wvl,logT=logT,Z=Z,ion=ion,jon=jon,$
      desig=desig,econf=econf,verbose=vv, _extra=e)
    ff=fold_ioneq(ff,z,jon,logt=logt,tmax=tmax,trng=trng,verbose=vv,_extra=e)
  endelse				;STOR)
endif else begin				;)(decode STOR
  ff=stor.line_int & logT=stor.logt & wvl=stor.wvl & z=stor.z
  ion=stor.ion & desig=stor.desig & econf=stor.config
  if min(abs(wvl)) gt w0 or max(abs(wvl)) lt w1 then begin
    print,'Emissivities cover the range:',min(abs(wvl)),' - ',max(abs(wvl))
    print,'current wavelength range:',w0,' - ',w1
    message,'WARNING: input emissivities do not cover the',/info
    message,'         wavelength range of interest',/info
    print,string(7B)
    print,'hit ctrl-C within 1 sec to stop'
    wait,1
  endif
endelse						;STORe)
;
;	get Tmax and Trange
nwvl=n_elements(wvl) & nt=n_elements(logT)
tmax=fltarr(nwvl) & trng=fltarr(2,nwvl)
for iw=0L,nwvl-1L do begin
  tmp=reform(ff[*,iw]) & tmpmax=max(tmp,imax) & tmax[iw]=logT[imax]
  oo=where(tmp ge 0.1*tmpmax,moo)
  if moo gt 0 then trng[*,iw]=[logT[oo[0]],logT[oo[moo-1]]] else $
    trng[*,iw]=[4.,8.]
endfor
;
;	get line fluxes
f=lineflx(ff,logT,wvl,Z,_extra=e)		;[ph/s(/cm^2)]

;	summarize the level descriptions
nw=n_elements(wvl) & labl=strarr(2,nw)
if n_elements(desig) eq 2*nw then labl=desig
if n_elements(econf) eq 2*nw then labl='('+econf+') '+labl

minf=min(f[[where(f ne 0)]])
if minf eq 0 then begin
  message,'	No lines were found; exiting',/info & return,idhlp
endif

;	in batch mode, get local maxima
ilm=lindgen(nx)
if keyword_set(locmax) then ilm=getlocmax(x,y,width=locmax)

;	in batch mode, figure out which wavelengths to ID
xx=x & yy=y
if keyword_set(batch) then begin
  xx=xx[ilm] & yy=yy[ilm]
endif
nwav=n_elements(xx)

;	too many to ID in batch mode?
IDLversion=float(!VERSION.RELEASE)
if keyword_set(batch) and nwav ge 126 then begin
  if IDLversion lt 5 then begin
    c1='Cannot handle more than 126 matches at a time!' & message,c1,/info
    c1='Will return only the first 126' & message,c1,/info
  endif else begin
    if vv gt 0 then begin
      c1='There will be '+strtrim(nwav,2)+' IDs!  Is this really what you want?'
      message,c1,/info
      if vv gt 10 then begin
        C1='hit q to quit, f to return first 126, any other key to continue'
        print,c1 & c1=get_kbrd(1)
        if strlowcase(c1) eq 'q' then return,idhlp
        if strlowcase(c1) eq 'f' then IDLversion=4	;yes, I know this is a hack.
      endif else begin
        print,'hit ctrl-C within 1 sec to stop'
        wait,1
      endelse
    endif
  endelse
endif

c1='	Choose input and matching wavelengths interactively'
if not keyword_set(batch) then message,c1,/info

matchwvl=fltarr(nw)-1.	;to keep track of what's been matched with what
;	if OLDID is supplied, fill up this array as best as possible
if keyword_set(xw) then begin
  ww=oldid.WVL
  for i=0L,n_elements(ww)-1L do begin
    tmpid=oldid.(i+1) & wwx=tmpid.WVL & zzx=tmpid.Z & iix=tmpid.ION
    for j=0L,n_elements(wwx)-1L do begin
      oo=where(abs(wvl-wwx[j]) lt 1e-4 and z eq zzx[j] and ion eq iix[j],moo)
      if moo gt 0 then matchwvl[oo]=ww[i]
    endfor
  endfor
endif

;	loop to select wavelength to id
while iwav lt nwav do begin		;{do, a loop, a simple loop

  if batch eq 0 then begin		;(interactively pick a wavelength to ID
    oo=-1L
    while oo[0] eq -1 do begin
      ;	get default window size
      if !d.window ge 0 then begin & xs=!d.x_size & ys=!d.y_size & endif
      ;	use PICKRANGE to select a range
      if keyword_set(outwvl) then begin
	wmark=[outwvl] & nwmark=n_elements(wmark)
	fmark=strarr(nwmark)
	for i=0,nwmark-1 do begin
	  tmpid=idwvl.(i) & zz=tmpid.Z & jon=tmpid.ION
	  nzz=n_elements(zz)
	  fmark[i]='  '+atom(zz[0])+rom(jon[0])
	  for j=1,nzz-1 do fmark[i]=fmark[i]+' ; '+atom[zz[j]]+rom[jon[j]]
	endfor
      endif else fmark=''
      if keyword_set(markw) then begin
	nmarkw=n_elements(markw) & nfeatr=n_elements(feature)
	ffmark=strarr(nmarkw)
	if nfeatr gt 0 then begin
	  if nfeatr le nmarkw then ffmark[0L:nfeatr-1L]=feature else $
	    ffmark=feature[0L:nmarkw-1L]
	endif
	if keyword_set(outwvl) then begin
	  wmark=[markw,outwvl]
	endif else begin
	  wmark=markw
	endelse
	fmark=[ffmark,fmark]
      endif
      oo=pickrange(xx,yy,xsize=xs,ysize=ys,oxr=oxr,oyr=oyr,$
	markx=wmark,markl=fmark,/legg,marko=90, _extra=e)
      if oo[0] eq -1 then begin
	kaboom,'?',/flash
	c1='no data selected: type q to quit, x/z to halt, any key to continue'
	message,c1,/info & c1=get_kbrd(1)
	if c1 eq 'x' or c1 eq 'z' then stop,'halting (.CON to continue)'
	if c1 eq 'q' then begin
	  if jwav eq 0 then return,idhlp else $
	    return,create_struct('WVL',outwvl,idwvl)
	endif
      endif
    endwhile
    ;	choose maximum in selected range as wavelength to ID
    i_wvl=(where(yy[oo] eq max(yy[oo])))[0] & x_wvl=(xx[oo])[i_wvl]
    ;	force to match with previous?
    if n_tags(oldid) gt 0 and keyword_set(omatch) then begin
      owvl=oldid.wvl & dx=min(abs(x_wvl-owvl),iwvl)
      if dx gt 1e-4*x_wvl then begin
	;force match if within abs(OMATCH), force no match if beyond
	if omatch lt 0 then begin
	  if dx le abs(omatch) then x_wvl=owvl[iwvl]
	endif
	if omatch gt 0 then begin		;(ask
	  print,string("7b)
	  print,'selected wavelength ('+strtrim(x_wvl,2)+') is close to '+$
	  	'previously IDd '+strtrim(owvl[iwvl],2)
	  c1='make the switch? [n/y]' & print,c1 & c1=get_kbrd(1)
	  if strlowcase(c1) eq 'y' then x_wvl=owvl[iwvl]
	endif					;OMATCH=1)
      endif
      print,'IDing @... '+strtrim(x_wvl,2)
    endif
  endif else begin			;)(batch mode
    x_wvl=xx[iwav] & iwav=iwav+1
  endelse				;end batch mode)
  ywvl=min(abs(xx-x_wvl),iywvl) & y_wvl=yy[iywvl]

  ;	look for matches in what range?  fractional/incremental
  if not keyword_set(cord) then cord='1'
  PARS:	;the one and only goto, just to reacquire the defaults (see CODE)
  ;	reacquire WRANGE
  if frang eq 1 then begin
    w0=x_wvl-x_wvl*abs(wr[0]) & w1=x_wvl+x_wvl*abs(wr[1])
  endif
  if frang eq 0 then begin
    w0=x_wvl-abs(wr[0]) & w1=x_wvl+abs(wr[1])
  endif
  w0=w0+ws & w1=w1+ws			;shift wavelength scale
  if topX lt 1 then begin
    tmp=abs(topX) & thresh=tmp-fix(tmp)		;topX acts as a threshold
  endif
  ;	to include multiple grating orders..
  ords=str_2_arr(cord,/squish) & nord=n_elements(ords)
  gwvl=wvl & gf=f	;default is first order
  for i=0,nord-1 do begin
    if i eq 0 then begin
      gwvl=wvl*ords[i]
      if ords[i] ne 0 then gf=f/abs(ords[i]^2) else gf=f
    endif else begin
      gwvl=[gwvl,wvl*ords[i]]
      if ords[i] ne 0 then gf=[gf,f/abs(ords[i]^2)] else gf=[gf,f]
    endelse
  endfor

  if ords[0] ne 1 or nord gt 1 then message,$
    'fluxes from grating orders other than the 1st are rough estimates',/info

  ;	now find the IDs for this chosen wavelength
  hh=where(abs(gwvl) ge w0 and abs(gwvl) le w1 and gf gt 0,moo)
  if moo gt 0 then begin
    hhh = hh mod nwvl
    fmax=max(f[hhh])
    hh=where(abs(gwvl) ge w0 and abs(gwvl) le w1 and gf gt thresh*fmax,moo)
    hhh = hh mod nwvl
  endif

  tipX=fix(abs(topX))
  if tipX eq 0 then tipX=20
  ;if moo gt tmp then moo=tmp	;upper limit to number of choices

  ;	add to collection
  if moo gt 0 then begin		;{something to add...
    outgw=gwvl[hh] & outw=wvl[hhh] & outz=z[hhh] & outi=ion[hhh]
    outl=labl[*,hhh] & outf=f[hhh] & outt=tmax[hhh]
    outm=trng[*,hhh] & outff=ff[*,hhh] & outmw=matchwvl[hhh]
    ;	sort in order of "best" or nearest
    tmp=abs(x_wvl-abs(outgw))
    if best ge 0 then tmp=(tmp)^(best)/(outf>minf)
    	;note that this appears (but only appears) to be the INVERSE of what
	;the documentation says: however, "best" == minimum (TMP), i.e., the
	;first element of SORT(TMP)
    oo=sort(tmp)
    if moo gt tipX and tipX gt 1 then moo=tipX	;upper limit to #choices
    oo=oo[0:moo-1]
    outf=outf[oo] & outgw=outgw[oo]
    outw=outw[oo] & outz=outz[oo] & outi=outi[oo] & outl=outl[*,oo]
    outt=outt[oo] & outm=outm[*,oo] & outff=outff[*,oo]
    outmw=outmw[oo]

    ;	pick some, any some
    if batch eq 0 then begin		;{if interactive, then pick
      tmp=strarr(moo)
      slen0=max(strlen(outl[0,*])) & slen1=max(strlen(outl[1,*]))
      print,'' & print,'Chosen wavelength='+strtrim(x_wvl,2)+' A'
      print,'ID by clicking with mouse'
      if not keyword_set(fmult) then fmult=1.
      if keyword_set(renorm) then begin
	if renorm lt 0 then begin
	  if renorm ne -1 then fmult=abs(renorm) else begin
	    if fmax gt 0 then begin
	      fmult=abs(y_wvl/fmax)
	      renorm=-fmult
	      if fmult eq 1 then renorm=renorm-1e-5*randomu(seed)
	    endif
	  endelse
	endif
	if renorm gt 0 and fmax gt 0 then fmult=abs(float(renorm)/fmax)
      endif
      tmp=string(outw,'(f10.4)')+$
	string(atom[outz],'(a3)')+string(outi,'(i4)')+' '+$
	string(fmult*outf,'(g12.5)')+' '+$
	string(outt,'(f6.3)')+' <'+$
	string(outm[0,0:moo-1],'(f6.3)')+'-'+$
	string(outm[1,0:moo-1],'(f6.3)')+'> '+$
	string(outl[0,0:moo-1],'(a'+strtrim(slen0,2)+')')+' '+$
	string(outl[1,0:moo-1],'(a'+strtrim(slen1,2)+')')
      for imw=0,n_elements(tmp)-1 do if outmw[imw] ne -1 then tmp[imw]='*'+$
		tmp[imw]+' ['+string(outmw[imw],'(f10.4)')+']'
      oh=lineid_menu(tmp,x_wvl,abs(outw),outf,ws,topx,wr,cord,code,scratch,$
		     wvls=outwvl,idwvls=idwvl,gwvl=abs(outgw),renorm=renorm)
      if code eq 3 then oh=-2
    endif else begin			;end interactive picking
      oh=lindgen(moo)
      if best ge 0 then oh=[0L]		;return only the best
    endelse				;end picking}

    if oh[0] ne -1 then begin		;(selected any?
      jwav=jwav+1
      if oh[0] eq -2 then begin
	unlabl=strarr(2,1) & unlabl[0]='Unknown'
	w=create_struct('WVL',[-x_wvl],'Z',[0],'ION',[0],'LABL',unlabl,$
	  'FLUX',[0.],'FLUXERR',[0.],'LOGT',logt,'EMIS',fltarr(nt,1))
	if keyword_set(scratch) then w=create_struct(w,'NOTES',scratch)
      endif else begin
        w=create_struct('WVL',[outw[oh]],'Z',[outz[oh]],'ION',[outi[oh]],$
	  'LABL',[outl[*,oh]],'FLUX',[outf[oh]],'FLUXERR',0*[outf[oh]],$
	  'LOGT',logt,'EMIS',outff[*,oh])
	if keyword_set(scratch) then w=create_struct(w,'NOTES',scratch)
	matchwvl[hhh[oo[oh]]]=x_wvl		;these are matched
      endelse
      if not keyword_set(outwvl) then begin	;first match
        outwvl=[x_wvl]
        idwvl=create_struct('ID'+strtrim(jwav,2),w)
      endif else begin				;subsequent matches
	kwvl=(where(outwvl eq x_wvl))[0]	;if a "repeat", overwrite!
	if kwvl eq -1 then begin		;(new
          outwvl=[outwvl,x_wvl]
          idwvl=create_struct(idwvl,'ID'+strtrim(jwav,2),w)
	endif else begin			;)(old
	  message,'overwriting @ '+strtrim(kwvl[0]+1,2),/info
	  jwav=jwav-1
	  tmp=idwvl & tt=tag_names(tmp)
	  case kwvl of
	    0: begin
	      idwvl=create_struct(tt[0],w)
	      for i=1,jwav-1 do idwvl=create_struct(idwvl,tt[i],tmp.(i))
	    end
	    else: begin
	      idwvl=create_struct(tt[0],tmp.(0))
	      for i=1,kwvl-1 do idwvl=create_struct(idwvl,tt[i],tmp.(i))
	      idwvl=create_struct(idwvl,tt[kwvl],w)
	      for i=kwvl+1,jwav-1 do idwvl=create_struct(idwvl,tt[i],tmp.(i))
	    endelse
	  endcase
	endelse					;KWVL)
      endelse
    endif				;oh)

  endif else begin			;}{nothing to add..
    c1='No matches for '+strtrim(x_wvl,2)+' [A] in ['+strtrim(w0,2)+$
	','+strtrim(w1,2)+'] [A]'
    message,'	'+c1,/info
    if (where([ords[*]] eq 1))[0] lt 0 then begin
      message,'resetting grating order to first only',/info & cord='1'
    endif
  endelse				;moo=0}

  ;	are we done?
  if not keyword_set(batch) then begin
    if code eq 1 then iwav=nwav			;quit
    if code eq 2 then goto, PARS		;redo
  endif

  ;	an embarassment of riches? (an IDL4 limitation)
  if jwav eq 125 and IDLversion lt 5 then begin
    c1='Reached maximum allowed number of matches; exiting'
    message,c1,/info & iwav=nwav
  endif

endwhile				;iwav < nwav}

;	this is the output
;	you know there were no matches if (IDWVL.WVL)[0] is 0
if jwav gt 0 then idwvl=create_struct('WVL',outwvl,idwvl) else idwvl=idhlp

return,idwvl
end
