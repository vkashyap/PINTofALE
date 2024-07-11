function fitlines,x,y,ysig=ysig,funcs=funcs,intens=intens,dchi=dchi,$
	pos=pos,perrp=perrp,perrm=perrm,perrc=perrc,$
	flx=flx,ferrp=ferrp,ferrm=ferrm,ferrc=ferrc,$
	wdt=wdt,werrp=werrp,werrm=werrm,werrc=werrc,$
	thaw=thaw,type=type,epithet=epithet,ties=ties,$
	conlev=conlev,consig=consig,comment=comment,$
	histerix=histerix,oldstr=oldstr, _extra=e
;+
;function	fitlines
;	widget-based procedure to fit Gaussians and Lorentzians (or
;	combinations, or other 3-parameter functions) to lines in a
;	spectrum and determine line fluxes.  returns a structure of
;	the form
;	  {POS,PERRP,PERRM,PERRC,$
;	   FLX,FERRP,FERRM,FERRC,$
;	   WDT,WERRP,WERRM,WERRC,$
;	   THAW,TYPE,TIES,EPITHET,CONLEVX,CONLEV,CONSIGX,CONSIG,COMMENT}
;	All these are also available as output individually as keywords.
;	See keyword descriptions for what the variables mean.
;
;syntax
;	fxstr=fitlines(x,y,ysig=ysig,funcs=funcs,/intens,dchi=dchi,$
;	pos=pos,perrp=perrp,perrm=perrm,perrc=perrc,$
;	flx=flx,ferrp=ferrp,ferrm=ferrm,ferrc=ferrc,$
;	wdt=wdt,werrp=werrp,werrm=werrm,werrc=werrc,$
;	thaw=thaw,type=type,epithet=epithet,ties=ties,$
;	conlev=conlev,consig=consig,comment=comment,$
;	histerix=histerix,$
;	/dumb,verbose=verbose,xsize=xsize,ysize=ysize,$
;	wid=wid,dynrng=dynrng,/posve,/negve,itmax=itmax,chithr=chithr,$
;	jumpup=jumpup,jumpdn=jumpdn,svdthr=svdthr,missing=missing,$
;	/noday,/notime,/nouser,/nopack,stacol=stacol,stacol=stacol,$
;	stathk=stathk,/nuthin)
;
;parameters
;	x	[INPUT; required] absissa, e.g., wavelength or energy
;	y	[INPUT; required] count spectrum Y(X)
;		* size of Y must match that of X
;
;keywords
;	ysig	[INPUT] errors on Y
;		* default is sqrt(abs(Y)+0.75)+1.
;	funcs	[INPUT] name of user defined function (actually a procedure)
;		that takes as input X and A (the model parameters), and returns
;		Y(X;A) and the partial derivatives dY/dA
;		* default is set to X3MODEL
;		* NOTE: if MPFIT is chosen as the optimization algorithm,
;		  then "_f" is automatically appended to the name.  This
;		  is because MPFIT requires a _function_, not a procedure.
;		  The program corresponding to X3MODEL is X3MODEL_F() and
;		  that corresponding to LIBMODEL is LIBMODEL_F()
;	intens	[INPUT] if set, Y(X) is assumed to be intensity, i.e.,
;		(counts|ergs|.../bin_unit/...)
;		* default is to take Y(X) to be (counts|ergs|.../bin/...)
;		* passed straight to FITLINES_EVENT
;		* WARNING: this has >not< been road-tested!
;	dchi	[INPUT] delta-chisq threshold for projection errors
;		* default is 1.0
;		* may be changed within FITLINES_EVENT, but those changes
;		  will not be propagated back.
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
;	histerix	[OUTPUT] a structure containing the state at each step
;			of the fitting process
;			* if explicitly set to 0, will not contain any output
;			* passed w/o comment to FITLINES_EVENT
;       oldstr  [INPUT] old fitlines structure from which to start with 
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
;
;subroutines
;	FITLINES_EVENT
;	FMTSTUFF
;	PICKRANGE
;	FUNCS
;	(X3MODEL/LIBMODEL/MPFITFUN)
;	(MK_3MODEL, MK_GAUSS, MK_LORENTZ, MK_SLANT, MK_ROGAUSS, MK_POWLAM, MK_POLY1D)
;	LINEREM
;	SETCONT
;	ERORS
;	FIT_LEVMAR
;	ADJUSTIE
;	LEVMARQ
;	LMCOEFF
;	CURVE_FIT
;	SPLAC
;	WHICH
;	KILROY
;	WHEE
;	STAMPLE
;	PEASECOLR
;	IS_KEYWORD_SET
;
;known bugs
;	cannot handle spectra with varying bin sizes
;	INTENS has not been tested
;	SPLAC keywords are not gracefully handled
;	help is available only in UNIX
;
;history
;	vinay kashyap (Oct98)
;	added renormalization option (VK; Nov98)
;	added keywords INTENS,DCHI,PERRL,WERRL,FERRL,PERRC,WERRC,FERRC;
;	  changed output structure format (VK; FebMM)
;	changed keywords PERR,WERR,FERR to PERRU,WERRU,FERRU; also the
;	  corresponding output structure fields (VK; MarMM)
;	changed ?ERRU to ?ERRP and ?ERRL to ?ERRM to avoid confusion
;	  with ERR[U,L] of ERORS (VK; MarMM)
;	moved QUIT to end of row, far from FIT; added ability to save
;	  to disk; changed ?ERR[M,C] to be I/O; added TIES to output
;	  structure (VK; JulMM)
;	added labeling capability via keyword EPITHET; reorganized
;	  to combine 6th and 7th rows (VK; AugMM)
;	changed location of doc file to ardb ; allowed call to STAMPLE
;	  (VK; DecMM)
;	changed default of ?ERRC to 0.0 (VK; JanMMI)
;	per Antonio Maggio's suggestion, removed keyword MULTIPLE from
;	  call to WIDGET_LIST (VK; FebMMI)
;	changed suffix of help file from ".doc" to ".hlp" (VK; FebMMI)
;	added code to make it easier to switch between color tables
;	  (VK; Oct02)
;	added monte-carlo errors as option (LL; Aug03) 
;	added keyword HISTERIX (VK; Feb04)
;       added keyword OLDSTR (LL; Jul05) 
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;-

;	usage
nx=n_elements(x) & ny=n_elements(y)
ostr={POS:[-1.], PERRP:[0.], PERRM:[0.], PERRC:[0.],$
      FLX:[-1.], FERRP:[0.], FERRM:[0.], FERRC:[0.],$
      WDT:[-1.], WERRP:[0.], WERRM:[0.], WERRC:[0.],$
      THAW:[0], TYPE:['NONE'], TIES:[''], EPITHET:[''],$
      CONLEV:[-1.], CONSIG:[0.], COMMENT:'help'}
if nx eq 0 or nx ne ny then begin
  print,'Usage: fxstr=fitlines(x,y,ysig=ysig,funcs=funcs,/intens,dchi=dchi,$'
  print,'	pos=pos,perrp=perrp,perrm=perrm,perrc=perrc,$'
  print,'	flx=flx,ferrp=ferrp,ferrm=ferrm,ferrc=ferrc,$'
  print,'	wdt=wdt,werrp=werrp,werrm=werrm,werrc=werrc,$'
  print,'	thaw=thaw,type=type,epithet=epithet,ties=ties,$'
  print,'       conlev=conlev,consig=consig,comment=comment,$'
  print,'       histerix=histerix,$'
  print,'       /dumb,verbose=verbose,xsize=xsize,ysize=ysize,wid=wid,$'
  print,'	dynrng=dynrng,/posve,/negve,itmax=itmax,chithr=chithr,$'
  print,'	jumpup=jumpup,jumpdn=jumpdn,svdthr=svdthr,missing=missing'
  print,'	/noday,/notime,/nouser,/nopack,stacol=stacol,/nuthin) '
  return,ostr
endif

;       check OLDSTR 
if n_tags(oldstr) gt 0 then begin 
  pos = oldstr.pos & perrp = oldstr.perrp & perrm = oldstr.perrm 
  wdt = oldstr.wdt & werrp = oldstr.werrp & werrm = oldstr.werrm 
  flx = oldstr.flx & ferrp = oldstr.ferrp & ferrm = oldstr.ferrm 
  epithet = oldstr.epithet & flx = oldstr.flx 
  ties = oldstr.ties & type = oldstr.type & conlevx = oldstr.conlevx 
  consigx = oldstr.consigx & consig = oldstr.consig  & conlev = oldstr.conlev
endif 

;	check YSIG
sigy=sqrt(abs(y)+0.75)+1. & ns=n_elements(ysig)
ms = ns < nx
if ns gt 0 then sigy(0:ms-1)=ysig(0:ms-1)

;	check user-defined function name
if not keyword_set(funcs) then funcs='x3model'

;	check if dCHI is given
dchisq=1.0 & if keyword_set(dchi) then dchisq=float(dchi(0))

;	check model inputs
np=n_elements(pos) & nf=n_elements(flx) & nw=n_elements(wdt)
npe=n_elements(perrp) & nfe=n_elements(ferrp) & nwe=n_elements(werrp)
if not is_keyword_set(perrp) then perrp=-1
if not is_keyword_set(werrp) then werrp=-1
if not is_keyword_set(ferrp) then ferrp=-1
npem=n_elements(perrm) & nfem=n_elements(ferrm) & nwem=n_elements(werrm)
npec=n_elements(perrc) & nfec=n_elements(ferrc) & nwec=n_elements(werrc)
nt=n_elements(thaw) & nm=n_elements(type) & nepi=n_elements(epithet)
if np eq 0 then begin
  message,'Forcing at least one component!',/info
  tmp=max(y,imx) & pos=x(imx) & np=1L
endif

if np ne 0 then begin		;(consistency checks for input model pars
  t=fltarr(np)
  mf = nf < np
  mw = nw < np
  mpe = npe < np
  mfe = nfe < np
  mwe = nwe < np
  mpem = npem < np
  mfem = nfem < np
  mwem = nwem < np
  mpec = npec < np
  mfec = nfec < np
  mwec = nwec < np
  mt = nt < 3*np
  mm = nm < np
  mepi = nepi < np
  ;
  t(*)=1.
  if nf gt 0 then begin t(*)=flx(0) & t(0:mf-1)=flx(0:mf-1) & endif
  flx=t
  t(*)=1.
  if nw gt 0 then begin t(*)=wdt(0) & t(0:mw-1)=wdt(0:mw-1) & endif
  wdt=t
  t(*)=-1.
  if npe gt 0 then begin t(*)=perrp(0) & t(0:mpe-1)=perrp(0:mpe-1) & endif
  perrp=t
  t(*)=-1.
  if nfe gt 0 then begin t(*)=ferrp(0) & t(0:mfe-1)=ferrp(0:mfe-1) & endif
  ferrp=t
  t(*)=-1.
  if nwe gt 0 then begin t(*)=werrp(0) & t(0:mwe-1)=werrp(0:mwe-1) & endif
  werrp=t
  t(*)=perrp
  if npem gt 0 then begin t(*)=perrm(0) & t(0:mpem-1)=perrm(0:mpem-1) & endif
  perrm=t
  t(*)=ferrp
  if nfem gt 0 then begin t(*)=ferrm(0) & t(0:mfem-1)=ferrm(0:mfem-1) & endif
  ferrm=t
  t(*)=werrp
  if nwem gt 0 then begin t(*)=werrm(0) & t(0:mwem-1)=werrm(0:mwem-1) & endif
  werrm=t
  t(*)=0.0
  if npec gt 0 then begin t(*)=perrc(0) & t(0:mpec-1)=perrc(0:mpec-1) & endif
  perrc=t
  t(*)=0.0
  if nfec gt 0 then begin t(*)=ferrc(0) & t(0:mfec-1)=ferrc(0:mfec-1) & endif
  ferrc=t
  t(*)=0.0
  if nwec gt 0 then begin t(*)=werrc(0) & t(0:mwec-1)=werrc(0:mwec-1) & endif
  werrc=t
  t=intarr(3*np)
  if nt gt 0 then t(0:mt-1)=thaw(0:mt-1)
  thaw=t
  t=strarr(np)+'gauss'
  if nm gt 0 then begin t(*)=type(0) & t(0:mm-1)=type(0:mm-1) & endif
  type=t
  t='['+strtrim(pos,2)+']'
  if nepi gt 0 then t(0:mepi-1)=epithet(0:mepi-1)
  epithet=t
endif				;NP>0)

;;	define optional output arrays
;perrm=perrp & werrm=werrp & ferrm=ferrp
;perrc=0.*perrp & werrc=0.*werrp & ferrc=0.*ferrp

nc=n_elements(conlev) & nce=n_elements(consig)
if nc ne nx then begin		;(consistency check for continuum
  if nc eq 1 then begin
    t=fltarr(nx) & t(*)=conlev(0) & conlev=t
    t=fltarr(nx) & t(*)=consig(0) & consig=t
  endif else begin
    if nc gt 1 then message,$
	'continuum does not match input spectrum.. ignoring',/info
    conlev=fltarr(nx) & consig=conlev
  endelse
endif				;NC .NE. NP)
nc=n_elements(conlev) & nce=n_elements(consig)

if nce ne nc then begin		;(consistency check for error on continuum
  t=sqrt(abs(conlev)+0.75)+1.
  mce = nce < nc
  if nce gt 0 then t(0:mce-1)=consig(0:mce-1) & consig=t
endif				;NCE .NE. NC)

;	this is the current model
ncom=n_elements(pos)
if ncom gt 0 then begin
  plist=strarr(4L*ncom)+string(' ','(a100)') & qlist=intarr(4L*ncom) & k=-1L
endif else begin
  plist=[string(' ','(a100)')] & qlist=[0]
endelse

;	initialze some global parameters
csize=3 & nsig=3.			;for LINEREM
dpos=0.1 & dwdt=0.1 & dflx=5.		;width of priors on POS,WDT,FLX
					;dFLX is in units of FERRP@POS
label_top='[Thaw?] :  Component : Parameter : Tag : Value [+- Error [@dCHI]]'
for i=0,ncom-1 do begin
  k=k+1 & k0=k
  plist(k)=string(i,'(i5)')+string(type(i),'(a18)')+' '+$
	string(epithet(i),'(a35)')
  for j=0,2 do begin
    jpar=3*i+j & k=k+1 & plist(k)='     '
    ;if thaw(jpar) eq 0 then plist(k)=' fix ' else plist(k)=' var '
    if j eq 0 then plist(k)=plist(k)+'     '+string('position','(a10)')
    if j eq 1 then plist(k)=plist(k)+'     '+string('width','(a10)')
    if j eq 2 then plist(k)=plist(k)+'     '+string('flux','(a10)')
    c='a'+strtrim(jpar,2)
    plist(k)=plist(k)+string(c,'(a7)')
    if j eq 0 then val=pos(i)
    if j eq 1 then val=wdt(i)
    if j eq 2 then val=flx(i)
    plist(k)=plist(k)+string(val,'(g14.4)')+string(' ','(a30)')
    if thaw(jpar) ne 0 then begin
      qlist[k]=1 & qlist[k0]=1
    endif
  endfor
endfor

;	define the widget base heirarchy to control the procedure
title='MEASURE LINE FLUXES'
base=widget_base(frame=2,/column,/map,title=title,uvalue='base')

;	first row contains:
;		button for FIT
;		drop down menu for VIEW
;		drop down menu for XRANGE, YRANGE
;		drop down menu for MODEL, CONTINUUM
;		entry fields for CELLSIZE, THRESHOLD
;		button for QUIT
base_menu=widget_base(base,frame=2,/row,uvalue='base_menu')
;
menu_fits=widget_button(base_menu,value='FIT',uvalue='menu_fits')
;
menu_view=widget_button(base_menu,value='VIEW',/menu,uvalue='menu_view')
view_show=widget_button(menu_view,value='SHOW',uvalue='view_show')
view_hcpy=widget_button(menu_view,value='Hardcopy',uvalue='view_hcpy')
view_colr=widget_button(menu_view,value='ColorTables',/menu,uvalue='view_colr')
  colr_ctpb=widget_button(view_colr,value='BlackPease',uvalue='colr_ctpb')
  colr_ctpw=widget_button(view_colr,value='WhitePease',uvalue='colr_ctpw')
  colr_idlt=widget_button(view_colr,value='LOADCT,#',uvalue='colr_idlt')
view_defs=widget_button(menu_view,value='Set Defaults',uvalue='view_defs')
;
menu_xrng=widget_button(base_menu,value='XRANGE',/menu,uvalue='menu_xrng')
xrng_crsr=widget_button(menu_xrng,value='Cursor',uvalue='xrng_crsr')
xrng_kbrd=widget_button(menu_xrng,value='Keyboard',uvalue='xrng_kbrd')
xrng_auto=widget_button(menu_xrng,value='Guess',uvalue='xrng_auto')
xrng_zero=widget_button(menu_xrng,value='RESET',uvalue='xrng_zero')
;
menu_yrng=widget_button(base_menu,value='YRANGE',/menu,uvalue='menu_yrng')
yrng_crnt=widget_button(menu_yrng,value='Current',uvalue='yrng_crnt')
yrng_kbrd=widget_button(menu_yrng,value='Keyboard',uvalue='yrng_kbrd')
yrng_zero=widget_button(menu_yrng,value='RESET',uvalue='yrng_zero')
;
mods_modl=widget_button(base_menu,value='MODEL',/menu,uvalue='mods_modl')
modl_friz=widget_button(mods_modl,value='Freeze All',uvalue='modl_friz')
modl_thaw=widget_button(mods_modl,value='Thaw All',uvalue='modl_thaw')
modl_fixx=widget_button(mods_modl,value='Freeze Unimportant',uvalue='modl_fixx')
modl_addc=widget_button(mods_modl,value='ADD Component',uvalue='modl_addc')
modl_delc=widget_button(mods_modl,value='DELETE Component',uvalue='modl_delc')
;
mods_cont=widget_button(base_menu,value='CONTINUUM',/menu,uvalue='mods_cont')
cont_fcmp=widget_button(mods_cont,value='Remove Lines',uvalue='cont_fcmp')
cont_rcmp=widget_button(mods_cont,value='LineRem (sub)',uvalue='cont_rcmp')
cont_pcws=widget_button(mods_cont,value='Piecewise',uvalue='cont_pcws')
cont_norm=widget_button(mods_cont,value='Adjust Norm',uvalue='cont_norm')
cont_acpt=widget_button(mods_cont,value='ACCEPT',uvalue='cont_acpt')
cont_orig=widget_button(mods_cont,value='RESET',uvalue='cont_orig')
;
cont_cell=cw_field(base_menu,/column,/long,/return_events,$
	title='width [pix]',xsize=5,value=csize,$
	uvalue='cont_cell')
cont_nsig=cw_field(base_menu,/column,/float,/return_events,$
	title='threshold',xsize=5,value=nsig,uvalue='cont_nsig')
;
menu_quit=widget_button(base_menu,value='QUIT',uvalue='menu_quit')

b_menu = [base_menu,menu_fits,menu_view,view_show,view_hcpy,$
	view_colr,colr_ctpb,colr_ctpw,colr_idlt,$
	view_defs,$
	menu_xrng,xrng_crsr,xrng_kbrd,xrng_auto,xrng_zero,$
	menu_yrng,yrng_crnt,yrng_kbrd,yrng_zero,$
	mods_modl,modl_friz,modl_thaw,modl_fixx,modl_addc,modl_delc,$
	mods_cont,cont_rcmp,cont_fcmp,cont_pcws,cont_norm,cont_acpt,$
	  cont_orig,cont_cell,cont_nsig,$
	menu_quit]

;	second row contains:
;		text STATUS
base_txts=widget_base(base,frame=2,/row,uvalue='base_txts',/align_left)
;
txts_stts=widget_text(base_txts,font='10x20',frame=2,/scroll,xsize=75,$
	value='Widget-based line fitting IDL program',uvalue='txts_stts')

b_txts = [base_txts,txts_stts]

;	third row contains:
;		button for HELP
;		button for undo fit, reset original, model renorm
;		droplist for fitting algorithm
;		button for projection-error calcs
;		entry field for delta-chisq threshold
base_fits=widget_base(base,frame=2,/row,uvalue='base_fits')
;
algval=['LevMarq','Curvefit','MPFIT','tifsrhg']
fits_help=widget_button(base_fits,value='HELP',uvalue='fits_help')
fits_norm=widget_button(base_fits,value='RENORM',uvalue='fits_norm')
fits_undo=widget_button(base_fits,value='UNDO FIT',uvalue='fits_undo')
fits_zero=widget_button(base_fits,value='SAVE',uvalue='fits_zero')
fits_orig=widget_button(base_fits,value='Reset Original',uvalue='fits_orig')
fits_dump=widget_button(base_fits,value='DUMP',uvalue='fits_dump')
menu_eror=widget_button(base_fits,value='ERRORS',/menu, uvalue='menu_eror')
  fts_steppar=widget_button(menu_eror, value='StepPar',uvalue='fts_steppar')
  fts_mcerror=widget_button(menu_eror, value='Monte-Carlo',uvalue='fts_mcerror')
fits_dchi=cw_field(base_fits,/row,/float,/return_events,$
	title='dCHI',xsize=5,value=dchisq,uvalue='fits_dchi')
fits_algo=widget_droplist(base_fits,value=algval,uvalue='fits_algo')

b_fits=[base_fits,fits_help,fits_norm,fits_zero,fits_undo,fits_orig,$
	fits_dump,menu_eror,fits_dchi,fits_algo]

;	fourth row contains:
;		entry fields to edit parameter values of specified component
;		list of the model parameters including checkboxes to freeze/thaw
base_pars=widget_base(base,frame=2,/row,uvalue='base_pars')
;
pars_comp=widget_base(base_pars,/column,uvalue='pars_comp')
comp_indx=cw_field(pars_comp,/column,/long,/all_events,title='Component',$
	value=-1L,uvalue='comp_indx')
comp_posv=cw_field(pars_comp,/column,/float,/return_events,title='Position',$
	value=0.,uvalue='comp_posv')
comp_wdtv=cw_field(pars_comp,/column,/float,/return_events,title='Width',$
	value=0.,uvalue='comp_wdtv')
comp_flxv=cw_field(pars_comp,/column,/float,/return_events,title='Flux',$
	value=0.,uvalue='comp_flxv')
comp_type=cw_field(pars_comp,/column,/string,/return_events,title='Type',$
	value=type(0),uvalue='comp_type')
comp_sign=cw_field(pars_comp,/column,/string,/return_events,title='Label',$
	value=epithet(0),uvalue='comp_sign')
;
bgp_ids=1
pars_list=cw_bgroup(base_pars,plist,/column,/nonexclusive,uvalue='pars_list',$
	/return_index,/scroll,y_scroll_size=340,x_scroll_size=540,$
	set_value=qlist,font='8x13',label_top=label_top,ids=list_ids)

b_pars = [base_pars,pars_list,pars_comp,comp_indx,comp_posv,comp_wdtv,$
	comp_flxv,comp_type,comp_sign,list_ids]

;	fifth row contains:
;		drop down menu for TIES, DELETE
;		entry fields for EDIT, ADD
base_ties=widget_base(base,frame=2,/row,uvalue='base_ties')
;
ties_acts=widget_button(base_ties,value='TIES',/menu,uvalue='ties_acts')
acts_frcp=widget_button(ties_acts,value='Positions',uvalue='acts_frcp')
acts_frcw=widget_button(ties_acts,value='Widths',uvalue='acts_frcw')
acts_frcf=widget_button(ties_acts,value='Fluxes',uvalue='acts_frcf')
;
ties_dele=widget_button(base_ties,value='DELETE',/menu,uvalue='ties_dele')
dele_some=widget_button(ties_dele,value='Selected',uvalue='dele_some')
dele_alls=widget_button(ties_dele,value='Clear All',uvalue='dele_alls')
;
flds_ties=widget_base(base_ties,/column,uvalue='flds_ties')
ties_edit=cw_field(flds_ties,/row,/string,/return_events,title='EDIT',$
	xsize=60,value='',uvalue='ties_edit')
ties_adds=cw_field(flds_ties,/row,/string,/return_events,title='ADD',$
	xsize=60,value='',uvalue='ties_adds')

b_ties = [base_ties,ties_acts,acts_frcp,acts_frcw,acts_frcf,ties_dele,$
	dele_some,dele_alls,flds_ties,ties_edit,ties_adds]

;	sixth row contains:
;		entry fields for dPOS, dWDT, dFLX in left column
;		list of all the existing ties in right column
base_dtls=widget_base(base,frame=2,/row,uvalue='base_dtls')
dtls_ddel=widget_base(base_dtls,frame=1,/column,uvalue='dtls_ddel')
ties_dpos=cw_field(dtls_ddel,/row,/float,/return_events,title='dPOS',$
	value=dpos,uvalue='ties_dpos')
ties_dwdt=cw_field(dtls_ddel,/row,/float,/return_events,title='dWDT',$
	value=dwdt,uvalue='ties_dwdt')
ties_dflx=cw_field(dtls_ddel,/row,/float,/return_events,title='dFLX',$
	value=dflx,uvalue='ties_dflx')
;
if not keyword_set(ties) then begin
  ties=0 & tt=''
endif else begin
  if strtrim(ties(0),2) eq '0' then begin
    ties=0 & tt=''
  endif else tt=ties
endelse
dtls_tlst=widget_base(base_dtls,frame=2,/row,uvalue='dtls_tlst')
;tlst_list=widget_list(dtls_tlst,/multiple,scr_xsize=500,scr_ysize=115,$
;	value=tt,uvalue='tlst_list')
;Antonio Maggio points out that keyword MULTIPLE is not allowed in
;older versions.  removing it causes little pain.
tlst_list=widget_list(dtls_tlst,scr_xsize=500,scr_ysize=115,$
	value=tt,uvalue='tlst_list')

b_dtls=[base_dtls,dtls_ddel,ties_dpos,ties_dwdt,ties_dflx,$
	dtls_tlst,tlst_list]

;	the help section
fit_help='fitlines.hlp'

;	this is how the widget IDs are transferred to FITLINES_EVENT
widg={base:base, base_menu:b_menu, base_txts:b_txts, base_fits: b_fits,$
	base_pars:b_pars, base_ties:b_ties, base_dtls:b_dtls,$
	fit_help:fit_help }

;	start up the widgets
widget_control,base,/realize
widget_control,/hourglass

;register with XMANAGER
;	this is the officially recommended way, but we're not going to do it
;	because there are too many parameters and keywords (not to mention
;	the _EXTRA keywords) to pass to the event handling routine.
;		xmanager,'fitlines',base
;therefore,
;	we handle the widget events in a separate subourtine, partly also
;	to keep things easy to read

fitlines_event,x,y,pos,wdt,flx,perrp,werrp,ferrp,thaw,type,ties,epithet,$
	conlev,consig,chisq,funcs,sigy,widg,$
	perrm=perrm,werrm=werrm,ferrm=ferrm,rmfstr=rmfstr,$
	perrc=perrc,werrc=werrc,ferrc=ferrc, intens=intens,$
	histerix=histerix,$
	_extra=e

;	clean up
widget_control,base,/destroy

;	output
;	compress CONLEV and CONSIG
;	ok, I know this is incomplete, and YTOL,NOBAND,DISCON need
;	to be properly handled, especially the last.
if not keyword_set(ytol) then ytol=-0.01	;correct to 1%
if n_elements(uniq(conlev)) gt 1 then begin
  splac,x,conlev,ytol,u1,v1,w1,/verbose
  splac,x,consig,ytol,u2,v2,w2,/verbose
  u1=x & v1=conlev
  u2=x & v2=consig
endif else begin
  umin=min(x,max=umax) & u1=[umin,umax] & u2=u1
  v1=conlev(0)*[1,1] & v2=consig(0)*[1,1]
endelse
;
comment='<chisq='+strtrim(chisq,2)+'>'
ntag=n_tags(e) & if ntag gt 0 then t=tag_names(e) else t=''
for i=0,ntag-1 do comment=comment+' <'+t(i)+'='+strtrim((e.(i))(0),2)+'>'
ntie=n_elements(ties)
if ntie gt 0 then begin
  if strtrim(ties(0),2) ne '0' then for i=0,ntie-1 do comment=$
	comment+' {'+ties(i)+'}'
endif
;
if n_elements(pos) gt 0 then ostr=$
	{POS:pos, PERRP:perrp, PERRM:perrm, PERRC:perrc,$
	 FLX:flx, FERRP:ferrp, FERRM:ferrm, FERRC:ferrc,$
	 WDT:wdt, WERRP:werrp, WERRM:werrm, WERRC:werrc,$
	 THAW:thaw, TYPE:type, TIES:ties, EPITHET:epithet,$
	 CONLEVX:u1, CONLEV:v1, CONSIGX:u2, CONSIG:v2, COMMENT:comment}

return,ostr
end
