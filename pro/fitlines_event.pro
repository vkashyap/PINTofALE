;+
;procedure	fitlines_event
;	widget event handler subroutine for FITLINES.  see that routine
;	for a description of variables, etc.  only rudimentary consistency
;	checks are carried out here.  uses subroutine FMTSTUFF, which is
;	included in this file.
;
;syntax
;	fitlines_event,x,y,pos,wdt,flx,perr,werr,ferr,thaw,type,ties,epithet,$
;	conlev,consig,chisq,funcs,sigy,widg, perrm=perrm,werrm=werrm,$
;	ferrm=ferrm,perrc=perrc,werrc=werrc,ferrc=ferrc,/intens,rmfstr=rmfstr,$
;	histerix=histerix
;
;subroutines
;	FMTSTUFF
;	PICKRANGE
;	FUNCS
;	(X3MODEL/LIBMODEL)
;	(MK_3MODEL, MK_GAUSS, MK_LORENTZ, MK_SLANT, MK_ROGAUSS, MK_POWLAM, MK_POLY1D)
;	LINEREM
;	SETCONT
;	ERORS
;	FIT_LEVMAR , CURVE_FIT , MPFITFUN
;	ADJUSTIE
;	LEVMARQ
;	LMCOEFF
;	STAMPLE
;	KILROY
;	WHEE
;	PEASECOLR
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Oct98)
;	added renormalization option (VK; Nov98)
;	various bug fixes and reorganizations (VK; Dec98)
;	now allows fitting a single parameter (VK; Aug99)
;	enhancements to include call to ERORS, setting YRANGEs,
;	  INTENSity units, etc. (VK; FebMM)
;	avoid repeat calls to CW_BGROUP, saving tremendous headaches;
;	  bug correction CONTINUUM->AdjustNorm (VK; MarMM)
;	bug correction -- undefined STYPE for new component; rudimentary
;	  changing parameter labels (VK; JunMM)
;	moved QUIT to end of row, far away from FIT; added DUMP to disk
;	  (VK; JulMM)
;	added error-bar plotting capability via VIEW->SetDefaults;
;	  pass x,y-range to setcont; added labeling capability via
;	  keyword EPITHET; merged 6th and 7th rows; bug correction
;	  with xcont in adjust norm (VK; AugMM)
;	changed location of help file to ardb; added calls to STAMPLE
;	  (VK; DecMM)
;	bug fixes: invisible error output, adding model deleted existing
;	  projected errors, streamlined saves (VK; JanMMI)
;	per Antonio Maggio's suggestion, removed keyword MULTIPLE from
;	  call to WIDGET_LIST; improved color-scale setting for 24-bit
;	  consoles; replaced call to WHICH by call to ROUTINE_INFO
;	  (VK; FebMMI)
;	various extra info messages (VK; AprMMI)
;	changed y to z on lines 1012 and 1285, because.. y was
;	  inappropriate for some reason? (VK; Jun01)
;	was crashing due to missing ZSIG in cont_acpt (VK; Aug01)
;	bug was not passing predefined CONLEV to SETCONT (piecewise) (VK; Apr02)
;	handle color tables in 24-bit displays (VK; Jun02)
;	bug correction -- if all ties were deleted, then TIES was turning into
;	  and integer array (VK; Jul02)
;	made changes to plotting so that hardcopy won't come out with annoying
;	  blank sheets (VK; Sep02)
;	changed default colors to be compatible with PEASECOLOR (VK);
;	  added hooks into MPFIT (Liwei Lin/VK; Oct02)
;	tied the left column of parameter values to changes in the main parameter
;	  list window on the right (VK; Apr03)
;	added monte-carlo errors as option and added keyword RMFSTR (LL; Aug03)
;	added keyword HISTERIX (VK; Jul03)
;	bug correction: algo_type undefined unless fit is run first (LL; Aug04)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	modified some widget instructions to be slightly clearer (VK; Aug06)
;
;*********************************************************************
;
;subroutine	fmtstuff
;	takes in the model as described, and returns it properly formatted
;	for display in the widget.  parameters are named as before.  only
;	minimal consistency checks are carried out.
;usage
;	fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
;	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
;	  perrc=perrc,werrc=werrc,ferrc=ferrc,$
;	  plist=plist,list=list,$
;	  tiep=tiep,tiew=tiew,tief=tief,ties=ties,delcmp=delcmp
;
;keywords
;	plist 	[OUTPUT] appropriately formatted model parameter list
;	list	[OUTPUT] denoting frozen/thawed components
;	tiep	[I/O] ties specially setting the range on positions
;	tiew	[I/O] ties specially setting the range on widths
;	tief	[I/O] ties specially setting the range on fluxes
;	ties	[I/O] rest of ties
;	delcmp	[ACTION] if set to component number, deletes component DELCOMP
;		and all TIES that contain references to parameters in that
;		component
;	_extra	[INPUT] pass defined keywords to STAMPLE (NODAY,NOTIME,
;		NOUSER,NOPACK,STACOL)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;{
pro fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	perrm=perrm,werrm=werrm,ferrm=ferrm,$
	perrc=perrc,werrc=werrc,ferrc=ferrc,$
	plist=plist,list=list,$
	tiep=tiep,tiew=tiew,tief=tief,ties=ties,delcmp=delcmp,$
	_extra=e
;(plus)
;procedure	fmtstuff
;		see above for description
;(minus)

if n_elements(delcmp) eq 1 then begin	;(delete specified component
  idel=delcmp(0)
  ncom=n_elements(pos)
  if ncom gt idel and idel ge 0 then begin	;(something exists to delete
    ii=lindgen(ncom) & oi=where(ii ne idel,moi)
    if moi eq 0 then begin
      print,string("7b)
      message,'	Cannot delete the last component',/info
      message,'	First add another component and then delete this',/info
      return
    endif
    ;	delete parameters
    pos=pos(oi) & perr=perr(oi)
    wdt=wdt(oi) & werr=werr(oi)
    flx=flx(oi) & ferr=ferr(oi)
    type=type(oi) & thaw=thaw([3*oi,3*oi+1,3*oi+2])
    epithet=epithet(oi)
    ;	change the numbers on all special ties past the current one
    ntp=n_elements(tiep) & ntw=n_elements(tiew) & ntf=n_elements(tief)
    for i=idel+1,ncom-1 do begin
      c_old='a'+strtrim(i,2) & c_new='a'+strtrim(i-1,2)
      while strlen(c_new) lt strlen(c_old) do c_new=c_new+' '
      for j=0,2 do begin
	cc=''
	if j eq 0 and ntp eq ncom then cc=tiep(i)
	if j eq 0 and ntw eq ncom then cc=tiew(i)
	if j eq 0 and ntf eq ncom then cc=tief(i)
        i0=-1L & j0=0
        while j0 ge 0 do begin
	  j0=strpos(cc,c_old,i0+1)
	  if j0 ge 0 then strput,cc,c_new,j0
        endwhile
      endfor
    endfor
    ;	delete from special ties
    if n_elements(tiep) eq ncom then tiep=tiep(oi)
    if n_elements(tiew) eq ncom then tiew=tiew(oi)
    if n_elements(tief) eq ncom then tief=tief(oi)
    ;	and look through the rest of the ties to check for potential trouble
    nt=n_elements(ties)
    if nt gt 0 then begin
      ;	possible things to look for..
      cs=[' ','=','(',')','<','>','+','-','*','/'] & ncs=n_elements(cs)
      cc='a'+strtrim(idel,2)+[' ','=','(',')','<','>','+','-','*','/']
      cd='a'+strtrim(idel-1,2)
      if strlen(cc(0)) gt strlen(cd(0)) then cd=cd+' '+$
	[' ','=','(',')','<','>','+','-','*','/']
      ncc=n_elements(cc)
      aiee=intarr(nt)
      for i=0,nt-1 do begin		;{edit each tie
	ct=strlowcase(ties(i)) & ctt=ct
	for j=idel,ncom-1 do begin	;{for each component
	  ;	all parameters after IDEL drop in number by 3
	  for jj=0,2 do begin		;{3 parameters per component
	    cc='a'+strtrim(3*j+jj,2)
	    if j gt idel then begin	;(replace with new A#
	      cd='a'+strtrim(3*j+jj-3,2)
	      if strlen(cc) gt strlen(cd) then cd=cd+' '
	      for k=0,ncs-1 do begin	;{check each possibility
	        ik=0L
	        while ik ge 0 do begin	;{find all occurances and replace
	          ikold=ik
	          ik=strpos(ctt,cc+cs(k),ikold)
	          if ik ge 0 then strput,ctt,cd+cs(k),ik
	        endwhile			;IK>0}
	      endfor			;K=0,NCS-1}
	    endif else begin		;)(throw away
	      for k=0,ncs-1 do begin	;{check each possibility
	        ik=strpos(ctt,cc+cs(k),0)	;need only one hit
	        if ik ge 0 then aiee(i)=1
	      endfor			;K=0,NCS-1}
	      if aiee(i) eq 1 then ctt=''
	    endelse			;J)
	  endfor			;JJ=0,2}
	endfor				;J=IDEL,NCOM-1}
	if ct ne ctt then begin
	  print,'OLD: '+ct & print,'NEW: '+ctt
	  ties(i)=ctt
	endif
      endfor				;I=0,NT-1}
      oy=where(aiee eq 0,moy)
      if moy gt 0 then ties=ties(oy)
    endif
  endif						;NCOM>IDEL)
endif					;N(DELCMP)=1)


ncom=n_elements(pos)
if ncom gt 0 then begin
  plist=strarr(4L*ncom)+string(' ','(a100)') & list=intarr(4L*ncom) & k=-1L
  if n_elements(perr) ne ncom then perr=fltarr(ncom)
  if n_elements(perrm) ne ncom then perrm=perr
  if n_elements(perrc) ne ncom then perrc=fltarr(ncom)
  if n_elements(werr) ne ncom then werr=fltarr(ncom)
  if n_elements(werrm) ne ncom then werrm=werr
  if n_elements(werrc) ne ncom then werrc=fltarr(ncom)
  if n_elements(ferr) ne ncom then ferr=fltarr(ncom)
  if n_elements(ferrm) ne ncom then ferrm=ferr
  if n_elements(ferrc) ne ncom then ferrc=fltarr(ncom)
endif else begin
  plist=[string(' ','(a100)')] & list=[0]
endelse
for i=0,ncom-1 do begin
  k=k+1 & k0=k
  plist(k)=string(i,'(i5)')+string(type(i),'(a18)')+' '+$
	string(epithet(i),'(a35)')
  if strlowcase(strmid(strtrim(type(i),2),0,1)) eq 'p' then begin
    pnam1='index1' & pnam2='index2' & pnam3='flux'
  endif else begin
    pnam1='position' & pnam2='width' & pnam3='flux'
  endelse
  for j=0,2 do begin
    jpar=3*i+j & k=k+1 & plist(k)='     ' & c='a'+strtrim(jpar,2)
    if j eq 0 then begin
      plist(k)=plist(k)+'     '+string(pnam1,'(a10)')
      val=pos(i) & valu=perr(i) & vall=perrm(i) & valc=perrc(i)
    endif
    if j eq 1 then begin
      plist(k)=plist(k)+'     '+string(pnam2,'(a10)')
      val=wdt(i) & valu=werr(i) & vall=werrm(i) & valc=werrc(i)
    endif
    if j eq 2 then begin
      plist(k)=plist(k)+'     '+string(pnam3,'(a10)')
      val=flx(i) & valu=ferr(i) & vall=ferrm(i) & valc=ferrc(i)
    endif
    plist(k)=plist(k)+string(c,'(a7)')
    plist(k)=plist(k)+string(val,'(g11.5)')
    if valc eq 0 then begin
      plist(k)=plist(k)+' +- '+string(valu,'(g10.5)')+string(' ','(a20)')
    endif else begin
      plist(k)=plist(k)+' +'+string(valu,'(g10.5)')+$
	' -'+string(vall,'(g10.5)')+' @'+string(valc,'(f5.3)')
    endelse
    if thaw(jpar) ne 0 then begin
      list(k)=1 & list(k0)=1
    endif
  endfor
endfor

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;}
pro fitlines_event,x,y,pos,wdt,flx,perr,werr,ferr,thaw,type,ties,epithet,$
	conlev,consig,chisq,funcs,sigy,widg,$
	histerix=histerix,$
	perrm=perrm,werrm=werrm,ferrm=ferrm, rmfstr=rmfstr,$
	perrc=perrc,werrc=werrc,ferrc=ferrc, intens=intens, _extra=e
;(plus)
;procedure	fitlines_event
;		see above for description
;(minus)

forward_function mpfitfun

;	some initializations
nx=n_elements(x) & oo=lindgen(nx)		;default selection
oxr=[min(x),max(x)] & oyr=[min(y),max(y)]
chisq=-1 & icomp=-1L & ktie=0L

;create keyword x_seg if rmfstr is set 
 if keyword_set(rmfstr) then begin 
    xx = x(oo)
    if not keyword_set(rmfgrid) then rmfgrid = (rmfstr.elo+rmfstr.ehi)/2 
    jnk1 = min( abs(max(xx)-rmfgrid), mx) 
    jnk2 = min( abs(min(xx)-rmfgrid), mn)
    x_seg = rmfgrid(mn<mx:mn>mx) 
 endif 

;
;ccon=100	;color for continuum
;cmod=150	;color for total model
;ccmp=200	;color for model components
;cspc=250	;color for spectrum
ccon=1	;color for continuum
cmod=2	;color for total model
ccmp=3	;color for model components
cspc=6	;color for spectrum
dncolors=256. > !D.N_COLORS	;24-bit color screen temporary fix
iccon=ccon ;	iccon=fix(ccon*float(!d.n_colors)/dncolors+1)
icmod=cmod ;	icmod=fix(cmod*float(!d.n_colors)/dncolors+1)
iccmp=ccmp ;	iccmp=fix(ccmp*float(!d.n_colors)/dncolors+1)
icspc=cspc ;	icspc=fix(cspc*float(!d.n_colors)/dncolors+1)
tcon=1. & tmod=1. & tcmp=1. & tspc=1.
			;thickness for continuum, model, components, spectrum
xtitle='Wavelength ['+string(byte(197))+']'
ytitle='Counts' & if keyword_set(intens) then ytitle='Intensity'
title=''
dwvl=1.

if keyword_set(intens) then znorm=0 else znorm=1
zerr=1		;to plot error bars, ZERR sets the width in units of sigma
zlab=1		;size of labels. to not plot labels, use ZLAB=0

;	save original values
ncom=n_elements(pos)
if n_elements(perrm) ne ncom then perrm=perr
if n_elements(werrm) ne ncom then werrm=werr
if n_elements(ferrm) ne ncom then ferrm=ferr
if n_elements(perrc) ne ncom then perrc=fltarr(ncom)
if n_elements(werrc) ne ncom then werrc=fltarr(ncom)
if n_elements(ferrc) ne ncom then ferrc=fltarr(ncom)
if n_elements(epithet) ne ncom then epithet=strarr(ncom)
if ncom gt 0 then begin
  opos=pos & owdt=wdt & oflx=flx
  operr=perr & owerr=werr & oferr=ferr
  operrm=perrm & owerrm=werrm & oferrm=ferrm
  operrc=perrc & owerrc=werrc & oferrc=ferrc
  othaw=thaw & otype=type & oepithet=epithet
endif
oconlev=conlev & oconsig=consig
if n_elements(ties) gt 1 then oties=ties else oties=['']

;	disentangle WIDG
base=widg.base
k=0+0 & base_menu=widg.base_menu(k)	;first row
k=k+1 & menu_fits=widg.base_menu(k)
k=k+1 & menu_view=widg.base_menu(k)
k=k+1 & view_show=widg.base_menu(k)
k=k+1 & view_hcpy=widg.base_menu(k)
k=k+1 & view_colr=widg.base_menu(k)
k=k+1 & colr_ctpb=widg.base_menu(k)
k=k+1 & colr_ctpw=widg.base_menu(k)
k=k+1 & colr_idlt=widg.base_menu(k)
k=k+1 & view_defs=widg.base_menu(k)
k=k+1 & menu_xrng=widg.base_menu(k)
k=k+1 & xrng_crsr=widg.base_menu(k)
k=k+1 & xrng_kbrd=widg.base_menu(k)
k=k+1 & xrng_auto=widg.base_menu(k)
k=k+1 & xrng_zero=widg.base_menu(k)
k=k+1 & menu_yrng=widg.base_menu(k)
k=k+1 & yrng_crnt=widg.base_menu(k)
k=k+1 & yrng_kbrd=widg.base_menu(k)
k=k+1 & yrng_zero=widg.base_menu(k)
k=k+1 & mods_modl=widg.base_menu(k)
k=k+1 & modl_friz=widg.base_menu(k)
k=k+1 & modl_thaw=widg.base_menu(k)
k=k+1 & modl_fixx=widg.base_menu(k)
k=k+1 & modl_addc=widg.base_menu(k)
k=k+1 & modl_delc=widg.base_menu(k)
k=k+1 & mods_cont=widg.base_menu(k)
k=k+1 & cont_rcmp=widg.base_menu(k)
k=k+1 & cont_fcmp=widg.base_menu(k)
k=k+1 & cont_pcws=widg.base_menu(k)
k=k+1 & cont_norm=widg.base_menu(k)
k=k+1 & cont_acpt=widg.base_menu(k)
k=k+1 & cont_orig=widg.base_menu(k)
k=k+1 & cont_cell=widg.base_menu(k)
k=k+1 & cont_nsig=widg.base_menu(k)
k=k+1 & menu_quit=widg.base_menu(k)
k=0+0 & base_txts=widg.base_txts(k)	;second row
k=k+1 & txts_stts=widg.base_txts(k)
k=0+0 & base_fits=widg.base_fits(k)	;third row
k=k+1 & fits_help=widg.base_fits(k)
k=k+1 & fits_norm=widg.base_fits(k)
k=k+1 & fits_zero=widg.base_fits(k)
k=k+1 & fits_undo=widg.base_fits(k)
k=k+1 & fits_orig=widg.base_fits(k)
k=k+1 & fits_dump=widg.base_fits(k)
k=k+1 & menu_eror=widg.base_fits(k)
k=k+1 & fits_dchi=widg.base_fits(k)
k=k+1 & fits_algo=widg.base_fits(k)
k=0+0 & base_pars=widg.base_pars(k)	;fourth row
k=k+1 & pars_list=widg.base_pars(k)
k=k+1 & pars_comp=widg.base_pars(k)
k=k+1 & comp_indx=widg.base_pars(k)
k=k+1 & comp_posv=widg.base_pars(k)
k=k+1 & comp_wdtv=widg.base_pars(k)
k=k+1 & comp_flxv=widg.base_pars(k)
k=k+1 & comp_type=widg.base_pars(k)
k=k+1 & comp_sign=widg.base_pars(k)
k=k+1 & list_ids =widg.base_pars(k:*)
k=0+0 & base_ties=widg.base_ties(k)	;fifth row
k=k+1 & ties_acts=widg.base_ties(k)
k=k+1 & acts_frcp=widg.base_ties(k)
k=k+1 & acts_frcw=widg.base_ties(k)
k=k+1 & acts_frcf=widg.base_ties(k)
k=k+1 & ties_dele=widg.base_ties(k)
k=k+1 & dele_some=widg.base_ties(k)
k=k+1 & dele_alls=widg.base_ties(k)
k=k+1 & flds_ties=widg.base_ties(k)
k=k+1 & ties_edit=widg.base_ties(k)
k=k+1 & ties_adds=widg.base_ties(k)
k=0+0 & base_dtls=widg.base_dtls(k)	;sixth row
k=k+1 & dtls_ddel=widg.base_dtls(k)
k=k+1 & ties_dpos=widg.base_dtls(k)
k=k+1 & ties_dwdt=widg.base_dtls(k)
k=k+1 & ties_dflx=widg.base_dtls(k)
k=k+1 & dtls_tlst=widg.base_dtls(k)
k=k+1 & tlst_list=widg.base_dtls(k)
;
fit_help=widg.fit_help & nhelp=n_elements(fit_help)	;help

;to begin with, call FMTSTUFF just to make sure everything's OK
fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	perrm=perrm,werrm=werrm,ferrm=ferrm,$
	perrc=perrc,werrc=werrc,ferrc=ferrc,$
	plist=plist,list=list & updatelist=1

timestamp=systime()
histr=create_struct('FITLINES',timestamp)
if n_elements(histerix) gt 0 then begin
  if not keyword_set(histerix) then histr=0
endif
ihist=0L

ok=''
while ok ne 'quit' do begin			;{endless loop
  
  ;catch a falling event
  event=widget_event(base) & widget_control,event.id,get_uvalue=ev

  case ev of					;{handle the events

    ;							first row
    'menu_fits': begin				;{fit
      ;		define line spectrum
      widget_control,txts_stts,set_value='Fit 3-component model to data..'
      z=y-conlev & xx=x(oo) & zz=z(oo) & zsig=sqrt(sigy(oo)^2+consig(oo)^2)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      widget_control,pars_list,get_value=list
      ncom=n_elements(pos)
      if ncom gt 0 then begin		;(if there are models to fit..
        thaw=intarr(3*ncom) & ii=lindgen(ncom)
        thaw(3*ii+0)=list(4*ii+1)		;positions
        thaw(3*ii+1)=list(4*ii+2)		;widths
        thaw(3*ii+2)=list(4*ii+3)		;fluxes
	freeze=where(thaw eq 0,nfrozen)
	aa=fltarr(3L*ncom)		;define the parameter array
	aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	if not keyword_set(alltie) then begin	;define prevailing constraints
	  if keyword_set(ties) then alltie=ties else alltie=''
	endif
	;
	nfree=3*ncom & nfree=nfree-nfrozen
	;if nfree gt 0 then begin
	  ;	save the old values
	  spos=pos & swdt=wdt & sflx=flx*dwvl
	  sperr=perr & swerr=werr & sferr=ferr*dwvl
	  sperrm=perrm & swerrm=werrm & sferrm=ferrm
	  sperrc=perrc & swerrc=werrc & sferrc=ferrc
	  stype=type & sthaw=thaw & sepithet=epithet
	  ;
	  ;	call fitting program
	  zwt=zsig & ozw=where(zwt gt 0,mozw)
	  if mozw gt 0 then zwt(ozw)=1./zwt(ozw)
	  if not keyword_set(algo_type) then algo_type='LevMarq+SVD'	;default
	  if algo_type eq 'LevMarq+SVD' then fit_levmar,$
		xx,zz,aa,freeze=freeze,erra=erra,chisq=chisq,ties=alltie,$
		funcs=funcs,sig=zsig,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr, _extra=e
	  if algo_type eq 'IDL-Curvefit' then tmp=curve_fit(xx,zz,zwt,aa,erra,$
		chisq=chisq,function_name=funcs,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr,$
		_extra=e)
	  if algo_type eq 'MPFIT' then begin
	    ;freeze is where thaw=0
	    PARINFO=replicate({fixed:0,limited:[0,0],limits:[0,0]},n_elements(aa))
	    ;freeze it
	    for j=0,nfrozen-1 do parinfo(freeze(j)).fixed(0)=1
	    userfuncs=funcs
	    functarg=create_struct(e,'normflx',znorm,'type',type)
	    if strpos(userfuncs,'_f',0) lt 0 then userfuncs=userfuncs+'_f'
	    results=MPFITFUN(userfuncs,xx,zz,zsig,aa,PERROR=erra,yfit=yfit,$
	    	parinfo=parinfo,functarg=functarg, _extra=e)
	    ;chisq=results(n_elements(results)-1L)
	    ;aa=results[0:n_elements(results)-2L]
	    aa=results
	    chisq=total( ((yfit-zz)^2)/((zsig)^2) )
	  endif
	  ;
	  ;	update the errors as necessary
	  for i=0L,ncom-1L do begin	;{for each component
	    jjp=0 & jjw=0 & jjf=0	;do nothing
	    if thaw(3*i+0) gt 0 then jjp=1 else if erra(3*i+0) ne 0 then jjp=1
	    if thaw(3*i+1) gt 0 then jjw=1 else if erra(3*i+1) ne 0 then jjw=1
	    if thaw(3*i+2) gt 0 then jjf=1 else if erra(3*i+2) ne 0 then jjf=1
	    pos(i)=aa(3L*i+0)
	    if jjp gt 0 then begin
	      if perrc(i) lt 0 then perrc(i)=0
	      if perrc(i) eq 0 then begin
	        perr(i)=erra(3L*i+0)
	        perrm(i)=erra(3L*i+0)
	      endif
	    endif
	    wdt(i)=aa(3L*i+1)
	    if jjw gt 0 then begin
	      if werrc(i) lt 0 then werrc(i)=0
	      if werrc(i) eq 0 then begin
	        werr(i)=erra(3L*i+1)
	        werrm(i)=erra(3L*i+1)
	      endif
	    endif
	    flx(i)=aa(3L*i+2)/dwvl
	    if jjf gt 0 then begin
	      if ferrc(i) lt 0 then ferrc(i)=0
	      if ferrc(i) eq 0 then begin
	        ferr(i)=erra(3L*i+2)
	        ferrm(i)=erra(3L*i+2)
	      endif
	    endif
	  endfor			;I=0,NCOM-1}
	  ;
          ;	update the widget state
	  npt=n_elements(xx)-nfree
	  c1='ChiSq = '+strtrim(chisq,2)+' / '+strtrim(npt,2)+$
		' (reduced: '+strtrim(chisq/float(npt),2)+')'
	  widget_control,txts_stts,set_value=c1
	  print,'' & print,c1 & print,''
          fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	    perrm=perrm,werrm=werrm,ferrm=ferrm,$
	    perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
          ;
	  updatelist=1
	  ;
	  ;	plot
	  ;plot,xx,zz,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
	  plot,xx,zz,xr=oxr,yr=oyr,/nodata,/xs,/ys
	  oplot,xx,zz,psym=10,col=icspc,thick=tspc
	  call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr, _extra=e
	  oplot,xx,ymod,col=icmod,thick=tmod
	;endif else begin
	;  print,string("7b)
	;  widget_control,txts_stts,set_value='Insufficient free parameters'
	;endelse
      endif				;NCOM>0)
      ;	make a note for historical purposes
      if n_tags(histr) gt 0 then begin		;(HISTR
	tagnam=tag_names(histr) & ntagnam=n_elements(tagnam) & nfit=0L
	if strpos(tagnam[ntagnam-1L],'FITTING',0) ge 0 then begin
	  tmphistr=histr.(ntagnam-1L)
	  oldthaw=tmphistr.THAW & oldnfit=tmphistr.NFIT & oldnfree=tmphistr.NFREE
	  oldtype=tmphistr.TYPE & oldties=tmphistr.TIES
	  oldalgo=tmphistr.ALGO_TYPE
	  if n_elements(thaw) ne n_elements(oldthaw) or $
		total(abs(thaw-oldthaw)) ne 0 or $
		total(abs(byte(type)-byte(oldtype))) ne 0 or $
		total(abs(byte(alltie)-byte(oldties))) ne 0 or $
		oldnfree ne npt or $
		algo_type ne oldalgo then begin
	    ihist=ihist+1L
	    tmphistr=create_struct('NFIT',Nfit,$
		'POS',pos,'WDT',wdt,'FLX',flx,'CHISQ',chisq,$
		'THAW',thaw,'NFREE',npt,'TYPE',type,'ALGO_TYPE',algo_type,$
		'TIES',alltie,'XRANGE',oxr)
	    histr=create_struct(histr,'FITTING'+strtrim(ihist,2),tmphistr)
	  endif else begin
	    histr.(ntagnam-1L).NFIT = oldnfit+1L
	    histr.(ntagnam-1L).POS = pos
	    histr.(ntagnam-1L).WDT = wdt
	    histr.(ntagnam-1L).FLX = flx
	    histr.(ntagnam-1L).CHISQ = chisq
	    histr.(ntagnam-1L).TIES = alltie
	  endelse
	endif else begin
	  ihist=ihist+1L
	  tmphistr=create_struct('NFIT',Nfit,$
	    'POS',pos,'WDT',wdt,'FLX',flx,'CHISQ',chisq,$
	    'THAW',thaw,'NFREE',npt,'TYPE',type,'ALGO_TYPE',algo_type,$
	    'TIES',alltie,'XRANGE',oxr)
	  histr=create_struct(histr,'FITTING'+strtrim(ihist,2),tmphistr)
	endelse
      endif					;HISTR)
    end						;fit}
    'menu_quit': begin				;{quit
      widget_control,txts_stts,set_value='EXITING...'
      wait,0.1
      ok='quit'
    end						;quit}

    'menu_view': 			;plot model (drop down menu)
    'view_show': begin				;{plot
      ;	plot the original spectrum, overplot continuum, overplot lines
      widget_control,txts_stts,set_value='Use VIEW->SetDefaults to '+$
	'change titles, colors, etc.'

      xx=x(oo) & yy=y(oo) & yysig=sigy(oo)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      ;plot,xx,yy,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,$
      ;	xtitle=xtitle,ytitle=ytitle,title=title,/xs,/ys
      plot,xx,yy,xr=oxr,yr=oyr,xtitle=xtitle,ytitle=ytitle,title=title,/xs,/ys,/nodata
      oplot,xx,yy,psym=10,col=icspc,thick=tspc
      if is_keyword_set(zerr) then begin
	moo=n_elements(oo)
	for ier=0L,moo-1L do oplot,xx[ier]*[1,1],yy(ier)+$
		yysig(ier)*zerr*[-1.,1.],col=icspc,thick=tspc
      endif
      if is_keyword_set(zlab) then begin
	olx=where(pos ge oxr(0) and pos le oxr(1),molx)
	for ilx=0L,molx-1L do xyouts,pos(olx(ilx)),oyr(1),epithet(olx(ilx)),$
		orient=90,align=1,charsize=zlab
      endif
      oplot,xx,conlev(oo),col=iccon,thick=tcon
      oplot,xx,conlev(oo)+consig(oo),col=iccon,line=1,thick=tcon
      oplot,xx,conlev(oo)-consig(oo),col=iccon,line=1,thick=tcon
      ncom=n_elements(pos)
      if ncom gt 0 then begin		;(if there are models to show
	ii=lindgen(ncom) & aa=fltarr(3L*ncom)
	aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr, _extra=e
	oplot,xx,ymod+conlev(oo),col=icmod,thick=tmod
	for i=0,ncom-1 do begin
	  aa=[pos(i),wdt(i),flx(i)*dwvl]
	  call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type(i),x_seg=x_seg, rmfstr=rmfstr, _extra=e
	  oplot,xx,ymod,col=iccmp,thick=tcmp
	endfor
      endif				;NCOM>0)
      stample, _extra=e

    end						;plot}
    'view_hcpy': begin				;{plot to postscript
      if not keyword_set(psfil) then psfil='idl.ps'
      if not keyword_set(iland) then iland=0
      if not keyword_set(iencap) then iencap=0
      if not keyword_set(icol) then icol=0
      widget_control,txts_stts,set_value='Dumping plot to file '+psfil
      set_plot,'ps'
      ccmd='device,file="'+psfil+'"'
      if keyword_set(iland) then ccmd=ccmd+',/landscape'
      if keyword_set(iencap) then ccmd=ccmd+',/encapsulated'
      if keyword_set(icol) then ccmd=ccmd+',/color'
      print,ccmd & jnk=execute(ccmd)
      ;device,file=psfil,landscape=iland,encapsulated=iencap,color=icol
      ;(same as VIEW->SHOW, except COL=Cxxx instead of COL=ICxxx
      xx=x(oo) & yy=y(oo)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      ;plot,xx,yy,xr=oxr,yr=oyr,psym=10,col=cspc,thick=tspc,$
      ;	xtitle=xtitle,ytitle=ytitle,title=title,/xs,/ys
      plot,xx,yy,xr=oxr,yr=oyr,xtitle=xtitle,ytitle=ytitle,title=title,/xs,/ys,/nodata
      oplot,xx,yy,psym=10,col=cspc,thick=tspc
      if is_keyword_set(zerr) then begin
	moo=n_elements(oo)
	for ier=0L,moo-1L do oplot,xx[ier]*[1,1],yy(ier)+$
		yysig(ier)*zerr*[-1.,1.],col=icspc,thick=tspc
      endif
      if is_keyword_set(zlab) then begin
	olx=where(pos ge oxr(0) and pos le oxr(1),molx)
	for ilx=0L,molx-1L do xyouts,pos(olx(ilx)),oyr(1),epithet(olx(ilx)),$
		orient=90,align=1,charsize=zlab
      endif
      oplot,xx,conlev(oo),col=ccon,thick=tcon
      oplot,xx,conlev(oo)+consig(oo),col=ccon,line=1,thick=tcon
      oplot,xx,conlev(oo)-consig(oo),col=ccon,line=1,thick=tcon
      ncom=n_elements(pos)
      if ncom gt 0 then begin		;(if there are models to show
	ii=lindgen(ncom) & aa=fltarr(3L*ncom)
	aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr, _extra=e
	oplot,xx,ymod+conlev(oo),col=cmod,thick=tmod
	for i=0,ncom-1 do begin
	  aa=[pos(i),wdt(i),flx(i)*dwvl]
	  call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type(i),x_seg=x_seg, rmfstr=rmfstr, _extra=e
	  oplot,xx,ymod,col=ccmp,thick=tcmp
	endfor
      endif				;NCOM>0)
      stample, _extra=e
      ;end same part as VIEW->SHOW)
      device,/close & set_plot,'x'
      if !D.N_COLORS gt 256 then device,decomposed=0
    end						;plot to postscript}
    'view_colr': 		;choosing color tables
    'colr_ctpb': begin				;{peasecolr, black bg
      plot,[0],/nodata
      ;setkolor,'black',0
      peasecolr,verbose=10
      cc='Loading Pease color table for black background'
      widget_control,txts_stts,set_value=cc
    end
    'colr_ctpw': begin
      plot,[0],/nodata
      ;setkolor,'white',0
      peasecolr,verbose=10,/white
      cc='Loading Pease color table for white background'
      widget_control,txts_stts,set_value=cc
    end
    'colr_idlt': begin
      widget_control,txts_stts,set_value='Type in color table number via keyboard at terminal window'
      if not keyword_set(iloadct) then iloadct=3
      cloadct=''
      read,prompt='IDL color table number ['+strtrim(iloadct,2)+']> ',cloadct
      if keyword_set(cloadct) then iloadct=fix(cloadct)
      cc='loadct,'+strtrim(iloadct,2) & jnk=execute(cc)
      widget_control,txts_stts,set_value='loading color table '+$
	strtrim(iloadct,2)
    end
    'view_defs': begin				;{set plot defaults
      widget_control,txts_stts,set_value='Change plot attributes from keyboard (at terminal window)'
      print,string("7b)
      print,'*** KEYBOARD INPUT : type Q when done to return to GUI ***'
      print,'     X,Y,M to set x-,y-,main-title'
      print,'     ~ to toggle error bars'
      print,'     @ to toggle labels'
      print,'     E,L,C to toggle encapsulated, landscape, color PS output'
      print,'     F to set output filename'
      print,'     S to change color attributes'
      print,'     T to set line thicknesses'
      print,'*** KEYBOARD INPUT : type Q when done to return to GUI ***'
      cc=''
      while cc ne 'q' do begin		;{a smaller endless loop
        kilroy
	cc=strlowcase(get_kbrd(1))
	case cc of			;{set plot attribute
	  'x': begin
	    tmp='' & read,prompt='X-title['+xtitle+']> ',tmp
	    if strtrim(tmp,2) ne '' then xtitle=tmp
            widget_control,txts_stts,set_value='XTITLE='+xtitle
	  end
	  'y': begin
	    tmp='' & read,prompt='Y-title['+ytitle+']> ',tmp
	    if strtrim(tmp,2) ne '' then ytitle=tmp
            widget_control,txts_stts,set_value='YTITLE='+ytitle
	  end
	  'm': begin
	    tmp='' & read,prompt='plot-title['+title+']> ',tmp
	    if strtrim(tmp,2) ne '' then title=tmp
            widget_control,txts_stts,set_value='TITLE='+title
	  end
	  '~': begin
	    tmp='' & szerr=strtrim(zerr,2)
	    read,prompt='error-bar width, in units of sigma ['+szerr+']> ',tmp
	    if strtrim(tmp,2) ne '' then zerr=float(tmp)
	    widget_control,txts_stts,set_value='Error bars plotted @ '+$
		strtrim(zerr,2)+'-sigma'
	  end
	  '@': begin
	    tmp='' & szlab=strtrim(zlab,2) ;& szlab2=strtrim(1-zlab,2)
	    read,prompt='charsize of labels: ['+szlab+'] >',tmp
	    if strtrim(tmp,2) ne '' then zlab=float(tmp)
	    ;if strtrim(tmp,2) eq '' then tmp=szlab2 & zlab=fix(tmp)
	    if not is_keyword_set(zlab) then widget_control,txts_stts,$
		set_value='drawing labels with size '+strtrim(zlab,2) else $
		widget_control,txts_stts,set_value='no labeling'
	  end
	  's': begin
	    print,'currently assigned colors are:'
	    print,'   #1. spectrum   .. '+strtrim(cspc,2)
	    print,'   #2. model      .. '+strtrim(cmod,2)
	    print,'   #3. component  .. '+strtrim(ccmp,2)
	    print,'   #4. continuum  .. '+strtrim(ccon,2)
	    jc=1 & jcol=1
	    while jc gt 0 do begin
	      read,prompt='type NUMBER COLORINDEX (0 0 when done)> ',jc,jcol
	      if jc eq 1 then cspc=jcol
	      if jc eq 2 then cmod=jcol
	      if jc eq 3 then ccmp=jcol
	      if jc eq 4 then ccon=jcol
	      if jc eq 1 then c1='Spectrum color index = '+strtrim(jcol,2)
	      if jc eq 2 then c1='Model color index = '+strtrim(jcol,2)
	      if jc eq 3 then c1='Component color index = '+strtrim(jcol,2)
	      if jc eq 4 then c1='Continuum color index = '+strtrim(jcol,2)
              iccon=ccon ;	iccon=fix(ccon*float(!d.n_colors)/dncolors+1)
              icmod=cmod ;	icmod=fix(cmod*float(!d.n_colors)/dncolors+1)
              iccmp=ccmp ;	iccmp=fix(ccmp*float(!d.n_colors)/dncolors+1)
              icspc=cspc ;	icspc=fix(cspc*float(!d.n_colors)/dncolors+1)
	      print,'setting '+c1
              widget_control,txts_stts,set_value=c1
	    endwhile
	  end
	  't': begin
	    print,'currently assigned line thicknesses are:'
	    print,'   #1. spectrum   .. '+strtrim(tspc,2)
	    print,'   #2. model      .. '+strtrim(tmod,2)
	    print,'   #3. component  .. '+strtrim(tcmp,2)
	    print,'   #4. continuum  .. '+strtrim(tcon,2)
	    jc=1 & lthk=1
	    while jc gt 0 do begin
	      read,prompt='type NUMBER THICKNESS (0 0 when done)> ',jc,lthk
	      if jc eq 1 then tspc=float(lthk)
	      if jc eq 2 then tmod=float(lthk)
	      if jc eq 3 then tcmp=float(lthk)
	      if jc eq 4 then tcon=float(lthk)
	      if jc eq 1 then c1='Spectrum line thickness = '+strtrim(lthk,2)
	      if jc eq 2 then c1='Model line thickness = '+strtrim(lthk,2)
	      if jc eq 3 then c1='Component line thickness = '+strtrim(lthk,2)
	      if jc eq 4 then c1='Continuum line thickness = '+strtrim(lthk,2)
	      print,'setting '+c1
              widget_control,txts_stts,set_value=c1
	    endwhile
	  end
	  'e': begin
	    if not keyword_set(iencap) then iencap=1 else iencap=0
	    print,'DEVICE: ENCAPSULATED='+strtrim(iencap,2)
            widget_control,txts_stts,set_value='ENCAPSULATED POSTSCRIPT =?= '+$
		strtrim(iencap,2)
	  end
	  'l': begin
	    if not keyword_set(iland) then iland=1 else iland=0
	    print,'DEVICE: LANDSCAPE='+strtrim(iland,2)
            widget_control,txts_stts,set_value='LANDSCAPE =?= '+$
		strtrim(iland,2)
	  end
	  'c': begin
	    if not keyword_set(icol) then icol=1 else icol=0
	    print,'DEVICE: COLOR='+strtrim(icol,2)
            widget_control,txts_stts,set_value='COLOR POSTSCRIPT =?= '+$
		strtrim(icol,2)
	  end
	  'f': begin
	    if not keyword_set(psfil) then psfil='idl.ps'
	    fil='' & read,prompt='type output postscript filename ['+psfil+']> ',fil
	    if strtrim(fil,2) ne '' then psfil=fil
	    c1='Postscript to file -> '+psfil
	    print,c1 & widget_control,txts_stts,set_value=c1
	  end
	  'q': begin
	    print,'Returning to GUI'
            widget_control,txts_stts,set_value='Done'
	  end
	  else: begin
            print,'*** @KEYBOARD : type Q when done to return to GUI ***'
            print,'     X,Y,M to set x-,y-,main-title'
            print,'     ~ to toggle error bars'
            print,'     @ to toggle labels'
            print,'     E,L,C to toggle encapsulated, landscape, color PS output'
            print,'     F to set output filename'
            print,'     S to change color attributes'
            print,'     T to set line thicknesses'
            print,'*** @KEYBOARD : type Q when done to return to GUI ***'
	    widget_control,txts_stts,set_value='type Q if done'
	  end
	endcase				;plot attributes}
      endwhile					;cc v/s q}
    end						;set plot defaults}

    'menu_xrng': ;set range -- drop down menu
    'xrng_crsr': begin				;{call PICKRANGE
      widget_control,txts_stts,set_value='Right button click when done'
      oldOO=oo & oldOXR=oxr & oldOYR=oyr & xx=x(oo)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      z=y-conlev & ncom=n_elements(pos) & mx=0 & my=0
      ;if ncom gt 0 then begin & mx=pos & my=flx/wdt/sqrt(2*!pi) & endif
      if ncom gt 0 then begin & mx=pos & my=flx & endif
	;	either the selection defined by left-cursor..
      oo=pickrange(x,z,markx=mx,marky=my,/legg,oxr=oxr,oyr=oyr,$
	_extra=e)
	;	.. or the selection defined by right-cursor zoom
      if oo(0) eq -1 then oo=where(x ge oxr(0) and x le oxr(1)) else $
	oxr=[min(x(oo)),max(x(oo))]
      if oo(0) eq -1 then begin
	message,'range too small; making no changes',/info
	oo=oldOO & oxr=oldOXR & oyr=oldOYR
      endif
      ;plot,x,z,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
      plot,x,z,xr=oxr,yr=oyr,/xs,/ys,/nodata
      oplot,x,z,psym=10,col=icspc,thick=tspc
      widget_control,txts_stts,set_value='Chosen range = ['+$
	strtrim(oxr(0),2)+','+strtrim(oxr(1),2)+']'
    end						;call PICKRANGE}
    'xrng_kbrd': begin				;{set range from keyboard
      oldOO=oo & oldOXR=oxr & oldOYR=oyr
      widget_control,txts_stts,set_value='type Xmin, Xmax via keyboard in terminal window'
      read,prompt='type X-range: ',x0,x1
      x0=min([x0,x1],max=x1) & oo=where(x ge x0 and x le x1)
      if oo(0) eq -1 then begin
	message,'range too small; making no changes',/info
	oo=oldOO & oxr=oldOXR & oyr=oldOYR
      endif else oxr=[x0,x1]
      ;	plot
      z=y-conlev & xx=x(oo) & zz=z(oo)
      ;plot,xx,zz,psym=10,xr=oxr,/xs,col=icspc,thick=tspc
      plot,xx,zz,xr=oxr,/xs,/nodata
      oplot,xx,zz,psym=10,col=icspc,thick=tspc
      widget_control,txts_stts,set_value='Chosen range = ['+$
	strtrim(oxr(0),2)+','+strtrim(oxr(1),2)+']'
    end						;set range from keyboard}
    'xrng_auto': begin				;{set range automagically
      oldOO=oo & oldOXR=oxr & oldOYR=oyr
      ncom=n_elements(pos)
      if ncom gt 0 then begin		;(if there are models to fit..
        thaw=intarr(3*ncom) & ii=lindgen(ncom)
        widget_control,pars_list,get_value=list
        thaw(3*ii+0)=list(4*ii+1)		;positions
        thaw(3*ii+1)=list(4*ii+2)		;widths
        thaw(3*ii+2)=list(4*ii+3)		;fluxes
	xmin=oxr(1) & xmax=oxr(0)
	for i=0L,ncom-1L do begin
	  if thaw(3*i+0) gt 0 then xmin=pos(i)-5*wdt(i) < xmin
	  if thaw(3*i+0) gt 0 then xmax=pos(i)+5*wdt(i) > xmax
	endfor
        ;	plot
        x0=min([xmin,xmax],max=x1) & oo=where(x ge x0 and x le x1,moo)
        if moo eq 0 then begin
	  widget_control,txts_stts,set_value='range too small; making no changes'
	  oo=oldOO & oxr=oldOXR
        endif else begin
	  oxr=[x0,x1]
          z=y-conlev & xx=x(oo) & zz=z(oo)
          ;plot,xx,zz,psym=10,xr=oxr,/xs,col=icspc,thick=tspc
          plot,xx,zz,xr=oxr,/xs,/nodata
          oplot,xx,zz,psym=10,col=icspc,thick=tspc
          widget_control,txts_stts,set_value='Chosen X range = ['+$
	    strtrim(oxr(0),2)+','+strtrim(oxr(1),2)+']'
	endelse
      endif else widget_control,txts_stts,set_value='Gimme a hint, please!'
    end						;set range automagically}
    'xrng_zero': begin				;{reset full range
      oo=lindgen(n_elements(x))
      z=y-conlev & xx=x(oo) & zz=z(oo) & zsig=sqrt(sigy(oo)^2+consig(oo)^2)
      oxr=[min(x),max(x)] & oyr=[min(y),max(y)]
      ;	plot
      ;plot,x,z,psym=10,/xs,col=icspc,thick=tspc
      plot,x,z,/xs,/nodata
      oplot,x,z,psym=10,col=icspc,thick=tspc
      widget_control,txts_stts,set_value='Chosen range = ['+$
	strtrim(oxr(0),2)+','+strtrim(oxr(1),2)+']'
    end						;reset full range}

    'menu_yrng': ;set range -- drop down menu
    'yrng_crnt': begin				;{reset to current selection
      z=y-conlev & xx=x(oo) & yy=y(oo) & zz=z(oo)
      delz=max(yy)-min(zz) & oyr=[min(zz)-0.05*delZ,max(yy)+0.05*delZ]
      if delZ eq 0 then oyr=zz(0)+0.5*[-1,1]
      ;plot,x,z,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
      plot,x,z,xr=oxr,yr=oyr,/xs,/ys,/nodata
      oplot,x,z,psym=10,col=icspc,thick=tspc
      widget_control,txts_stts,set_value='Chosen range = ['+$
	strtrim(oyr(0),2)+','+strtrim(oyr(1),2)+']'
    end						;reset to current selection}
    'yrng_kbrd': begin				;{set range from keyboard
      widget_control,txts_stts,set_value='type Ymin, Ymax via keyboard in terminal window'
      read,prompt='type Y-range: ',y0,y1
      y0=min([y0,y1],max=y1)
      if y1 gt y0 then begin
	oyr=[y0,y1]
        ;	plot
        z=y-conlev & xx=x(oo) & zz=z(oo)
        ;plot,xx,zz,psym=10,xr=oxr,/xs,col=icspc,thick=tspc
        plot,xx,zz,xr=oxr,/xs,/nodata
        oplot,xx,zz,psym=10,col=icspc,thick=tspc
        widget_control,txts_stts,set_value='Chosen Y range = ['+$
	  strtrim(oyr(0),2)+','+strtrim(oyr(1),2)+']'
      endif else begin
        widget_control,txts_stts,set_value='Incomprehensible input.  Making no changes'
      endelse
    end						;set range from keyboard}
    'yrng_zero': begin				;{reset full range
      z=y-conlev & oyr=[min(z),max(z)]
      ;	plot
      xx=x(oo) & zz=z(oo)
      ;plot,xx,zz,psym=10,/xs,col=icspc,thick=tspc
      plot,xx,zz,/xs,/nodata
      oplot,xx,zz,psym=10,col=icspc,thick=tspc
      widget_control,txts_stts,set_value='Resetting Y range to: ['+$
	strtrim(oyr(0),2)+','+strtrim(oyr(1),2)+']'
    end						;reset full range}

    'mods_modl': ;modify model -- this is a drop down menu
    'modl_friz': begin				;{freeze all
      ncom=n_elements(pos)
      if ncom gt 0 then begin
	widget_control,txts_stts,set_value='All parameters frozen'
        widget_control,pars_list,get_value=list
	list(*)=0
	ncom=n_elements(pos) & if ncom gt 0 then thaw=intarr(3*ncom)
	widget_control,pars_list,set_value=list
      endif else widget_control,txts_stts,$
	set_value='No models defined, nothing to freeze!'
    end						;freeze all}
    'modl_thaw': begin				;{thaw all
      ncom=n_elements(pos)
      if ncom gt 0 then begin
	widget_control,txts_stts,set_value='All parameters free'
        widget_control,pars_list,get_value=list
	list(*)=1
	ncom=n_elements(pos) & if ncom gt 0 then thaw=intarr(3*ncom)+1
	widget_control,pars_list,set_value=list
      endif else widget_control,txts_stts,$
	set_value='No models defined, nothing to thaw!'
    end						;thaw all}
    'modl_fixx': begin				;{freeze unimportant
      ncom=n_elements(pos)
      if ncom gt 0 then begin
	oc=where(pos+2*abs(wdt) lt oxr(0) or pos-2*abs(wdt) gt oxr(1),moc)
	if moc gt 0 then begin
	  widget_control,txts_stts,set_value='Freezing all components '+$
		'that lie far beyond chosen range'
          widget_control,pars_list,get_value=list
	  list(4*oc+0)=0	;freeze components
	  list(4*oc+1)=0	;freeze position parameters
	  list(4*oc+2)=0	;freeze width parameters
	  list(4*oc+3)=0	;freeze flux parameters
	  widget_control,pars_list,set_value=list
	  thaw=intarr(3*ncom) & ii=lindgen(ncom)
	  thaw(3*ii+0)=list(4*ii+1)		;positions
	  thaw(3*ii+1)=list(4*ii+2)		;widths
	  thaw(3*ii+2)=list(4*ii+3)		;fluxes
	endif
      endif else widget_control,txts_stts,$
	set_value='No models defined, nothing to freeze!',/info
    end						;freeze unimportant}
    'modl_addc': begin				;{add component
      widget_control,txts_stts,set_value='Click & Drag to mark, Click when done'
      z=y-conlev & xx=x(oo) & zz=z(oo) & ncom=n_elements(pos)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      ;plot,xx,zz,xr=oxr,/xs,psym=10,col=icspc,thick=tspc
      plot,xx,zz,xr=oxr,/xs,/nodata
      oplot,xx,zz,psym=10,col=icspc,thick=tspc
      print,string("7b)
      print,'CLICK & DRAG mouse to denote location, peak, and width'
      print,'CLICK to exit'
      get_comp=1
      while get_comp eq 1 do begin	;{interactively define new components
	cursor,x0,y0,/down,/data & cursor,x1,y1,/up,/data
	peak=y0 & width=abs(x1-x0)
	draglim=1e-3*sqrt((oxr(1)-oxr(0))^2+(oyr(1)-oyr(0))^2)
	draglim=1e-3*abs(oxr(1)-oxr(0))
	if width gt draglim then begin			;(click+drag
	  ;convert PEAK to NORMalization assuming (for now) gaussian form
	  z0=peak & if not keyword_set(intens) then z0=z0/width/sqrt(2.*!pi)
	  if ncom gt 0 then begin
	    flxerr=sqrt(abs(z0/dwvl)+0.75)+1.
	    pos=[pos,x0] & perr=[perr,0.] & perrm=[perrm,0.]
	    flx=[flx,z0/dwvl] & ferr=[ferr,flxerr] & ferrm=[ferrm,flxerr]
	    wdt=[wdt,width] & werr=[werr,0.] & werrm=[werrm,0.]
	    perrc=[perrc,0.] & werrc=[werrc,0.] & ferrc=[ferrc,0.]
	    thaw=[thaw,1,1,1] & type=[type,type(0)]
	    epithet=[epithet,strtrim(x0,2)]
	    ncom=ncom+1
	    ;
	    ;	when original is reset, don't lose the added components
	    opos=[opos,x0] & owdt=[owdt,width] & oflx=[oflx,z0/dwvl]
	    operr=[operr,0.] & owerr=[owerr,0.] & oferr=[oferr,flxerr]
	    operrm=[operrm,0.] & owerrm=[owerrm,0.] & oferrm=[oferrm,flxerr]
	    operrc=[operrc,0.] & owerrc=[owerrc,0.] & oferrc=[oferrc,0.]
	    otype=[otype,type(0)] & othaw=[othaw,1,1,1]
	    oepithet=[oepithet,strtrim(x0,2)]
	    ;
	    ;	nor when previous fit is reset
	    if n_elements(spos) gt 0 then begin
	      spos=[spos,x0] & swdt=[swdt,width] & sflx=[sflx,z0/dwvl]
	      sperr=[sperr,0.] & swerr=[swerr,0.] & sferr=[sferr,flxerr]
	      sperrm=[sperrm,0.] & swerrm=[swerrm,0.] & sferrm=[sferrm,flxerr]
	      sperrc=[sperrc,0.] & swerrc=[swerrc,0.] & sferrc=[sferrc,0.]
	      stype=[stype,type(0)] & sthaw=[sthaw,1,1,1]
	      sepithet=[sepithet,epithet(0)] & sthaw=[sthaw,1,1,1]
	    endif
	    ;
	    ;	and make sure the global error storage doesn't screw up
	  endif else begin
	    pos=[x0] & perr=[0.] & perrm=[0.]
	    flx=[z0/dwvl] & ferr=[sqrt(abs(z0/dwvl)+0.75)+1.] & ferrm=ferr
	    wdt=[width] & werr=[0.] & werrm=[0.]
	    perrc=[0.] & werrc=[0.] & ferrc=[0.]
	    thaw=[1,1,1] & type=['gauss'] & epithet=[strtrim(x0,2)]
	    ncom=1
	    ;
	    ;	as above, save..
	    opos=[x0] & owdt=[width] & oflx=[z0/dwvl]
	    operr=[0.] & owerr=[0.] & oferr=[0.]
	    operrm=operr & owerrm=owerr & oferrm=oferr
	    operrc=[0.] & owerrc=[0.] & oferrc=[0.]
	    otype=['gauss'] & othaw=[1,1,1]
	    oepithet=[strtrim(x0,2)]
	    spos=[x0] & swdt=[width] & sflx=[z0/dwvl]
	    sperr=[0.] & swerr=[0.] & sferr=[0.]
	    sperrm=sperr & swerrm=swerr & sferrm=sferr
	    sperrc=[0.] & swerrc=[0.] & sferrc=[0.]
	    stype=['gauss'] & sthaw=[1,1,1] & sepithet=[strtrim(x0,2)]
	  endelse
	  ;	overplot the model component
	  aaa=[x0,width,peak]
	  call_procedure,funcs,xx,aaa,ycomp,type='gauss', _extra=e
	  oplot,xx,ycomp,col=icmod,thick=tmod
	endif else get_comp=0				;click)
      endwhile				;get new components}
      ;
      ;	update the widget state
      fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
        perrm=perrm,werrm=werrm,ferrm=ferrm,$
        perrc=perrc,werrc=werrc,ferrc=ferrc,$
        plist=plist,list=list, tiep=tiep,tiew=tiew,tief=tief,ties=ties
      ;
      updatelist=1
      widget_control,txts_stts,set_value='Done'
    end						;add component}
    'modl_delc': begin				;{delete component
      widget_control,txts_stts,set_value='Choose component to delete '+$
	'(type component number via keyboard in terminal window)'
      print,string("7b) & idel=-1L
      read,prompt='type component # to delete (-1 to not)> ',idel
      fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
        perrm=perrm,werrm=werrm,ferrm=ferrm,$
        perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list,$
	tiep=tiep,tiew=tiew,tief=tief,ties=ties,delcmp=idel
      ;
      ;	update the parameter list and the list of ties
      updatelist=1
      updateties=1
      widget_control,txts_stts,set_value='Done'
    end						;delete component}

    'mods_cont': ;recompute continuum -- this is a drop down menu
    'cont_rcmp': begin				;{recompute over selected range
      cnt_histr = 'Remove_Lines' 
      widget_control,txts_stts,set_value='Removing lines '+$
	'over choosen range'
      xx=x(oo) & yy=y(oo) & ysig=sigy(oo)
      yc=linerem(xx,yy,sig=ysig,cell=csize,nsigma=nsig,$
	bkgval=cony,bkgerr=scony,/quiet, _extra=e)
      ;	plot
      ;plot,xx,yy,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
      plot,xx,yy,xr=oxr,yr=oyr,/xs,/ys,/nodata
      oplot,xx,yy,psym=10,col=icspc,thick=tspc
      oplot,xx,cony,col=iccon,thick=tspc
      oplot,xx,cony+scony,line=1,col=iccon,thick=tspc
      oplot,xx,cony-scony,line=1,col=iccon,thick=tspc
      widget_control,txts_stts,set_value='Done.  Remember to ACCEPT the new continuum!'
    end						;recompute}
    'cont_fcmp': begin				;{recompute over full range
      cnt_histr = 'LineRem(sub)'
      widget_control,txts_stts,set_value='Recomputing the continuum '+$
	'over full range'
      xx=x & yy=y & ysig=sigy
      yc=linerem(xx,yy,sig=ysig,cell=csize,nsigma=nsig,$
	bkgval=cony,bkgerr=scony,/quiet, _extra=e)
      ;	plot
      ;plot,xx,yy,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
      plot,xx,yy,xr=oxr,yr=oyr,/xs,/ys,/nodata
      oplot,xx,yy,psym=10,col=icspc,thick=tspc
      oplot,xx,cony,col=iccon,thick=tspc
      oplot,xx,cony+scony,line=1,col=iccon,thick=tspc
      oplot,xx,cony-scony,line=1,col=iccon,thick=tspc
      widget_control,txts_stts,set_value='Done.  Remember to ACCEPT the new continuum!'
    end						;recompute}
    'cont_pcws': begin				;{define it piecewise
      cnt_histr = 'Piecewise'
      widget_control,txts_stts,set_value='Define continuum with cursor'
      xx=x(oo) & yy=y(oo) & ycont=conlev(oo)
      oxrng=oxr & oyrng=oyr
      cony=setcont(xx,yy,scony,ycont=ycont,xcont=xcont,const=const,$
	xrange=oxrng,yrange=oyrng)
      widget_control,txts_stts,set_value='Done.  Choose CONTINUUM->ACCEPT to accept the continuum'
    end						;piecewise model}
    'cont_norm': begin				;{global renorm continuum
      cnt_histr = 'Adjust_Norm'
      widget_control,txts_stts,set_value='change global normalization of '+$
	'continuum using the cursor'
      xx=x(oo) & yy=y(oo) & ycont=conlev(oo)
      oxrng=oxr & oyrng=oyr
      cony=setcont(xx,yy,ycont=ycont,xcont=xcont,/lock,$
	xrange=oxrng,yrange=oyrng)
      widget_control,txts_stts,set_value='Done.  Choose CONTINUUM->ACCEPT to accept the continuum'
    end						;CONT_NORM}
    'cont_acpt': begin				;{accept new continuum
      if not keyword_set(cnt_histr) then cnt_histr = 'ACCEPT'
      widget_control,txts_stts,set_value='Replacing old continuum with '+$
	'the newly minted'
      nn=n_elements(cony) & moo=n_elements(oo)
      if n_elements(scony) ne nn then scony=0.*cony
      if nn eq moo and nn ne nx then begin
	conlev(oo)=cony(*) & consig(oo)=scony(*)
      endif
      if nn eq nx then begin
	conlev=cony & consig=scony
      endif
      ;	and plot the results
      xx=x(oo) & z=y-conlev & zz=z(oo) & zsig=sqrt(sigy(oo)^2+consig(oo)^2)
      ;plot,xx,zz,psym=10,/xs,/ys,col=icspc,thick=tspc
      plot,xx,zz,/xs,/ys,/nodata
      oplot,xx,zz,psym=10,col=icspc,thick=tspc
      if is_keyword_set(zerr) then begin
	moo=n_elements(oo)
	for ier=0L,moo-1L do oplot,zz(ier)+zsig(ier)*zerr*[-1.,1.]
      endif
      ncom=n_elements(pos)
      if ncom gt 0 then begin		;(if there are models to show
	ii=lindgen(ncom) & aa=fltarr(3L*ncom)
	aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr, _extra=e
	oplot,xx,ymod,col=icmod,thick=tmod
      endif				;NCOM>0)
      if n_tags(histr) gt 0 then begin 
         tagnam=tag_names(histr) & ntagnam=n_elements(tagnam) 
         ihist = ihist+1L 
         if keyword_set(cnt_histr) then $
         tmphistr = create_struct('conlev',conlev,'consig',consig,'button',cnt_histr,'XRANGE',oxr) 
         histr=create_struct(histr,'CONTINUUM'+strtrim(ihist,2),tmphistr)
      endif
    end						;accept new continuum}
    'cont_orig': begin				;{reset to original continuum
      cnt_histr = 'RESET'
      widget_control,txts_stts,set_value='Resetting continuum to initial values'
      print,string("7b)
      print,'WARNING: This will reset the continuum over the ENTIRE range!'
      print,'Type any key to continue, N if you changed your mind'
      c1='' & c1=get_kbrd(1)
      if strlowcase(c1) ne 'n' then begin
        conlev=oconlev & consig=oconsig
      endif else widget_control,txts_stts,set_value='continuum NOT reset'
    end						;reset to original continuum}
		;entry fields
    'cont_cell': begin				;{cell size for LINEREM
      cnt_histr = 'cont_cell' 
      widget_control,txts_stts,set_value='Smoothing size for removing lines '+$
	'(a bit bigger than line width in pixels)'
      widget_control,cont_cell,get_value=csize
    end						;cell size for LINEREM}
    'cont_nsig': begin				;{NSIGMA for LINEREM
      widget_control,txts_stts,set_value='Threshold for removing lines'
      widget_control,cont_nsig,get_value=nsig
    end						;NSIGMA for LINEREM}

    ;							third row
    'fits_help': begin				;{display help file
      widget_control,txts_stts,set_value='working...'
      if not keyword_set(TOPDIR) then begin
	;	where would the helpfile be?
	tmp=routine_info('FITLINES_EVENT',/source)
	scardir=tmp.PATH
	jvar=rstrpos(scardir,'pro/')
	TOPDIR=strmid(scardir,0L,jvar-1L)
      endif 
      ;if not keyword_set(TOPDIR) then which,'initale.pro',$
      ;	TOPDIR,/dironly,verbose=1
      ;ipath=strpos(TOPDIR[0],'pro/scrypt',0)
      ;TOPDIR=strmid(TOPDIR[0],0,ipath)
      if fix(strmid(!version.release,0,1)) lt 5 then begin
	openr,uhelp,TOPDIR(0)+'/ardb/'+fit_help(0),/get_lun
	while not eof(uhelp) do begin
	  clin='' & readf,uhelp,clin & print,clin
	endwhile
	close,uhelp & free_uhelp
      endif else begin
	xdisplayfile,TOPDIR(0)+'/ardb/'+fit_help(0),$
		font='10x20',title='FITLINES HELP',width=80
		;font='10x20',/modal,title='FITLINES HELP',width=80
      endelse
      widget_control,txts_stts,set_value='Done.'
    end						;FITS_HELP}

    'fits_zero': begin				;{update the "Original"
      if n_elements(pos) gt 0 then begin
        widget_control,txts_stts,set_value='Updating the "Original" set'
	opos=pos & owdt=wdt & oflx=flx
	operr=perr & owerr=werr & oferr=ferr
	operrm=perrm & owerrm=werrm & oferrm=ferrm
	operrc=perrc & owerrc=werrc & oferrc=ferrc
	otype=type & oties=ties
	oepithet=epithet
      endif else widget_control,txts_stts,set_value='whachutakinaboot, willis?'
    end						;FITS_ZERO}

    'fits_undo': begin				;{undo previous fit
      if n_elements(spos) gt 0 then begin
        widget_control,txts_stts,set_value='Parameters reset to previous'
	pos=spos & wdt=swdt & flx=sflx
	perr=sperr & werr=swerr & ferr=sferr
	perrm=sperrm & werrm=swerrm & ferrm=sferrm
	perrc=sperrc & werrc=swerrc & ferrc=sferrc
	type=stype & epithet=sepithet
        fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
	  perrc=perrc,werrc=werrc,ferrc=ferrc,$
	  plist=plist,list=list,tiep=tiep,tiew=tiew,tief=tief,ties=ties
	updatelist=1
      endif else widget_control,txts_stts,set_value='whachutakinaboot, willis?'
    end						;FITS_UNDO}

    'fits_orig': begin				;{Reset *everything*
      ncom=n_elements(opos)
      if ncom gt 0 then begin
        widget_control,txts_stts,set_value=$
		'Type any key to restore original parameters'
        print,string("7b)
        print,'This will reset ALL parameters to their original values!'
        print,'Type any key to continue, N if you changed your mind'
        c1='' & c1=get_kbrd(1)
        if strlowcase(c1) ne 'n' then begin
	  pos=opos & wdt=owdt & flx=oflx
	  perr=operr & werr=owerr & ferr=oferr
	  perrm=operrm & werrm=owerrm & ferrm=oferrm
	  perrc=operrc & werrc=owerrc & ferrc=oferrc
	  type=otype & ties=oties & thaw=othaw & epithet=oepithet
	  if keyword_set(odchisq) then begin
	    dchisq=odchisq
	    widget_control,fits_dchi,set_value=dchisq
	  endif
          widget_control,txts_stts,set_value='Parameters reset to original'
          fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	    perrm=perrm,werrm=werrm,ferrm=ferrm,$
	    perrc=perrc,werrc=werrc,ferrc=ferrc,$
	    plist=plist,list=list,tiep=tiep,tiew=tiew,tief=tief,ties=ties
	  updatelist=1
	  updateties=1
        endif else widget_control,txts_stts,set_value='Parameters NOT reset'
      endif else begin
	widget_control,txts_stts,set_value='Parameters NOT reset'
      endelse
    end						;FITS_ORIG}

    'fits_dump': begin				;{dump to disk
      if not keyword_set(savfil) then savfil='fitlines.save'
      c1='' & read,prompt='Save to file? ['+savfil+']> ',c1
      cc=strlowcase(strtrim(c1,2))
      if cc ne '' and cc ne 'y' and cc ne 'yes' then savfil=c1
      widget_control,txts_stts,set_value=savfil+':: use FITLINES_UNDUMP'+$
	' to return to present state'
      save,file=savfil,/compress
    end						;FITS_DUMP}

    'fits_norm': begin				;{renormalize model amplitudes
      ;	renormalize models to match the total counts in observed spectrum,
      ;	taking into account such things as frozen amplitudes, etc.
      widget_control,txts_stts,set_value='Renormalize the model amplitudes '+$
	'to match the overall observed counts'

      z=y-conlev & xx=x(oo) & zz=z(oo) & tot_obs=total(zz)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      if tot_obs ne 0 then begin
	;	get the model predictions
        ncom=n_elements(pos)
        if ncom gt 0 then begin		;(if there are models to show
	  ii=lindgen(ncom) & aa=fltarr(3L*ncom)
	  aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	  call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr, _extra=e
	  oz=where(ymod ge 0.01*max(ymod),moz)
	  if moz eq 0 then stop		;debugging step
	  tot_obs_mod=total(zz(oz))
	  tot_mod=total(ymod)
	  if tot_obs_mod eq 0 then tot_obs_mod=tot_mod
	  if tot_mod ne 0 then begin		;(model fluxes non-zero
	    widget_control,pars_list,get_value=list
	    ii=lindgen(ncom) & ampthaw=list(4*ii+3)
	    oamp=where(ampthaw gt 0,moamp)
	    if moamp gt 0 then begin		;(thawed amplitudes exist
	      ;	set T(Mfrozen)+f*T(Mthaw)=T(Obs), given that
	      ;	T(Mfroz)+T(Mthaw)=T(M) and T(Mthaw)=(T(Mampthaw)/T(Mamp))*T(M)
	      ;	i.e., (f-1)*T(Mthaw)=T(Obs)-T(M)
	      tot_mthaw=0.
	      for i=0,moamp-1 do begin
		aa=[pos(oamp(i)),wdt(oamp(i)),flx(oamp(i))*dwvl]
		call_procedure,funcs,xx,aa,ymod,normflx=znorm,type=type,x_seg=x_seg, rmfstr=rmfstr,$
		_extra=e
		tot_mthaw=tot_mthaw+total(ymod)
	      endfor
	      fracamp=1.
	      if tot_mthaw ne 0 then fracamp=(tot_obs_mod-tot_mod)/tot_mthaw+1.
	      newflx=flx
	      ;if tmamp ne 0 then begin
	        ;tmthaw=(total(newflx(oamp))/tmamp)*tot_mod
	        ;if tmthaw ne 0 then fracamp=(tot_obs_mod-tot_mod)/tmthaw + 1.
	      ;endif
	      widget_control,txts_stts,set_value='Multiplying by '+$
		strtrim(fracamp,2)
	      if fracamp le 0 then begin
		wait,1
		widget_control,txts_stts,set_value='Uh oh.  Problem.  See terminal window for details.'
		message,'Renormalization factor is coming out as a negative number',/informational
		help,FRACAMP
		message,'Reset the variable FRACAMP manually and type .CON to continue',/informational
		stop
	      endif
	      newflx(oamp)=newflx(oamp)*fracamp
	      ;
	      spos=pos & swdt=wdt & sflx=flx
	      sperr=perr & swerr=werr & sferr=ferr
	      sperrm=perrm & swerrm=werrm & sferrm=ferrm
	      sperrc=perrc & swerrc=werrc & sferrc=ferrc
	      stype=type & sthaw=thaw
	      flx=newflx 	;renorm the amplitudes!
	      ;
	      fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	        perrm=perrm,werrm=werrm,ferrm=ferrm,$
	        perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist
	      updatelist=1
	    endif else $			;MOAMP>0)
	      widget_control,txts_stts,set_value='All amplitudes frozen'
	  endif else $				;TOT_MOD > 0)
	    widget_control,txts_stts,set_value='Unable to renormalize'
        endif else $			;NCOM>0)
	  widget_control,txts_stts,set_value='No models to renormalize'
      endif else widget_control,txts_stts,set_value=$
	'Normalize to 0?!  Surely ye jest?'

    end						;renormalize}

   'menu_eror': ;compute assymetric errors -- this is a drop down menu
    'fts_steppar': begin ; use erors() StepPar
      widget_control,txts_stts,set_value='Working...'
      ;	initialize stuff
      algo=widget_info(fits_algo,/droplist_select)
      widget_control,fits_dchi,get_value=delchi & dchisq=delchi(0)
      z=y-conlev & xx=x(oo) & zz=z(oo) & zsig=sqrt(sigy(oo)^2+consig(oo)^2)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      widget_control,pars_list,get_value=list
      ncom=n_elements(pos)
      ;
      if ncom gt 0 then begin		;(if there are models to fit..
	thaw=intarr(3*ncom) & ii=lindgen(ncom)
	thaw(3*ii+0)=list(4*ii+1)		;positions
	thaw(3*ii+1)=list(4*ii+2)		;widths
	thaw(3*ii+2)=list(4*ii+3)		;fluxes
	freeze=where(thaw eq 0,nfrozen)
	aa=fltarr(3L*ncom)		;define the parameter array
	aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	if not keyword_set(alltie) then begin	;define prevailing constraints
	  if keyword_set(ties) then alltie=ties else alltie=''
	endif
	;
	nfree=3*ncom & nfree=nfree-nfrozen
	if nfree ge 1 then begin	;(any free components?
	  ;	save the old values
	  spos=pos & swdt=wdt & sflx=flx*dwvl
	  sperr=perr & swerr=werr & sferr=ferr*dwvl
	  sperrm=perrm & swerrm=werrm & sferrm=ferrm
	  sperrc=perrc & swerrc=werrc & sferrc=ferrc
	  stype=type & sthaw=thaw & sepithet=epithet
	  ;
	  ;	call ERORS 
	  erors,xx,zz,aa,erru,errl,ysig=zsig,freeze=freeze,dchi=dchisq,$
	    algo=algo,yfunc=ymod,erra=erra, funcs=funcs,function_name=funcs,$
	    ties=alltie,type=type,normflx=znorm,rmfstr=rmfstr,x_seg=x_seg,/dumb, _extra=e            
	  ;	update the errors as necessary
	  for i=0L,ncom-1L do begin		;{for each component
	    jjp=0 & jjw=0 & jjf=0	;do nothing
	    if thaw(3*i+0) gt 0 then jjp=1 else if erra(3*i+0) ne 0 then jjp=1
	    if thaw(3*i+1) gt 0 then jjw=1 else if erra(3*i+1) ne 0 then jjw=1
	    if thaw(3*i+2) gt 0 then jjf=1 else if erra(3*i+2) ne 0 then jjf=1
	    pos(i)=aa(3L*i+0)
	    if jjp gt 0 then begin
	      if erru(3L*i+0) ne 0 or errl(3L*i+0) ne 0 then begin
	        perr(i)=erru(3L*i+0)-pos(i)
		perrm(i)=pos(i)-errl(3L*i+0)
		perrc(i)=dchisq
	      endif
	    endif
	    wdt(i)=aa(3L*i+1)
	    if jjw gt 0 then begin
	      if erru(3L*i+1) ne 0 or errl(3L*i+1) ne 0 then begin
	        werr(i)=erru(3L*i+1)-wdt(i)
		werrm(i)=wdt(i)-errl(3L*i+1)
		werrc(i)=dchisq
	      endif
	    endif
	    flx(i)=aa(3L*i+2)/dwvl
	    if jjf gt 0 then begin
	      if erru(3L*i+2) ne 0 or errl(3L*i+2) ne 0 then begin
	        ferr(i)=erru(3L*i+2)/dwvl-flx(i)
		ferrm(i)=flx(i)-errl(3L*i+2)/dwvl
		ferrc(i)=dchisq
	      endif
	    endif
	  endfor				;I=0,NCOM-1}
	  ;
	  ;	update the widget state
	  widget_control,txts_stts,set_value='Done finding projected errors'
	  fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	    perrm=perrm,werrm=werrm,ferrm=ferrm,$
	    perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	  ;
	  updatelist=1
	  ;
	  ;	plot
	  ;plot,xx,zz,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
	  plot,xx,zz,xr=oxr,yr=oyr,/xs,/ys,/nodata
	  oplot,xx,zz,psym=10,col=icspc,thick=tspc
	  oplot,xx,ymod,col=icmod,thick=tmod
          ;
          ;      make a note for historical purposes      
          ;
          if n_tags(histr) gt 0 then begin                ;(HISTR
            tagnam=tag_names(histr) & ntagnam=n_elements(tagnam) 
            ihist = ihist+1L 
            tmphistr = create_struct('NFIT',nfit,'POS',pos,'WDT',wdt,'FLX',flx,$
            'CHISQ',total(((ymod-zz)^2)/(zsig)^2),'THAW',thaw,'NFREE',npt,'TYPE',$
            type,'ALGO_TYPE',algo_type,'TIES',alltie,'XRANGE',oxr,'BESTPARM',aa,$
            'ERRL',errl,'ERRU',erru,'freeze',freeze,'algo',algo,'erra',erra)  
            histr=create_struct(histr,'ERROR_fitpar'+strtrim(ihist,2),tmphistr)
         endif
	  if is_keyword_set(zerr) then begin
	    moo=n_elements(oo)
	    for ier=0L,moo-1L do oplot,zz(ier)+zsig(ier)*zerr*[-1.,1.]
	  endif
	endif else widget_control,txts_stts,set_value=$		;NFREE>1)
	  'Cannot project errors without at least 1 thawed parameter'
      endif else widget_control,txts_stts,set_value=$		;NCOM>0)
	'No parameters to find errors on'
      end					;use StepPar}
 

   'fts_mcerror':begin ; mcerror() Monte-Carlo
      widget_control,txts_stts,set_value='Working...'
      ;	initialize stuff
      algo=widget_info(fits_algo,/droplist_select)
      widget_control,fits_dchi,get_value=delchi & dchisq=delchi(0)
      z=y-conlev & xx=x(oo) & zz=z(oo) & zsig=sqrt(sigy(oo)^2+consig(oo)^2)
      dwvl=median(abs(xx(1:*)-xx)) & if dwvl eq 0 then dwvl=1.
      if keyword_set(intens) then dwvl=1.
      widget_control,pars_list,get_value=list
      ncom=n_elements(pos)
      ;
      if ncom gt 0 then begin		;(if there are models to fit..
	thaw=intarr(3*ncom) & ii=lindgen(ncom)
	thaw(3*ii+0)=list(4*ii+1)		;positions
	thaw(3*ii+1)=list(4*ii+2)		;widths
	thaw(3*ii+2)=list(4*ii+3)		;fluxes
	freeze=where(thaw eq 0,nfrozen)
	aa=fltarr(3L*ncom)		;define the parameter array
	aa(3L*ii+0)=pos(ii) & aa(3L*ii+1)=wdt(ii) & aa(3L*ii+2)=flx(ii)*dwvl
	if not keyword_set(alltie) then begin	;define prevailing constraints
	  if keyword_set(ties) then alltie=ties else alltie=''
	endif
	nfree=3*ncom & nfree=nfree-nfrozen
	if nfree ge 1 then begin	;(any free components?
	  ;	save the old values
	  spos=pos & swdt=wdt & sflx=flx*dwvl
	  sperr=perr & swerr=werr & sferr=ferr*dwvl
	  sperrm=perrm & swerrm=werrm & sferrm=ferrm
	  sperrc=perrc & swerrc=werrc & sferrc=ferrc
	  stype=type & sthaw=thaw & sepithet=epithet
	  ;	call MCERROR but check if MPFIT used first
          if algo eq 2 then begin 
             PARINFO=replicate({fixed:0,limited:[0,0],limits:[0,0]},n_elements(aa))
             for j=0,nfrozen-1 do parinfo(freeze(j)).fixed(0)=1
	        userfuncs=funcs
	        functarg=create_struct(e,'normflx',znorm,'type',type)
	     if strpos(userfuncs,'_f',0) lt 0 then userfuncs=userfuncs+'_f'
          endif  
          if n_tags(histr) gt 0 then mcsims=1
          mcerror,xx,zz,aa,erru,errl,ysig=zsig,freeze=freeze,dchi=dchisq,$
	    algo=algo,yfunc=ymod,erra=erra, funcs=funcs,function_name=funcs,$
	    ties=alltie,type=type,normflx=znorm,mpfunc=userfuncs,x_seg=x_seg, rmfstr=rmfstr,$
            parinfo=parinfo,functarg=functarg,mcsims=mcsims,/dumb, _extra=e	 
	  ;	update the errors as necessary
	  for i=0L,ncom-1L do begin		;{for each component
	    jjp=0 & jjw=0 & jjf=0	;do nothing
	    if thaw(3*i+0) gt 0 then jjp=1 else if erra(3*i+0) ne 0 then jjp=1
	    if thaw(3*i+1) gt 0 then jjw=1 else if erra(3*i+1) ne 0 then jjw=1
	    if thaw(3*i+2) gt 0 then jjf=1 else if erra(3*i+2) ne 0 then jjf=1
	    pos(i)=aa(3L*i+0)
	    if jjp gt 0 then begin
	      if erru(3L*i+0) ne 0 or errl(3L*i+0) ne 0 then begin
	        perr(i)=erru(3L*i+0)-pos(i)
		perrm(i)=pos(i)-errl(3L*i+0)
		perrc(i)=dchisq
	      endif
	    endif
	    wdt(i)=aa(3L*i+1)
	    if jjw gt 0 then begin
	      if erru(3L*i+1) ne 0 or errl(3L*i+1) ne 0 then begin
	        werr(i)=erru(3L*i+1)-wdt(i)
		werrm(i)=wdt(i)-errl(3L*i+1)
		werrc(i)=dchisq
	      endif
	    endif
	    flx(i)=aa(3L*i+2)/dwvl
	    if jjf gt 0 then begin
	      if erru(3L*i+2) ne 0 or errl(3L*i+2) ne 0 then begin
	        ferr(i)=erru(3L*i+2)/dwvl-flx(i)
		ferrm(i)=flx(i)-errl(3L*i+2)/dwvl
		ferrc(i)=dchisq
	      endif
	    endif
	  endfor				;I=0,NCOM-1}
	  ;
	  ;	update the widget state
	  widget_control,txts_stts,set_value='Done finding projected errors'
	  fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	    perrm=perrm,werrm=werrm,ferrm=ferrm,$
	    perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	  ;
	  updatelist=1
	  ;
	  ;	plot
	  ;plot,xx,zz,xr=oxr,yr=oyr,psym=10,col=icspc,thick=tspc,/xs,/ys
	  plot,xx,zz,xr=oxr,yr=oyr,/xs,/ys,/nodata
	  oplot,xx,zz,psym=10,col=icspc,thick=tspc
	  oplot,xx,ymod,col=icmod,thick=tmod
          ;
          ;      make a note for historical purposes      
          ;
          if n_tags(histr) gt 0 then begin                ;(HISTR
            tagnam=tag_names(histr) & ntagnam=n_elements(tagnam) 
            ihist = ihist+1L 
         tmphistr = create_struct('NFIT',nfit,'POS',pos,'WDT',wdt,'FLX',flx,$
		'CHISQ',total(((ymod-zz)^2)/(zsig)^2),'THAW',thaw,'NFREE',npt,$
                'TYPE',type,'ALGO_TYPE',algo_type,'TIES',alltie,'XRANGE',oxr,'BESTPARM',aa,$
                'ERRL',errl,'ERRU',erru,'freeze',freeze,'algo',algo,'erra',$
                erra,'MCSIMS',mcsims)  
             histr=create_struct(histr,'ERROR_MC'+strtrim(ihist,2),tmphistr)
         endif
	  if is_keyword_set(zerr) then begin
	    moo=n_elements(oo)
	    for ier=0L,moo-1L do oplot,zz(ier)+zsig(ier)*zerr*[-1.,1.]
	  endif
	endif else widget_control,txts_stts,set_value=$		;NFREE>1)
	  'Cannot project errors without at least 1 thawed parameter'
      endif else widget_control,txts_stts,set_value=$		;NCOM>0)
	'No parameters to find errors on'
      end

    'fits_dchi': begin				;{read dCHISQ threshold
      widget_control,fits_dchi,get_value=delchi & dchisq=delchi(0)
      if keyword_set(odchisq) then begin
        widget_control,txts_stts,set_value='Old value of dCHI: '+$
	  strtrim(odchisq,2)
      endif else widget_control,txts_stts,set_value='dCHI='+strtrim(dchisq,2)
      odchisq=dchisq
    end						;read dCHISQ threshold}

    'fits_algo': begin				;{change fitting algorithm
      algo=widget_info(fits_algo,/droplist_select)
      algo_type='LevMarq+SVD'
      case algo of				;{which algorithm to use?
	0: begin
          algo_type='LevMarq+SVD'
          c1='Using Levenberg-Marquardt with SVD (FIT_LEVMAR)'
	end
	1: begin
          algo_type='IDL-Curvefit'
          c1='gradient expansion non-linear least-squares (CURVE_FIT)'
	end
	2: begin
	  algo_type='MPFIT'
	  c1="hooks into Craig Markwardt's MPFIT suite"
	end
	else: begin
          algo_type='LevMarq+SVD'
          c1='Sorry, this method not implemented; using '+algo_type
	end
      endcase					;ALGO type}
      widget_control,txts_stts,set_value=c1
    end						;FITS_ALGO}

    ;							fourth row
    'pars_list': begin				;{freeze/thaw parameters
      widget_control,txts_stts,set_value=''
      widget_control,pars_list,get_value=list
      ibtn=event.value
      ncom=n_elements(list)/4
      if ncom gt 0 then begin
	for i=0,ncom-1 do begin
	  i4=4L*i & i3=i4+[1,2,3]
	  if ibtn eq i4 and list(i4) eq 0 then list(i3)=0
	  if ibtn eq i4 and list(i4) eq 1 then list(i3)=1
	  if ibtn eq i4+1 and list(ibtn) eq 1 then list(i4)=1
	  if ibtn eq i4+2 and list(ibtn) eq 1 then list(i4)=1
	  if ibtn eq i4+3 and list(ibtn) eq 1 then list(i4)=1
	  if (ibtn eq i4+1 or ibtn eq i4+2 or ibtn eq i4+3) and total(list(i3)) eq 0 then list(i4)=0
	endfor
        widget_control,pars_list,set_value=list
	thaw=intarr(3*ncom) & ii=lindgen(ncom)
	thaw(3*ii+0)=list(4*ii+1)		;positions
	thaw(3*ii+1)=list(4*ii+2)		;widths
	thaw(3*ii+2)=list(4*ii+3)		;fluxes
      endif
      widget_control,txts_stts,set_value='Component '+strtrim(ibtn/4,2)
      ;	update the left column on parameters..
      widget_control,comp_indx,set_value=ibtn/4
      widget_control,comp_posv,set_value=pos(ibtn/4)
      widget_control,comp_wdtv,set_value=wdt(ibtn/4)
      widget_control,comp_flxv,set_value=flx(ibtn/4)
      widget_control,comp_type,set_value=type(ibtn/4)
      widget_control,comp_sign,set_value=epithet(ibtn/4)
    end						;freeze/thaw parameters}
    'comp_indx': begin				;{edit component..
      widget_control,txts_stts,set_value='edit component'
      widget_control,pars_list,get_value=list
      widget_control,comp_indx,get_value=icomp
      ncom=n_elements(pos)
      if icomp ge 0 and icomp lt ncom then begin
	widget_control,comp_posv,set_value=pos(icomp)
        widget_control,comp_wdtv,set_value=wdt(icomp)
        widget_control,comp_flxv,set_value=flx(icomp)
        widget_control,comp_type,set_value=type(icomp)
        widget_control,comp_sign,set_value=epithet(icomp)
      endif else begin
	widget_control,comp_posv,set_value=-999
        widget_control,comp_wdtv,set_value=-999
        widget_control,comp_flxv,set_value=-999
        widget_control,comp_type,set_value='***'
        widget_control,comp_sign,set_value='***'
      endelse
    end						;edit component}
    'comp_posv': begin				;{edit position value
      widget_control,txts_stts,set_value=''
      widget_control,pars_list,get_value=list
      widget_control,comp_indx,get_value=icomp
      widget_control,comp_posv,get_value=val
      if icomp ge 0 and icomp lt n_elements(pos) then begin
        widget_control,txts_stts,set_value='Old value = '+strtrim(pos(icomp),2)
	pos(icomp)=val
        fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
	  perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	updatelist=1
      endif
    end						;edit position value}
    'comp_wdtv': begin				;{edit width value
      widget_control,txts_stts,set_value=''
      widget_control,pars_list,get_value=list
      widget_control,comp_indx,get_value=icomp
      widget_control,comp_wdtv,get_value=val
      if icomp ge 0 and icomp lt n_elements(wdt) then begin
        widget_control,txts_stts,set_value='Old value = '+strtrim(wdt(icomp),2)
	wdt(icomp)=val
        fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
	  perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	updatelist=1
      endif
    end						;edit width value}
    'comp_flxv': begin				;{edit flux value
      widget_control,txts_stts,set_value=''
      widget_control,pars_list,get_value=list
      widget_control,comp_indx,get_value=icomp
      widget_control,comp_flxv,get_value=val
      if icomp ge 0 and icomp lt n_elements(flx) then begin
        widget_control,txts_stts,set_value='Old value = '+strtrim(flx(icomp),2)
	flx(icomp)=val
        fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
	  perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	updatelist=1
      endif
    end						;edit flux value}
    'comp_type': begin				;{edit component type
      widget_control,txts_stts,set_value=''
      widget_control,pars_list,get_value=list
      widget_control,comp_indx,get_value=icomp
      widget_control,comp_type,get_value=val
      if icomp ge 0 and icomp le n_elements(flx) then begin
        widget_control,txts_stts,set_value='Old value = '+type(icomp)
	type(icomp)=val
        fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
	  perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	updatelist=1
      endif
    end						;edit component type}
    'comp_sign': begin				;{edit label
      widget_control,txts_stts,set_value=''
      widget_control,pars_list,get_value=list
      widget_control,comp_indx,get_value=icomp
      widget_control,comp_sign,get_value=val
      if icomp ge 0 and icomp le n_elements(flx) then begin
        widget_control,txts_stts,set_value='Old value = '+epithet(icomp)
	epithet(icomp)=val(0)
        fmtstuff,pos,wdt,flx,thaw,type,epithet,perr,werr,ferr,$
	  perrm=perrm,werrm=werrm,ferrm=ferrm,$
	  perrc=perrc,werrc=werrc,ferrc=ferrc, plist=plist,list=list
	updatelist=1
      endif
    end						;edit label}

    ;							fourth row
    'ties_acts': 	;handle ties among parameters -- drop-down menu
    'acts_frcp': begin				;{constrain positions
      widget_control,txts_stts,set_value='Limit range for line positions'
      widget_control,ties_dpos,get_value=dpos
      ncom=n_elements(pos) & if ncom gt 0 then tiep=strarr(ncom)
      for i=0,ncom-1 do begin
	ic=strtrim(3*i,2)
	tiep(i)='a'+ic+' = (( a'+ic+' < ('+strtrim(pos(i)+dpos,2)+') ) > ('+$
		strtrim(pos(i)-dpos,2)+') )'
      endfor
      updateties=1
    end						;constrain positions}
    'acts_frcw': begin				;{constrain widths
      widget_control,txts_stts,set_value='Limit range for line widths'
      widget_control,ties_dwdt,get_value=dwdt
      ncom=n_elements(wdt) & if ncom gt 0 then tiew=strarr(ncom)
      for i=0,ncom-1 do begin
	ic=strtrim(3*i+1,2)
	tiew(i)='a'+ic+' = (( a'+ic+' < ('+strtrim(wdt(i)+dwdt,2)+') ) > ('+$
		strtrim(wdt(i)-dwdt,2)+') )'
      endfor
      updateties=1
    end						;constrain widths}
    'acts_frcf': begin				;{constrain fluxes
      widget_control,txts_stts,set_value='Limit range for line fluxes'+$
	' (fraction of SIGMA)'
      widget_control,ties_dflx,get_value=dflx
      ncom=n_elements(flx) & if ncom gt 0 then tief=strarr(ncom)
      zsig=sqrt(sigy^2+consig^2)
      for i=0,ncom-1 do begin
	ic=strtrim(3*i+2,2)
	tmp=min(abs(pos-x),j)
	fup=flx(i)+dflx*zsig(j) & fdn=flx(i)-dflx*zsig(j) > 0
	if not keyword_set(dwvl) then dwvl=1.
	tief(i)='a'+ic+' = (( a'+ic+' < ('+strtrim(fup*dwvl,2)+') ) > ('+$
		strtrim(fdn*dwvl,2)+') )'
      endfor
      updateties=1
    end						;constrain fluxes}
    'ties_dele': ;delete ties -- drop down menu
    'dele_some': begin				;{delete selected
      widget_control,txts_stts,set_value='Deleting selected constraints'
      ntie=n_elements(itie)
      for i=0,ntie-1 do begin
	if tietyp(i) eq 0 then begin
	  iit=lindgen(n_elements(ties)) & oi=where(iit ne itie(i),moi)
	  if moi eq 0 then ties='' else ties=ties(oi)
	endif
	if tietyp(i) eq 1 then tiep(itie(i)-nts)=''
	if tietyp(i) eq 2 then tiew(itie(i)-nts-ntp)=''
	if tietyp(i) eq 3 then tief(itie(i)-nts-ntp-ntw)=''
      endfor
      updateties=1
    end						;delete selected}
    'dele_alls': begin				;{delete all
      widget_control,txts_stts,set_value='Roam free, ye parameter brethren!'
      ties=''
      if n_elements(tiep) gt 0 then tiep(*)=''
      if n_elements(tiew) gt 0 then tiew(*)=''
      if n_elements(tief) gt 0 then tief(*)=''
      updateties=1
    end						;delete all}
    'ties_edit': begin				;{edit chosen constraint
      nts=n_elements(ties) & if strtrim(ties(0),2) eq '0' then nts=0
      ntp=n_elements(tiep) & ntw=n_elements(tiew) & ntf=n_elements(tief)
      if keyword_set(ktie) then begin
	widget_control,ties_edit,get_value=val
        if tietyp(0) eq 0 then ties(ktie-1)=val
        if tietyp(0) eq 1 then tiep(ktie-nts-1)=val
        if tietyp(0) eq 2 then tiew(ktie-nts-ntp-1)=val
        if tietyp(0) eq 3 then tief(ktie-nts-ntp-ntw-1)=val
      endif
      updateties=1
    end						;edit chosen constraint}
    'ties_adds': begin				;{add new constraint
      widget_control,txts_stts,set_value='Remember to freeze the appropriate parameter if necessary'
      widget_control,ties_adds,get_value=val
      nts=n_elements(ties)
      if nts eq 0 then ties=[val] else begin
        if nts eq 1 and strtrim(ties(0),2) eq '0' then ties=val else $
		ties=[ties,val]
      endelse
      updateties=1
    end						;add new constraint}
    'ties_dpos': begin				;{POS-dPOS < POS < POS+dPOS
      widget_control,txts_stts,set_value='Constrain positions to lie in this range'
      widget_control,ties_dpos,get_value=dpos
      updateties=1
    end						;POS-dPOS < POS < POS+dPOS}
    'ties_dwdt': begin				;{WDT-dWDT < WDT < WDT+dWDT
      widget_control,txts_stts,set_value='Constrain widths to lie in this range'
      widget_control,ties_dwdt,get_value=dwdt
      updateties=1
    end						;WDT-dWDT < WDT < WDT+dWDT}
    'ties_dflx': begin				;{FLX-dFLX < FLX < FLX+dFLX
      widget_control,ties_dflx,get_value=dflx
      widget_control,txts_stts,set_value='Constrain fluxes to lie within '+$
      	strtrim(dflx,2)+'*SIGMA_INTENSITY'
      updateties=1
    end						;FLX-dFLX < FLX < FLX+dFLX}

    ;							fifth row
    'tlst_list': begin				;{list of constraints
      widget_control,txts_stts,set_value=''
      itie=widget_info(tlst_list,/list_select)
      j=itie(0) & ktie=j+1 & set_list_top=j
      nts=n_elements(ties) & if strtrim(ties(0),2) eq '0' then nts=0
      ntp=n_elements(tiep) & ntw=n_elements(tiew) & ntf=n_elements(tief)
      tietyp=0*itie-1
      ot=where(itie le nts-1,mot) & if mot gt 0 then tietyp(ot)=0
      ot=where(itie gt nts-1,mot) & if mot gt 0 then tietyp(ot)=1
      ot=where(itie gt nts+ntp-1,mot) & if mot gt 0 then tietyp(ot)=2
      ot=where(itie gt nts+ntp+ntw-1,mot) & if mot gt 0 then tietyp(ot)=3
      jpar=tietyp(0)
      if jpar eq 0 then widget_control,ties_edit,set_value=ties(j)
      if jpar eq 1 then widget_control,ties_edit,set_value=tiep(j-nts)
      if jpar eq 2 then widget_control,ties_edit,set_value=tiew(j-nts-ntp)
      if jpar eq 3 then widget_control,ties_edit,set_value=tief(j-nts-ntp-ntw)
    end						;list of constraints}

    else: ;nothing
  endcase					;handle EV}

  if keyword_set(updatelist) then begin
    ;repeat call to CW_BGROUP if a model has been added or deleted
    if ev eq 'modl_addc' or ev eq 'modl_delc' then begin
      widget_control,pars_list,/destroy
      widget_control,event.top,update=0
      label_top='[Thaw?] :  Component : Parameter : Tag : Value +- Error [@dCHI]'
      pars_list=cw_bgroup(base_pars,plist,/column,/nonexclusive,$
        uvalue='pars_list',/return_index,/scroll,y_scroll_size=340,$
        x_scroll_size=540,set_value=list,font='8x13',label_top=label_top,$
        ids=list_ids)
      widget_control,event.top,/update
      widget_control,pars_list,/realize,/map
    endif
    for i=0L,n_elements(list_ids)-1L do $
	widget_control,list_ids(i),set_value=plist(i)
    updatelist=0
  endif

  if keyword_set(updateties) then begin
    ;		figure out the ties
    nts=n_elements(ties) & if strtrim(ties(0),2) eq '0' then nts=0
    ntp=n_elements(tiep) & ntw=n_elements(tiew) & ntf=n_elements(tief)
    if ntp+ntw+ntf+nts gt 0 then begin
      alltie=['']
      if nts gt 0 then alltie=[alltie,ties]
      if ntp gt 0 then alltie=[alltie,tiep]
      if ntw gt 0 then alltie=[alltie,tiew]
      if ntf gt 0 then alltie=[alltie,tief]
      alltie=alltie(1:*)
    endif else alltie=''
    widget_control,tlst_list,/destroy
    widget_control,event.top,update=0
    if not keyword_set(alltie) then tt='' else tt=alltie
    if not keyword_set(set_list_top) then set_list_top=0
    ;tlst_list=widget_list(dtls_tlst,/multiple,scr_xsize=500,scr_ysize=115,$
    ;  value=tt,uvalue='tlst_list')
    ;Antonio Maggio points out that keyword MULTIPLE is not allowed in
    ;older versions.  removing it causes little pain.
    tlst_list=widget_list(dtls_tlst,scr_xsize=500,scr_ysize=115,$
      value=tt,uvalue='tlst_list')
    widget_control,tlst_list,set_list_top=set_list_top
    widget_control,event.top,/update
    widget_control,tlst_list,/realize,/map,/hour
    updateties=0
  endif

  if n_tags(histr) gt 0 then histerix=histr
endwhile				;end of endless loop}

;	output
if keyword_set(alltie) then ties=alltie
timestamp=systime()
if n_tags(histr) gt 0 then histerix=create_struct(histr,'DONE',timestamp)

return
end
