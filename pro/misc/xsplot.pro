pro xsplot,data,model,ee,dsig,lcol,perbin=perbin,minct=minct,$
	resid=resid,fracres=fracres,delchi=delchi,ratio=ratio,$
	data0=data0,mod0=mod0,ee0=ee0,dsig0=dsig0,lcol0=lcol0,$
	data1=data1,mod1=mod1,ee1=ee1,dsig1=dsig1,lcol1=lcol1,$
	data2=data2,mod2=mod2,ee2=ee2,dsig2=dsig2,lcol2=lcol2,$
	data3=data3,mod3=mod3,ee3=ee3,dsig3=dsig3,lcol3=lcol3,$
	data4=data4,mod4=mod4,ee4=ee4,dsig4=dsig4,lcol4=lcol4,$
	data5=data5,mod5=mod5,ee5=ee5,dsig5=dsig5,lcol5=lcol5,$
	data6=data6,mod6=mod6,ee6=ee6,dsig6=dsig6,lcol6=lcol6,$
	data7=data7,mod7=mod7,ee7=ee7,dsig7=dsig7,lcol7=lcol7,$
	data8=data8,mod8=mod8,ee8=ee8,dsig8=dsig8,lcol8=lcol8,$
	data9=data9,mod9=mod9,ee9=ee9,dsig9=dsig9,lcol9=lcol9,$
	ylog=ylog,xtitle=xtitle,ytitle=ytitle,xthick=xthick,ythick=ythick,$
	title=title,subtitle=subtitle,charsize=charsize,charthick=charthick,$
	yrange=yrange,zrange=zrange,xstyle=xstyle,psym=psym,color=color,$
	wsplit=wsplit, verbose=verbose,$
	_extra=e
;+
;procedure	xsplot
;	make XSPEC like double-panel plots for given counts spectrum
;	and associated model
;
;syntax
;	xsplot,data,model,egrid,dataerror,linecolor,/perbin,minct=minct,$
;	dataN=dataN,modN=modN,eeN=eeN,dsigN=dsigN,lcolN=lcolN,$	;N=0..9
;	/resid,/fracres,/delchi,/ratio,verbose=vv,wsplit=wsplit,$
;	/ylog,xtitle=xtitle,ytitle=ytitle,xthick=xthick,ythick=ythick,$
;	title=title,subtitle=subtitle,charsize=charsize,charthick=charthick,$
;	yrange=yrange,zrange=zrange, /xlog,xrange=xrange,thick=thick,$
;	/halfit,eps=eps, /nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol
;
;parameters
;	data	[INPUT; required] array of observed counts
;	model	[INPUT] array of model counts
;		* must match size of DATA
;		- extra bins are ignored
;		- missing bins are assumed to be 0
;	ee	[INPUT] x-axis grid defining DATA and MODEL
;		* if size is N(DATA), assumed to be mid-bin values
;		* if size is N(DATA)+1, assumed to be bin boundaries
;		* if not given, 1-based integer indices are used
;	dsig	[INPUT] error bars on DATA
;		* if not given, assumed to be Gehrels' approximation
;		  of Poisson counts error, sqrt(DATA+0.75)+1
;		* if set to a scalar, then
;		  - if -ve, errors are assumed constant = abs(DSIG)
;		  - if integer, assumed to be in units of SIGMAs, i.e.,
;		    DSIG*default are plotted
;		  - if not integer and 0<DSIG<0.5, then assumed to be a
;		    fractional error DSIG*abs(DATA)
;		  - if not integer and 0.5<DSIG<10, then assumed to be
;		    in units of SIGMAs, i.e., DSIG*default are plotted
;		  - if not integer and DSIG>10, assumed to be a percentage
;		    error, (DSIG/100)*abs(DATA)
;	lcol	[INPUT] 3-element byte array specifying the line colors
;		  for MODEL, DATA, DSIG in that order
;		* if incompletely given, missing elements are assumed
;		  to be [0,10,20] or [255,245,235] depending on whether
;		  !p.background is 255 or 0.
;
;keywords
;	perbin	[INPUT] if set, will treat the input as counts/bin
;		and plot accordingly.  if not set, will compute the
;		bin widths and plot the density instead
;		* if EE is not given, the bin widths will be 1 and
;		  this has no relevance except to the YTITLE
;	minct	[INPUT] if given, regroups the input DATA on the fly
;		to have a minimum number of counts in each bin, and
;		then regroups all the other arrays to the same grid
;		* forces PERBIN to be ignored
;	verbose	[INPUT] controls chatter
;	dataN	[INPUT] for N=0..9, additional DATA arrays
;		* size must match DATA
;	modN	[INPUT] for N=0..9, additional MODEL arrays
;		* size must match DATA
;	eeN	[INPUT] for N=0..9, additional EE arrays
;		* if not given, default is to use EE
;	dsigN	[INPUT] for N=0..9, additional DSIG arrays
;		* size must match DATA
;	lcolN	[INPUT] for N=0..9, line colors, by default, N*25+[25,35,45]
;	resid	[INPUT] if set, plots (DATA-MODEL) in lower panel
;	fracres	[INPUT] if set, plots (DATA-MODEL)/MODEL in lower panel
;	delchi	[INPUT] if set, plots (DATA-MODEL)/DSIG in lower panel
;	ratio	[INPUT] if set, plots DATA/MODEL in lower panel
;		* if more than one of RESID,FRACRES,DELCHI,RATIO are set,
;		  the priority is in that order.  i.e., the later ones
;		  in the list override anything given before it.
;		* if none are set, assumes /DELCHI
;	wsplit	[INPUT] fraction by which to split the upper and lower windows
;		* default is 0.35
;		* ignored if <0.05 or > 0.95
;	ylog	[INPUT] catch and release to upper PLOT
;	xtitle	[INPUT] catch and release to lower PLOT
;	ytitle	[INPUT] catch and release to upper PLOT
;	xthick	[INPUT] catch and release to PLOT
;	ythick	[INPUT] catch and release to PLOT
;	title	[INPUT] catch and release to upper PLOT
;	charsize [INPUT] catch and release to PLOT
;	charthick [INPUT] catch and release to PLOT
;	subtitle [INPUT] catch and release to lower OPLOT
;	yrange	[INPUT] catch and release to upper PLOT
;	zrange	[INPUT] yrange catch and release to lower PLOT
;	xstyle	[INPUT] catch and discard
;	psym	[INPUT] catch and discard
;	color	[INPUT] catch and discard
;	_extra	[INPUT ONLY] pass defined keywords to subroutines:
;		PLOT: XLOG, XRANGE, THICK
;		STAMPLE: NUTHIN, NODAY, NOTIME, NOUSER, NOPACK, STACOL
;		MID2BOUND: HALFIT, EPS
;
;examples
;	xsplot,data,model
;	xsplot,data,model,/verbose	;to stamp creator name and time
;	xsplot,data,model,wvl		;to include physical x-axis grid
;	xsplot,data,model,wvl,0		;to not plot error bars
;	xsplot,counts,model,minct=25	;adaptively group the counts on the fly
;	xsplot,data,model,wvl,1,150	;to change color of model curve
;	xsplot,data,model,data3=data3,mod3=mod3	;overplot another set
;	xsplot,data,model,mod0=mod0,mod1=mod1	;overplot alternate models
;	xsplot,data,model,mod9=mod9,lcol9=5	;alt model with specific color
;	xsplot,counts,model,keV,/perbin	;plot counts/bin rather than counts/keV
;	xsplot,data,model,/ratio	;plot data/model in lower panel
;	xsplot,data,model,/delchi	;plot delta chi in lower panel
;	xsplot,data,model,/resid	;plot data-model in lower panel
;	xsplot,data,model,/fracres	;plot (data-model)/model in lower panel
;	xsplot,data,model,wvl,/xlog,/ylog,xrange=xr,yrange=yr,zrange=zr
;
;	x=findgen(100)+1 & y=50*(sin(x/10)+1.1) & d=0*y & plot,x,y
;	for i=0,99 do d[i]=randomu(seed,poisson=y[i]) & xsplot,d,y,x
;
;subroutines
;	STAMPLE [SETSYSVAL, LEGALVAR()]
;	MID2BOUND()
;	REGROUP()
;
;history
;	vinay kashyap (Oct09)
;-

;	usage
ok='ok' & np=n_params() & nd=n_elements(data)
if np eq 0 then ok='Insufficient parameters' else $
 if nd eq 0 then ok='DATA undefined' else $
  if nd eq 1 then ok='DATA cannot be plotted'
if ok ne 'ok' then begin
  print,'Usage: xsplot,data,model,egrid,dataerror,linecolor,/perbin,minct=minct,$'
  print,'       dataN=dataN,modN=modN,eeN=eeN,dsigN=dsigN,lcolN=lcolN,$	;{N=0..9}'
  print,'       /resid,/fracres,/delchi,/ratio,verbose=vv,wsplit=wsplit,$'
  print,'       /ylog,xtitle=xtitle,ytitle=ytitle,xthick=xthick,ythick=ythick,$'
  print,'       title=title,subtitle=subtitle,charsize=charsize,charthick=charthick,$'
  print,'       yrange=yrange,zrange=zrange, /xlog,xrange=xrange,thick=thick,$'
  print,'       /halfit,eps=eps, /nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol'
  print,'  make XSPEC-like double panel plots'
  if np ne 0 then message,ok,/informational
  return
endif

;	optional parameters

yy=data
;
nm=n_elements(model) & nx=n_elements(ee) & ns=n_elements(dsig) & nc=n_elements(lcol)
;
mm=0.*yy & km=nm<(nd-1L) & if nm gt 0 then mm[0L:km]=model[0L:km]
;
defxt='index' & xx=findgen(nd+1)+1
if nx eq nm then begin & xx=mid2bound(ee,_extra=e) & defxt='X' & endif
if nx eq nm+1L then begin & xx=ee & defxt='X' & endif
;
ysig=sqrt(abs(yy)+0.75)+1. & ks=ns<(nd-1L)
if ns gt 0 then begin
  if ns gt 1 then ysig[0L:ks]=dsig[0L:ks] else begin
    if dsig[0] lt 0 then ysig[*]=abs(dsig[0])
    if size(dsig,/type) ge 4 then begin
      if dsig[0] ge 0 and dsig[0] lt 0.5 then ysig=dsig[0]*abs(yy)
      if dsig[0] ge 0.5 and dsig[0] le 10 then ysig=dsig[0]*ysig
      if dsig[0] gt 10 then ysig=(dsig[0]/100.)*abs(yy)
    endif else ysig=dsig[0]*ysig
  endelse
endif
;
defcol=[255,245,235]-!p.background & lcolarr=abs(defcol)
if nc eq 1 then lcolarr[0]=lcol[0] else $
 if nc eq 2 then lcolarr[0:1]=lcol else $
  if nc eq 3 then lcolarr=lcol else $
   if nc gt 3 then lcolarr=lcol[0:2]

;	keywords

vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
ctthr=0 & if keyword_set(minct) then ctthr=minct[0]
if keyword_set(ctthr) then begin
  if total(yy,/nan) gt ctthr then begin
    yyr=regroup(yy,ctthr,grid,/newgrid,iout=iout)
    mmr=regroup(mm,ctthr,grid)
    ssr=sqrt(regroup(ysig^2,ctthr,grid))
    xxr=xx[grid]
    if vv gt 0 then $
    	message,'regrouping the input arrays according to DATA',/informational
    yy=yyr & mm=mmr & ysig=ssr & xx=xxr
  endif else begin
    message,'unable to regroup the input arrays; check DATA',/informational
    ctthr=0
  endelse
endif
ny=n_elements(yy)
;
if keyword_set(perbin) and ctthr eq 0 then begin	;(if PERBIN is set _and_ arrays are not regrouped
  defyt='Y' & dxx=fltarr(nd)+1.
endif else begin					;/PERBIN&&CTTHR=0)(not PERBIN, or arrays are grouped
  defyt='intensity' & dxx=abs(xx[1:*]-xx)
endelse							;PERBIN=0||CTTHR>0)
;
defytb='(D-M)/SIG' & iresid=0 & ifrac=0 & idelchi=1 & iratio=0
if keyword_set(resid) then begin   & defytb='D-M'     & iresid=1 & ifrac=0 & idelchi=0 & iratio=0 & endif
if keyword_set(fracres) then begin & defytb='(D-M)/M' & iresid=0 & ifrac=1 & idelchi=0 & iratio=0 & endif
if keyword_set(delchi) then begin  & defytb='d!4v!X'  & iresid=0 & ifrac=0 & idelchi=1 & iratio=0 & endif
if keyword_set(ratio) then begin   & defytb='D/M'     & iresid=0 & ifrac=0 & idelchi=0 & iratio=1 & endif
;
winsplit=0.35
if keyword_set(wsplit) then begin
  if wsplit[0] ge 0.05 or wsplit[0] le 0.95 then winsplit=wsplit[0]
endif
;
if not keyword_set(ylog) then ylog=0
if keyword_set(xtitle) then defxt=xtitle
if keyword_set(ytitle) then defyt=ytitle

;	how many alternate sets?
nalt=0
;
nd0=n_elements(data0) & nm0=n_elements(mod0) & ns0=n_elements(dsig0)
if nd0 gt 0 or nm0 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd0<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data0[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm0<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod0[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns0<(nd-1L) & sig'+j+'[0L:k]=dsig0[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns0<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns0 gt 0 then dtmp=dsig0')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig0[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol0) then jnk=execute('col'+j+'=lcol0') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd1=n_elements(data1) & nm1=n_elements(mod1) & ns1=n_elements(dsig1)
if nd1 gt 0 or nm1 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd1<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data1[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm1<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod1[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns1<(nd-1L) & sig'+j+'[0L:k]=dsig1[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns1<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns1 gt 0 then dtmp=dsig1')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig1[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol1) then jnk=execute('col'+j+'=lcol1') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd2=n_elements(data2) & nm2=n_elements(mod2) & ns2=n_elements(dsig2)
if nd2 gt 0 or nm2 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd2<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data2[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm2<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod2[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns2<(nd-1L) & sig'+j+'[0L:k]=dsig2[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns2<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns2 gt 0 then dtmp=dsig2')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig2[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol2) then jnk=execute('col'+j+'=lcol2') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd3=n_elements(data3) & nm3=n_elements(mod3) & ns3=n_elements(dsig3)
if nd3 gt 0 or nm3 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd3<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data3[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm3<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod3[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns3<(nd-1L) & sig'+j+'[0L:k]=dsig3[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns3<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns3 gt 0 then dtmp=dsig3')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig3[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol3) then jnk=execute('col'+j+'=lcol3') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd4=n_elements(data4) & nm4=n_elements(mod4) & ns4=n_elements(dsig4)
if nd4 gt 0 or nm4 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd4<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data4[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm4<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod4[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns4<(nd-1L) & sig'+j+'[0L:k]=dsig4[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns4<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns4 gt 0 then dtmp=dsig4')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig4[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol4) then jnk=execute('col'+j+'=lcol4') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd5=n_elements(data5) & nm5=n_elements(mod5) & ns5=n_elements(dsig5)
if nd5 gt 0 or nm5 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd5<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data5[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm5<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod5[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns5<(nd-1L) & sig'+j+'[0L:k]=dsig5[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns5<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns5 gt 0 then dtmp=dsig5')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig5[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol5) then jnk=execute('col'+j+'=lcol5') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd6=n_elements(data6) & nm6=n_elements(mod6) & ns6=n_elements(dsig6)
if nd6 gt 0 or nm6 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd6<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data6[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm6<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod6[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns6<(nd-1L) & sig'+j+'[0L:k]=dsig6[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns6<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns6 gt 0 then dtmp=dsig6')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig6[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol6) then jnk=execute('col'+j+'=lcol6') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd7=n_elements(data7) & nm7=n_elements(mod7) & ns7=n_elements(dsig7)
if nd7 gt 0 or nm7 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd7<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data7[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm7<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod7[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns7<(nd-1L) & sig'+j+'[0L:k]=dsig7[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns7<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns7 gt 0 then dtmp=dsig7')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig7[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol7) then jnk=execute('col'+j+'=lcol7') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd8=n_elements(data8) & nm8=n_elements(mod8) & ns8=n_elements(dsig8)
if nd8 gt 0 or nm8 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd8<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data8[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm8<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod8[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns8<(nd-1L) & sig'+j+'[0L:k]=dsig8[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns8<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns8 gt 0 then dtmp=dsig8')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig8[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol8) then jnk=execute('col'+j+'=lcol8') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
nd9=n_elements(data9) & nm9=n_elements(mod9) & ns9=n_elements(dsig9)
if nd9 gt 0 or nm9 gt 0 then begin
  j=strtrim(nalt,2) & defcol=nalt*25+[25,35,45]
  jnk=execute('d'+j+'=0*yy & k=nd9<(nd-1L) & if k gt 0 then d'+j+'[0L:k]=data9[0L:k]')
  jnk=execute('m'+j+'=0*yy & k=nm9<(nd-1L) & if k gt 0 then m'+j+'[0L:k]=mod9[0L:k]')
  ;jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & k=ns9<(nd-1L) & sig'+j+'[0L:k]=dsig9[0L:k]')
  jnk=execute('sig'+j+'=sqrt(d'+j+'+0.75)+1 & ms=ns9<(nd-1L) & stmp=sig'+j+' & ytmp=d'+j+'& if ns9 gt 0 then dtmp=dsig9')
  if ms gt 0 then begin
    if ms gt 1 then jnk=execute('sig'+j+'[0L:k]=dsig9[0L:k]') else begin
      if dtmp[0] lt 0 then stmp[*]=abs(dtmp[0])
      if size(dtmp,/type) ge 4 then begin
        if dtmp[0] ge 0 and dtmp[0] lt 0.5 then stmp=dsig[0]*abs(ytmp)
	if dtmp[0] ge 0.5 and dtmp[0] le 10 then stmp=dsig[0]*stmp
	if dtmp[0] gt 10 then stmp=(dtmp[0]/100.)*abs(ytmp)
      endif else stmp=dtmp[0]*stmp
    endelse
  endif
  jnk=execute('sig'+j+'=stmp')
  if keyword_set(lcol9) then jnk=execute('col'+j+'=lcol9') else jnk=execute('col'+j+'=defcol')
  nalt=nalt+1
endif
;
if vv gt 0 then message,'number of extra curves to plot: '+strtrim(nalt,2),/informational

;	figure out what to plot
ydynrng=1e-5 & if keyword_set(drange) then ydnrng=abs(drange[0])	;remember to add DRANGE to calling sequence!
ymin=min(yy/dxx,max=ymax,/nan)
if n_elements(yrange) eq 2 then yuprange=yrange else $
	yuprange=ymax*[ydynrng,1.5]
xmid=0.5*(xx[1:*]+xx)

;	save current plot state
pmulti=!p.multi & pposition=!p.position & ymargin=!y.margin & pregion=!p.region

;	set up new plot state
!p.multi=[0,1,2] & !p.position=0

!p.region=[0,winsplit,1,1]
!y.margin=[0,3]
plot,xx,yy/dxx,/nodata,xstyle=5,yrange=yuprange,ytitle=defyt,$
	ylog=ylog,xthick=xthick,ythick=ythick,title=title,$
	charsize=charsize,charthick=charthick,$
	_extra=e
axis,xaxis=1,xstyle=1,xticks=1,xtickname=[' ',' ']
if total(abs(yy),/nan) gt 0 then begin
  for i=0L,ny-1L do oplot,$
  	[xx[i],xx[i+1L]],(yy[i]/dxx[i])*[1,1],color=lcolarr[1],_extra=e
  for i=0L,ny-1L do oplot,$
  	xmid[i]*[1,1],(yy[i]+ysig[i]*[-1,1])/dxx[i],color=lcolarr[2], _extra=e
  oplot,xmid,yy/dxx,psym=10,color=lcolarr[1],_extra=e
endif
if total(abs(mm),/nan) gt 0 then $
	oplot,xmid,mm/dxx,psym=10,color=lcolarr[0], _extra=e
for i=0,nalt-1 do begin
  j=strtrim(i,2)
  jnk=execute('yN=d'+j)
  jnk=execute('mN=m'+j)
  jnk=execute('ysigN=sig'+j)
  jnk=execute('colN=col'+j)
  if ctthr gt 0 then begin
    yNr=0*yy & mNr=0*mm & ysigNr=0*ysig
    if total(abs(yN),/nan) gt 0 then begin
      yNr=regroup(yN,ctthr,grid)
      ysigNr=sqrt(regroup(ysigN^2,ctthr,grid))
    endif
    if total(abs(mN),/nan) gt 0 then mNr=regroup(mN,ctthr,grid)
    yN=yNr & ysigN=ysigNr & mN=mNr
  endif
  ncN=n_elements(colN) & lcolNarr=nalt*25+[25,35,45]
  if ncN eq 1 then lcolNarr[0]=colN[0] else $
   if ncN eq 2 then lcolNarr[0:1]=colN else $
    if ncN eq 3 then lcolNarr=colN else $
     if ncN gt 3 then lcolNarr=colN[0:2]
  if total(abs(yN),/nan) gt 0 then begin
    for k=0L,ny-1L do oplot,$
    	[xx[k],xx[k+1L]],(yN[k]/dxx[k])*[1,1],color=lcolNarr[1], _extra=e
    for k=0L,ny-1L do oplot,$
    	xmid[k]*[1,1],(yN[k]+ysigN[k]*[-1,1])/dxx[k],color=lcolNarr[2], _extra=e
    oplot,xmid,yN/dxx,psym=10,color=lcolNarr[1], _extra=e
  endif
  if total(abs(mN),/nan) gt 0 then $
  	oplot,xmid,mN/dxx,psym=10,color=lcolNarr[0], _extra=e
endfor

if keyword_set(iresid) then zz=(yy-mm)/dxx
if keyword_set(ifrac) then begin
  zz=0.*yy & o0=where(mm gt 0,mo0)
  if mo0 gt 0 then zz[o0]=(yy[o0]-mm[o0])/float(mm[o0])
endif
if keyword_set(idelchi) then zz=(yy-mm)/ysig
if keyword_set(iratio) then begin
  zz=0.*yy+1 & o0=where(mm gt 0,mo0)
  if mo0 gt 0 then zz[o0]=yy[o0]/float(mm[o0])
endif
zmin=min(zz,max=zmax,/nan)
if n_elements(zrange) eq 2 then ydnrange=zrange else ydnrange=[zmin,zmax]

!y.margin=[5,0]
!p.region=[0,0,1,winsplit]
if total(abs(zz),/nan) gt 0 then begin
  plot,xx,zz,/nodata,xstyle=1,yrange=ydnrange,ystyle=1,xtitle=defxt,ytitle=defytb,$
	xthick=xthick,ythick=ythick,subtitle=subtitle,$
	charsize=charsize,charthick=charthick,$
	_extra=e
  oplot,xx,zz,psym=10,color=lcolarr[0], _extra=e
endif else begin
  plot,xx,0*xx,/nodata,xstyle=1,yrange=ydnrange,ystyle=1,xtitle=defxt,ytitle=defytb,$
  	xthick=xthick,ythick=ythick,subtitle=subtitle,$
	charsize=charsize,charthick=charthick,$
	_extra=e
endelse
for i=0,nalt-1 do begin
  j=strtrim(i,2)
  jnk=execute('yN=d'+j)
  jnk=execute('mN=m'+j)
  jnk=execute('ysigN=sig'+j)
  jnk=execute('colN=col'+j)
  if ctthr gt 0 then begin
    yNr=0*yy & mNr=0*mm & ysigNr=0*ysig
    if total(abs(yN),/nan) gt 0 then begin
      yNr=regroup(yN,ctthr,grid)
      ysigNr=sqrt(regroup(ysigN^2,ctthr,grid))
    endif
    if total(abs(mN),/nan) gt 0 then mNr=regroup(mN,ctthr,grid)
    yN=yNr & ysigN=ysigNr & mN=mNr
  endif
  ncN=n_elements(colN) & lcolNarr=nalt*25+[25,35,45]
  if ncN eq 1 then lcolNarr[0]=colN[0] else $
   if ncN eq 2 then lcolNarr[0:1]=colN else $
    if ncN eq 3 then lcolNarr=colN else $
     if ncN gt 3 then lcolNarr=colN[0:2]
  if keyword_set(iresid) then zN=(yN-mN)/dxx
  if keyword_set(ifrac) then begin
    zN=0.*yN & o0=where(mN gt 0,mo0)
    if mo0 gt 0 then zN[o0]=(yN[o0]-mN[o0])/float(mN[o0])
  endif
  if keyword_set(idelchi) then zN=(yN-mN)/ysigN
  if keyword_set(iratio) then begin
    zN=0.*yN+1 & o0=where(mN gt 0,mo0)
    if mo0 gt 0 then zN[o0]=yN[o0]/float(mN[o0])
  endif
  if total(abs(yN),/nan) gt 0 and $
     total(abs(mN),/nan) gt 0 and $
     total(abs(zN),/nan) gt 0 then $
       oplot,xx,zN,psym=10,color=lcolNarr[0], _extra=e
endfor

if vv gt 0 then stample, _extra=e

;	restore original plot state
!p.multi=pmulti
!p.position=pposition
!y.margin=ymargin
!p.region=pregion

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return
end
