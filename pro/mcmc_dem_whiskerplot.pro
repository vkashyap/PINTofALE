pro mcmc_dem_whiskerplot,savedem,epsroot=epsroot,thrfrac=thrfrac,thrnum=thrnum,dvdt=dvdt,$
	logT=logT,simdem=simdem,storidx=storidx,$
	cbox=cbox,cwhisk=cwhisk,cdots=cdots,cbest=cbest,cmode=cmode,cgrey=cgrey,cnum=cnum,$
	labul=labul,labur=labur,labll=labll,lablr=lablr,yposul=yposul,yposur=yposur,yposll=yposll,yposlr=yposlr,$
	xrange=xrange,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title,cthick=cthick,csize=csize,$
	verbose=verbose, _extra=e
;+
;procedure	mcmc_dem_whiskerplot
;	make a nice whisker plot of the DEM solutions from MCMC_DEM.
;	The central box covers 50% of the solutions at each temperature,
;	and the whiskers cover the full range.  Individual points are
;	also overplotted, along with the mode at each temperature, and
;	the full DEM(T) that corresponds to the mode at each T, as well
;	as the overall best-fit solution.
;
;syntax
;	mcmc_dem_whiskerplot,savedem,epsroot=epsroot,thrfrac=thrfrac,thrnum=thrnum,/dvdt,$
;	logT=logT,simdem=simdem,storidx=storidx,$
;	cbox=cbox,cwhisk=cwhisk,cdots=cdots,cbest=cbest,cmode=cmode,cgrey=cgrey,cnum=cnum,$
;	labul=labul,labur=labur,labll=labll,lablr=lablr,yul=yul,yur=yur,yll=yll,ylr=ylr,$
;	xrange=xrange,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title,cthick=cthick,csize=csize,$
;	verbose=verbose,$
;	other PLOT keywords
;
;parameters
;	savedem	[INPUT; required] this is the name of the file where
;		MCMC_DEM saves its state
;		* if set to 'NONE' or 'none', looks for keywords
;		  LOGT,SIMDEM,STORIDX for input
;
;keywords
;	epsroot	[INPUT] root name for encapsulated postscript file
;		* if set, puts the plot in EPSROOT.eps
;	thrfrac	[INPUT] the fraction of the maximum of the histogram
;		of selected temperatures in MCMC iteration to use
;		* this looks at how many times a given temperature index
;		  is present in the MCMC trace
;		* default is 0.05
;		* if 0, assumes default
;		* if >1, the reciprocal is used
;		* if <0, abs value fraction is taken of the _total_ number of
;		  of iterations, not of the max amongst different T
;		* NOTE: the greater of THRFRAC and THRNUM is used
;	thrnum	[INPUT] the minimum _number_ of points of a given T
;		that should be present
;		* default is 100
;		* if 0, assumes default
;		* if <0, abs value is used
;		* there is a hardcoded minimum of 5
;		* NOTE: the greater of THRFRAC and THRNUM is used
;	dvdt	[INPUT] if set, converts the units of DEM to [.../degK]
;		before plotting
;		* MCMC_DEM computes DEMs as dV/dlogT == dV/dlnT/ln(10) == (T/ln(10))*dV/dT
;		  so if this is set, divides SIMDEM by T/ln(10)
;
;	logT	[INPUT] logT grid, used iff SAVEDEM='NONE'
;	simdem	[INPUT] simulated DEMs, used iff SAVEDEM='NONE'
;	storidx	[INPUT] logT index selections during MCMC iterations, used iff SAVEDEM='NONE'
;
;	cbox	[INPUT; default=none] color index for box
;		* a color is applied only if set
;	cwhisk	[INPUT; default=3] color index for whiskers
;	cdots	[INPUT; default=4] color index for dots on whiskers
;		* recommend setting this to 0 for postscript output
;	cbest	[INPUT; default=8] color index for best-fit curve
;	cmode	[INPUT; default=1] color index for mode marker
;	cgrey	[INPUT; default=95] color index for example curves
;		* set CWHISK,CDOTS,CBEST,CMODE,CGREY explicitly to 0 to prevent plotting
;	cnum	[INPUT; default=9] color index for marking the number of samples
;
;	labul	[INPUT] LABel to put on Upper Left of plot
;	labur	[INPUT] LABel to put on Upper Right of plot
;	labll	[INPUT] LABel to put on Lower Left of plot
;	lablr	[INPUT] LABel to put on Lower Right of plot
;	yposul	[INPUT] y-location of LABUL
;	yposur	[INPUT] y-location of LABUR
;	yposll	[INPUT] y-location of LABLL
;	yposlr	[INPUT] y-location of LABLR
;		* there are no equivalent XPOS because they can always be adjusted with
;		  judicious use of spaces within the label strings
;
;	xrange	[INPUT] passed on to PLOT without comment
;	yrange	[INPUT] passed on to PLOT without comment
;	xtitle	[INPUT] passed on to PLOT without comment
;	ytitle	[INPUT] passed on to PLOT without comment
;	title	[INPUT] passed on to PLOT without comment
;	cthick	[INPUT] line thickness
;		* default is 2, unless EPSROOT is set, when it is 5
;	csize	[INPUT] character and symbol sizes
;		* default is 2, unless EPSROOT is set, when it is 1.2
;
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;
;requires subroutines
;	PEASECOLR
;	MODALPOINT
;
;history
;	Vinay Kashyap (2015sep)
;	bug fix with lower bound of whisker box not including mode; ensure clean quit
;	  when NSIM is small (VK; 2015nov)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(savedem) & szs=size(savedem,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SAVEDEM is not specified' else $
  if ns gt 1 then ok='SAVEDEM must be a scalar' else $
   if szs ne 7 then ok='SAVEDEM must be character string' else begin
     if strupcase(savedem[0]) ne 'NONE' then begin	;(SAVEDEM is given
       sfil=file_search(savedem[0],count=nsav)
       if nsav eq 0 then ok=savedem[0]+': file not found' else $
        if nsav gt 1 then ok=savedem[0]+': too many files found'
     endif else begin					;)(SAVEDEM is none
       nT=n_elements(logT) & nd=n_elements(simdem) & ni=n_elements(storidx)
       szd=size(simdem) & szi=size(storidx)
       if nT eq 0 or nd eq 0 or ni eq 0 then ok='LOGT,SIMDEM,STORIDX must be specified' else $
        if szd[0] ne 2 then ok='SIMDEM must be 2D' else $
	 if szd[1]*(szd[2]-1L) ne szi[1] then ok='SIMDEM and STORIDX are incompatible' else $
	  nsim=szd[2]-1L
     endelse						;SAVEDEM=NONE)
   endelse
if ok ne 'ok' then begin
  print,'Usage: mcmc_dem_whiskerplot,savedem,epsroot=epsroot,thrfrac=thrfrac,thrnum=thrnum,/dvdt,$'
  print,'       logT=logT,simdem=simdem,storidx=storidx,$'
  print,'       cbox=cbox,cwhisk=cwhisk,cdots=cdots,cbest=cbest,cmode=cmode,cgrey=cgrey,cnum=cnum,$'
  print,'       labul=labul,labur=labur,labll=labll,lablr=lablr,yposul=yposul,yposur=yposur,yposll=yposll,yposlr=yposlr,$'
  print,'       xrange=xrange,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title,cthick=cthick,csize=csize,$'
  print,'       verbose=verbose, other PLOT keywords'
  print,'  make whisker plots for output of MCMC_DEM'
  if np ne 0 then message,ok,/informational
  return
endif

;	read in MCMC_DEM output and initialize
peasecolr & tvlct,rr,gg,bb,/get & tvlct,rr,gg,bb & peasecolr
;	set default labels for plot
xttl='log!d10!n(T [K])'
yttl='DEM [cm!u-5!n logK!u-1!n]'
ttl=' '
if strupcase(savedem[0]) ne 'NONE' then begin
  restore,sfil[0]	;SIMDEM, LOGT, STORIDX, NT, NSIM
  spawn,'basename '+sfil[0],tmp & ttl=tmp[0]
endif
hidx=histogram(storidx,min=0,max=nT-1L,binsize=1,reverse_indices=ri)
maxhidx=max(hidx)
;	convert units if necessary
demarr=simdem
if keyword_set(dvdt) then begin
  for i=0L,nT-1L do demarr[i,*]=simdem[i,*]*alog(10.)/10.D^(logT[i])	;convert dV/dlogT to dV/dT
  yttl='DEM [cm!u-5!n K!u-1!n]'
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
nps=n_elements(epsroot)
th=2 & cs=2
if nps ne 0 then begin
  psfil=strtrim(epsroot[0],2)+'.eps'
  dname=!d.NAME
  set_plot,'ps' & device,file=psfil,/color,/encapsulated
endif else begin
  if !d.window lt 0 and !d.name eq 'X' then window,0,xsize=1000,ysize=800,title=ttl	;if window is already open, plot into that, otherwise make a nice big labeled canvas
endelse
if !d.name eq 'PS' then begin & th=5 & cs=1.2 & endif	;thickness/charsize for postscript
if keyword_set(cthick) then th=float(abs(cthick[0]))
if keyword_set(csize) then cs=float(abs(csize[0]))
poaintsym,'circle',/pfil,psiz=cs/4.
;
fthr=0.05 & nmax=maxhidx
if keyword_set(thrfrac) then begin
  if thrfrac[0] lt 0 then nmax=total(hidx)
  fthr=abs(thrfrac[0])
  if abs(thrfrac[0]) gt 1 then fthr=1./abs(thrfrac[0])
endif
nthr=200L
if keyword_set(thrnum) then nthr=abs(thrnum[0]) > 5L
thr=(fthr*nmax) > nthr
;
ibox=255 & if !d.name eq 'PS' then ibox=0 & if n_elements(cbox) gt 0 then ibox=byte(cbox[0])
iwhisk=3 & if !d.name eq 'PS' then iwhisk=34 & if n_elements(cwhisk) gt 0 then iwhisk=byte(cwhisk[0])
idots=4 & if n_elements(cdots) gt 0 then idots=byte(cdots[0])
ibest=8 & if n_elements(cbest) gt 0 then ibest=byte(cbest[0])
imode=1 & if n_elements(cmode) gt 0 then imode=byte(cmode[0])
inum=9 & if n_elements(cnum) gt 0 then inum=byte(cnum[0])
igrey=95 & if n_elements(cgrey) gt 0 then igrey=byte(cgrey[0])
;
dx=median(logT[1:*]-logT)
ot=where(hidx gt thr,mot)
if mot eq 0 then begin
  message,'THRESHOLD for selecting T bins to plot ('+strtrim(thr,2)+'; cf. NSIM='+strtrim(NSIM,2)+') seems to be set too high.',/informational
  message,'There is nothing to plot.  Quitting.',/informational
  return
endif
xmin=min(logT[ot],max=xmax) & xrng=[xmin-dx,xmax+dx]
if n_elements(xrange) eq 2 then xrng=xrange
ymin=min(demarr,max=ymax) & yrng=[ymin/2.,ymax*2.]	;DEMARR comes from SIMDEM comes from SAVEDEM
if n_elements(yrange) eq 2 then yrng=yrange
if keyword_set(xtitle) then xttl=strjoin(xtitle)
if keyword_set(ytitle) then yttl=strjoin(ytitle)
if keyword_set(title) then ttl=strjoin(title)
;
ullab='' & if keyword_set(labul) then ullab=strtrim(labul,2)
urlab='' & if keyword_set(labur) then urlab=strtrim(labur,2)
lllab='' & if keyword_set(labll) then lllab=strtrim(labll,2)
lrlab='' & if keyword_set(lablr) then lrlab=strtrim(lablr,2)

;	draw the axes
plot,[0],/nodata,xrange=xrng,yrange=yrng,/ylog,/xs,/ys,$
	xtitle=xttl,ytitle=yttl,title=ttl,$
	thick=th,xthick=th,ythick=th,charthick=th,charsize=cs,$
	_extra=e

;	step through the temperatures
ileft=nT & iright=-1L
for i=0L,nT-1L do begin
  ;this is not needed, but here just to show how it could be done: ok=ri[ri[i]:ri[i+1L]-1L]
  nidx=ri[i+1L]-ri[i]
  if nidx gt thr then begin
    sdem=demarr[i,0L:nsim-1L]
    os=sort(sdem) & ssdem=sdem[os]				;sorted DEMARR == SIMDEM
    lsdem=alog10(ssdem)						;log sorted DEMARR == SIMDEM
    bdem=demarr[i,nsim]	;best DEM
    lmdem=modalpoint(lsdem) & jnk=min(abs(lsdem-lmdem),jmode)	;mode of log DEMARR == SIMDEM and index of mode in sorted DEMARR == SIMDEM
    mindem=min(lsdem,max=maxdem)				;range of log DEMARR == SIMDEM
    w50=fltarr(nsim)+(maxdem-mindem)
    j0=(jmode-nsim/2+1L) > 0L
    ;j1=(jmode+nsim/2-1L) < (nsim/2-1L)	;old code -- incorrect
    j1=jmode < (nsim/2-1L)	;this stops looking when the lower bound goes above the mode
    for j=j0,j1 do w50[j]=lsdem[j+nsim/2-1L]-lsdem[j]
    wmin=min(w50,imin)						;smallest width that includes mode
    jnk=min(abs(alog10(ssdem[imin])-alog10(sdem)),jsim)		;iteration index that corresponds to this mode
    kmin=imin
    kmax=imin+nsim/2-1L
    ;if kmin lt 0 then begin & kmin=0L & kmax=nsim/2 & endif		;this is an unnecessary check
    ;if kmax gt nsim-1L then begin & kmax=nsim-1L & kmin=nsim/2 & endif	;this is an unnecessary check
    demlo=ssdem[kmin] & demhi=ssdem[kmax]		;50% range of DEMARR == SIMDEM around mode of log DEMARR == SIMDEM

    ileft = ileft < i
    iright = iright > i

    ;	draw the dots
    oy=where(ssdem le demlo or ssdem ge demhi,moy)
    if idots gt 0 and moy gt 0 then oplot,logT[i]+fltarr(moy),ssdem[oy],psym=8,color=idots
    ;	draw the whiskers
    if iwhisk gt 0 then begin
      oplot,logT[i]*[1,1],[demhi,10.D^(maxdem)],color=iwhisk
      oplot,logT[i]*[1,1],[demlo,10.D^(mindem)],color=iwhisk
      oplot,logT[i]+dx*0.25*[-1,1],10.D^(maxdem)*[1,1],color=iwhisk
      oplot,logT[i]+dx*0.25*[-1,1],10.D^(mindem)*[1,1],color=iwhisk
    endif
    ;	draw the solution that goes through the mode
    if igrey ne 0 then oplot,logT,demarr[*,jsim],color=igrey
    ;	mark the number of samples in this bin
    if inum ne 0 then xyouts,logT[i],10.D^(mindem)*0.7,strtrim(nidx,2),align=0.5,charsize=0.5,color=inum
    ;	(always) draw the 50% box
    oplot,logT[i]+dx*0.5*[-1,1,1,-1,-1],[demlo,demlo,demhi,demhi,demlo],color=ibox,thick=th
    ;	draw the mode
    if imode gt 0 then oplot,logT[i]+dx*0.5*[-1,1],10.D^(lmdem)*[1,1],color=imode,thick=th

    if lmdem lt alog10(demlo) or lmdem gt alog10(demhi) then stop,'BUG!'

  endif else begin
    if nidx gt 1 then begin
      sdem=demarr[i,0L:nsim-1L] & os=sort(sdem) & ssdem=sdem[os]
      lsdem=alog10(ssdem)
      ;	draw the dots for the other cases
      if idots gt 0 then oplot,logT[i]+fltarr(nsim),ssdem,psym=8,color=idots,symsize=0.2
      ;	draw the whiskers (Note: design decision -- no T-stops at the ends)
      if iwhisk gt 0 then oplot,logT[i]*[1,1],minmax(sdem),color=iwhisk
    endif
  endelse

endfor

;	plot the best-fit
if ibest gt 0 then oplot,logT,demarr[*,nsim],color=ibest,thick=th,line=2,psym=10

;	plot the labels
;	(first figure out default y-locations
yul=max(demarr[ileft,*])*2. > (yrng[1]/alog10(yrng[1]/yrng[0]))
yur=max(demarr[iright,*])*2. > (yrng[1]/alog10(yrng[1]/yrng[0]))
yul = yul > yur
yur = yur > yul
yll=min(demarr[ileft,*])*0.5 < (yrng[0]*alog10(yrng[1]/yrng[0]))
ylr=min(demarr[iright,*])*0.5 < (yrng[0]*alog10(yrng[1]/yrng[0]))
yll = yll < ylr
ylr = ylr < yll
;	Y[UL][LR])(override if specified during call
if keyword_set(yposul) then yul=yposul[0]
if keyword_set(yposur) then yur=yposur[0]
if keyword_set(yposll) then yll=yposll[0]
if keyword_set(yposlr) then ylr=yposlr[0]
;	YPOS[UL][LR])
if keyword_set(ullab) then xyouts,xrng[0]+dx,yul,ullab,align=-1,charsize=cs
if keyword_set(urlab) then xyouts,xrng[1]-dx,yur,urlab,align=1,charsize=cs
if keyword_set(lllab) then xyouts,xrng[0]+dx,yll,lllab,align=-1,charsize=cs
if keyword_set(lrlab) then xyouts,xrng[1]-dx,ylr,lrlab,align=1,charsize=cs

;	close out postscript output
if nps ne 0 then begin
  device,/close & set_plot,dname
  print,'' & print,'$gv '+psfil+' &' & print,''
  spawn,'ls -l '+psfil
endif

if vv gt 1000 then stop,'halting; type .CON to continue'

return
end
