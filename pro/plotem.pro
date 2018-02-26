pro plotem,logT,EM,idstr=idstr,markmin=markmin,multid=multid,$
	lcol=lcol,lsiz=lsiz,tlog=tlog,wdth=wdth,$
	xtitle=xtitle,ytitle=ytitle,title=title,xrange=xrange,yrange=yrange,$
	xstyle=xstyle,ystyle=ystyle,xlog=xlog,ylog=ylog, _extra=e
;+
;procedure	plotem
;	plots emission measures derived from the output of POTEM_TOOL
;	(or other means)
;
;syntax
;	plotem,logT,EM,idstr=idstr,markmin=markmin,multid=multid,lcol=lcol,$
;	lsiz=lsiz,tlog=tlog,wdth=wdth,xtitle=xtitle,ytitle=ytitle,title=title,$
;	xrange=xrange,yrange=yrange,xstyle=xstyle,ystyle=ystyle,/xlog,/ylog,$
;	/nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$
;	stathk=stathk
;
;parameters
;	logT	[INPUT] temperatures at which to plot the emission measure
;		estimates
;	EM	[INPUT] array of emission measures at the various temperatures
;		for each line
;		* if EM is not given, it is assumed that EM <- "LOGT" and
;		  LOGT = default
;		* EM is an array of the form EM(LOGT,WVL); if the size of
;		  the first dimension does not match input LOGT, the plots
;		  are interpolated
;
;keywords
;	idstr	[INPUT] ID structure, essentially to figure out labels for
;		the curves.
;		* if size does not match that of EM, cest la vie.
;	markmin	[INPUT] if set, plots only a single point for each curve,
;		corresponding to the minimum in EM
;		* the actual value of MARKMIN sets the PSYM for the points
;	multid	[INPUT] an integer array saying how many IDs there are
;		to each component of the ID structure (comes in handy
;		if IDSTR has all the emissivities collapsed into one)
;	lcol	[INPUT; default=200] color for the labels
;	lsiz	[INPUT; default=0.4] size of the labels
;	tlog	[I/O] if LOGT does not match size of EM, >and< EM does
;		not match the size of the temperature grid in emissivities
;		in ID being considered, then look for TLOG to supply the
;		actual temperature grid. (this latter feature is not yet
;		implemented)
;		* default is to go from log(T)=4..8 in as equidistant steps
;		* beware -- this will get overwritten if incorrect!!
;	wdth	[INPUT] the widths of each of the EM curves, from POTTASCH(),
;		via POTEM_TOOL()
;		* the size of this array must match the 2nd dimension of
;		  EM, or else it will be ignored
;
;	xtitle	[INPUT] passed w/o comment to PLOT
;	ytitle	[INPUT] passed w/o comment to PLOT
;	title	[INPUT] passed w/o comment to PLOT
;	xrange	[INPUT] overrides locally determined XRANGE
;	yrange	[INPUT] overrides locally determined YRANGE
;	xstyle	[INPUT] passed w/o comment to PLOT
;	ystyle	[INPUT] passed w/o comment to PLOT
;	xlog	[INPUT] passed w/o comment to PLOT, but default is 0
;	ylog	[INPUT] passed w/o comment to PLOT, but default is 1
;	_extra	[INPUT] pass defined keywords to subroutines
;		STAMPLE: NUTHIN,NODAY,NOTIME,NOUSER,NOPACK,STACOL,STASIZ,STATHK
;
;history
;	vinay kashyap (Nov98)
;	added keyword MULTID (VK; Dec98)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	added call to STAMPLE (VK; JanMMI)
;	improved color-scale setting for 24-bit consoles (VK; FebMMI)
;	added keyword WDTH (VK; SepMMVII)
;-

;	usage
np=n_params()
if np lt 1 then begin
  print,'Usage: plotem,logT,EM,idstr=idstr,markmin=markmin,multid=multid,$'
  print,'       lcol=lcol,lsiz=lsiz,tlog=tlog,wdth=wdth,xtitle=xtitle,$'
  print,'       ytitle=ytitle,title=title,xrange=xrange,yrange=yrange,$'
  print,'       xstyle=xstyle,ystyle=ystyle,/xlog,/ylog,/nuthin,/noday,$'
  print,'       /notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$'
  print,'       stathk=stathk'
  print,"  plots predicted emission measures for ID'd lines"
  return
endif

;	check input
nT=n_elements(logT) & szEM=size(EM) & nszEM=n_elements(szEM)
ok='ok'
if nT eq 0 then begin			;(logT not given
  nT=81L & logT=findgen(nT)*0.05+4.
  if szEM(nszEM-2) eq 0 then ok='Emission Measure not given' else $
   if szEM(nszEM-1) lt 2 then ok='Emission Measure must be an array'
endif					;NT=0->81)
if szEM(nszEM-2) eq 0 then begin	;(EM not given?
  szT=size(logT)
  if szT(0) eq 2 then begin
    EM=logT & szEM=szT & nszEM=n_elements(szEM)
    nT=81L & logT=findgen(nT)*0.05+4.
  endif else ok='confused.. where is the Emission Measure?'
endif					;EM<-logT)
if szEM(0) eq 1 then begin		;(EM is 1D?
  mT=szEM(1) & szEM=reform(EM,mT,1)
endif					;EM is 1D)
mT=szEM(1) & mW=szEM(2)
if mT lt 2 then ok='Emission Measure must be an array'
if ok ne 'ok' then begin
  message,ok,/info & return
endif

;	check keywords
dncolors=256. > !D.N_COLORS	;24-bit color screen temporary fix
lc=200 & if keyword_set(lcol) then lc=fix(lcol)
lc=fix(lc*(!d.n_colors/dncolors))
ls=0.5 & if keyword_set(lsiz) then ls=float(lsiz)
defT=findgen(mT)*(4./(mT-1.))+4. & moT=n_elements(tlog)
if moT ne mT then tlog=defT
xt='log!d10!n(T)' & if keyword_set(xtitle) then xt=strtrim(xtitle,2)
yt='EM' & if keyword_set(ytitle) then yt=strtrim(ytitle,2)
tt='Emission Measures' & if keyword_set(title) then tt=strtrim(title,2)
xs=1 & if keyword_set(xstyle) then xs=strtrim(xstyle,2)
ys=1 & if keyword_set(ystyle) then ys=strtrim(ystyle,2)
;
Tmin=min(logT,max=Tmax)
if Tmax eq Tmin then begin
  Tmin=Tmin-0.5 & Tmax=Tmax+0.5
endif
xr=[Tmin,Tmax]
if n_elements(xrange) eq 2 then xr=xrange	;override default XRange
;
EMin=min(EM,max=EMax)
if EMin eq -1 then begin	;(EMin=-1 ==> standard format for ignore
  oo=where(EM gt 0,moo)
  if moo gt 0 then EMin=min(EM(oo),max=EMax)
endif				;EMin=-1)
if EMax le EMin then EMax=EMin+!pi
if EMax gt 100. then begin		;(EM in normal numbers
  Emin=EMin/5. & EMax=EMax*5.			;stretch a little
endif else begin			;)(log(EM)
  EMin=EMin-0.5 & EMax=EMax+0.5		;stretch a little
endelse				;logEM)
yr=[EMin,EMax]
if n_elements(yrange) eq 2 then yr=yrange	;override default YRange
;
logx=0 & if keyword_set(xlog) then logx=1
logy=1 & if n_elements(ylog) ne 0 and not keyword_set(ylog) then logy=0
;
twdth=fltarr(mW)+median(logT(1:*)-logT)
if n_elements(wdth) eq mW then twdth=wdth

;	are there labels to be had?
labl=strarr(mW) & itlog=intarr(mW) & nid=n_tags(idstr)
if nid gt 0 then begin
  ;	initialize stuff
  ok='ok'
  atom=1 & rom=1 & inicon,atom=atom,roman=rom & atom=['?',atom] & rom=['',rom]
  ;	make sure the structure is in right format
  idnam=tag_names(idstr)
  ok='ok'
  if idnam(0) ne 'WVL' then ok='not ID structure' else $
   if idnam(1) eq 'WVL_COMMENT' then ok='no data in ID structure' else $
    if nid ne n_elements(idstr.(0))+1 then ok='no IDs?'
  if ok ne 'ok' then goto,noID		;{skip the label maker
  ;	if the number of WVLs in IDSTR match MW, excellent..
  nw=n_elements(idstr.WVL) & multID=0
  if n_elements(multid) eq nw then idmult=multid else idmult=intarr(nw)+1
  ;	..otherwise, see how many WVLs there are anyway..
  if nw ne mW then begin
    kID=0L
    for i=0,nW-1 do kID=kID+n_elements(idstr.(i+1).WVL)
    if kID eq mW then multID=1 else begin
      ;	..and as last resort, do the best ye can..
      multID=2
    endelse
  endif
  ;
  j=0L
  for i=0,nW-1 do begin				;{step through IDSTR
    tmp=idstr.(i+1) & zz=tmp.Z & jon=tmp.ion & ww=tmp.WVL
    n=n_elements(ww)
    ;
    em_added=0
    if idmult(i) gt 1 then em_added=1
    szz=size(zz) & nszz=n_elements(szz) & if szz(nszz-2) eq 4 then em_added=1
    ;
    for k=0,n-1 do begin			;{step through IDs
      if multID eq 0 then begin			;(IDs collapsed into 1
	labl(j)=labl(j)+'!C'+'!3'+atom(zz(k))+' '+rom(jon(k))+' !4k!3'+$
		strtrim(string(ww(k),'(g10.5)'),2)
	if k eq n-1 then j=j+1
      endif else begin				;)(each ID separate
	if j lt mW then begin
	  labl(j)='!3'+atom(zz(k))+' '+rom(jon(k))+' !4k!3'+$
		strtrim(string(ww(k),'(g10.5)'),2)
	  if keyword_set(em_added) then labl(j)=labl(j)+' ETC.'
	  j=j+1
	endif
      endelse					;MULTID)
    endfor					;K=0,N-1}
  endfor					;I=0,NW-1}
  noID:					;get here from GOTO}
endif

;	now go forth and set up the plot
plot,logT,0*logT,/nodata,xtitle=xt,ytitle=yt,title=tt,xrange=xr,yrange=yr,$
	xstyle=xs,ystyle=ys,xlog=logx,ylog=logy
stample, _extra=e
;	and plot!
for i=0,mW-1 do begin
  EMplot=reform(EM(*,i))
  if not keyword_set(markmin) then oplot,tlog,EMplot,col=(lc([i]))(0)
  oo=where(EMplot gt EMin,moo)
  if moo gt 0 then begin
    yl=min(EMplot(oo),ixl) & xl=tlog(oo(ixl))
    xyouts,xl,yl,labl(i),color=(lc([i]))(0),charsize=(ls([i]))(0),align=0.5
    if n_elements(wdth) eq mW then $
     oplot,tlog(i)+twdth(i)*0.5*[-1,1],yl*[1,1],color=(lc([i]))(0)
    if keyword_set(markmin) then begin
      ;	find the locations of all the minima (needed because emissivities
      ;	may be multi-bumped.  then mark all these points with PSYM=MARKMIN.
      xx=tlog & yy=-EMplot & nx=n_elements(xx) & ilm=lonarr(nx)-1L
      for j=0L,nx-1L do $
        if yy(j) ne -1 and (where(yy(j) le yy([j-1,j+1])))(0) eq -1 then $
	  ilm(j)=j
      olm=where(ilm ge 0,molm)
      if molm gt 0 then for j=0,molm-1 do $
	oplot,[tlog(olm(j))],[EMplot(olm(j))],psym=markmin;,col=(lc([i]))(0)
      if molm gt 0 then oplot,tlog(olm),EMplot(olm),line=1,col=(lc([i]))(0)
    endif
  endif
endfor

return
end
