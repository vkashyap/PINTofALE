pro ecfcontour,xin,yin,nbg,ecflev=ecflev,collev=collev,outroot=outroot,$
	xout=xout,yout=yout,fout=fout,xcen=xcen,ycen=ycen,$
	polycon=polycon,verbose=verbose,$
	landscape=landscape,encapsulated=encapsulated, _extra=e
;-
;procedure	ecfcontour
;	given a list of 2D locations, constructs the enclosed counts fraction
;	using HULLSORT(), makes contours at specified values and plots them.
;
;syntax
;	ecfcontour,xin,yin,nbg,ecflev=ecflev,collev=collev,outroot=outroot,$
;	xout=xout,yout=yout,fout=fout,xcen=xcen,ycen=ycen,polycon=polycon,$
;	verbose=verbose,/landscape,/encapsulated, twait=twait
;
;parameters
;	xin	[INPUT; required] 1D array of x-coordinate points
;	yin	[INPUT; required] 1D array of y-coordinate points
;		* sizes must match
;		* sent straight to HULLSORT()
;	nbg	[INPUT] expected number of background events
;		* it is assumed that the background is uniform
;		* the cdfs are corrected at each level based on the
;		  expected number of the points that could be
;		  attributed to the background, scaled by the
;		  area at each point.
;		* note that NBG does not have to be an integer, it is
;		  the _expected_ background counts
;
;keywords
;	ecflev	[INPUT] array of enclosed counts fractions at which to
;		draw contours
;		* default is [0.95, 0.90, 0.50, 0.39]
;		* values >1 or <0 are pegged to the bounds
;		* can be given in any order, but will always get resorted
;		  from outside to in
;	collev	[INPUT] colors to fill the contours within each level
;		* default is byte((findgen(N(ECFLEV))+1)*(256/N(ECFLEV)))
;		* if size of input COLLEV is smaller than that of ECFLEV,
;		  the remainder is filled out with the default
;	outroot	[INPUT] root name of postscript plots
;		* default is 'hull'
;		* two plots are generated --
;		  - OUTROOTecfcontour.ps, which makes colored polyfills
;		  - OUTROOTecfpoints.ps, which colors the points
;		* if OUTROOT='none' files are not written out
;		* if OUTROOT='display', 'x', 'plot', then
;		  if !d.name='X' then plots are drawn to screen else ignored
;		* file extensions change to .eps if ENCAPSULATED is set
;	xout	[OUTPUT] reordered XIN, straight out of HULLSORT()
;	yout	[OUTPUT] reordered XIN, straight out of HULLSORT()
;	fout	[OUTPUT] enclosed fraction corresponding to (XOUT,YOUT)
;	xcen	[OUTPUT] centroid x-coordinate computed at each step in HULLSORT()
;	ycen	[OUTPUT] centroid y-coordinate computed at each step in HULLSORT()
;	polycon	[OUTPUT] structure containing the polygon to define the polyfill
;		for each ECFLEV
;	verbose	[INPUT] controls chatter
;	landscape	[INPUT] set for landscape postscript output
;	encapsulated	[INPUT] set for encapsulated postscript output
;
;	_extra	[INPUT ONLY] pass defined variables to subroutines
;		HULLSORT: TWAIT
;
;subroutines
;	HULLSORT()
;
;history
;	Vinay Kashyap (2013sep)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(xin) & ny=n_elements(yin)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='XIN is undefined' else $
  if ny eq 0 then ok='YIN is undefined' else $
   if nx ne ny then ok='XIN and YIN are incompatible' else $
    if nx lt 4 then ok='too few points in the input'
if ok ne 'ok' then begin
  print,'Usage: ecfcontour,xin,yin,nbg,ecflev=ecflev,collev=collev,outroot=outroot,$'
  print,'       xout=xout,yout=yout,fout=fout,xcen=xcen,ycen=ycen,polycon=polycon,$'
  print,'       verbose=verbose,/landscape,/encapsulated, twait=twait'
  print,'  constructs enclosed count contours'
  if np ne 0 then message,ok,/informational
  return
endif
if not keyword_set(nbg) then nbkg=0L else nbkg=float(nbg[0])>0
if nbkg ge nx then begin
  message,'background is set too high. returning',/informational
  return
endif

;	keywords
;
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
iland=0 & if keyword_set(landscape) then iland=1
;
ieps=0 & if keyword_set(encapsulated) then ieps=1
;
oroot='hull' & if keyword_set(outroot) then begin
  if size(outroot,/type) eq 7 then oroot=strtrim(outroot[0],2)
endif
;
levecf=[0.95,0.9,0.5,0.39] & necf=n_elements(ecflev)
if necf gt 0 then levecf=ecflev
necf=n_elements(levecf)
;
levcol=byte((findgen(necf)+1)*(255./float(NECF))) & ncol=n_elements(collev)
if ncol eq necf then levcol=byte(collev) else begin
  if ncol gt 0 then levcol[0L:(ncol<necf)-1L]=collev[0L:(ncol<necf)-1L]
endelse
;
os=reverse(sort(levecf)) & levecf=levecf[os] & levcol=levcol[os]
;
file1=oroot+'ecfcontour' & file2=oroot+'ecfpoints'
if keyword_set(ieps) then begin
  file1=file1+'.eps' & file2=file2+'.eps'
endif else begin
  file1=file1+'.ps' & file2=file2+'.ps'
endelse

;	get the enclosed counts fractions
jnk=hullsort(xin,yin,xout=xout,yout=yout,xcen=xcen,ycen=ycen,areap=areap,verbose=vv<5, _extra=e)
fout=reverse((findgen(nx)+1.)/float(nx))
cout=bytarr(nx)

;	instead of fancy stuff promised below, do the background correction
;	by correcting fout for expected number of events at each level.
if nbkg gt 0 then gout=(fout*nx-(areap/max(areap))*nbkg)/float(nx) else gout=fout
gmax=max(gout,imx)
if gmax le 0 then begin
  message,'background is out of control; ignoring completely',/informational
  gout=fout
endif
o0=where(gout lt 0,mo0)
if mo0 gt 0 then begin
  message,'background is too high beyond '+strtrim(imx,2)+'; ignoring all points past this level',/informational
  gout[imx:*]=gmax
endif

;	for each of the ecf levels, find the convex hull and plot it up
if vv gt 10 then plot,xin,yin,psym=3,/xs,/ys
for i=0L,necf-1L do begin
  oo=where(fout le levecf[i],moo)
  ;	background correction
  ;	FOUT<LEVECF says what fraction of the points are below LEVECF
  ;	GOUT at this cut says what fraction of those belong to the source
  ;	so to get the enclosed fraction of the source alone, must look at the LEVECF*GOUT'th point of FOUT
  if nbkg gt 0 then begin
    if moo gt 0 then gmax=max(gout[oo]) else gmax=levecf[i]
    oo=where(fout le gmax,moo)
  endif
  if moo gt 3 then begin
    cout[oo]=i+1
    xx=xout[oo] & yy=yout[oo] & arp=max(areap[oo])
    qhull,xx,yy,tr & jj=reform(tr[0,*]) & xp=xx[jj] & yp=yy[jj]
    avx=mean(xp,/double) & avy=mean(yp,/double)
    tht=atan(yp-avy,xp-avx) & os=sort(tht)
    xp=xp[os] & yp=yp[os] & tht=tht[os]
    xp=[xp,xp[0]] & yp=[yp,yp[0]] & tht=[tht,tht[0]]
    if vv gt 10 then oplot,xp,yp,psym=-1,color=levcol[i]
    if vv gt 10 then oplot,xx,yy,psym=3,color=levcol[i]
    tmp=create_struct('CLEV',levecf[i],'X',xp,'Y',yp,'area',arp,'NPT',moo)
    if n_tags(contourstr) eq 0 then contourstr=create_struct('Level'+strtrim(i+1,2),tmp) else $
    	contourstr=create_struct(contourstr,'Level'+strtrim(i+1,2),tmp)
  endif
  ;stop,i,moo,levecf[i]
endfor

;{
;	this block set aside for background corrections
;	basically, repeat the above block, with adjustments for number of points.
;	e.g., with N points total, the C^th level used to be at N*C, but if M of
;	these N are due to background, and the area at the C^th level is A_C, then
;	the C^th level will be at N*C-M*A_C/A_1
;}

;	output
polycon=contourstr

;	now make plots, if necessary
if strlowcase(oroot) eq 'none' then return
if !d.name ne 'X' and (strlowcase(oroot) eq 'display' or strlowcase(oroot) eq 'x' or strlowcase(oroot) eq 'plot') then return
hardcopy=1
if !d.name eq 'X' and (strlowcase(oroot) eq 'display' or strlowcase(oroot) eq 'x' or strlowcase(oroot) eq 'plot') then hardcopy=0
ywindow=1
if strlowcase(oroot) eq 'plain' then begin
  hardcopy=0 & ywindow=0
endif

;	polyfill plot
th=2
if keyword_set(hardcopy) then begin
  oldname=!d.name
  th=5 & set_plot,'ps' & device,file=file1,encapsulated=ieps,landscape=iland,color=1
endif else begin
  if keyword_set(ywindow) then if !d.name eq 'X' then window,0
endelse

plot,/nodata,xin,yin,/xs,/ys,xtitle='X',ytitle='Y',$
	xthick=th,ythick=th,thick=th,charthick=th,charsize=1.4,$
	_extra=e
delx=(!x.crange[1]-!x.crange[0])/float(necf+2)
dely=(!y.crange[1]-!y.crange[0])/float(necf+2)
for i=0L,necf-1L do begin
  tmp=contourstr.(i)
  xp=tmp.X & yp=tmp.Y
  polyfill,xp,yp,color=levcol[i]
  xyouts,!x.crange[0]+(i+1)*delx,!y.crange[1]-dely/2.,string(levecf[i],'(f4.2)'),charthick=th,charsize=1.4,color=levcol[i]
  xyouts,!x.crange[0]+(i+1)*delx,!y.crange[0]+dely/2.,string(tmp.NPT,'(i5)'),charthick=th,charsize=1.4,color=levcol[i]
endfor

if keyword_set(hardcopy) then begin
  device,/close & device,file=file2,encapsulated=ieps,landscape=iland,color=1
endif else begin
  if keyword_set(ywindow) then if !d.name eq 'X' then window,1
endelse

plot,/nodata,xin,yin,/xs,/ys,xtitle='X',ytitle='Y',$
	xthick=th,ythick=th,thick=th,charthick=th,charsize=1.4,$
	_extra=e
delx=(!x.crange[1]-!x.crange[0])/float(necf+2)
dely=(!y.crange[1]-!y.crange[0])/float(necf+2)
for i=0L,necf-1L do begin
  tmp=contourstr.(i)
  oo=where(cout eq i+1,moo)
  if moo gt 0 then oplot,xout[oo],yout[oo],psym=1,thick=th,color=levcol[i]
  xyouts,!x.crange[0]+(i+1)*delx,!y.crange[1]-dely/2.,string(levecf[i],'(f4.2)'),charthick=th,charsize=1.4,color=levcol[i]
  xyouts,!x.crange[0]+(i+1)*delx,!y.crange[0]+dely/2.,string(tmp.NPT,'(i5)'),charthick=th,charsize=1.4,color=levcol[i]
endfor

if keyword_set(hardcopy) then begin
  device,/close
  set_plot,oldname
  spawn,'ls -l '+file1+' '+file2
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return
end
