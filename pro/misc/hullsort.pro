function hullsort,xin,yin,xout=xout,yout=yout,xcen=xcen,ycen=ycen,$
	areap=areap,bycen=bycen,twait=twait,verbose=verbose, _extra=e
;+
;function	hullsort
;	reorder a 2D list from outside in and return the indices of
;	the reorderd list.  At each stage, compute the convex hull
;	and remove the "farthest" point from the hull centroid.
;
;syntax
;	iout=hullsort(xin,yin,xout=xout,yout=yout,xcen=xcen,ycen=ycen,$
;	areap=areap,/bycen,twait=twait,verbose=verbose)
;
;parameters
;	xin	[INPUT; required] 1D array of x-coordinate points
;	yin	[INPUT; required] 1D array of y-coordinate points
;		* sizes must match
;		* only useful if there are more than 4 points
;
;keywords
;	xout	[OUTPUT] reorderd XIN
;	yout	[OUTPUT] reorderd YIN
;	xcen	[OUTPUT] centroid x-coordinate computed at each step
;	ycen	[OUTPUT] centroid y-coordinate computed at each step
;	areap	[OUTPUT] area covered by the convex hull at each stage
;	bycen	[INPUT] if set, eliminates points from convex hull based
;		on distance to center
;		* default is to eliminate those farthest from all other
;		  points on the convex hull
;	twait	[INPUT] waiting time for replotting on screen [seconds]
;		* only applies if VERBOSE>10
;		* default is 0.001
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run hullsort
;
;subroutines
;	QHULL (IDL built-in)
;	KILROY
;	AREAPOLY
;
;history
;	vinay kashyap (MMXII.VII)
;	added call to AREAPOLY() (VK; MMXII.X)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(xin) & ny=n_elements(yin)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='XIN is undefined' else $
  if ny eq 0 then ok='YIN is undefined' else $
   if nx ne ny then ok='XIN and YIN are incompatible' else $
    if nx lt 4 then ok='too few points in the input'
if ok ne 'ok' then begin
  print,'Usage: iout=hullsort(xin,yin,xout=xout,yout=yout,xcen=xcen,ycen=ycen,$'
  print,'       areap=areap,/bycen,twait=twait,verbose=verbose)'
  print,'  reorder a 2D list from the outside in and return indices of reordered list'
  if np ne 0 then message,ok,/informational
  if nx gt 0 then return,lindgen(nx) else return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
waitsec=0.001 & if keyword_set(twait) then waitsec=double(twait[0])>0

;	outputs
iout=lonarr(nx) & jout=lindgen(nx)
xx=xin[*] & yy=yin[*] & xout=xx & yout=yy
xcen=fltarr(nx)+mean(xx,/double,/nan) & ycen=fltarr(ny)+mean(yy,/double,/nan)
areap=fltarr(nx)

for i=0L,nx-3L-1L do begin
  qhull,xx,yy,tr & jj=reform(tr[0,*]) & mm=n_elements(xx)
  ;avx=mean(xx[jj],/double,/nan) & avy=mean(yy[jj],/double,/nan)
  avx=mean(xx,/double,/nan) & avy=mean(yy,/double,/nan)
  dd=sqrt((xx[jj]-avx)^2+(yy[jj]-avy)^2)
  xcen[i:*]=avx & ycen[i:*]=avy
  areap[i:*]=areapoly(xx[jj],yy[jj])
  if keyword_set(bycen) then begin
    jnk=max(dd,imx)
    kk=jj[imx] & iout[i]=jout[kk] & xout[i]=xx[kk] & yout[i]=yy[kk]
  endif else begin
    ndd=n_elements(dd) & dd2=dblarr(ndd,ndd)
    for j=0L,ndd-1L do dd2[j,*]=sqrt((xx[jj]-xx[jj[j]])^2+(yy[jj]-yy[jj[j]])^2)
    dd22=total(dd2,2)/double(ndd)
    ;dd2j=dblarr(ndd) & for j=0L,ndd-1L do begin & oj=where(dd2[j,*] gt 0) & dd2j[j]=min(dd2[j,oj]) & endfor
    ;dd2j=dblarr(ndd) & for j=0L,ndd-1L do begin & os=sort(dd2[j,*]) & dd2j[j]=mean(dd2[j,os[0:ndd/3]],/double) & endfor
    dd2j=dblarr(ndd) & for j=0L,ndd-1L do begin & os=sort(dd2[j,*]) & dd2j[j]=mean(dd2[j,os[0:1]],/double) & endfor
    ;for j=0L,ndd-1L do dd[j]=dd[j]+dd22[j]+dd2j[j]
    for j=0L,ndd-1L do dd[j]=dd2j[j]
    ;for j=0L,ndd-1L do dd[j]=dd[j]+mean(sqrt((xx[jj]-xx[jj[j]])^2+(yy[jj]-yy[jj[j]])^2),/double)
    ;for j=0L,ndd-1L do dd[j]=dd[j]+mean(sqrt((xx-xx[jj[j]])^2+(yy-yy[jj[j]])^2),/double)
    jnk=max(dd,imx)
    kk=jj[imx] & iout[i]=jout[kk] & xout[i]=xx[kk] & yout[i]=yy[kk]
  endelse
  if vv gt 10 then begin
    plot,xx,yy,psym=1,/xs,/ys,title=i,/nodata
    oplot,xx,yy,psym=1,col=100
    if vv gt 50 then oplot,xx[jj],yy[jj],psym=4,col=200
    if vv gt 100 then oplot,[xout[i]],[yout[i]],psym=4,col=150,symsize=2
    wait,waitsec
  endif
  ;ii=lindgen(mm) & oo=where(ii ne jj[imx])
  jout[kk]=-1 & o0=where(jout ge 0) & jout=jout[o0]
  xx=xx[o0] & yy=yy[o0]
  if vv gt 0 and vv lt 100 then $
   if i*vv eq 1000*long((i*vv)/1000) then kilroy
endfor
xout[i:*]=xx & yout[i:*]=yy & iout[i:*]=jout[o0]

if vv gt 1000 then stop,xout[i],yout[i]
  
return,iout
end

;	example run
if not keyword_set(npt) then npt=500L
if not keyword_set(verbose) then verbose=101
if not keyword_set(twait) then twait=0.01

;xin=randomn(seed,npt)*3+2 & yin=randomn(seed,npt)*6+4
xin=randomn(seed,npt/2)*3+2 & yin=randomn(seed,npt/2)*6+2
xin=[xin,randomn(seed,npt-npt/2)*1+4]
yin=[yin,randomn(seed,npt-npt/2)*2+4]

iout=hullsort(xin,yin,xout=xout,yout=yout,xcen=xcen,ycen=ycen,verbose=verbose,twait=twait)
icol=byte((lindgen(npt)/float(npt))*192.)+63

if !d.name eq 'X' then window,0 & plot,xin,yin,xtitle='X',ytitle='Y',/xs,/ys,psym=3,title='all points and centroid evolution'
for i=1L,npt-1L do oplot,[xcen[i-1L],xcen[i]],[ycen[i-1L],ycen[i]],col=icol[i]

if !d.name eq 'X' then window,2 & plot,xin,yin,/nodata,xtitle='X',ytitle='Y',/xs,/ys,psym=3,title='points color coded by EE'
for i=0L,npt-1L do plots,xout[i],yout[i],col=icol[i],psym=1

;	calling sequence
print,''
jnk=hullsort()
print,''

end
