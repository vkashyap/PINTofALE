function untailgate,time,chipx,chipy,iflag,delt=delt,delr=delr,$
	hemis=hemis,tdiff=tdiff,rdiff=rdiff,sdiff=sdiff,imult=imult,$
	verbose=verbose, _extra=e
;+
;function	untailgate
;	find and return a set of indices of event arrival times that are unaffected by tailgating
;
;syntax
;	idx=untailgate(time,chipx,chipy,iflag,delt=delt,delr=delr,$
;	/hemis,tdiff=tdiff,rdiff=rdiff,sdiff=sdiff,imult=imult,$
;	verbose=verbose)
;
;parameters
;	time	[INPUT; required] event arrival times
;	chipx	[INPUT; required] chip X positions of events
;	chipy	[INPUT; required] chip Y positions of events
;	iflag	[OUTPUT] indices of events that are flagged as having been tailgated
;
;keywords
;	delt	[INPUT; default=0.05 sec] height of pillbox
;	delr	[INPUT; default=20 pix] radius of pillbox
;	hemis	[INPUT] if set, flags events as tailgating if they are within a
;		3-D hemispherical volume with the times rescaled as (TDIFF/DELT)*DELR
;	tdiff	[OUTPUT] time offset of events from preceding event [sec]
;	rdiff	[OUTPUT] pixel offset of events from preceding event [pix]
;	sdiff	[OUTPUT] space-time offset of events from preceding event [pix]
;	imult	[OUTPUT] count of how many times a given event got flagged
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2018-nov)
;-

;	usage
ok='ok' & np=n_params() & nt=n_elements(time) & nx=n_elements(chipx) & ny=n_elements(chipy)
if np lt 3 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='TIME is not given' else $
  if nx eq 0 then ok='chipX is not given' else $
   if ny eq 0 then ok='chipY is not given' else $
    if nt ne nx then ok='Time and chipX should be of the same size' else $
     if nt ne ny then ok='Time and chipY should be of the same size' else $
      if nx ne ny then ok='chipX and chipY should be of the same size' else $
       if nt eq 1 then ok='No point in running this if input is not an array'
if ok ne 'ok' then begin
  print,'Usage: idx=untailgate(time,chipx,chipy,iflag,delt=delt,delr=delr,$'
  print,'           /hemis,tdiff=tdiff,rdiff=rdiff,verbose=verbose)'
  print,'  find and return a set of indices of event arrival times that are'
  print,'  unaffected by tailgating'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
tdel=0.05 & if keyword_set(delt) then tdel=delt[0]
rdel=20.0 & if keyword_set(delr) then rdel=delr[0]

;	keep track of each event -- 1 for tailgated, 0 for not tailgated
tarr=(time-time[0])
iarr=lindgen(nt) & karr=interpol(iarr,tarr,tarr+tdel) < (nt-1L)
ok=where(long(karr) gt iarr,mok)
ii=lonarr(nt) & rdiff=0.*chipx+10.*rdel & tdiff=0.*tarr+10.*tdel & sdiff=rdiff & imult=intarr(nt)
if keyword_set(hemis) then timexy=((time-time[0])/tdel)*rdel else timexy=0.*time

for i=0L,mok-1L do begin
  if vv gt 0 and i eq long(i/5000L)*5000L then kilroy
  j=ok[i] & k=long(karr[j])
  ztime=time[j+1L:k]-time[j]
  zt=timexy[j+1L:k]-timexy[j]
  zx=float(chipx[j+1L:k]-chipx[j]) & zy=float(chipy[j+1L:k]-chipy[j])
  zdist=sqrt(zx^2+zy^2)
  zd=sqrt(zt^2+zx^2+zy^2) & oy=where(zd le rdel,moy)
  if moy gt 0 then begin
    imult[oy+j+1L]=imult[oy+j+1L]+1
    ii[oy+j+1L]=1
    sdiff[oy+j+1L]=zd[oy] < (sdiff[oy+j+1L])
    tdiff[oy+j+1L]=ztime[oy] < (tdiff[oy+j+1L])
    rdiff[oy+j+1L]=zdist[oy] < (rdiff[oy+j+1L])
  endif
  ;tdiff[j+1L:k]=ztime < (tdiff[j+1L:k])
  ;rdiff[j+1L:k]=zdist < (rdiff[j+1L:k])
endfor

idx=where(ii eq 0,complement=iflag)

if vv gt 1000 then stop,'halting; type .CON to continue'

return,idx
end
;	example run

if not keyword_set(delt) then delt=0.05	;[sec]
if not keyword_set(delr) then delr=20.0	;[pix]
if not keyword_set(verbose) then verbose=1
if n_elements(hemis) eq 0 then hemis=0	;hemispheric pix-time (1) or cylindrical pillbox (0)

if n_tags(t) eq 0 then t=mrdfits('/snafu/largedata/HR1099/21945/secondary/hrcf21945_000N001_evt1.fits.gz',1,h) 
time=t.TIME & cx=t.CHIPX & cy=t.CHIPY & tzero=min(time) & time=time-tzero
xx=t.X & yy=t.Y & x0=16400. & y0=16200. & dd=sqrt((xx-x0)^2+(yy-y0)^2)
oo=where(dd lt 500) & x0=mean(xx[oo],/double) & y0=mean(yy[oo],/double) & dd=sqrt((xx-x0)^2+(yy-y0)^2)
oo=where(dd lt 100) & x0=mean(xx[oo],/double) & y0=mean(yy[oo],/double) & dd=sqrt((xx-x0)^2+(yy-y0)^2)
oo=where(dd lt 50) & x0=mean(xx[oo],/double) & y0=mean(yy[oo],/double) & dd=sqrt((xx-x0)^2+(yy-y0)^2)
ixtail=untailgate(time,cx,cy,iytail,delt=delt,delr=delr,hemis=hemis,tdiff=tdiff,rdiff=rdiff,sdiff=sdiff,imult=imult,verbose=verbose)

r85=fltarr(10) & per85=r85 & num85=lonarr(10)	;all, good, bad, imult>1, imult>2, dt=0-0.01,0.01-0.02,0.02-0.03,0.03-0.04,0.04-0.05
oo=where(dd lt 20,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[0]=zd[long(0.85*moo)] & per85[0]=sqrt(0.85*0.15/moo) & num85[0]=moo
oo=where(dd[ixtail] lt 20,moo) & zd=(dd[ixtail])[oo] & os=sort(zd) & zd=zd[os] & r85[1]=zd[long(0.85*moo)] & per85[1]=sqrt(0.85*0.15/moo) & num85[1]=moo
oo=where(dd[iytail] lt 20,moo) & zd=(dd[iytail])[oo] & os=sort(zd) & zd=zd[os] & r85[2]=zd[long(0.85*moo)] & per85[2]=sqrt(0.85*0.15/moo) & num85[2]=moo
oo=where(imult gt 1 and dd lt 20,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[3]=zd[long(0.85*moo)] & per85[3]=sqrt(0.85*0.15/moo) & num85[3]=moo
oo=where(imult gt 2 and dd lt 20,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[4]=zd[long(0.85*moo)] & per85[4]=sqrt(0.85*0.15/moo) & num85[4]=moo
dmin=0.00 & dmax=0.01 & oo=where(imult gt 0 and tdiff gt dmin and tdiff le dmax and dd lt 30,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[5]=zd[long(0.85*moo)] & per85[5]=sqrt(0.85*0.15/moo) & num85[5]=moo
dmin=0.01 & dmax=0.02 & oo=where(imult gt 0 and tdiff gt dmin and tdiff le dmax and dd lt 30,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[6]=zd[long(0.85*moo)] & per85[6]=sqrt(0.85*0.15/moo) & num85[6]=moo
dmin=0.02 & dmax=0.03 & oo=where(imult gt 0 and tdiff gt dmin and tdiff le dmax and dd lt 30,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[7]=zd[long(0.85*moo)] & per85[7]=sqrt(0.85*0.15/moo) & num85[7]=moo
dmin=0.03 & dmax=0.04 & oo=where(imult gt 0 and tdiff gt dmin and tdiff le dmax and dd lt 30,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[8]=zd[long(0.85*moo)] & per85[8]=sqrt(0.85*0.15/moo) & num85[8]=moo
dmin=0.04 & dmax=0.05 & oo=where(imult gt 0 and tdiff gt dmin and tdiff le dmax and dd lt 30,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & r85[9]=zd[long(0.85*moo)] & per85[9]=sqrt(0.85*0.15/moo) & num85[9]=moo
print,r85,num85,per85*100,form='(/,5(f6.2,1x),5x,5(f6.2,1x),/,5(i6,1x),5x,5(i6,1x),/,5(f6.3,1x),5x,5(f6.3,1x))'

;oo=where(imult gt 0 and tdiff gt 0.00220310688018799 and tdiff le 0.00439059734344482 and dd lt 100./0.13175,moo) & zd=dd[oo] & os=sort(zd) & zd=zd[os] & print,zd[0.85*moo]/0.13175,moo

end
