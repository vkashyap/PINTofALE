pro hrcs_process_evt,chip_id,au1,au2,au3,av1,av2,av3,amp_sf,crsu,crsv,$
  pha,statbit, ampscl,fineu,finev,fracu,fracv,u,v,au3cor=au3cor,av3cor=av3cor,$
  rawu=rawu,rawv=rawv,xdet=xdet,ydet=ydet,Hfilt=Hfilt,HXfilt=HXfilt,$
  caldb=caldb, degaps=degaps,degapf=degapf,tringf=tringf,$
  dampcor=dampcor,ampgain_s=ampgain_s,ampratio_s=ampratio_s,$
  pha_rng_min=pha_rng_min,pha_rng_mid=pha_rng_mid,pha_rng_max=pha_rng_max,$
  hevt=hevt,rangelev=rangelev,$
  pha_1to2=pha_1to2,pha_2to3=pha_2to3,w_1to2=w_1to2,w_2to3=w_2to3,$
  va=va,vb=vb,vc=vc,vd=vd,vbeta=vbeta,vgamma=vgamma,vthr=vthr,voff=voff,$
  ua=ua,ub=ub,uc=uc,ud=ud,ubeta=ubeta,ugamma=ugamma,uthr=uthr,uoff=uoff,$
  vBh=vBh,vHh=vHh,vAh=vAh,vdHh=vdHh,$
  uBh=uBh,uHh=uHh,uAh=uAh,udHh=udHh,$
  verbose=verbose,$
  _extra=e
;+
;procedure	hrcs_process_evt
;	process HRC-S amplifier data to produce positions
;
;syntax
;	hrcs_process_evt,chip_id,au1,au2,au3,av1,av2,av3,amp_sf,crsu,crsv,$
;	pha,statbit, ampscl,fineu,finev,fracu,fracv,u,v,au3cor=au3cor,av3cor=av3cor,$
;	rawu=rawu,rawv=rawv,xdet=xdet,ydet=ydet,Hfilt=Hfilt,HXfilt=HXfilt,$
;	caldb=caldb,degaps=degaps,degapf=degapf,tringf=tringf,$
;	dampcor=dampcor,ampgain_s=ampgain_s,ampratio_s=ampratio_s,$
;	pha_rng_min=pha_rng_min,pha_rng_mid=pha_rng_mid,pha_rng_max=pha_rng_max,$
;	hevt=hevt,rangelev=rangelev,$
;	pha_1to2=pha_1to2,pha_2to3=pha_2to3,w_1to2=w_1to2,w_2to3=w_2to3,$
;	va=va,vb=vb,vc=vc,vd=vd,vbeta=vbeta,vgamma=vgamma,vthr=vthr,$
;	ua=ua,ub=ub,uc=uc,ud=ud,ubeta=ubeta,ugamma=ugamma,uthr=uthr,$
;	vBh=vBh,vHh=vHh,vAh=vAh,vdHh=vdHh,$
;	uBh=uBh,uHh=uHh,uAh=uAh,udHh=udHh,$
;	verbose=verbose
;
;parameters
;	chip_id	[INPUT; required] on which plate did the count register?
;	au1	[INPUT; required] U tap signal from amplifier 1
;	au2	[INPUT; required] U tap signal from amplifier 2
;	au3	[INPUT; required] U tap signal from amplifier 3
;	av1	[INPUT; required] V tap signal from amplifier 1
;	av2	[INPUT; required] V tap signal from amplifier 2
;	av3	[INPUT; required] V tap signal from amplifier 3
;	amp_sf	[INPUT; required] amplifier scale
;	crsu	[INPUT; required] U tap (coarse position)
;	crsv	[INPUT; required] V tap (coarse position)
;	pha	[INPUT; required] event pulse height amplitudes
;	statbit	[INPUT; required] status bits
;	ampscl	[OUTPUT] rescaled AMP_SF, according to PHA range
;	fineu	[OUTPUT] U fine position based on AU1,AU2,AU3
;	finev	[OUTPUT] V fine position based on AU1,AU2,AU3
;	fracu	[OUTPUT] fraction of center-tap events in U
;	fracv	[OUTPUT] fraction of center-tap events in V
;	u	[OUTPUT] degapped U position
;	v	[OUTPUT] degapped V position
;
;keywords
;	au3cor	[OUTPUT] tap ringing corrected AU3
;	av3cor	[OUTPUT] tap ringing corrected AV3
;	rawu	[OUTPUT] V converted to raw detector coordinates
;	rawv	[OUTPUT] U converted to raw detector coordinates
;	xdet	[OUTPUT] RAWV converted to detector X coordinates *** NOT IMPLEMENTED ***
;	ydet	[OUTPUT] RAWU converted to detector Y coordinates *** NOT IMPLEMENTED ***
;	Hfilt	[OUTPUT] indices of all photons that pass the H-test
;	HXfilt	[OUTPUT] indices of all photons that fail the H-test
;	caldb	[INPUT] the location of the Chandra CALDB
;		* if not set, assumed to be /soft/ciao/CALDB/
;	degapf	[INPUT] CALDB file containing the degap coefficients
;		* if not set (and DEGAPS is not set), read from file
;		$CALDB/data/chandra/hrc/bcf/degap/hrcsD1999-07-22gapN0002.fits
;		* if string, read from specified file
;		* ignored if DEGAPS is well-defined
;	degaps	[I/O] the degap coefficients in a structure of the
;		format defined by the file DEGAPF
;		* if not set, read in from DEGAPF
;		* if set and passes rudimentary format checks, ignores
;		  DEGAPF and uses these values
;		* if set and does not pass rudimentary format checks,
;		  will be overwritten by the contents of DEGAPF
;		* if set to 0, will NOT remove gaps
;	tringf	[INPUT] CALDB file containing the tap ringing coefficients
;		* must be pathname relative to CALDB, e.g.,
;		  '/data/chandra/hrc/bcf/tapring/hrcsD1999-07-22tapringN0001.fits'
;		* if not defined, then looks at keywords
;		    VA,VB,VC,VD,VBETA,VGAMMA,VTHR,VOFF,
;		    UA,UB,UC,UD,UBETA,UGAMMA,UTHR,UOFF
;		  for the values, which btw default to the values in
;		  $CALDB/data/chandra/hrc/bcf/tapring/hrcsD1999-07-22tapringN0001.fits
;	dampcor	[INPUT] if set, dynamically determines the PHA range
;		over which the SUMAMPS scale changes
;	hevt	[INPUT] the header from the EVT1 file, to be used to
;		determine the value of RANGELEV
;		* see http://cxc.harvard.edu/ciao/threads/hrc_ampsf/
;		* if RANGELEV or DATE-OBS is present, overrides keyword input
;	ampgain_s,ampratio_s,pha_rng_min,pha_rng_mid,pha_rng_max,rangelev,pha_1to2,pha_2to3,w_1to2,w_2to3
;		[INPUT] AMP_SF scale correction parameters
;	va,vb,vc,vd,vbeta,vgamma,vthr	[INPUT] tap-ringing coefficients
;	ua,ub,uc,ud,ubeta,ugamma,uthr	[INPUT] tap-ringing coefficients
;	vBh,vHh,vAh,vdHh	[INPUT] H-test coefficients
;	uBh,uHh,uAh,udHh	[INPUT] H-test coefficients
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
;
;history
;	vinay kashyap (Apr02)
;	added keyword HEVT,RANGELEV (VK; Aug03)
;	added keywords TRINGF,UOFF,VOFF, and changed algorithm on which
;	  events to apply tap ringing correction (VK; Sep03)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(chip_id)
nu1=n_elements(au1) & nu2=n_elements(au2) & nu3=n_elements(au3)
nv1=n_elements(av1) & nv2=n_elements(av2) & nv3=n_elements(av3)
nas=n_elements(amp_sf) & ncu=n_elements(crsu) & ncv=n_elements(crsv)
npha=n_elements(pha) & nst=n_elements(statbit)
if np lt 10 then ok='Insufficient input parameters' else $
 if ni eq 0 then ok='CHIP_ID: missing' else $
  if nu1 eq 0 then ok='AU1: missing' else $
   if nu2 eq 0 then ok='AU2: missing' else $
    if nu3 eq 0 then ok='AU3: missing' else $
     if nv1 eq 0 then ok='AV1: missing' else $
      if nv2 eq 0 then ok='AV2: missing' else $
       if nv3 eq 0 then ok='AV3: missing' else $
	if nas eq 0 then ok='AMP_SF: missing' else $
	 if ncu eq 0 then ok='CRSU: missing' else $
	  if ncv eq 0 then ok='CRSV: missing' else $
	   if npha eq 0 then ok='PHA: missing' else $
	    if nst eq 0 then ok='STATBIT: missing' else $
	     if ni ne nu1 then ok='CHIP_ID incompatible with AU1' else $
	      if ni ne nu2 then ok='CHIP_ID incompatible with AU2' else $
	       if ni ne nu2 then ok='CHIP_ID incompatible with AU2' else $
	        if ni ne nv1 then ok='CHIP_ID incompatible with AV1' else $
	         if ni ne nv2 then ok='CHIP_ID incompatible with AV2' else $
	          if ni ne nv2 then ok='CHIP_ID incompatible with AV2' else $
	           if ni ne nas then ok='CHIP_ID incompatible with AMP_SF' else $
	            if ni ne ncu then ok='CHIP_ID incompatible with CRSU' else $
	             if ni ne ncv then ok='CHIP_ID incompatible with CRSV' else $
	              if ni ne npha then ok='CHIP_ID incompatible with PHA' else $
		       if (ni ne nst and ni ne nst/4) then ok='CHIP_ID incompatible with STATBIT'
if ok ne 'ok' then begin
  print,'Usage: hrcs_process_evt,chip_id,au1,au2,au3,av1,av2,av3,amp_sf,crsu,crsv,$'
  print,'       pha,statbit, ampscl,fineu,finev,fracu,fracv,u,v,au3cor=au3cor,av3cor=av3cor,$'
  print,'       rawu=rawu,rawv=rawv,xdet=xdet,ydet=ydet,Hfilt=Hfilt,HXfilt=HXfilt,$'
  print,'       degaps=degaps,degapf=degapf,caldb=caldb,verbose=v,$'
  print,'       dampcor=dampcor,ampgain_s=ampgain_s,ampratio_s=ampratio_s,$'
  print,'       pha_rng_min=pha_rng_min,pha_rng_mid=pha_rng_mid,$'
  print,'       pha_rng_max=pha_rng_max,hevt=hevt,rangelev=rangelev,$'
  print,'       pha_1to2=pha_1to2,pha_2to3=pha_2to3,w_1to2=w_1to2,w_2to3=w_2to3,$'
  print,'       va=va,vb=vb,vc=vc,vd=vd,vbeta=vbeta,vgamma=vgamma,vthr=vthr,voff=voff,$'
  print,'       ua=ua,ub=ub,uc=uc,ud=ud,ubeta=ubeta,ugamma=ugamma,uthr=uthr,uoff=uoff,$'
  print,'       vBh=vBh,vHh=vHh,vAh=vAh,vdHh=vdHh,$'
  print,'       uBh=uBh,uHh=uHh,uAh=uAh,udHh=udHh'
  print,'  process HRC-S amplifier data to produce positions'
  if np ne 0 then message,ok,/info
  return
endif

;	check keywords
;VERBOSE
vv=0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;CALDB
if not keyword_set(caldb) then begin
  caldb=getenv('CALDB')
  if caldb eq '' then caldb='/soft/ciao/CALDB/'
endif
;DEGAPF
ndf=n_elements(degapf) & nds=n_elements(degaps) & mds=n_tags(degaps)
degapfil=CALDB+'/data/chandra/hrc/bcf/degap/hrcsD1999-07-22gapN0002.fits'
if ndf ne 0 then begin
  szd=size(degapf) & nszd=n_elements(szd)
  if szd[nszd-2] eq 7 then degapfil=degapf[0]
endif
;DEGAPS
nodegap=0 & rddegap=0
if nds eq 0 then rddegap=1 else begin
  if mds eq 0 then begin	;(DEGAPS is not a structure
    if degaps[0] eq 0 then nodegap=1 else rddegap=1
  endif else begin		;not structure)(structure
    sds=tag_names(degaps) & k=0
    for i=0,mds-1 do begin
      if sds[i] eq 'CHIP_ID' then k=k+1
      if sds[i] eq 'AXIS' then k=k+1
      if sds[i] eq 'TAPNUM' then k=k+1
      if sds[i] eq 'LA' then k=k+1
      if sds[i] eq 'RA' then k=k+1
      if sds[i] eq 'LB' then k=k+1
      if sds[i] eq 'RB' then k=k+1
    endfor
    if k lt 7 then begin
      if vv gt 0 then message,'DEGAPS is in unknown format',/info
      rddegap=1
    endif else degapstr=DEGAPS
  endelse			;DEGAPS is a structure)
endelse
if keyword_set(rddegap) then begin
  if vv gt 1 then message,'reading degap file: '+degapfil,/info
  degapstr=mrdfits(degapfil,1,h)
  degaps=degapstr
endif
;TRINGF and Tap-ringing coefficients
;	see Table 1 of
;	http://hea-www.harvard.edu/~juda/memos/spie2000_tap_correction.html
;	also 
;va=-40.0 & vb=-655.0 & vc=0.279 & vd=655.0 & vbeta=0.0 & vgamma=1.0 & vthr=0.2 & voff=0.
;ua=-24.0 & ub=-690.0 & uc=0.279 & ud=530.0 & ubeta=6.7 & ugamma=1.75 & uthr=0.2 & uoff=0.
if not keyword_set(va) then	va=-40.0
if not keyword_set(vb) then	vb=-655.0
if not keyword_set(vc) then	vc=0.279
if not keyword_set(vd) then	vd=655.0
if not keyword_set(vbeta) then	vbeta=0.0
if not keyword_set(vgamma) then	vgamma=1.0
if not keyword_set(vthr) then	vthr=0.2
if not keyword_set(voff) then	voff=0.0
;
if not keyword_set(ua) then	ua=-24.0
if not keyword_set(ub) then	ub=-690.0
if not keyword_set(uc) then	uc=0.279
if not keyword_set(ud) then	ud=530.0
if not keyword_set(ubeta) then	ubeta=6.7
if not keyword_set(ugamma) then	ugamma=1.75
if not keyword_set(uthr) then	uthr=0.2
if not keyword_set(uoff) then	uoff=0.0
ntf=n_elements(tringf)
if ntf ne 0 then begin
  sztf=size(tringf,/type)
  if sztf eq 7 then begin
    trings=mrdfits(CALDB[0]+tringf[0],1,htring) & nts=n_tags(trings)
    if nts gt 0 then begin
      ivaxis=0 & iuaxis=1
      if trings[0].AXIS eq 'U' then begin & ivaxis=1 & iuaxis=0 & endif
      ua=trings[iuaxis].TRING_A          & va=trings[ivaxis].TRING_A
      ub=trings[iuaxis].TRING_B          & vb=trings[ivaxis].TRING_B
      uc=trings[iuaxis].TRING_C          & vc=trings[ivaxis].TRING_C
      ud=trings[iuaxis].TRING_D          & vd=trings[ivaxis].TRING_D
      ubeta=trings[iuaxis].TRING_BETA    & vbeta=trings[ivaxis].TRING_BETA
      ugamma=trings[iuaxis].TRING_GAMMA  & vgamma=trings[ivaxis].TRING_GAMMA
      uthr=trings[iuaxis].TRING_THRESH12 & vthr=trings[ivaxis].TRING_THRESH12
      if nts gt 8 then uoff=trings[iuaxis].TRING_OFFSET
      if nts gt 8 then voff=trings[ivaxis].TRING_OFFSET
    endif
  endif
endif

;make amplifier scale corrections, based on Mike Juda's memo at
;http://hea-www.harvard.edu/~juda/memos/amp_scale/improved_scale_correction.html
;either with a dynamical determination of the range switch points and range
;over which they switch, or prespecified ones
	;	original code, based on
	;	http://hea-www.harvard.edu/~juda/memos/amp_scale/scale_correction.html
	;	is here, but note that it has a bug anyway
	;ok=where(pha gt 0 and pha lt 255,mok)
	;if mok gt 0 then ampscl[ok]=3-fix(sumamps[ok]/(ampratio_s/pha[ok]))
	;ok=where(pha ge 255,mok) & if mok gt 0 then ampscl[ok]=3
sumav=av1+av2+av3 & sumau=au1+au2+au3 & sumamps=sumav+sumau
if not keyword_set(ampratio_s) then ampratio_s=38.
if not keyword_set(ampgain_s) then ampgain_s=52.9*2.
	;default value of AMPGAIN_S of comes from Steve Murray memo quoted in
	;http://hea-www.harvard.edu/~juda/memos/amp_scale/compare.html
if not keyword_set(pha_rng_min) then pha_rng_min=30
if not keyword_set(pha_rng_mid) then pha_rng_mid=100
if not keyword_set(pha_rng_max) then pha_rng_max=180
	;default values of pha_rng_min,pha_rng_mid,pha_rng_max are set
	;based on PKS2155 coadded data
	;default values of rangelev,pha_1to2,pha_2to3,w_1to2,w_2to3 from
	;http://cxc.harvard.edu/contrib/juda/memos/amp_scale/improved_scale_correction.html
if not keyword_set(rangelev) then rangelev=125
if keyword_set(hevt) then begin
  rnglev=sxpar(hevt,'RANGELEV') & WDTHAST=1
  date_obs=sxpar(hevt,'DATE-OBS')
  if strlen(date_obs) gt 10 then begin
    yyyy=long(strmid(date_obs,0,4))
    mm=fix(strmid(date_obs,5,2))
    dd=fix(strmid(date_obs,8,2))
    if not keyword_set(rnglev) then begin
      if yyyy lt 1999 then rangelev=90
      if yyyy eq 1999 and mm lt 12 then rangelev=90
      if yyyy eq 1999 and mm eq 12 and dd lt 6 then rangelev=90
    endif else rangelev=rnglev
    ;2000:279:12:57
    if yyyy lt 2000 then WDTHAST=0
    if yyyy eq 2000 and mm lt 10 then WDTHAST=0
    if yyyy eq 2000 and mm eq 10 and dd lt 5 then WDTHAST=0
  endif else rangelev=rnglev
endif
if rangelev eq 90 then begin
  if not keyword_set(pha_1to2) then pha_1to2=51.0
  if not keyword_set(pha_2to3) then pha_2to3=99.5
endif else begin
  if not keyword_set(pha_1to2) then pha_1to2=70.5
  if not keyword_set(pha_2to3) then pha_2to3=137.5
endelse
if not keyword_set(w_1to2) then w_1to2=5.
if not keyword_set(w_2to3) then w_2to3=5.
if keyword_set(dampcor) then begin
  if vv gt 10 then plot,pha,sumamps,/xs,/ys,/nodata
  ampscl=0*amp_sf+1 & s0=0. & s1=fltarr(256) & iamp=1
  oo=where(pha lt pha_rng_min,moo) & if moo gt 0 then ampscl[oo]=1
  if moo gt 0 and vv gt 10 then oplot,pha[oo],sumamps[oo],psym=3
  oo=where(pha gt pha_rng_max,moo) & if moo gt 0 then ampscl[oo]=3
  if moo gt 0 and vv gt 10 then oplot,pha[oo],sumamps[oo],psym=3
  for i=pha_rng_min,pha_rng_max,1 do begin    	;{
    oo=where(pha eq i,moo)    
    if moo gt 0 then begin    			;(
      ss=total(sumamps[oo])/moo & s1[i]=ss    
      if vv gt 10 then oplot,pha[oo],sumamps[oo],psym=3,col=iamp    
      if ss lt s0 then begin    			;((
        d1=abs(pha[oo]-1*sumamps[oo]/ampgain_s)    
        d2=abs(pha[oo]-2*sumamps[oo]/ampgain_s)    
        d3=abs(pha[oo]-4*sumamps[oo]/ampgain_s)    
        if i lt pha_rng_mid then begin    	;(((
	  if iamp eq 1 then iamp=2    
	  ampscl[oo]=1 & oi=where(d2 lt d1,moi) & if moi gt 0 then ampscl[oo[oi]]=2    
        endif else begin    			;)))(((
	  if iamp eq 2 then iamp=3    
	  ampscl[oo]=2 & oi=where(d3 lt d2,moi) & if moi gt 0 then ampscl[oo[oi]]=3    
        endelse    				;)))
      endif else ampscl[oo]=iamp & s0=ss    	;))
    endif    					;)
  endfor						;}
  sumamp=long(sumamps)*2^(ampscl-1)
endif else begin
  if vv gt 10 then plot,pha,sumamps,/xs,/ys,/nodata
  ampscl=0*amp_sf+1
  oo=where(pha lt pha_1to2-w_1to2,moo) & if moo gt 0 then ampscl[oo]=1
  oo=where(pha ge pha_1to2-w_1to2 and pha lt pha_1to2+w_1to2,moo)
  if moo gt 0 then begin
    d1=abs(pha[oo]-1*sumamps[oo]/ampgain_s)
    d2=abs(pha[oo]-2*sumamps[oo]/ampgain_s)
    oi=where(d2 lt d1,moi) & if moi gt 0 then ampscl[oo[oi]]=2
    if vv gt 10 and moi gt 0 then oplot,pha[oo[oi]],sumamps[oo[oi]],psym=3,col=2
  endif
  oo=where(pha ge pha_1to2+w_1to2 and pha lt pha_2to3-w_2to3,moo)
  if vv gt 10 and moo gt 0 then oplot,pha[oo],sumamps[oo],psym=3,col=2
  if moo gt 0 then ampscl[oo]=2
  oo=where(pha ge pha_2to3-w_2to3 and pha lt pha_2to3+w_2to3,moo)
  if moo gt 0 then begin
    d2=abs(pha[oo]-2*sumamps[oo]/ampgain_s)
    d3=abs(pha[oo]-4*sumamps[oo]/ampgain_s)
    oi=where(d3 lt d2,moi) & if moi gt 0 then ampscl[oo[oi]]=3
    if vv gt 10 and moi gt 0 then oplot,pha[oo[oi]],sumamps[oo[oi]],psym=3,col=3
  endif
  oo=where(pha ge pha_2to3+w_2to3,moo) & if moo gt 0 then ampscl[oo]=3
  if vv gt 10 and moo gt 0 then oplot,pha[oo],sumamps[oo],psym=3,col=3
endelse
;o1=where(ampscl eq 1,mo1) & o2=where(ampscl eq 2,mo2) & o3=where(ampscl eq 3,mo3)
;oplot,pha[o1],sumamps[o1],psym=3,col=1
;oplot,pha[o2],sumamps[o2],psym=3,col=2
;oplot,pha[o3],sumamps[o3],psym=3,col=3
;plot,pha,long(sumamps)*2^(ampscl-1),psym=3

;	make tap-ringing corrections for AU3 and AV3
av3cor=av3 & au3cor=au3
;	Original code below, based on http://cxc.harvard.edu/contrib/juda/memos/tap_correction.html
;	corrected Sep03, based on http://cxc.harvard.edu/contrib/juda/memos/ring_evt_id/event_id.html
;;{original
;;oo=where(av2 gt 0,moo) & va1a2=0.*av2 & av3cor=av3
;;if moo gt 0 then va1a2[oo]=float(av1[oo])/float(av2[oo])
;;ov=where(ampscl eq 3 and av1 gt av3 and va1a2 ge vthr,mov)
;;;
;;oo=where(au2 gt 0,moo) & ua1a2=0.*au2 & au3cor=au3
;;if moo gt 0 then ua1a2[oo]=float(au1[oo])/float(au2[oo])
;;ou=where(ampscl eq 3 and au1 gt au3 and ua1a2 ge uthr,mou)
;;original}
;{Sep03
ostatus=b11001001(statbit[1,*],/otto)
ostatus11=reform(ostatus[5,*]) & ostatus12=reform(ostatus[4,*])
ov=where(ampscl eq 3 and WDTHAST*ostatus11 eq 0 and av1 gt av3 and av1 gt vthr*av2+voff,mov)
ou=where(ampscl eq 3 and WDTHAST*ostatus12 eq 0 and au1 gt au3 and au1 gt uthr*au2+uoff,mou)
;Sep03}
if mov gt 0 then av3cor[ov]=av3[ov]-((av2[ov]+vb)/va)*sin(2.*!pi*(av2[ov]-$
	vbeta*((float(av2[ov])/av1[ov])^(vgamma))-1.)/(vc*av2[ov]+vd))
if mou gt 0 then au3cor[ou]=au3[ou]-((au2[ou]+ub)/ua)*sin(2.*!pi*(au2[ou]-$
	ubeta*((float(au2[ou])/au1[ou])^(ugamma))-1.)/(uc*au2[ou]+ud))

;	get fine positions for H-test
fineu=0.*au1 & finev=0.*av1 & fracu=fineu & fracv=finev
ou=where(sumau gt 0,mou) & ov=where(sumav gt 0,mov)
;	fine positions
if mou gt 0 then fineu[ou]=(au3cor[ou]-au1[ou])/float(sumau[ou])
if mov gt 0 then finev[ov]=(av3cor[ov]-av1[ov])/float(sumav[ov])
;	fraction in center tap
if mou gt 0 then fracu[ou]=au2[ou]/float(sumau[ou])
if mov gt 0 then fracv[ov]=av2[ov]/float(sumav[ov])

;	the H-test: counts that have
;	((f1-(H+dH))/A)^2-(f2/B) < 1 OR ((f1-(H-dH))/A)^2-(f2/B) > 1
;	FAIL the H-test, where
;	f1=A2/SUM{Ai}==fr, f2=abs(A3-A1)/SUM{Ai}==fp, hf=B*sqrt(((f1-H)/A)^2-1),
;	and B,H,A,dH are defined separately for U and V,
;	see status bit definitions at
;	http://asc.harvard.edu/ciao/data_products_guide/hrc_status_bits.html
;	vBh=0.248 & vHh=1.0710 & vAh=0.2706 & vdHh=0.0350
;	uBh=0.262 & uHh=1.0180 & uAh=0.2706 & udHh=0.0350
if not keyword_set(vBh) then	vBh=0.248
if not keyword_set(vHh) then	vHh=1.0710
if not keyword_set(vAh) then	vAh=0.2706
if not keyword_set(vdHh) then	vdHh=0.0350
;
if not keyword_set(uBh) then	uBh=0.262
if not keyword_set(uHh) then	uHh=1.0180
if not keyword_set(uAh) then	uAh=0.2706
if not keyword_set(udHh) then	udHh=0.0350
;
lt1v=((fracv-(vHh+vdHh))/vAh)^2-(abs(finev)/vBh)^2
lt1u=((fracu-(uHh+udHh))/uAh)^2-(abs(fineu)/uBh)^2
gt1v=((fracv-(vHh-vdHh))/vAh)^2-(abs(finev)/vBh)^2
gt1u=((fracu-(uHh-udHh))/uAh)^2-(abs(fineu)/uBh)^2
HXfilt=where((lt1u lt 1 or gt1u gt 1) or (lt1v lt 1 or gt1v gt 1),moHX)
Hfilt=where((lt1u ge 1 and gt1u le 1) and (lt1v ge 1 and gt1v le 1),moH)
;
if vv ge 5 and !d.name eq 'X' then begin
  pmulti=!p.multi & !p.multi=[0,1,2]
  plot,fineu,fracu,/xs,/ys,psym=3,$
	xtitle='U fine positions',ytitle='U center tap fraction'
  oplot,fineu[Hfilt],fracu[Hfilt],psym=3,col=1
  if vv ge 10 then oplot,fineu[HXfilt],fracu[HXfilt],psym=3,col=2
  plot,finev,fracv,/xs,/ys,psym=3,$
	xtitle='V fine positions',ytitle='V center tap fraction'
  oplot,finev[Hfilt],fracv[Hfilt],psym=3,col=1
  if vv ge 10 then oplot,finev[HXfilt],fracv[HXfilt],psym=3,col=2
  !p.multi=pmulti
endif

;	apply degap
u=crsu+fineu & v=crsv+finev	;initialize
if not keyword_set(nodegap) then begin
  nrow=n_elements(degapstr.CHIP_ID) & mds=n_tags(degapstr)
  ncoeff=(mds-3)/2 & al=fltarr(ncoeff) & ar=al
  for i=0L,nrow-1L do begin
    d=degapstr[i]
    if vv ge 1 then kilroy,dot=strtrim(i,2)+'.. '
    for k=0,ncoeff-1 do al[k]=d.(k+3)
    for k=0,ncoeff-1 do ar[k]=d.(k+3+ncoeff)

    if d.AXIS eq 'U' then begin
      ou=where(chip_ID eq d.CHIP_ID and crsu eq d.TAPNUM,mou)
      if mou gt 0 then begin
	if vv gt 1 then kilroy,dot=strtrim(mou,2)
	xx=fineu[ou] & xxx=0.*xx & s=intarr(mou)+1
	olt=where(xx lt 0,molt) & ort=where(xx ge 0,mort)
	if molt gt 0 then s[olt]=-1
	xx=abs(xx)
	if molt gt 0 then for k=0,ncoeff-1 do $
		xxx[olt]=xxx[olt]+al[k]*(xx[olt])^(k+1)*s[olt]
	if mort gt 0 then for k=0,ncoeff-1 do $
		xxx[ort]=xxx[ort]+ar[k]*(xx[ort])^(k+1)*s[ort]
	u[ou]=crsu[ou]+xxx
	plot,fineu[ou],xxx,psym=3,/xs,/ys,$
		xtitle='U fine',ytitle='degapped U fine',title=d.TAPNUM
      endif
    endif

    if d.AXIS eq 'V' then begin
      ov=where(chip_ID eq d.CHIP_ID and crsv eq d.TAPNUM,mov)
      if mov gt 0 then begin
	if vv gt 2 then kilroy,dot=strtrim(mov,2)
	xx=finev[ov] & xxx=0.*xx & s=intarr(mov)+1
	olt=where(xx lt 0,molt) & ort=where(xx ge 0,mort)
	if molt gt 0 then s[olt]=-1
	xx=abs(xx)
	if molt gt 0 then for k=0,ncoeff-1 do $
		xxx[olt]=xxx[olt]+al[k]*(xx[olt])^(k+1)*s[olt]
	if mort gt 0 then for k=0,ncoeff-1 do $
		xxx[ort]=xxx[ort]+ar[k]*(xx[ort])^(k+1)*s[ort]
	v[ov]=crsv[ov]+xxx
	plot,finev[ov],xxx,psym=3,/xs,/ys,$
		xtitle='V fine',ytitle='degapped V fine',title=d.TAPNUM
      endif
    endif

    if vv ge 5 then print,d

  endfor
endif
rawu=u*256.+128. & rawv=v*256.+128.
;xdet=58223.056-rawv

end
