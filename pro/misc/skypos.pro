pro skypos,ra,dec,out1,out2,verbose=verbose,i24=i24,osxg=osxg,$
	hrs=hrs,idnum=idnum,idpre=idpre,idform=idform,PM=PM,$
	chop=chop,osep=osep, _extra=e
;+
;procedure	skypos
;	process input RA and Dec and produce the appropriate output 
;
;syntax
;	skypos,ra,dec,out1,out2,/i24,/osxg,/hrs,/idnum,idpre=idpre,$
;	idform=idform,/PM,chop=chop,osep=osep,verbose=verbose
;
;parameters
;	ra	[INPUT; required] Right Ascension
;		* may be a real or string array
;		* if real, assumed to be decimal degrees unless
;		  the keyword I24 is set
;		* if string, will check whether decimal or sexagesimal
;		  -- sexagesimal if ':' or ' ' exist, decimal otherwise
;	dec	[INPUT; required] Declination
;		* size _must_ match RA
;		* may be real or string array
;		* if string, will check whether decimal or sexagesimal
;		  -- sexagesimal if ':' or ' ' or '	' exist, decimal otherwise
;		NOTE: all calculations are performed in double precision
;	out1	[OUTPUT; required] output, depends on keyword (see below)
;	out2	[OUTPUT; optional] output, depends on keyword (see below)
;		* by default,
;		  OUT1 = RA [decimal deg] and OUT2 = Dec [decimal deg]
;		* if keyword IDNUM is set, then
;		  OUT1 = IAU style ID number = OUT2
;		  -- in addition, keyword IDFORM overrides keywords
;		     OSXG and HRS (see below)
;		* if keyword OSXG is set, then
;		  OUT1 = RA [h:m:s] and OUT2 = Dec [d:m:s]
;		* if keyword HRS is set, then
;		  OUT1 = RA [decimal hours] and OUT2 = Dec [decimal deg]
;
;keywords
;	i24	[INPUT] set this to indicate whether the input RA is
;		in decimal hours, not degrees
;	osxg	[INPUT] if set, output RA and Dec will be in sexagesimal
;		units, not decimal
;	hrs	[INPUT] if set, output RA will be in hours
;		(in sexagesimal if OSXG is set, in decimal otherwise)
;		* if neither OSXG nor HRS is set, then RA will be
;		  in decimal degrees
;	idnum	[INPUT] if set, an IAU style name is constructed and
;		returned as a string in OUT1
;	idpre	[INPUT] a uniform prefix to prepend to IDNUM
;		* default is "J"
;	idform	[INPUT] a string denoting the format for how IDNUM
;		is written out
;		* default is 'HMS.sXDMS'
;		* ignored if IDNUM is not set
;		* format rules:
;		  -- always has RA first and Dec next
;		  -- X denotes the sign of the declination (optional)
;		  -- RA is specified either as
;			"R.r" [deg] or "HMS.s" [sexagesimal]
;		  -- Dec is specified either as
;			"D.d" [deg] or "DMS.s" [sexagesimal]
;		  -- repetions of small letters after the period denote the
;		     number of digits, e.g., "R.rrr" is written as f7.3
;		  -- repetitions of capital letters before the period has
;		     no effect unless they are greater than the minimum
;		     required, e.g.,
;		     "R.r" == f5.1 == "RR.r" == "RRR.r", but "RRRR.r" == f6.1
;	PM	[INPUT] if set, replaces "+/-" with "P/M" in the sign of Dec
;	chop	[INPUT] chops off the specified number of digits off the
;		right hand end of the IDFORMatted output
;	osep	[INPUT] set the separator between fields in a sexagesimal
;		output (see keyword OSXG) to OSEP rather than ":"
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	DEG2HRS
;	STR_2_ARR
;	SYZE
;
;history
;	vinay kashyap (Jul03)
;	bug correction when output was R.rrr, and added keyword PM (VK; Dec03)
;	added keyword OSEP (VK; May04)
;	updated RA and Dec outputs to ss.ssss by default and changed
;	  all internal calcs to double precision (VK; Aug15)
;-

;	usage
ok='ok' & np=n_params() & nra=n_elements(ra) & ndec=n_elements(dec)
if np lt 3 then ok='Insufficient parameters' else $
 if nra eq 0 then ok='RA: undefined' else $
  if ndec eq 0 then ok='Dec: undefined' else $
   if nra ne ndec then ok='RA and Dec are incompatible'
if ok ne 'ok' then begin
  print,'Usage: skypos,ra,dec,out1,out2,/i24,/osxg,/hrs,$'
  print,'       /idnum,idpre=idpre,idform=idform,/pm,chop=chop,$'
  print,'       osep=osep,verbose=verbose'
  print,'  convert input (RA,Dec) to [degrees,degrees] or:'
  print,'  if OSXG:	[h:m:s,d:m:s]'
  print,'  if HRS:	[hours,degrees]'
  print,'  if IDNUM:	IDPRE+RA+Dec'
  if np ne 0 then message,ok,/informational
  return
endif

;	keywords
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
;
ss=':' & if keyword_set(osep) then ss=string(osep[0])

;	what is the input like?
szra=size(ra,/type)
if szra ne 7 then begin			;(all numbers
  xra = [double(ra)]
endif else begin			;RA)(string
  xra = dblarr(nra)
  for i=0L,nra-1L do begin	;{for each entry
    rsep=''	;default is to assume degrees
    icol=where([strpos(ra[i],':',0)] ge 0,micol)
    if micol eq 0 then begin		;(no ":" present
      isp=where([strpos(ra[i],' ',0)] ge 0,misp)
      if misp eq 0 then begin		;(no " " present
        itab=where([strpos(ra[i],'	',0)] ge 0,mitab)
        if mitab ne 0 then rsep='	'	;"<tab>" present
      endif else rsep=' '		;MISP)
    endif else rsep=':'			;MICOL)
    if keyword_set(rsep) then begin
      cc=strsplit(ra[i],rsep,/extract) & ncc=n_elements(cc)
      rh=cc[0] & rm=0. & rs=0.
      if ncc gt 1 then rm=cc[1]
      if ncc gt 2 then rs=cc[2]
      xra[i]=rh+(rm+(rs/60.D))/60.D
      if rsep eq ':' then i24=1	;assume RA is in hours
    endif else begin
      xra[i]=float(ra[i])
    endelse
  endfor			;I=0,NRA-1}
endelse					;SZRA=7)
if keyword_set(i24) then begin
  xra = xra*15.D
  oo=where(xra gt 360.,moo)
  if moo gt 0 then begin
    message,'Some RA were already input in [degrees]; keyword I24',/informational
    message,'appears to have been set, incorrectly.  Ignoring.',/informational
    xra = xra/15.D
  endif
endif
;
szdec=size(dec,/type)
if szdec ne 7 then begin			;(all numbers
  xdec = [double(dec)]
endif else begin			;DEC)(string
  xdec = dblarr(ndec)
  for i=0L,ndec-1L do begin	;{for each entry
    dsep=''	;default is to assume degrees
    icol=where([strpos(dec[i],':',0)] ge 0,micol)
    if micol eq 0 then begin		;(no ":" present
      isp=where([strpos(dec[i],' ',0)] ge 0,misp)
      if misp eq 0 then begin		;(no " " present
        itab=where([strpos(dec[i],'	',0)] ge 0,mitab)
        if mitab ne 0 then dsep='	'	;"<tab>" present
      endif else dsep=' '		;MISP)
    endif else dsep=':'			;MICOL)
    if keyword_set(dsep) then begin
      cc=strsplit(dec[i],dsep,/extract) & ncc=n_elements(cc)
      sd=1 & if long(cc[0]) lt 0 then sd=-1
      dd=abs(long(cc[0])) & dm=0.D & ds=0.D
      if ncc gt 1 then dm=long(cc[1])
      if ncc gt 2 then ds=double(cc[2])
      xdec[i]=sd*(dd+(dm+(ds/60.D))/60.D)
    endif else begin
      xdec[i]=double(dec[i])
    endelse
  endfor			;I=0,NDEC-1}
endelse					;SZDEC=7)

;	decode IDFORM
runit='deg' & dunit='deg'	;the output formats
rrform='f6.2' & rhform='i2.2' & rmform='i2.2' & rsform='f8.4'
ddform='f6.2' & dmform='i2.2' & dsform='f8.4'
if keyword_set(osxg) then ddform='i2.2'
if keyword_set(idnum) then begin
  prefixID = 'J' & if keyword_set(idpre) then prefixID=strtrim(idpre[0],2)
  idf='HMS.sXDMS.s' & if keyword_set(idform) then idf=strtrim(idform[0],2)
  if strpos(idf,'R',0) ge 0 then runit='deg' else $
   if strpos(idf,'H',0) ge 0 then runit='hms' else begin
     message,'IDFORM cannot be understood',/informational
     runit='junk'
   endelse
  if strpos(idf,'D',0) ge 0 then begin
    if strpos(idf,'M',0) ge 0 then dunit='dms' else dunit='deg'
  endif else begin
    message,'IDFORM cannot be understood',/informational
    dunit='junk' 
  endelse

  if runit ne 'junk' and dunit ne 'junk' then begin	;(go ahead and decode

    ;	possible formats
    if runit eq 'deg' then rexpr='(R+)(\.*)(r*)' else rexpr='(H+)(M*)(S*)(\.*)(s*)'
    if dunit eq 'deg' then dexpr='(X*)(D+)(\.*)(d*)' else dexpr='(X*)(D+)(M*)(S*)(\.*)(s*)'
    pos=stregex(idf,rexpr+dexpr,length=len,/subex)

    kk=0L
    if runit eq 'deg' then begin	;(R.r
      Ir0=pos[kk+1] & Ir1=pos[kk+2] & Ir2=pos[kk+3]
      Lr0=len[kk+1] & Lr1=len[kk+2] & Lr2=len[kk+3]
      kk=3
      rlen=Lr0 > 3
      rrform='i'+strtrim(rlen,2)+'.'+strtrim(rlen,2)	;RRR
      if Lr1 gt 0 then begin
	rlen=rlen+1
	rrform='f'+strtrim(rlen,2)+'.0'			;RRR.
      endif
      if Lr2 gt 0 then begin
	rlen=rlen+Lr2
	rrform='f'+strtrim(rlen,2)+'.'+strtrim(Lr2,2)	;RRR.rrr
      endif
      rform = rrform
    endif else begin			;R.r)(HMS.s
      Ir0=pos[kk+1L] & Ir1=pos[kk+2] & Ir2=pos[kk+3] & Ir3=pos[kk+4] & Ir4=pos[kk+5]
      Lr0=len[kk+1L] & Lr1=len[kk+2] & Lr2=len[kk+3] & Lr3=len[kk+4] & Lr4=len[kk+5]
      kk=5
      rhlen=Lr0 > 2
      rhform='i'+strtrim(rhlen,2)+'.'+strtrim(rhlen,2)		;HH
      if Lr1 gt 0 then begin
	rmlen=Lr1 > 2
	rmform='i'+strtrim(rmlen,2)+'.'+strtrim(rmlen,2)	;MM
	if Lr2 gt 0 then begin
	  rslen=Lr2 > 2
	  rsform='i'+strtrim(rslen,2)+'.'+strtrim(rslen,2)	;SS
	endif
	if Lr3 gt 0 then begin
	  rslen=rslen+1
	  rsform='f'+strtrim(rslen,2)+'.0'			;SS.
	endif
	if Lr4 gt 0 then begin
	  rslen=rslen+Lr4
	  rsform='f'+strtrim(rslen,2)+'.'+strtrim(Lr4,2)	;SS.ss
	endif
	rform = rhform+','+rmform+','+rsform
      endif
    endelse				;HMS.s)

    ;	the sign of the Dec
    Ix=pos[kk] & Lx=len[kk] & sform='a1'
    if Lx gt 1 then sform='a'+strtrim(Lx,2)
    kk=kk+1L

    ;	Dec
    if dunit eq 'deg' then begin	;(D.d
      Id0=pos[kk+1] & Id1=pos[kk+2] & Id2=pos[kk+3]
      Ld0=len[kk+1] & Ld1=len[kk+2] & Ld2=len[kk+3]
      dlen=Ld0 > 2
      ddform='i'+strtrim(dlen,2)+'.'+strtrim(dlen,2)	;DD
      if Ld1 gt 0 then begin
	dlen=dlen+1
	ddform='f'+strtrim(dlen,2)+'.0'			;DD.
      endif
      if Ld2 gt 0 then begin
	dlen=dlen+Ld2
	ddform='f'+strtrim(dlen,2)+'.'+strtrim(Ld2,2)
      endif
      dform = ddform
    endif else begin			;D.d)(DMS.s
      Id0=pos[kk+1L] & Id1=pos[kk+2] & Id2=pos[kk+3] & Id3=pos[kk+4] & Id4=pos[kk+5]
      Ld0=len[kk+1L] & Ld1=len[kk+2] & Ld2=len[kk+3] & Ld3=len[kk+4] & Ld4=len[kk+5]
      ddlen=Ld0 > 2
      ddform='i'+strtrim(ddlen,2)+'.'+strtrim(ddlen,2)
      if Ld1 gt 0 then begin
	dmlen=Ld1 > 2
	dmform='i'+strtrim(dmlen,2)+'.'+strtrim(dmlen,2)
	if Ld2 gt 0 then begin
	  dslen=Ld2 > 2
	  dsform='i'+strtrim(dslen,2)+'.'+strtrim(dslen,2)
	endif
	if Ld3 gt 0 then begin
	  dslen=dslen+1
	  dsform='f'+strtrim(dslen,2)+'.0'
	endif
	if Ld4 gt 0 then begin
	  dslen=dslen+Ld4
	  dsform='f'+strtrim(dslen,2)+'.'+strtrim(Ld4,2)
	endif
	dform = ddform+','+dmform+','+dsform
      endif
    endelse				;DMS.s)

    if vv gt 1 then begin
      print,'IDFORM : '+idf
      print,'format : ('+rform+','+sform+','+dform+')'
    endif

  endif							;no 'junk')

endif
if keyword_set(hrs) or keyword_set(osxg) then runit='hms'
if keyword_set(osxg) then dunit='dms'

;	figure out the output
out1 = xra & out2 = xdec	;the default, in [degrees]

;	RA in [hours]
if runit eq 'hms' then begin
  xra = xra/15.
  out1 = xra
endif

;	fancier output
if runit eq 'hms' then begin
  rh=fix(xra) & rm=fix((xra-rh)*60.) & rs=double( (xra-rh-(rm/60.))*3600. )
  crh=strtrim(string(rh,'('+rhform+')'),2)
  crm=strtrim(string(rm,'('+rmform+')'),2)
  crs=strtrim(string(rs,'('+rsform+')'),2)
  o1=where(rs lt 10,mo1)
  if mo1 gt 0 and strpos(rsform,'i',0) lt 0 then crs[o1]='0'+crs[o1]
  out1=crh+ss+crm+ss+crs
  if keyword_set(idnum) then out1=prefixID+crh+crm+crs
  if vv gt 100 then for i=0,nra-1 do print,xra[i],' ',out1[i],' ',(deg2hrs(15.*xra[i]))[0]
endif else begin
  crr=strtrim(string(xra,'('+rrform+')'),2)
  o1=where(xra lt 100,mo1) & if mo1 gt 0 then crr[o1]='0'+crr[o1]
  o2=where(xra lt 10,mo2) & if mo2 gt 0 then crr[o2]='0'+crr[o2]
  if keyword_set(idnum) then out1=prefixID+crr
  if keyword_set(chop) then begin
    out1=strmid(out1,0,strlen(out1)-chop[0])
  endif
endelse
;
sdec='+'+strarr(ndec)
om=where(xdec lt 0,mom) & if mom gt 0 then sdec[om]='-'
if keyword_set(PM) then begin
  sdec='P'+strarr(ndec) & if mom gt 0 then sdec[om]='M'
endif
if dunit eq 'dms' then begin
  xdec=abs(xdec) & dd=fix(xdec) & dm=fix((xdec-dd)*60.) & ds=double( (xdec-dd-(dm/60.))*3600. )
  cdd=strtrim(string(dd,'('+ddform+')'),2)
  cdm=strtrim(string(dm,'('+dmform+')'),2)
  cds=strtrim(string(ds,'('+dsform+')'),2)
  o1=where(ds lt 10,mo1)
  if mo1 gt 0 and strpos(dsform,'i',0) lt 0 then cds[o1]='0'+cds[o1]
  out2=sdec+cdd+ss+cdm+ss+cds
  if keyword_set(idnum) then out1=out1+sdec+cdd+cdm+cds
  if vv gt 100 then for i=0,ndec-1 do print,xdec[i],' ',out2[i],' ',(deg2hrs(xdec[i],/dec))[0]
endif else begin
  cdd=strtrim(string(abs(xdec),'('+ddform+')'),2)
  o1=where(abs(xdec) lt 10,mo1) & if mo1 gt 0 then cdd[o1]='0'+cdd[o1]
  ;o2=where(abs(xdec) lt 10,mo2) & if mo2 gt 0 then cdd[o2]='0'+cdd[o2]
  if keyword_set(idnum) then out1=out1+sdec+cdd
  if keyword_set(chop) then begin
    out1=strmid(out1,0,strlen(out1)-chop[0])
  endif
endelse

return
end
