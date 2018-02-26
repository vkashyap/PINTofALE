function id_to_emis,idstr,abund=abund,ldir=ldir,dwvl=dwvl,okmult=okmult,$
	outstr=outstr,Z=Z,ion=ion,wvl=wvl,logT=logT,multidx=multidx, _extra=e
;+
;function	id_to_emis
;	return emissivities EMIS(LOGT,WVL) for all lines that have been ID'd
;	(but see keyword OUTSTR)
;	
;	in case of multiple IDs for a single observed line, the emissivities
;	of the IDs are added together, and assumed for all practical purposes
;	to be the strongest line modified by the other IDs.
;	in case of IDs spanning different elements, the emissivities are
;	added as before, but weighted with their relative abundances.
;
;warning
;	unknown IDs are ignored, and the output won't refer to them.
;
;parameters
;	idstr	[INPUT; required] the ID structure that comes out of LINEID
;
;keywords
;	abund	[INPUT] abundances
;		* default is to use Anders & Grevesse
;	ldir	[INPUT] string array containing line database directory
;		names, one for each ID.
;		* if scalar or single element, look in this directory for all
;		* if size does not match number of IDs but matches number of
;		  observed lines, then for each line get the IDs from stated
;		  directory
;		* if size does not make sense, use just the first element
;	dwvl	[INPUT; default=5e-3] search the line database for the
;		wavelength of interest with this assumed roundoff error
;	okmult	[INPUT] if set, does NOT add up the emissivities but returns
;		a positional mapping from emissivity to ID in MULTIDX
;	outstr	[INPUT] if set, returns the output in the form of a RD_LINE
;		structure
;	Z	[OUTPUT] atomic numbers of dominant ID for each line
;	ion	[OUTPUT] ionic state of dominant ID for each line
;	wvl	[OUTPUT] peak-emissivity weighted average of wvls of all IDs
;	logT	[OUTPUT] log_10(Temperatures [K]) at which output is
;		defined
;	multidx	[OUTPUT] long array mapping emissivity location to ID'd line
;
;	_extra	[INPUT ONLY] use to pass defined keywords to subroutines:
;		RD_LIST: INCIEQ
;		RD_LINE: PRES,LOGP,N_E
;		FOLD_IONEQ: CHIFIL,CHIDIR,EQFILE
;
;history
;	vinay kashyap (Oct98)
;-

message,'OBSOLETE!',/informational

;	usage
nid=n_tags(idstr)
if nid eq 0 then begin
  print,'Usage: emis=id_to_emis(idstr,abund=abund,ldir=ldir,dwvl=dwvl,/okmult,$'
  print,'       /outstr,Z=Z,ion=ion,wvl=wvl,multidx=multidx;'
  print,'       RD_LIST:incieq; RD_LINE:pres,logp,n_e; FOLD_IONEQ:chifil,chidir,eqfile)'
  return,[-1L]
endif

;	check input
tid=tag_names(idstr) & owvl=idstr.(0) & nwvl=n_elements(owvl)
ok='ok'
if tid(0) ne 'WVL' then ok='input structure of unknown form' else $
 if nwvl ne nid-1 then ok='ID structure not in standard format' else $
  if tid(1) eq 'WVL_COMMENT' then ok='ID structure contains no data'
if ok ne 'ok' then begin
  message,ok,/info & return,[-1L]
endif

;	extract info from IDSTR
mZ=[0] & mIon=[0] & mwvl=[0.] & mIdx=[-1L] & mm=0L
for i=0,nwvl-1 do begin
  tmp=idstr.(i+1)	;structure of match
  mZ=[mZ,tmp.Z] & mIon=[mIon,tmp.ION]		;matching Z,ion
  mwvl=[mwvl,tmp.WVL]				;matching wvl
  m=n_elements(tmp.Z)				;how many matches
  mIdx=[mIdx,lonarr(m)+i]			;keep track
  mm=mm+m					;total matches
endfor
if mm gt 1 then begin
  mZ=mZ(1:*) & mIon=mIon(1:*) & mwvl=mwvl(1:*) & mIdx=mIdx(1:*)
endif else begin
  message,'no matches found?!?',/info & return,[-1L]
endelse

;	keywords
if n_elements(abund) lt 30 then abund=getabund('anders & grevesse')
;
if not keyword_set(dwvl) then dwvl=5e-3
;
ndir=n_elements(ldir)
;
initstuff,atom,rom

;	make up list for input to RD_LIST
lnlst=strarr(mm) & sep='	'
for i=0,mm-1 do begin
  if mZ(i) gt 0 then c=atom(mZ(i)-1) else c='X'		;if mZ=0, unknown
  if mIon(i) gt 0 then c=c+' '+rom(mIon(i)-1)		;if mIon=0, unknown
  ;c=atom(mZ(i)-1)+' '+rom(mIon(i)-1)
  c=c+sep+strtrim(mwvl(i),2)+' ? '+strtrim(dwvl,2)
  if ndir ne 0 then begin
    ldbdir=ldir(0)				;default
    if ndir eq mm then ldbdir=ldir(i)		;one for each ID
    if ndir eq nwvl then ldbdir=ldir(mIdx(i))	;one for each line
    c=c+sep+ldbdir
  endif
  lnlst(i)=c
endfor

;	read in line emissivities
lstr=rd_list(lnlst,sep=sep, _extra=e)

;	just get the appropriate emissivities and exit
if keyword_set(okmult) then begin
  multidx=mIdx & if keyword_set(outstr) then return,lstr
  Z=lstr.Z & ion=lstr.ION & wvl=lstr.WVL & logT=lstr.LOGT
  return,lstr.LINE_INT
endif

;	catch errors
if n_elements(lstr.wvl) ne mm then stop,'RD_LIST did not work?'

;	now cast the line emissivity structure to match the IDs
if nwvl eq mm then begin	;(oh sweet coincidence!
  multidx=lindgen(nwvl) & if keyword_set(outstr) then return,lstr
  Z=lstr.Z & ion=lstr.ION & wvl=lstr.WVL & logT=lstr.LOGT
  return,lstr.LINE_INT
endif else begin		;)(the more usual case
  ;	outputs
  multidx=lindgen(nwvl) & Z=intarr(nwvl) & ion=Z & wvl=fltarr(nwvl)
  logT=lstr.logT & nT=n_elements(logT) & emis=dblarr(nT,nwvl)

  idesig=0 & ieconf=0
  if n_elements(lstr.desig) ne 1 then idesig=1
  if n_elements(lstr.config) ne 1 then ieconf=1
  if idesig then desig=strarr(2,nwvl) else desig=['']
  if ieconf then econf=strarr(2,nwvl) else econf=['']
  src=lonarr(nwvl)

  for i=0,nwvl-1 do begin		;{for each line, accumulate output
    oo=where(mIdx eq i,moo) & if moo eq 0 then message,'bug!'
    zz=mZ(oo) & jon=mIon(oo) & ww=mwvl(oo) & xmx=0.*ww
    line=((lstr.LINE_INT)(*,oo))
    aa=abund(zz-1) & idom=0 & xdom=0. & tmp=reform(emis(*,i))
    if moo gt 1 then begin			;(multiple IDs
      for j=0,moo-1 do begin		;{for each ID
        xmax=max(line(*,j)) & xmx(j)=xmax
        if xmax gt xdom then begin & xdom=xmax & idom=j & endif
        tmp=tmp+aa(j)*line(*,j)
      endfor				;J=0,MOO-1}
      emis(*,i)=tmp(*)/aa(idom)
      Z(i)=zz(idom) & ion(i)=jon(idom)
      if total(xmx) eq 0 then xmx(*)=1. & wvl(i)=total(xmx*ww)/total(xmx)
      if idesig then desig(*,i)=(lstr.desig)(*,idom)
      if ieconf then econf(*,i)=(lstr.config)(*,idom)
      src(i)=lstr.src(idom)
    endif else begin				;)(one line, one ID
      Z(i)=zz(0) & ion(i)=jon(0) & wvl(i)=ww(0) & emis(*,i)=line(*,0)
      if idesig then desig(*,i)=(lstr.desig)(*,i)
      if ieconf then econf(*,i)=(lstr.config)(*,i)
      src(i)=lstr.src(i)
    endelse					;lines v/s IDs)
  endfor				;I=0,NWVL-1}
endelse				;number of lines v/s number of IDs)

if keyword_set(outstr) then return,create_struct('LINE_INT',emis,$
	'LOGT',logT,'WVL',wvl,'Z',z,'ION',ion,$
	'DESIG',desig,'CONFIG',econf,'SRC',src)

return,emis
end
