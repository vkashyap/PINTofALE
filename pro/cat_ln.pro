function cat_ln,AA,B,pick=pick,ask=ask,okeV=okeV,comm=comm,flst=flst,$
	reord=reord, _extra=e
;+
;function	cat_ln
;	returns a concatenated structure of Line emissivities
;
;syntax
;	C=cat_ln(A,B,pick=pick,ask=ask,comm=comm,reord=reord,/okev,$
;	abund=abund,/noph,effar=effar,wvlar=wvlar,dem=dem,/temp,/ikev,/regrid)
;
;parameters
;	A	[INPUT; required] output of RD_LINE
;	B	[INPUT] output of RD_LINE to be merged with A
;		* if any wavelengths are common, the IDs in A take
;		  precedence (but see keyword ASK)
;		* if the temperature grids do not match, the emissivities
;		  in B are interpolated to the same grid as in A
;		* if not given, then the program *prints* a summary
;		  of the lines to the screen
;
;keywords
;	pick	[INPUT; default: all] position indices of elements in A
;		to be selected (or not).
;		* if scalar, returns all elements EXCEPT the specified one
;		* has no effect on B
;		* use these two to delete entries, for e.g.
;	ask	[INPUT] specifies how to handle duplicate lines
;		* if not set, doesn't care about duplicates and keeps all
;		* if set (see exceptions below), will discard duplicates
;		  from B
;		* ASK='n' will reverse the normal by throwing away the
;		  line from A instead of the one from B
;		* ASK='r' will >replace< the emissivity, etc. in A with
;		  that of B
;		* ASK='k' will act just as though ASK has not been set
;	okeV	[INPUT] if set, converts wavelengths [Ang] to energy [keV]
;		in the screen listing
;	comm	[INPUT] appends a comment to each line while printing
;		* takes the comments from A.DESIG[1,i]+A.DESIG[0,i]
;		* truncates the comment to a maximum of 30 characters
;		  unless COMM is set to a higher number
;	flst	[INPUT] if set, prints out the line list in the form that
;		RD_LIST likes
;		* if FLST is a string, writes to file FLST
;	reord	[INPUT] if set, rearrange the order of the output
;		* +ve: descending order of fluxes
;		* -ve: ascending order of fluxes
;		* if abs(REORD) .NE. 1, limit output to first REORD entries
;	_extra	[INPUT ONLY] use to pass defined keywords to
;		LINEFLX: ABUND, NOPH, EFFAR, WVLAR, DEM, TEMP, REGRID
;
;restrictions
;	* works ONLY on the outputs of RD_LINE
;	* requires subroutines:
;	  -- LINEFLX [GETABUND, WHEE]
;	  -- RD_LINE
;	  -- CREATE_STRUCT
;	  -- INICON
;
;history
;	vinay kashyap (Jun98, based on CAT_ID)
;	changed keyword LNLST to FLST, added ABUND and REORD (VK; Jul98)
;	delete keyword ABUND, and added call to LINEFLX (VK; Aug98)
;	corrected bug that was skipping B because of roundoff error in
;	  first element; bug in handling ECONF (VK; Nov98)
;	changed keyword KEV to OKEV to avoid conflict with LINEFLX keyword
;	  of same name (VK; 99Apr)
;	added support for structure field JON (VK; 99May)
;	allowed merging of mismatched LOGT grids (VK; 99Jul)
;	converted to IDL5 array notation (VK; OctMM)
;	added full RD_LIST compatible printout for FLST; changed output
;	  file wavelength format from f8.4 to f10.5 (VK; JanMMI)
;	bug -- crashing on <empty> for FLST -- fixed (VK; Jul04)
;	added extra columns (FLUX, logTMAX) to output of FLST (VK; Feb08)
;	guard against index mismatches (VK; Oct16)
;-

;	Usage
sza=size(AA) & szb=size(B) & nsa=n_elements(sza) & nsb=n_elements(szb)
if sza[nsa-2] ne 8 then begin
  print,'Usage: C=cat_ln(A,B,pick=pick,ask=ask,/okeV,comm=comm,reord=+-r,$'
  print,'       abund=a,/noph,effar=ea,wvlar=wa,dem=dem,/temp,/ikev,/regrid,$
  print,'	flst=flst)'
  print,'  returns duplicate-less concatenation of RD_LINE structures A and B'
  print,'--------------------------------------------------------------------'
  return,create_struct('LINE_INT',[-1.],'LOGT',[-1.],'WVL',[-1.],$
	'Z',[-1],'ION',[-1],'DESIG',['X'],'CONFIG',['X'],'SRC',[-1],$
	'JON',[-1])
endif
A=AA

;	catch input errors
c1='ok'
;
if szb[0] ne 0 and szb[nsb-2] ne 8 then begin
  message,'Input not a structure',/info & return,A
endif
;
tagA=tag_names(A)
if tagA[0] ne 'LINE_INT' or tagA[1] ne 'LOGT' then c1=$
	'Structure not in right format'
;
if szb[nsb-2] eq 8 then begin
  tagB=tag_names(B)
  if tagB[0] ne 'LINE_INT' or tagB[1] ne 'LOGT' then c1=$
	'Structure in incorrect form'
endif
;
if c1 ne 'ok' then begin & message,c1,/info & return,A & endif

;	handle temperature grid mismatches
ntA=n_elements(A.LOGT)
if szb[nsb-2] eq 8 then begin
  ntB=n_elements(B.LOGT)
  if ntA ne ntB then begin
    message,'temperature grid mismatch: morphing B to suit A',/info
    fBB=B.LINE_INT & wBB=B.WVL & nwB=n_elements(wBB)
    fB0=0*fBB[0]+fltarr(ntA,nwB)
    if ntB gt 1 then begin
      for jw=0L,nwB-1L do fB0[*,jw]=interpol(fBB[*,jw],B.LOGT,A.LOGT)
    endif else for jw=0L,nwB-1L do fB0[*,jw]=fBB[jw]
  endif else fB0=B.LINE_INT
endif

;{	a small digression to "reset" A
fA=A.LINE_INT & wA=A.WVL & zA=A.Z & ionA=A.ION & jonA=A.JON & srcA=A.SRC
dsgA=A.desig & cfgA=A.config
nwA=n_elements(wA) & iA=lindgen(nwA)
if n_elements(pick) ne 0 then begin
  szo=size(pick)
  if szo[0] eq 0 then begin		;{scalar -- delete this element!
    oA=where(iA ne pick,moA)
    if moA eq 0 then begin		;nothing to return?
      if n_tags(B) ne 0 then return,B else return,A
    endif
  endif else begin			;}{vector -- pick only these!
    oA=pick
  endelse				;}
  fA=fA[*,oA]
  wA=wA[oA] & zA=zA[oA] & ionA=ionA[oA] & jonA=jonA[oA] & srcA=srcA[oA]
  if n_elements(dsgA) gt 1 then dsgA=dsgA[*,oA]
  if n_elements(cfgA) gt 1 then cfgA=cfgA[*,oA]
  A=create_struct('LINE_INT',fA,'LOGT',A.logT,'WVL',wA,'Z',zA,'ION',$
  	ionA,'DESIG',dsgA,'CONFIG',cfgA,'SRC',srcA,'JON',jonA)
endif
;reset A}

ww=[wA] & nwA=n_elements(wA)
if szb[0] ne 0 then begin
  wB=B.wvl & ww=[ww,wB] & nwB=n_elements(wB)
endif else nwB=0L
wC=ww[uniq(ww,sort(ww))] & nwC=n_elements(wC)
;if nwC gt 126 and float(!version.release) lt 5 then begin
;  c1='Too many fields; cannot concatenate'
;  message,c1,/info & return,A
;endif
ncomm=30 & if keyword_set(comm) then ncomm=fix(comm)>ncomm

;	initialize
k=0 & outw=ww
atom=1 & rom=1 & inicon,atom=atom,roman=rom
;atom=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al','Si',$
;	'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',$
;	'Cu','Zn']				;elements from 1-30
;rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
;rom=[rom,'X'+rom,'XX'+rom,'XXX'+rom]		;roman numerals from 1-40

;	merge
if nwB gt 0 then begin				;{there be sumpin' to merge
  ;	the output will include..
  fC=fA & wC=A.WVL & zC=A.Z & ionC=A.ION & jonC=A.JON
  dsgC=A.DESIG & cfgC=A.CONFIG & srcC=A.SRC
  ;	and possibly a lot of..
  fB=fB0 & wB=B.WVL & zB=B.Z & ionB=B.ION & jonB=B.JON
  dsgB=B.DESIG & cfgB=B.CONFIG & srcB=B.SRC

  ;	figure out whether DESIG/CONFIG are filled in either A or B
  ilvl=0 & iecf=0
  nlA=n_elements(A.DESIG) & nlB=n_elements(B.DESIG)
  neA=n_elements(A.CONFIG) & neB=n_elements(B.CONFIG)
  if nlA gt 1 or nlB gt 1 then ilvl=1
  if neA gt 1 or neB gt 1 then iecf=1
  ;
  if nlA eq 1 and ilvl eq 1 then dsgC=strarr(2,nwA)
  if neA eq 1 and iecf eq 1 then cfgC=strarr(2,nwA)

  ;	how to handle duplicates, if any?
  if keyword_set(ask) then begin
    c1=strlowcase(strmid(strtrim(ask,2),0,1))
    case c1 of
      'r': setask='r'
      'k': setask=''
      'n': setask='n'
      else: setask='d'
    endcase
  endif else setask=''

  ;	figure out where the duplicates lie
  oBC=lonarr(nwB)-1L
  for i=0L,nwB-1L do oBC[i]=(where((abs(wB[i])-abs(wA)) lt 1e-4 and zB[i] eq zA and ionB[i] eq ionA))[0]
  ok=where(oBC lt 0,mok) & odup=where(oBC ge 0,modup)

  ;	handle the duplicates here
  if setask eq 'd' then begin		;(delete the duplicates from B
    if mok gt 0 and modup gt 0 then begin
      fB=fB[*,ok]
      wB=wB[ok] & zB=zB[ok] & ionB=ionB[ok] & jonB=jonB[ok] & srcB=srcB[ok]
      if nlB gt 1 then dsgB=dsgB[2,ok]
      if neB gt 1 then cfgB=cfgB[2,ok]
    endif else fB=[-1.]
  endif					;SETASK='d')
  if setask eq 'r' then begin		;(delete from B, replace in A
    if mok gt 0 and modup gt 0 then begin
      fB=fB[*,ok]
      wB=wB[ok] & zB=zB[ok] & ionB=ionB[ok] & jonB=jonB[ok] & srcB=srcB[ok]
      if nlB gt 1 then dsgB=dsgB[2,ok]
      if neB gt 1 then cfgB=cfgB[2,ok]
      ;
      for i=0L,modup-1L do fC[*,oBC[odup[i]]]=fB[*,odup[i]]
      for i=0L,modup-1L do srcC[oBC[odup[i]]]=srcB[odup[i]]
    endif else fB=[-1.]
  endif					;SETASK='r')
  if setask eq 'n' then begin		;(delete from A
    if modup gt 0 then begin
      oCB=lonarr(nwA)-1L
      for i=0L,nwA-1L do oCB[i]=(where(i eq oBC))[0]
      oh=where(oCB lt 0,moh)
      if moh gt 0 then begin
        fC=fC[*,oh] & wC=wC[oh]
	zC=zC[oh] & ionC=ionC[oh] & jonC=jonC[oh] & srcC=srcC[oh]
        if ilvl eq 1 then dsgC=dsgC[2,oh]
	if iecf eq 1 then cfgC=cfgC[2,oh]
      endif else begin
	message,'All of A being deleted!',/info
	fC=fB & wC=wB & zC=zB & ionC=ionB & jonC=jonB & srcC=srcB
        if nlB eq 1 and ilvl eq 1 then dsgC=strarr(2,nwB) else dsgC=dsgB
        if neB eq 1 and iecf eq 1 then cfgC=strarr(2,nwB) else cfgC=cfgB
	fB=[-1.]
      endelse
    endif
  endif					;SETASK='n')

  ;	now append B to C
  ;if fB[0] gt -1 then begin
    fC=[fC[*],fB[*]]
    wC=[wC,wB] & zC=[zC,zB] & ionC=[ionC,ionB] & jonC=[jonC,jonB]
    srcC=[srcC,srcB]
    if nlB eq 1 and ilvl eq 1 then dsgB=strarr(2,nwB)
    if neB eq 1 and iecf eq 1 then cfgB=strarr(2,nwB)
    if ilvl eq 1 then dsgC=[dsgC[*],dsgB[*]]
    if iecf eq 1 then cfgC=[cfgC[*],cfgB[*]]
    ;	and now reformat..
    nwC=n_elements(wC)
    fC=reform(fC,ntA,nwC)
    if ilvl eq 1 then dsgC=reform(dsgC,2,nwC)
    if iecf eq 1 then cfgC=reform(cfgC,2,nwC)
  ;endif

  ;	form structure
  nwC=n_elements(wC)
  C=create_struct('LINE_INT',fC,'LOGT',A.logT,'WVL',wC,'Z',zC,'ION',$
  	ionC,'DESIG',dsgC,'CONFIG',cfgC,'SRC',srcC,'JON',jonC)

endif else begin				;}{just print, OK?

  C=A

  ;  determine fluxes?
  fC=C.LINE_INT & wC=C.WVL & zC=C.Z
  flx=dblarr(nwA) & tlog=fltarr(nwA)
  if ntA eq 1 then for iw=0,nwA-1 do flx[iw]=lineflx(fC[iw],A.LOGT,$
	wC[iw],zC[iw], _extra=e) else $
  	flx=lineflx(fC,A.LOGT,wC,zC, _extra=e)
  for i=0L,nwA-1L do begin
    tmp=max((A.line_int)[*,i],it) & tlog[i]=(A.logT)[it]
  endfor
  ;
  ;  reorder the output
  nwC=nwA & ok=lindgen(nwC)
  if keyword_set(reord) then begin
    ok=sort(flx) & if reord gt 0 then ok=reverse(ok)
    if abs(reord) gt 1 and reord lt nwC then begin
      rr=abs(reord) < nwA
      ok=ok[0:rr-1] & nwC=rr
    endif
  endif

  if not keyword_set(flst) then begin		;(stanadrd output
    c1=string('Z ION','(a10)')
    if keyword_set(okeV) then c1=c1+' '+string('[keV]','(a10)') else $
  	c1=c1+' '+string('WVL [Ang]','(a10)')
    c1=c1+'  '+string('FLUX   ','(a11)')+string('Tmax','(a6)') & print,c1
    for i=0,nwC-1 do begin
      j=ok[i]				;in case things have been reordered
      wvl=(A.wvl)[j] & z=(A.z)[j] & ion=(A.ion)[j]
      flux=flx[j] & tt=tlog[j]
      if n_elements(A.desig) gt 1 then labl=(A.desig)[*,j] else labl=['?','?']
      if labl[0] eq 'Unknown' then labl=['Un','known']
      print,strtrim(i,2)+$
	'---------------------------------------------'+strtrim(j+1,2)
      if flux ge 0 then begin
        if keyword_set(okeV) then wk=12.3985/(abs(wvl)>1e-3) else wk=abs(wvl)
        c1=string(atom[z-1]+' '+rom[ion-1],'(a10)')+string(wk,'(f10.4)')+$
	  ' '+string(flux,'(g11.5)')+' '+string(tt,'(f6.3)')
        c2='  '+labl[1]+' -> '+labl[0] & c2=strmid(c2,0,ncomm)
        if keyword_set(comm) then c1=c1+c2
      endif else c1='	_empty_'
      print,c1
    endfor
  endif else begin				;)(RD_LIST compatible output
    szl=size(flst) & nszl=n_elements(szl) & ulst=-1L
    if szl[nszl-2] eq 7 then openw,ulst,strtrim(flst[0],2),/get_lun
    for i=0,nwC-1 do begin
      j=ok[i]				;in case things have been reordered
      wvl=(A.wvl)[j] & Z=(A.Z)[j] & ion=(A.ion)[j]
      flux=flx[j] & tt=tlog[j]
      if z gt 0 then c1=string(atom[Z-1]) else c1='None'
      if ion gt 0 then c1=c1+' '+rom[ion-1]
      c1=c1+'	'+string(wvl,'(f11.5)')
      c1=c1+'	'+strtrim(flux,2)
      c1=c1+'	'+string(tt,'(f6.2)')
      ;
      ;this bit for compatibility reasons
      src=(A.src)[j] & szs=size(src) & nszs=n_elements(szs)
      if szs[nszs-2] eq 7 then c1=c1+'	'+src else begin
	if fix(src) eq 1 then c1=c1+'	$SPEX' else $
	 if fix(src) eq 2 then c1=c1+'	$CHIANTI' else $
	  if fix(src) eq 3 then c1=c1+'	$APED' else c1=c1+'	$CHIANTI'
      endelse
      ;
      if n_elements(A.desig) gt 1 then labl=(A.desig)[*,j] else labl=[' ',' ']
      if n_elements(A.config) gt 1 then labl=labl+' '+(A.config)[*,j]
      c1=c1+'	'+labl[0]+' -> '+labl[1]
      ;
      if Z gt 0 then printf,ulst,c1
      print,c1
    endfor
    if szl[nszl-2] eq 7 then begin & close,ulst,/all & free_lun,ulst & endif
  endelse					;FLST)

endelse						;nwB=0}

return,C
end
