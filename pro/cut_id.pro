function cut_id,idstr,idx,delet=delet,addon=addon,eps=eps,verbose=verbose,$
	_extra=e
;+
;function	cut_id
;	add or delete IDs to an existing ID structure
;
;syntax
;	newid=cut_id(idstr,idx,delet=delet,addon=addon,verbose=verbose,$
;	/incieq,dbdir=dbdir,sep=sep,pres=pres,logP=logP,n_e=n_e,$
;	chifil=chifil,chidir=chidir,eqfile=eqfile)
;
;warning
;	will not keep track of relative fluxes or their errors properly,
;	though the total flux will be conserved.
;
;parameters
;	idstr	[INPUT; required] an ID structure (see LINEID.PRO) containing
;		line identifications of spectral features
;	idx	[INPUT; required] a zero-based index of the feature
;		being (re)edited
;		* must be a scalar
;		* nothing happens if IDX is outside the legal range defined
;		  by IDSTR
;
;keywords
;	delet	[INPUT] set to a zero-based index to IDSTR.(IDX+1)
;		* if float, assumed to refer to IDSTR.(IDX+1).WVL
;		* if string, assumed to be in same format as understood by
;		  RD_LIST: "Z ION <sep> WAVE <sep> SOURCE <sep> DESCRIPTION"
;		  (RD_LIST will in fact be called in order to decipher it)
;		* may be an array
;	addon	[INPUT] string describing the line to be added as an ID
;		to IDSTR.(IDX+1)
;		* must be in same format as that understood by RD_LIST:
;		  "Z ION <sep> WAVE <sep> SOURCE <sep> DESCRIPTION"
;		* may be an array
;	eps	[INPUT] a small number
;		* default is 1e-5
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to
;		RD_LIST: INCIEQ,DBDIR,SEP,PRES,LOGP,N_E,CHIFIL,CHIDIR,EQFILE
;		LINEFLX: DEM,ABUND,NOPH,EFFAR,WVLAR
;
;restrictions
;	requires subroutines:
;	  RD_LIST, RD_LINE, FOLD_IONEQ, READ_IONEQ, RD_IONEQ, LINEFLX,
;	  GETABUND, WHEE, SYMB2ZION, LAT2ARAB, CAT_LN
;	requires IDL 5.3+ (because of use of STRSPLIT, STRJOIN, STRCMP)
;
;history
;	vinay kashyap (JanMMI)
;-

;	usage
ok='ok'
np=n_params() & n1=n_elements(idstr) & nid=n_tags(idstr) & n2=n_elements(idx)
if np lt 2 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='Input ID structure undefined' else $
  if n2 eq 0 then ok='Must specify which ID to edit' else $
   if nid eq 0 then ok='IDSTR: not a structure' else $
    if n2 gt 1 then ok='IDX: can only handle one feature at a time'
if ok ne 'ok' then begin
  print,'Usage: newid=cut_id(idstr,idx,delet=delet,addon=addon,verbose=v,$'
  print,'       /incieq,dbdir=dbdir,sep=sep,pres=pres,logP=logP,n_e=n_e,$'
  print,'       chifil=chifil,chidir=chidir,eqfile=eqfile)'
  if np gt 0 then message,ok,/info
  if nid eq 0 then return,-1L else return,idstr
endif

;	check that IDSTR contains useful stuff
ok='ok' & tnam=tag_names(idstr)
if tnam[0] ne 'WVL' then ok='ID structure not in correct format'
if tnam[0] eq 'WVL_COMMENT' then ok='ID structure contains no data'
if tnam[0] eq 'Z' then ok='ID structure has no wrapper'
if ok ne 'ok' then begin
  message,ok,/info
  return,-1L
endif

;	check that IDX is legal
wvl=IDSTR.WVL & nwvl=n_elements(wvl)
if idx[0] lt 0 or idx[0] ge nwvl then begin
  message,'Index to ID feature '+strtrim(idx[0],2)+$
	' outside valid range: [0,'+strtrim(nwvl-1)+']',/info
  if keyword_set(verbose) then message,'Returning with no changes',/info
  return,idstr
endif

;	other initializations
if not keyword_set(eps) then ee=1e-5 else ee=abs(eps)

;	now deal with IDSTR.(IDX+1)
idxtr=idstr.(idx[0]+1L) & ncomp=n_elements(idxtr.WVL)
tnam=tag_names(idxtr) & nflds=n_elements(tnam)
for i=0,nflds-1 do begin
  if tnam[i] eq 'WVL' then xWVL=idxtr.WVL
  if tnam[i] eq 'Z' then xZ=idxtr.Z
  if tnam[i] eq 'ION' then xION=idxtr.ION
  if tnam[i] eq 'LABL' then xLABL=idxtr.LABL
  if tnam[i] eq 'FLUX' then xFLUX=idxtr.FLUX
  if tnam[i] eq 'FLUXERR' then xFLUXERR=idxtr.FLUXERR
  if tnam[i] eq 'LOGT' then xLOGT=idxtr.LOGT
  if tnam[i] eq 'EMIS' then xEMIS=idxtr.EMIS
  if tnam[i] eq 'NOTES' then xNOTES=idxtr.NOTES
endfor
if n_elements(xLABL) eq 1 then begin
  xLABL=strarr(2,ncomp)+xLABL[0]
endif else begin
  if n_elements(xLABL) eq 2 and ncomp gt 2 then begin
    tmp=strarr(2,ncomp) & tmp[0,*]=xLABL[0] & tmp[1,*]=xLABL[1]
    xLABL=tmp
  endif
endelse
if n_elements(xNOTES) eq 0 then xNOTES=' '

;	figure out the fluxes and flux errors, in order to make sure that
;	total flux will stay the same.  be warned that relative fluxes will
;	NOT necessarily be correct
tflx=total(xFLUX) & tflxe=sqrt(total(xFLUXERR^2))

;	anything to delete?
if n_elements(delet) gt 0 then begin		;(DELET
  szd=size(delet) & nszd=n_elements(szd)
  if szd[nszd-2] eq 7 then begin	;(RD_LIST compatible description
    xNOTES[0]=xNOTES[0]+' CUT_ID: DELET='+strjoin(delet)
    lstr=rd_list(delet,verbose=verbose,/desig,/econf,brklst=brklst, _extra=e)
    dbdesc=strarr(n_elements(lstr.WVL))
    cd=lstr.DESIG & ce=lstr.CONFIG
    ;	unless DELET already contains a description, use the derived
    ;	e-configuration and level designation
    dbdesc=brklst.DESCRIPTION
    if strtrim(strjoin(dbdesc),2) eq '' then begin
      for j=0L,n_elements(lstr.WVL)-1L do dbdesc[j]=$
	cd[0,j]+' '+cd[1,j]+' '+ce[0,j]+' '+ce[1,j]
    endif
    for i=0,ncomp-1 do begin	;{check each component for a match
      oZ=where(xZ[i] eq lstr.Z,moZ)
      if moZ gt 0 then begin		;(Z matches
	oI=where(xION[i] eq lstr.ION,moI)
        if moI gt 0 then begin		;(Ion matches
          ow=where(abs(xWVL[i]-lstr.WVL) lt ee,mow)
          if mow gt 0 then begin		;(wvl matches
	    ;	now check for descriptions to match
	    cc=strsplit(xLABL[0,i]+' '+xLABL[1,i],/extract)
	    ncc=n_elements(cc)
	    for j=0L,ncc-1L do begin
	      ccc=strtrim(cc[j],2)
	      i0=strpos(ccc,'(',0) & i1=strpos(ccc,')',0)
	      if i0 ge 0 then strput,ccc,' ',i0
	      if i1 ge 0 then strput,ccc,' ',i1
	      ccc=strtrim(ccc,2) & ii=0
	      if ccc ne '' then ii=strmatch(dbdesc,'*'+ccc+'*',/fold_case)
	      if total(ii) gt 0 then begin
		if keyword_set(verbose) then begin
		  print,'Matching '+ccc+' to "',dbdesc,'"'
		  message,'Component '+strtrim(i,2)+' will be deleted',/info
		endif
		if n_elements(xidx) eq 0 then xidx=[i] else begin
		  oi=where(xidx eq i,moi) & if moi eq 0 then xidx=[xidx,i]
		endelse
	      endif
	    endfor
          endif					;wvl)
        endif				;Ion)
      endif				;Z)
    endfor			;I=0,NCOMP-1}
  endif else xidx=[long(delet)]			;indices)
  ;
  if n_elements(xidx) eq 0 then begin
    message,'specified criteria do not match any existing component; ignoring',/info
    xidx=[-1L]
  endif
  oo=where(xidx ge 0 and xidx lt ncomp,moo) & ox=lindgen(ncomp)
  if moo gt 0 then begin
    xidx=xidx[oo]
    ox[xidx]=-1L & oox=where(ox ge 0,moox)
    if moox gt 0 then begin
      xWVL=xWVL[ox[oox]] & xZ=xZ[ox[oox]] & xION=xION[ox[oox]]
      xLABL=xLABL[*,ox[oox]] & xFLUX=xFLUX[ox[oox]] & xFLUXERR=xFLUXERR[ox[oox]]
      xEMIS=xEMIS[*,ox[oox]]
    endif else begin
      message,'that would leave nothing. ignoring deletion',/info
    endelse
  endif
endif						;DELET)

;	anything to add on?
if keyword_set(addon) then begin
  lstr=rd_list(addon,verbose=verbose,/desig,/econf, _extra=e)
  if (lstr.LINE_INT)[0] gt -1 then begin
    mlin=n_elements(lstr.WVL)
    xWVL=[xWVL,lstr.WVL] & newnlin=n_elements(xWVL)
    xZ=[xZ,lstr.Z] & xION=[xION,lstr.ION]
    yLABL='('+lstr.CONFIG+') '+lstr.DESIG
    xLABL=reform([xLABL[*],yLABL[*]],2,newnlin)
    nT=n_elements(xLOGT) & mT=n_elements(lstr.LOGT) & yEMIS=lstr.LINE_INT
    if mT ne nT then begin
      if keyword_set(verbose) then message,'interpolating database T-grid'+$
	' to match existing grid',/info
      yEMIS=dblarr(nT,mlin)
      for i=0L,mlin-1L do yEMIS[*,i]=(interpol((lstr.LINE_INT)[*,i],$
	lstr.LOGT,xLOGT) > 0) < (max((lstr.LINE_INT)[*,i]))
    endif
    xEMIS=reform([xEMIS[*],yEMIS[*]],nT,newnlin)
    mFLUX=lineflx(lstr.LINE_INT,lstr.LOGT,lstr.WVL,lstr.Z, _extra=e)
    xFLUX=[xFLUX,mFLUX] & xFLUXERR=[xFLUXERR,fltarr(mlin)]
    xNOTES[0]=xNOTES[0]+' CUT_ID: ADDON='+strjoin(addon)
  endif else begin
    message,strjoin(addon,' ; ')+': No lines found: ',/info
  endelse
endif

;	renormalize the fluxes
tflx1=total(xFLUX) & tflx1e=sqrt(total(xFLUXERR^2))
if tflx1 ne 0 then xFLUX=xFLUX*tflx/tflx1
if tflx1e ne 0 then xFLUXERR=xFLUXERR*tflxe/tflx1e

;	reconstruct new ID structure
for i=0L,nwvl-1L do begin
  if i eq idx[0] then tmp=create_struct('WVL',xWVL,'Z',xZ,'ION',xION,$
    	'LABL',xLABL,'FLUX',xFLUX,'FLUXERR',xFLUXERR,$
	'LOGT',xLOGT,'EMIS',xEMIS,'NOTES',xNOTES) else tmp=idstr.(i+1)
  if i eq 0 then idwvl=create_struct('ID'+strtrim(i+1L,2),tmp) else $
	idwvl=create_struct(idwvl,'ID'+strtrim(i+1L,2),tmp)
endfor
newid=create_struct('WVL',wvl,idwvl)

return,newid
end
