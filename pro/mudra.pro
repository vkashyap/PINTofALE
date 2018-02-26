pro mudra,ststr,idstr,floc=floc,flab=flab,inclev=inclev,inclam=inclam,$
	filter=filter,oststr=oststr, _extra=e
;+
;procedure	mudra
;	imprint the contents on IDSTR onto the state structure describing
;	the plotting
;
;	"mudra" means "seal" or "imprint" (in the context of body paint --
;	your deity defines which mudra you wear)
;	pronounced "muh-dh-raw"
;
;parameters
;	ststr	[I/O; required] the state structure containing all
;		the necessary information on what labels to put where.
;		* see >>description<< in KALPANA.PRO
;	idstr	[INPUT] the output of LINEID, containing the wavelengths,
;		IDs, labels, etc. of identified features in the spectrum.
;
;keywords
;	floc	[INPUT] x-locations of features to add to STSTR
;	flab	[INPUT] labels for features to add to STSTR
;		* size of FLAB >must< match that of FLOC
;	inclev	[INPUT] if set, includes level designations and such in labels
;	inclam	[INPUT] if set, includes wavelength info
;	filter	[INPUT] if set to a long-integer array, deletes appropriate
;		labels from final STSTR, and if set to a scalar integer,
;		deletes everything *but* the appropriate integer
;	oststr	[OUTPUT] returns the original STSTR, in case modifications
;		are not acceptable
;	_extra	[JUNK] here only to prevent crashing the system
;
;history
;	vinay kashyap (SepMIM)
;-

;	usage
if n_params() lt 1 then begin
  print,'Usage: mudra,ststr,idstr,floc=floc,flab=flab,/inclev,/inclam,$'
  print,'       filter=filter,oststr=oststr'
  print,'  imprint ID info into plot state'
  return
endif

;	initialize
inicon,atom=atom,roman=roman

;	is state structure defined?
sstr=1		;not defined (0) or defined (1)
ns=n_elements(ststr) & ms=n_tags(ststr) & ok='ok'
if ns eq 0 then ok='State undefined.' else $
 if ms eq 0 then ok='State variable is junk.' else begin
   stnam=tag_names(ststr)
   if stnam(0) ne 'WINDOW' then ok='incorrect field.' else $
    if stnam(1) ne 'LOC' then ok='incorrect field.' else begin
      nl=n_elements(ststr.(1).(0))
      if ms ne nl+2L then ok='insufficient fields.'
    endelse
 endelse
if ok ne 'ok' then begin
  message,'Input state error: '+ok+'  Overwriting.',/info
  sstr=0
endif

;	if not defined, then define a skeleton
wstr=create_struct('MULTI',!P.multi,'TITLE','JUNK','SUBTITLE','',$
  'XTITLE','[Ang]','YTITLE','[Counts]','XSTYLE',1,'YSTYLE',1,$
  'XLOG',0,'YLOG',0,'CHARS',1.0,$
  'XMIN',0.0,'XMAX',1.0,'YMIN',0.0,'YMAX',0.0)
xstr=create_struct('X',[0.0],'Y',[0.0],'GROUP',[0L])
lstr=create_struct('POS',[0.0,0.0],'LABEL','','SIZE',1.0,$
  'XPATH',[0.0],'YPATH',[0.0],'ALIGN','LEFT','ORIENT',90.0,$
  'THICK',1.0,'LABCOLOR',100,'LINCOLOR',150,$
  'ARRANGE','ROW','UNDERLINE',1.0,'SIDELINE',0.0)
if sstr eq 0 then ststr=create_struct('WINDOW',wstr,'LOC',xstr,'L1',lstr)

;	verify input
istr=0		;ignore ID structure (0) or not (1)?
if n_params() eq 2 then begin
  ok='ok' & istr=1 & ni=n_tags(idstr)
  if ni gt 0 then begin
    idnam=tag_names(idstr)
    if idnam(0) ne 'WVL' then ok='ID structure in unknown format' else begin
      nw=n_elements(idstr.(0))
      if nw ne ni-1L then ok='ID structure is missing parts' else $
       if (idnam([1]))(0) eq 'WVL_COMMENT' then ok='no IDs worth the salt here'
    endelse
  endif else ok='ID structure not a structure'
  if ok ne 'ok' then begin
    message,ok+'. Ignoring.',/info & istr=0
  endif
endif

;	extract info from IDSTR
if istr eq 1 then begin
  locf=idstr.WVL & nx=n_elements(locf) & labf=strarr(nx)
  for i=1L,nx do begin			;{for each wavelength in IDSTR
    zz=idstr.(i).Z & jon=idstr.(i).ION & labl=idstr.(i).LABL
    nid=n_elements(zz) & cc=''
    ll=strtrim(string(abs(idstr.(i).WVL),'(f8.2)'),2)
    oo=where(idstr.(i).WVL lt 0,moo)
    if moo gt 0 then ll(oo)=ll(oo)+'*'	;theoretical wavelengths
    for j=0L,nid-1L do begin		;{for each ID

      if strlowcase(labl(0)) eq 'unknown' then begin
        cc='Unknown'
      endif else begin
	cc=cc+atom(zz(j)-1)+roman(jon(j)-1)	;ELEMENT IONIC_STATE
	if keyword_set(inclev) then begin	;(level designations
	  case inclev of
	    1: begin				;only level designations
		;for no reason whatsoever, I'm putting a "(" here
	      c1=strtrim((str_sep(labl(1,j),')'))(1),2) & lc1=strlen(c1)
	      cc=cc+' [!U'+strmid(c1,0,1)+'!N'+strmid(c1,1,1)+$
		 '!D'+strmid(c1,3,lc1-3)+'!N -> !U'
		;for no reason whatsoever, I'm putting a "(" here
	      c1=strtrim((str_sep(labl(0,j),')'))(1),2) & lc1=strlen(c1)
	      cc=cc+strmid(c1,0,1)+'!N'+strmid(c1,1,1)+'!D'+$
		 strmid(c1,3,lc1-3)+'!N]'
	    end
	    else: cc=cc+' ['+labl(1,j)+' -> '+labl(0,j)+']'	;everything
	  endcase
	endif					;INCLEV)
	if keyword_set(inclam) then begin		;(wavelengths
	  cc=cc+' !4k!3'+ll(j)
	endif					;INCLAM)
      endelse

      cc=cc+'|'					;field separator
    endfor				;J=0,NID-1}
    labf(i-1)=strmid(cc,0,strlen(cc)-1)		;remove the last "|"
  endfor				;I=1,NX}
endif

;	check for extra features to add
mx=n_elements(floc) & ml=n_elements(flab)
if mx gt 0 and mx eq ml then begin	;(worth a bother
  if n_elements(locf) gt 0 then locf=[locf,floc] else locf=floc
  if n_elements(labf) gt 0 then labf=[labf,flab] else labf=flab
endif					;MX>0 && MX==ML)
mx=n_elements(locf) & ml=n_elements(labf)

;	update ranges in STSTR.WINDOW.X(MIN|MAX)
if mx gt 0 then begin
  xmin=min(locf,max=xmax)
  ststr.WINDOW.XMIN = ststr.WINDOW.XMIN < xmin
  ststr.WINDOW.XMAX = ststr.WINDOW.XMAX > xmax
endif

oststr=ststr		;save old version
nlabel=0L
if mx gt 0 and sstr eq 0 then begin
  ;	start it off correctly
  ststr.(1).X=locf(0)
  ststr.(2).POS(0)=locf(0)
  ststr.(2).XPATH(0)=locf(0)
endif
if mx gt 0 and sstr ne 0 then begin	;appending to preexisting structure
  nlabel=n_elements(ststr.LOC.GROUP)
endif

;	break state structure into constituent pieces
stnam=tag_names(ststr)
wstr=ststr.(0)
xstr=ststr.(1) & xloc=xstr.(0) & yloc=xstr.(1) & iloc=xstr.(2)
nloc=n_elements(xloc)

;	update preexisting LOCs or append at end
k=0L	;appendage
for i=0L,mx-1L do begin			;{run through locations
  oo=where(xloc eq locf(i),moo)
  if moo gt 0 then begin		;(found a match
    ii=oo(0)
    ststr.(ii+2).LABEL=labf(i)	;update label
  endif else begin			;)(no matches, new label
    ;k=k+1L
    k=n_elements(iloc)
    lstr_i=lstr
    lstr_i.POS=[locf(i),0.0]
    lstr_i.LABEL=labf(i)
    lstr_i.XPATH=[locf(i)]
    lstr_i.YPATH=[0.0]
    xloc=[xloc,locf(i)] & yloc=[yloc,0.0] & iloc=[iloc,k]
    xstr=create_struct('X',xloc,'Y',yloc,'GROUP',iloc)
    tstr=create_struct('WINDOW',wstr,'LOC',xstr)
    for j=0L,k-1L do begin
      tstr=create_struct(tstr,'L'+strtrim(j+1L,2),ststr.(j+2L))
    endfor
    tstr=create_struct(tstr,'L'+strtrim(j+1L,2),lstr_i)
    ststr=tstr
  endelse				;MOO)
endfor					;I=0,NLOC-1}

;	now filter as needed
if n_elements(filter) gt 0 then begin		;(FILTER
  szf=size(filter) & nszf=n_elements(szf)
  if szf(nszf-2) eq 3 or szf(nszf-2) eq 2 then begin	;(bother?
    ;	return only FILTER(0) if scalar, and all but FILTER if vector
    exclude=szf(0) & nexcl=szf(nszf-1L)
    tstr=ststr & wstr=ststr.WINDOW & ststr=create_struct('WINDOW',wstr)
    nlabel=n_elements(tstr.LOC.X)
    if keyword_set(exclude) then begin		;(exclude FILTER
      ok=lindgen(nlabel)
      for i=0L,nexcl-1L do $
	if filter(i) ge 0 and filter(i) lt nlabel then ok(filter(i))=-1L
      ok=where(ok ge 0,mok)
    endif else begin				;)(include only FILTER(0)
      ok=[0L] & mok=1	;return the first element, by default
      if filter(0) ge 0 and filter(0) lt nlabel then ok=[filter(0)]
    endelse					;EXCLUDE)
    ;
    xloc=tstr.LOC.X & yloc=tstr.LOC.Y & grp=tstr.LOC.GROUP
    xloc=xloc(ok) & yloc=yloc(ok) & grp=grp(ok)
    oy=where(grp ge mok-1,moy) & if moy gt 0 then grp(oy)=oy ;reset lost groups
    lstr=create_struct('X',xloc,'Y',yloc,'GROUP',grp)
    ststr=create_struct(ststr,'LOC',lstr)
    ;
    for i=0L,mok-1L do begin
      j=ok(i)
      ststr=create_struct(ststr,'L'+strtrim(j+1L,2),tstr.(j+2))
    endfor
  endif							;FILTER is integer)
endif						;FILTER)

return
end
