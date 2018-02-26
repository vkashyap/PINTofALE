function id2emis2id,idstr,ldbdir,dWVL=dWVL,verbose=verbose, _extra=e
;+
;function	id2emis2id
;	takes an ID structure, strips it into component parts,
;	reads in new emissivities, and puts them back in.
;
;syntax
;	newidstr=id2emis2id(idstr,ldbdir,dWVL=dWVL,verbose=verbose,$
;	eps=eps,/incieq,mapping=mapping,pres=pres,logP=logP,n_e=n_e,$
;	chifil=chifil,chidir=chidir,eqfile=eqfile)
;
;warning
;	new emissivities will not necessarily reflect the fluxes
;	as distributed among the components.  also, if the database
;	is changed, there is no reason to believe a priori that
;	the correct lines are really read back in.  we strongly
;	recommend using a high (say 10) verbosity level.
;
;parameters
;	idstr	[INPUT; required] ID structure, see LINEID for description
;	ldbdir	[INPUT] directory in which to look for line database
;		* default is "$CHIANTI"
;		* if array size matches the number of features, is
;		  distributed among components appropriately; if size
;		  matches that of total number of lines, maps one-to-one
;		  onto line list; if size is incompatible with IDSTR,
;		  then uses only what is in the first element.
;
;keywords
;	dWVL	[INPUT] slack in wavelength search in which to find the
;		matching line (this comes in handy when databases are
;		being changed around)
;		* default is 0.005
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to
;		RD_LIST: EPS,INCIEQ,MAPPING,PRES,LOGP,N_E,CHIFIL,CHIDIR,EQFILE
;
;restrictions
;	requires subroutines:
;	  RD_LIST, RD_LINE, FOLD_IONEQ, RD_IONEQ, READ_IONEQ,
;	  SYMB2ZION, ZION2SYMB, LAT2ARAB, CAT_LN, LINEFLX,
;	  INICON, WHEE, GETABUND, RDABUND, SYZE
;	requires IDL 5.3+
;
;history
;	vinay kashyap (AugMM)
;	catch if correct emissivities not found; added keywords
;	  VERBOSE, DWVL; allowed LDBDIR to be array (VK; JanMMI)
;-

;	usage
ok='ok'
np=n_params() & ni=n_elements(idstr) & nid=n_tags(idstr)
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='IDSTR is undefined' else $
  if nid eq 0 then ok='IDSTR is not a structure' else $
   if nid eq 1 then ok='IDSTR too small to be an ID structure' else begin
    idnam=tag_names(idstr)
    if idnam[0] ne 'WVL' then ok='IDSTR in unknown format'
    if idnam[1] eq 'WVL_COMMENT' then ok='IDSTR is empty'
  endelse
if float(strmid(!version.RELEASE,0,3)) lt 5.3 then ok=$
  'Requires IDL v5.3 or higher'
if ok ne 'ok' then begin
  print,'Usage: newidstr=id2emis2id(idstr,ldbdir,dWVL=dWVL,verbose=verbose,$'
  print,'       eps=eps,/incieq,mapping=mapping,pres=pres,logP=logP,n_e=n_e,$'
  print,'       chifil=chifil,chidir=chidir,eqfile=eqfile)'
  print,"  reread new emissivities corresponding to ID'd lines and store"
  print,"  inside new ID structure"
  message,ok,/info
  return,-1L
endif

;	explode IDSTR
obswvl=idstr.WVL		;wavelengths of the observed features
ncomp=n_elements(obswvl)	;number of features ID'd
nwvl=0L				;total number of IDs
for ic=0L,ncomp-1L do nwvl=nwvl+n_elements(idstr.(ic+1L).WVL)
idx=lonarr(nwvl)		;index pointing ID wavelength to ID component
wvl=fltarr(nwvl)		;wavelength of each ID
Z=intarr(nwvl)			;atomic numbers of IDs
ion=Z				;ionic states of IDs
labl=strarr(2,nwvl)		;level designation/e-configuration of ID
	;??? NEED IMPROVEMENT IN HANDLING ???
elem=strarr(nwvl)		;atomic symbol for ID
notes=elem			;notes, if any *** NOT IMPLEMENTED ***
flux=fltarr(nwvl)		;fluxes (see UPDATID/SQUISHEM) of ID
fluxerr=flux			;errors on FLUX
logT=idstr.(1).LOGT		;temperature grid
nlogT=n_elements(logT)		;LOGT *must* be identical for all
emis=dblarr(nlogT,nwvl)		;emissivities, usually including ion balance

;	now step through the structure and extract the contents into arrays
k=0L
for ic=0L,ncomp-1L do begin
  idtmp=idstr.(ic+1L) & mw=n_elements(idtmp.WVL) & slabl=size(idtmp.LABL)
  idx[k:k+mw-1L]=ic+1L
  wvl[k:k+mw-1L]=idtmp.WVL
  z[k:k+mw-1L]=idtmp.Z
  ion[k:k+mw-1L]=idtmp.ION
  flux[k:k+mw-1L]=idtmp.FLUX
  fluxerr[k:k+mw-1L]=idtmp.FLUXERR
  emis[*,k:k+mw-1L]=(idtmp.EMIS)
  if slabl[0] eq 1 then begin
    if slabl[1] eq 1 then labl[*,k:k+mw-1L]=(idtmp.LABL)[0]
    if slabl[1] eq mw then labl[0,k:k+mw-1L]=idtmp.LABL
    if mw eq 1 and slabl[1] eq 2 then begin
      labl[0,k:k+mw-1L]=(idtmp.LABL)[0]
      labl[1,k:k+mw-1L]=(idtmp.LABL)[1]
    endif
  endif
  if slabl[0] eq 2 then begin
    if slabl[2] eq mw then begin
      labl[0,k:k+mw-1L]=(idtmp.LABL)[0,*]
      labl[1,k:k+mw-1L]=(idtmp.LABL)[1,*]
    endif
  endif
  k=k+mw
endfor
zion2symb,z,ion,elem,ziform='Z ION'

;	check input
dbdir=strarr(nwvl)+'$CHIANTI'
if np eq 1 then ldbdir=dbdir
if np eq 2 then begin
  ndb=n_elements(ldbdir)
  if ndb eq 0 then begin
    message,'LDBDIR is undefined: using default -- '+dbdir[0],/info
  endif
  if ndb eq 1 then dbdir[*]=ldbdir[0]
  if ndb gt 1 then begin
    if ndb eq nwvl then begin
      ;	the number of DBDIRs match the number of lines in ID structure
      dbdir=ldbdir
      for i=0L,nwvl-1L do if keyword_set(ldbdir[i]) then dbdir[i]=ldbdir[i]
    endif else begin
      if ndb eq ncomp then begin
	;	the number of DBDIRs match the number of features in IDSTR
	for ic=0L,ncomp-1L do begin
	  oc=where(idx eq ic+1L,moc)
	  if keyword_set(ldbdir[ic]) then dbdir[oc]=ldbdir[ic]
	endfor
      endif else begin
	message,'LDBDIR is incompatible with IDSTR.',/info
	message,'Using only LDBDIR[0]='+strtrim(ldbdir[0],2),/info
	dbdir=ldbdir[0]
      endelse
    endelse
  endif
endif

dW=0.005 & if keyword_set(dWVL) then dW=abs(dWVL[0]) > 1e-10

;	re read the emissivities
linelist=elem+'|'+strtrim(wvl,2)+'+-'+strtrim(dW,2)+'|'+dbdir+'|'+$
	strtrim(labl[0,*],2)+' '+strtrim(labl[1,*],2)
lstr=rd_list(linelist,sep='|',/desig,/econf,verbose=verbose, _extra=e)

;	well, did it work as expected?
ok='ok'
for i=0L,n_elements(wvl)-1L do begin
  ow=where(abs(wvl[i]-lstr.WVL) lt dW,mow)
  if mow eq 0 then begin
    ok='line missing : '+linelist[i]
    if keyword_set(verbose) then print,ok
  endif
endfor
if n_elements(lstr.WVL) ne k then ok='returned emissivities do not match number of input lines'
if ok ne 'ok' then begin
  print,''
  message,ok,/info
  message,'consider setting the keyword MAPPING',/info
  message,'returning with no changes',/info
  print,''
  return,idstr
endif

;	overwrite the old emissivities
emis=lstr.LINE_INT
wvl=lstr.WVL
labl=lstr.DESIG+' '+lstr.CONFIG

;	push em back into IDstr
newidstr=idstr
for i=0L,ncomp-1L do begin
  ok=where(idx eq i+1,mok) & if mok eq 0 then message,'BUG!'
  newidstr.(i+1).WVL = wvl[ok]
  newidstr.(i+1).EMIS = emis[*,ok]
  newidstr.(i+1).LABL = labl[*,ok]
endfor

return,newidstr
end
