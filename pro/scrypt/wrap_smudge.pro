;+
;WRAP_SMUDGE
;	wrapper program for SMUDGE, also does rudimentary plotting
;
;vinay kashyap
;-

;	initialize
if not keyword_set(keys) then begin
  keys=''
  print,'keywords to be looked for must be in KEYS'
endif
if not keyword_set(only) then begin
  only=''
  print,'keywords to be combined must be in ONLY'
endif
if not keyword_set(nisi) then begin
  nisi=''
  print,'keywords to be excluded must be in NISI'
endif
if not keyword_set(range) then begin
  range=[-1,1e7]
  print,'range in which to search for features must be in RANGE'
endif
if not keyword_set(infile) then begin
  infile=1
  print,'name of file that contains the features must be in INFILE'
endif
if n_elements(sep) eq 0 then begin
  sep=0
  print,'fields in INFILE are separated by SEP'
endif
if not keyword_set(prefix) then begin
  prefix=0
  print,'extra prefixes to be treated as comments are in PREFIX'
endif
if not keyword_set(keV) then begin
  keV=0
  print,'if output must be in keV, set KEV=1'
endif

;	call SMUDGE
smudge,locate,labels,keys=keys,only=only,nisi=nisi,range=range,$
	infile=infile,sep=sep,prefix=prefix,keV=keV

help,locate,labels

;	overplot
nloc=n_elements(locate)
if nloc gt 0 then begin
  if nloc gt 1 then range=[min(locate),max(locate)] else $
	range=locate(0)*[0.9,1.1]
  if !d.window lt 0 then plot,[1],xrange=range,/xs,xlog=xlog,ylog=ylog
  ;
  if not keyword_set(orient) then orient=90
  if not keyword_set(align) then align=-0.5
  if not keyword_set(col) then col=fix(200.*!d.n_colors/256.)
  ;
  yr=!y.crange & ylocate=0.*locate+0.5*(yr(0)+yr(1))
  for i=0,nloc-1 do oplot,locate(i)*[1,1],[yr(0),ylocate(i)],color=col
  for i=0,nloc-1 do oplot,[locate(i)],[ylocate(i)],psym=1,color=col
  xyouts,locate,ylocate,labels,orient=orient,align=align,color=col,clip=0
endif

end
