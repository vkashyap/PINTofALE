function regroup,counts,minct,grid,newgrid=newgrid,$
	iout=iout,verbose=verbose, _extra=e
;+
;function	regroup
;	group an input array of counts such that each group has
;	at least a minimum number of counts.  start accumulating
;	from one end, and keep going until the threshold is
;	exceeded, then start on the next group, etc.
;
;syntax
;	gct=regroup(counts,minct,grid,/newgrid,iout=iout,verbose=vv)
;
;parameters
;	counts	[INPUT; required] counts array to be regrouped
;	minct	[INPUT; required] minimum total counts in each group
;	grid	[I/O] bin indices indicating starting and ending index
;		of each group, such that group K runs from index
;		GRID[K] to GRID[K+1]-1
;		* if given as an array that contains indices from 0
;		  to N_ELEMENTS(COUNTS) on input, uses that gridding
;		  and does not calculate a new grid, unless NEWGRID
;		  is set
;		* if not well-defined, or if NEWGRID is set, will be
;		  overwritten on output
;
;keywords
;	newgrid	[INPUT] if set, ignores the input GRID and computes
;		the grouping and the resulting new grid
;		* if not set, and GRID is well-defined, GRID is used
;		* if not set, and GRID is not well-defined, a new
;		  GRID is computed and the input is overwritten
;	iout	[OUTPUT] an integer index array designating which
;		group each element of COUNTS has ended up in
;		* this illustrates the essential non-destructive nature
;		  of the regrouping -- it should always be possible to
;		  go back to the original array and rearrange it.
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	y=randomu(seed,20) & gy=regroup(y,2.5,x,/newgrid,iout=iout)
;	print,y,iout & print,gy,x
;	for i=0,n_elements(gy)-1 do print,i,x[i],x[i+1]-1,$
;	  gy[i],total(y[x[i]:x[i+1]-1]),total(y[where(iout eq i)])
;
;these are all the smoothing tools in PINTofALE
;	ALALOESS() makes a loess curve
;	CLTSMOOTH() removes -ves from background subtracted spectra or light curves
;	CONV_RMF convolves with specified RMF
;	HAARTRAN() smoothes by filtering on Haar wavelet coefficients
;	LINEREM() removes lines from spectrum by iterative S/N filtering
;	NOISMOOTH() does boxcar accumulation a la Ebeling's asmooth
;	REGROUP() accumulates from one end
;	SCRMF does fast convolution using optimized RMF
;	SMOOTHIE does peak-descent accumulation to build up S/N
;	SPLAC computes a segmented piecewise linear approximations
;	UNKINK() removes sharp kinks from a curve
;	VARSMOOTH() does local smoothing at specified scales
;	VOORSMOOTH() does fast local smoothing at specified scales
;
;history
;	vinay kashyap (Oct09)
;-

;	usage
ok='ok' & np=n_params() & nc=n_elements(counts) & nm=n_elements(minct)
if np eq 0 then ok='Insufficient parameters' else $
 if nc eq 0 then ok='COUNTS are undefined' else $
  if nm eq 0 then ok='MINCT is undefined' else $
   if nc eq 1 then ok='COUNTS already in a group of one' else $
    if nm gt 1 then ok='cannot deal with more than one MINCT'
ming=0L & maxg=0L & ng=n_elements(grid) & if ng gt 0 then ming=min(grid,max=maxg)
if ok ne 'ok' then begin
  print,'Usage: gct=regroup(counts,minct,grid,/newgrid,iout=iout,verbose=vv)
  print,'  group an input counts array to have minimum accumulated'
  print,'  counts in each group'
  if np ne 0 then message,ok,/informational
  if nc eq 1 then begin
    if keyword_set(newgrid) or ng eq 0 then grid=[0L]
    return,counts
  endif else return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
iout=lindgen(nc)	;this is the default grid index output

;	is GRID well defined?
compute_grid='yes'
if not keyword_set(newgrid) then begin
  if ming eq 0 and maxg eq nc then compute_grid='no'
endif

;	preexisting GRID
if compute_grid eq 'no' then begin
  if ng eq 1 then ct=[total(counts)] else begin
    ct=lonarr(ng-1L)+0*counts[0]
    for i=0L,ng-2L do begin
      i0=grid[i] & i1=grid[i+1]-1L
      ct[i]=total(counts[i0:i1])
      iout[i0:i1]=i
    endfor
  endelse
  return,ct
endif

;	make new GRID
tct=total(counts,/cumul) & totct=tct[nc-1L]
if totct le minct[0] then begin
  if vv gt 0 then message,'all counts bundled into a single bin',/informational
  grid=[0L] & return,[totct]
endif
go_on=1 & subct=tct & ct=0L & grid=[0L]
while go_on do begin
  oo=where(subct ge minct,moo)
  if moo ne 0 then begin
    j=oo[0]
    if not keyword_set(ct) then ct=subct[j] else ct=[ct,subct[j]]
    grid=[grid,j+1L]
    subct=subct-subct[j]
  endif else begin
    go_on=0
    ct=[ct,subct[nc-1L]]	;last bin
    grid=[grid,nc]
  endelse
  if j eq nc-1L then go_on=0
endwhile
for i=0L,n_elements(ct)-1L do iout[grid[i]:grid[i+1L]-1L]=i

return,ct
end
