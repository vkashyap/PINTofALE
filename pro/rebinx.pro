function rebinx,array,inX,outX,Xindex=Xindex,absmin=absmin,absmax=absmax,$
	norng=norng,verbose=verbose, _extra=e
;+
;function	rebinx
;	returns the input array suitably reformatted (resampled, stretched,
;	constricted, what have you) to the new range.
;
;syntax
;	f=rebinx(array,inX,outX,Xindex=Xindex,/absmin,/absmax,/norng,$
;	verbose=verbose)
;
;parameters
;	array	[INPUT; required] array to be reformatted
;	inX	[INPUT; required] abscissa values defining the input grid
;	outX	[INPUT] abscissa values defining the output grid
;		* if not given, assumed to be FINDGEN(81)*0.05+4.
;
;keywords
;	Xindex	[INPUT] specifies which dimension of ARRAY corresponds to inX
;		* if not given, picks that dimension that first matches to the
;		  number of elements in inX
;		* e.g., if ARRAY[1000,81] and inX[81], then Xindex=1
;	absmin	[INPUT] if set, force oARRAY to be > ABSMIN
;		* default is min(ARRAY)
;	absmax	[INPUT] if set, force oARRAY to be < ABSMAX
;		* default is max(ARRAY)
;	norng	[INPUT] if set, ignores ABSMIN and ABSMAX
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	1. rebin DEM in range logT=5..7 in 11 steps to DEM in range
;	   logT=4..8 in 81 steps:
;		fullDEM=rebinx(DEM,findgen(11)*0.2+5,findgen(81)*0.05+4)
;	2. resample emissivity from full range to smaller sample:
;		subEMIS=rebinx(EMIS,logT,sublogT,Xindex=0)
;
;restrictions
;	* can only handle at most 2-D arrays
;
;history
;	vinay kashyap (Oct98)
;	bug correction with Xindex, added keyword VERBOSE, converted
;	  to IDL5 format (VK; Jul02)
;-

;	usage
ok='ok' & np=n_params()
sza=size(array) & nsza=n_elements(sza) & na=sza(nsza-1) & ia=sza(nsza-2)
niX=n_elements(inX) & noX=n_elements(outX)
if np lt 2 then ok='Insufficient parameters' else $
 if ia eq 0 then ok='Input array undefined' else $
  if niX eq 0 then ok='input gridding undefined'
if ok ne 'ok' then begin
  print,'Usage: outArray=rebinx(Array,inX,outX,Xindex=Xindex,/absmin,/absmax,$'
  print,'       /norng,verbose=verbose)'
  print,'  reformat/resample/etc Array[inX,?] to outArray[outX,?]
  print,'' & if np ne 0 then message,ok,/info
  return,[-1L]
endif

;	check inputs
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
;	defaults for missing input
if noX eq 0 then begin
  if vv gt 2 then message,'Using default logT grid for output grid',/info
  outX=findgen(81)*0.05+4.
endif
noX=n_elements(outX)
;
if not keyword_set(absmin) then absmin=min(array)
if not keyword_set(absmax) then absmax=max(array)

;	which of the columns in ARRAY correspond to inX?
if n_elements(Xindex) eq 0 then begin
  Xindex=0L
  ;for i=sza[0]-1,0,-1 do if sza[i] eq nix then Xindex=i
  for i=1,sza[0] do if sza[i] eq nix then Xindex=i-1
  	;this essentially picks the first column that matches size
endif else if n_elements(Xindex) gt 1 then message,$
	'only the first element of XINDEX is being used',/info

;	consistency checks
ok='ok'
if Xindex[0] lt 0 then ok='non-existent column?' else $
 if Xindex[0] ge sza[0] then ok='non-existent column?' else $
  if na lt 2 then ok='Array size too small' else $
   if sza[Xindex[0]+1] ne niX then ok='ARRAY v/s inX mismatch'
if ok ne 'ok' then begin
  message,ok,/info & return,array
endif

;	define the output
szo=sza & szo[Xindex[0]+1]=noX & oArray=make_array(size=szo)

;	interpolate, stretch, etc.
case sza[0] of				;{handle each dimension separately
  1: begin
    oX=sort(inX) & oArray=Array(oX) & inY=inX[oX]
    oArray=interpol(oArray,inY,outX)
    if not keyword_set(norng) then oArray = (oArray > absmin) < absmax
  end
  2: begin
    n2=sza[1] & if Xindex[0] eq 0 then n2=sza[2]
    oX=sort(inX) & inY=inX[oX]
    for i=0L,n2-1L do begin
      if Xindex[0] eq 0 then tmp=reform(Array[*,i]) else tmp=reform(Array[i,*])
      tmp=tmp[oX]
      if vv gt 50 then plot,inY,tmp,title=i
      tmp=interpol(tmp,inY,outX)
      if vv gt 50 then oplot,outX,tmp,psym=-1 & if vv gt 50 then wait,vv/1000.
      if Xindex[0] eq 0 then oArray[*,i]=tmp[*] else oArray[i,*]=tmp[*]
    endfor
    if not keyword_set(norng) then oArray = (oArray > absmin) < absmax
  end
  else: begin
    message,'cannot handle arrays of dimension > 2',/info & return,Array
  end
endcase					;DIMENSION(ARRAY)}

return,oArray
end
