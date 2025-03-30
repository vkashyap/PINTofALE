function chandra_status_bit,inp,filter,bitstr=bitstr,stabit=stabit,$
	verbose=verbose, _extra=e
;+
;function	chandra_status_bit
;	filters a Chandra event file or status bit array on each status bit
;	and reports how many events are flagged, and returns the indices of
;	the filtered events
;
;syntax
;	idx=chandra_status_bit(inp,filter,bitstr=bitstr,verbose=verbose)
;
;parameters
;	inp	[INPUT; required] the input from which the ststus bits are
;		to be filtered
;		* if byte array, must be of size (4,NPH), and this is the
;		  IDL read-in of a 32-bit string of 0s and 1s.  Note that
;		  the columns are ordered last to first.
;		* if a structure, there must be a field in it named .STATUS
;		  and the byte array of above case is taken from that
;		* if string, assumed to be the full path to a filename
;		  that contains an events list, and it is read in and
;		  the STATUS column is extracted from it
;	filter	[INPUT; optional] the filter to be used, in the same
;		format as CIAO.  e.g., HRC-I status-bit filtering uses
;		'xxxxxx00xxxx0xxx00000000x0000000', i.e., 0 or 1 matches
;		exactly, and 'x' ignores it.
;
;keywords
;	bitstr	[OUTPUT] a structure that contains the bits that are
;		turned on, the number of events that are flagged, and
;		for each of those, the array of the flagged indices
;	stabit	[OUTPUT] ordered array of status bits, converted from
;		IDL's 4xNPH to 32xNPH
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;requires
;	B11001001
;	IDL-Astro library
;
;history
;	vinay kashyap (2014jul, inspired by CIAO's summarize_status_bits)
;	fixed printout to conform to bit order of 32 to 0 from l. to r. (VLK 2025mar)
;-

;	usage
ok='ok' & np=n_params() & ni=n_elements(inp) & si=size(inp)
if np eq 0 then ok='Insufficient parameters' else $
 if ni eq 0 then ok='INPUT not defined'
if ok ne 'ok' then begin
  print,'Usage: idx=chandra_status_bit(inp,filter,bitstr=bitstr,verbose=verbose)'
  print,'  filters on status bits in Chandra event lists'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

;	parse input
ti=n_tags(inp) & si=size(inp) & nsi=n_elements(si)
if ti eq 0 and si[nsi-2] eq 7 then begin
  evt=mrdfits(inp[0],1,hdr)
  if not keyword_set(evt) then begin
    message,INP[0]+': file does not appear to contain an event list',/informational
    return,-1L
  endif
  tnam=tag_names(evt)
  ok=where(strpos(tnam,'STATUS') eq 0,mok)
  if mok eq 0 then begin
    message,INP[0]+': does not contain the STATUS columns',/informational
    return,-1L
  endif
  starr=evt.STATUS
  si=size(starr)
endif
if ti gt 0 then begin
  tnam=tag_names(inp)
  ok=where(strpos(tnam,'STATUS') eq 0,mok)
  if mok eq 0 then begin
    message,'INPUT does not contain STATUS column',/informational
    return,-1L
  endif
  starr=inp.STATUS
  si=size(starr)
endif
if si[nsi-2] eq 1 then begin
  if si[1] ne 4 then begin
    message,'STATUS is not 32 bits',/informational
    return,-1L
  endif
  if n_elements(starr) eq 0 then starr=inp
endif

;	parse filter
nf=n_elements(filter) & sf=size(filter,/type)
filt=strjoin(replicate('x',32))
if nf gt 0 then begin
  if nf gt 1 and vv gt 0 then message,$
  	'FILTER must be scalar; only first element will be used',/informational
  tmp=filter[0]
  if sf ne 7 then begin
    if vv gt 0 then print,tmp
    message,'FILTER must be a string array in CIAO format; ignoring this',/informational
  endif else begin
    lf=strlen(tmp)
    if lf ne 32 then begin
      message,tmp+": is not 32 bits in length; ignoring',/informational
    endif else begin
      ;there is room here to check that it only contains 0,1,x
      filt=tmp
    endelse
  endelse
endif

;	at this stage, there should be two variables defined
;	STARR -- byte 4xN
;	FILT -- string of length 32 with 0s, 1s, and xs

;	convert STARR to a 32xN byte array in same order as FILT
si=size(starr) & nevt=si[2] & stabit=bytarr(32,nevt)
i0=b11001001(starr[0,*],/otto) & for i=0,7 do stabit[0*8+i,*]=i0[7-i,*]
i1=b11001001(starr[1,*],/otto) & for i=0,7 do stabit[1*8+i,*]=i1[7-i,*]
i2=b11001001(starr[2,*],/otto) & for i=0,7 do stabit[2*8+i,*]=i2[7-i,*]
i3=b11001001(starr[3,*],/otto) & for i=0,7 do stabit[3*8+i,*]=i3[7-i,*]

;	now start filtering
idx=ulonarr(nevt) & bitstr=0
for i=0,31 do begin	;{step through each status bit
  cf=strmid(filt,i,1) & sarr=reform(stabit[i,*])
  o0=where(sarr eq 0,mo0,complement=o1,ncomplement=mo1)
  if mo0 ne nevt then $
    print,string(strtrim(32-(i+1),2),'(a2)')+' of 32: '+$
    	'0s - '+strtrim(mo0,2)+', 1s - '+strtrim(mo1,2)+' events'
  if mo0 ne nevt and arg_present(bitstr) then begin
    tmp=create_struct('flagged',mo1,'idx',o1,'zip',mo0,'ixnay',o0)
    if n_tags(bitstr) eq 0 then begin
      bitstr=create_struct('s'+string(i+1,'(i2.2)'),tmp)
    endif else begin
      bitstr=create_struct(bitstr,'s'+string(i+1,'(i2.2)'),tmp)
    endelse
  endif
  if cf ne 'x' then begin
    if cf eq '0' then begin
      if mo1 gt 0 then idx[o1]=1
    endif else begin
      if mo0 gt 0 then idx[o0]=1
    endelse
  endif
endfor			;I=0,31}

;	debug
if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,where(idx eq 0)
end
