pro setcolor,colors,icols,rgbfil=rgbfil,oldr=oldr,oldg=oldg,oldb=oldb,$
	_ref_extra=ex
;+
;procedure	setcolor
;	loads specified colors into color table in specified indices
;
;warning
;	this routine is being renamed SETKOLOR to avoid conflict with
;	JHUAPL routine of the same name.  this will disappear very soon,
;	and is hanging around only to aid in the transition.
;
;syntax
;	setcolor,colors,icols,rgbfil=rgbfil,oldr=oldr,oldg=oldg,oldb=oldb
;
;parameters
;	colors	[INPUT; required] string array containing the names
;		of the colors to be loaded
;	icols	[INPUT] integer array specifying into which index
;		each color goes to
;		* if size does not match COLORS, ignored
;		* default is to go from 0:255
;
;keywords
;	rgbfil	[INPUT] name of file containing (r,g,b) color values
;		and corresponding names in same format as:-
;		* the default, /usr/lib/X11/rgb.txt
;	oldr	[OUTPUT] the old R colors
;	oldg	[OUTPUT] the old G colors
;	oldb	[OUTPUT] the old B colors
;		* NOTE: to undo the setting, just do
;		  	TVLCT,OLDR,OLDG,OLDB
;	_ref_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	* RGBFIL bit works only on UNIX
;	* requires subroutine NARY2DEC
;
;history
;	vinay kashyap (99Jul)
;	bug correction for ICOLS=!D.N_COLORS
;	STRMID must have 3 args for back-compatibility (VK FebMMI/A.Maggio)
;	the bell tolls (VK; FebMMI)
;-

message,'OBSOLETE!  use SETKOLOR instead',/info

;	usage
ncol=n_elements(colors)
if ncol eq 0 then begin
  print,'Usage: setcolor,colors,icols,rgbfil=rgbfil,oldr=oldr,oldg=oldg,oldb=oldb'
  print,'  loads specified colors into color table'
  return
endif

;	check inputs
nic=n_elements(icols) & jcols=lindgen(ncol)+1
if not keyword_set(rgbfil) then rgbfil='/usr/lib/X11/rgb.txt'

;	initialize color table
tvlct,rr,gg,bb,/get & mcol=n_elements(rr)

;	optional outputs
oldr=rr & oldg=gg & oldb=bb

;	stupid user tricks
szc=size(colors) & nszc=n_elements(szc)
if szc(nszc-2) ne 7 then begin
  message,'returning: input colors must be strings',/info & return
endif
if nic gt 0 and nic ne ncol then message,$
	'mismatch: ignoring color indices',/info else $
	if nic gt 0 then jcols=icols
if ncol ge mcol then begin
  message,'too many colors; ignoring past '+strtrim(mcol,2),/info
  cols=colors(0:mcol-1L)
  jcols=jcols(0:mcol-1L)
  ncol=mcol
endif else cols=colors

;	translate each color into RGB components using RGBFIL if needed
rc=rr & gc=gg & bc=bb
for i=1L,ncol do begin			;{for each color
  cc=strlowcase(strtrim(cols(i-1L),2)) & ii=jcols(i-1L)
  if ii ge mcol then begin
    message,'Color '+strtrim(ii,2)+' currently unavailable (max:'+$
	strtrim(mcol-1L,2)+')',/info
    goto,nextcol		;{yeah, a goto
  endif
  c1=strmid(cc,0,1)
  if c1 eq '#' then begin		;(decode color from hex
    lcol=strlen(cc)-1
    case lcol of
      3: begin				;{of the form #rgb
	cr0=strmid(cc,1,1) & cg0=strmid(cc,2,1) & cb0=strmid(cc,3,1)
	nary2dec,cr0,r0,/hex & nary2dec,cg0,g0,/hex & nary2dec,cb0,b0,/hex
	rc(ii)=r0*17 & gc(ii)=g0*17 & bc(ii)=b0*17
      end				;#rgb}
      6: begin				;{of the form #rrggbb
	cr0=strmid(cc,1,2) & cg0=strmid(cc,3,2) & cb0=strmid(cc,5,2)
	nary2dec,cr0,r0,/hex & nary2dec,cg0,g0,/hex & nary2dec,cb0,b0,/hex
	rc(ii)=r0 & gc(ii)=g0 & bc(ii)=b0
      end				;#rrggbb}
      else: message,'cannot understand color: '+cc,/info	;do nothing
    endcase
  endif else begin			;)(decode color from RGBFIL
    cmd='grep -i "'+cc+'" '+rgbfil
    spawn,cmd,cmtch & nmtch=n_elements(cmtch)
    if not keyword_set(cmtch) then nmtch=0
    case nmtch of
      0: message,'no matches found for color: '+cc,/info	;do nothing
      else: begin
	k=0L & norm0=total(fix(byte(cc))) & dnorm=norm0
	for j=0L,nmtch-1L do begin
	  c1=strtrim(strlowcase(strmid(cmtch(j),11,strlen(cmtch(j))-11)),2)
	  norm=total(fix(byte(c1)))
	  if abs(norm-norm0) lt dnorm then begin
	    k=j & dnorm=abs(norm-norm0)
	  endif
	endfor
	message,'choosing color: '+cmtch(k),/info
	rc(ii)=long(strmid(cmtch(k),0,3))
	gc(ii)=long(strmid(cmtch(k),4,3))
	bc(ii)=long(strmid(cmtch(k),8,3))
      end
    endcase
  endelse				;decode color)
  message,'setting COL='+strtrim(ii,2)+' to ['+$
	strtrim(long(rc(ii)),2)+','+$
  	strtrim(long(gc(ii)),2)+','+$
	strtrim(long(bc(ii)),2)+']',/info
  nextcol:	;the goto came here}
endfor					;I=1,NCOL}

;	load new colors
tvlct,rc,gc,bc

return
end
