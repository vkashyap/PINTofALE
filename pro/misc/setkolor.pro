pro setkolor,colors,icols,ocols,rgbfil=rgbfil,oldr=oldr,oldg=oldg,oldb=oldb,$
	quiet=quiet, _ref_extra=ex
;+
;procedure	setkolor
;	loads specified colors into color table in specified indices
;	(this used to be called "setcolor" but the name was changed to
;	avoid a clash with JHUAPL's procedure of the same name)
;
;syntax
;	setkolor,colors,icols,ocols,rgbfil=rgbfil,oldr=oldr,oldg=oldg,oldb=oldb,/quiet
;
;parameters
;	colors	[INPUT; required] string array containing the names
;		of the colors to be loaded
;	icols	[INPUT] integer array specifying into which index
;		each color goes to
;		* if size does not match COLORS, ignored
;		* default is to go from 0:255
;	ocols	[OUTPUT] byte array of size [3,N(COLORS)] containing
;		the hexadecimal values determined for each COLORS, as an
;		rgb triplet
;
;keywords
;	rgbfil	[INPUT] name of file containing (r,g,b) color values
;		and corresponding names in same format as:-
;		* the default, /usr/lib/X11/rgb.txt
;		* note: some new linux systems don't seem to have this file.
;		  those users can copy over this file from another system, or
;		  use the copy provided in !ARDB+'/rgb.txt'
;	oldr	[OUTPUT] the old R colors
;	oldg	[OUTPUT] the old G colors
;	oldb	[OUTPUT] the old B colors
;		* NOTE: to undo the setting, just do
;		  	TVLCT,OLDR,OLDG,OLDB
;	quiet	[INPUT] if set, does not echo spew stuff onto screen
;	_ref_extra	[JUNK] here only to prevent crashing the program
;
;example
;	plot,findgen(10)/10.,color=1 & setkolor,'green',1
;	oplot,sin(findgen(10.)),color=2 & setkolor,'#ff0000',2
;
;restrictions
;	* requires subroutine NARY2DEC
;	* RGBFIL bit works only on UNIX
;	* messes with the current color table
;	* feature, not a bug: if !D.N_COLORS < 256, and a SETKOLOR
;	  command is issued before a plot window is opened, a color
;	  table is loaded in which gets lopped off at !D.N_COLORS,
;	  and all subsequent plots appear dim.
;
;history
;	vinay kashyap (99Jul)
;	bug correction for ICOLS=!D.N_COLORS
;	STRMID must have 3 args for back-compatibility (VK FebMMI/A.Maggio)
;	changed name from SETCOLOR to SETKOLOR (VK FebMMI)
;	added optional output OCOLS and check for UNIX (VK; Jun02)
;	added checks to find rgb.txt in !ARDB if, as in some linux
;	  distributions, /usr/lib/X11/rgb.txt doesn't exist (VK; Mar07)
;	bug correction: local copy was trying to override RGBFIL (JJD; May08)
;	added switch to choose file_search() instead of findfile()
;	  for IDL 7+ (JJD/VK; Apr14)
;	added keyword QUIET (VK; Nov15)
;	if GDL, now does not use findfile (VK; Nov16)
;	added special colors 'CfA Red', 'CfA Violet', 'CfA Dark Blue' (VK; Nov18)
;-

forward_function findfile,file_search

;	usage
ncol=n_elements(colors)
if ncol eq 0 then begin
  print,'Usage: setkolor,colors,icols,ocols,rgbfil=rgbfil,$'
  print,'       oldr=oldr,oldg=oldg,oldb=oldb,/quiet'
  print,'  loads specified colors into color table'
  print,'  [see also: PEASECOLR]
  return
endif

;	check inputs
nic=n_elements(icols) & jcols=lindgen(ncol)+1
if not keyword_set(rgbfil) then rgbfil='/usr/lib/X11/rgb.txt'

;	should we use findfile or file_search?
idlversion=float(strmid(!version.release,0,1))
defsysv,'!GDL',exists=igdl

;	does RGBFIL exist?  if not, look in !ARDB
if (igdl eq 0 and idlversion le 6) then fil=findfile(rgbfil,count=nfil) else $
	fil=file_search(rgbfil,count=nfil)
if nfil eq 0 then begin
  zTOPDIR='/data/fubar/SCAR/'
  ivar=0 & defsysv,'!TOPDIR',exists=ivar	;if !TOPDIR exists
  if ivar ne 0 then setsysval,'TOPDIR',zTOPDIR,/getval else begin
    tmp=routine_info('SETKOLOR',/source) & scardir=tmp.PATH
    jvar=strpos(scardir,'pro/misc/setkolor.pro') ;we know where SETKOLOR is
    if jvar[0] lt 0 then $
     jvar=strpos(scardir,'pro/util/setkolor.pro') ;we DO know where it is too
    zTOPDIR=strmid(scardir,0,jvar-1)
  endelse
  rgbfil=filepath('rgb.txt',root_dir=zTOPDIR,subdirectory='ardb')
endif

;	initialize color table
tvlct,rr,gg,bb,/get & mcol=n_elements(rr)

;	optional outputs
oldr=rr & oldg=gg & oldb=bb
ocols=bytarr(3,ncol)

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
    if !version.OS_FAMILY eq 'unix' then begin		;(only on UNIX
      cmd='grep -i "'+cc+'" '+rgbfil
      spawn,cmd,cmtch & nmtch=n_elements(cmtch)
      if not keyword_set(cmtch) then nmtch=0
      case nmtch of
        0: begin
	  kcfa=''
	  ;	special cases
	  if strpos(strupcase(cc),'CFAR') ge 0 or strpos(strupcase(cc),'CFA R') ge 0 then kcfa='cfared'
	  if strpos(strupcase(cc),'CFAV') ge 0 or strpos(strupcase(cc),'CFA V') ge 0 then kcfa='cfaviolet'
	  if strpos(strupcase(cc),'CFAB') ge 0 or strpos(strupcase(cc),'CFAD') ge 0 or strpos(strupcase(cc),'CFA B') ge 0 or strpos(strupcase(cc),'CFA D') ge 0 then kcfa='cfadarkblue'
	  case kcfa of
	    'cfared': begin & rc(ii)=141 & gc(ii)=0 & bc(ii)=52 & end
	    'cfaviolet': begin & rc(ii)=43 & gc(ii)=53 & bc(ii)=117 & end
	    'cfadarkblue': begin & rc(ii)=19 & gc(i)=26 & bc(ii)=60 & end
	    else: message,'no matches found for color: '+cc,/info	;do nothing
	  endcase
	end
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
    endif else begin					;UNIX)(not UNIX
      message,'cannot decode color name on non-UNIX platform',/info
    endelse						;not UNIX)
  endelse				;decode color)
  if not keyword_set(quiet) then $
    message,'setting COL='+strtrim(ii,2)+' to ['+$
	strtrim(long(rc(ii)),2)+','+$
  	strtrim(long(gc(ii)),2)+','+$
	strtrim(long(bc(ii)),2)+']',/info
  ocols[0,i-1]=rc[ii] & ocols[1,i-1]=gc[ii] & ocols[2,i-1]=bc[ii]
  nextcol:	;the goto came here}
endfor					;I=1,NCOL}

;	load new colors
tvlct,rc,gc,bc

return
end
