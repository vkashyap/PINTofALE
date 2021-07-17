pro blindcolr,white=white,verbose=verbose,help=help,$
	retain=retain,oldr=oldr,oldg=oldg,oldb=oldb, _extra=e
;+
;procedure	blindcolr
;	set up a color table by loading in some nice colors
;
;syntax
;	blindcolr,/white,/brewer,verbose=verbose,/help,$
;	retain=retain,oldr=oldr,oldg=oldg,oldb=oldb, XYOUTS keywords
;
;parameters	NONE
;
;keywords
;	white	[INPUT] if set, chooses colors appropriate for a
;		white background, as is the case for a postscript
;		device
;		* mind you, setting this does NOT force the background
;		  color to be white.  it is just that colors appropriate
;		  for a white background are loaded in.
;		* if not set, will get automatically set if the device
;		  is postscript
;	brewer	[INPUT] if set, uses single-hue 9-class color sequences
;		from Color Brewer to fill in the 11-19, 21-29, 31-39,
;		51-59, 81-89, slots
;		(Cynthia Brewer et al., http://colorbrewer2.org)
;	verbose	[INPUT] controls chatter
;	help	[INPUT] if set, prints out the calling sequence
;		and notes which color is loaded into which index
;		* automatically set if VERBOSE > 50
;	retain	[INPUT] set to 0, 1, or 2 to specify how backing store
;		gets handled for the window (see RETAIN keyword for WINDOW
;		and DEVICE procedures):
;		0 = no backing store
;		1 = server or window system provides backing store
;		2 = IDL provides backing store directly
;		* default is 2
;		* WARNING: RETAIN is not remembered if plotting device changes
;	oldr	[OUTPUT] old R colors
;	oldg	[OUTPUT] old G colors
;	oldb	[OUTPUT] old B colors
;	_extra	[INPUT ONLY] pass defined keywords such as CHARTHICK,
;		CHARSIZE, and ALIGN to XYOUTS
;		* matters only if VERBOSE.ge.10, when a plot
;		  showing all the colors and names is displayed.
;
;description
;	this is based on PEASECOLR, which uses Deron Pease's scheme of
;	assigning specific colors to color table indices and more to the
;	point, being able to remember which color goes where.  Deron's
;	idea was to have 1=BLUE, 2=RED, 3=GREEN and then have 10's to be
;	different shades of blue, 20's to be shades of red, and 30's to be
;	shades of greens, etc.  Further, within each dex, the colors run as
;	X0=base, X1=dark, X9=light.  This basic scheme is enhanced by using
;	4=YELLOW, 5=PINK, 6=CYAN, 7=BROWN, 8=SAFFRON, 9=GREY.
;	For a white background, all shades are reversed so X1=light/X9=dark
;	and also (4) <-> (7) are swapped.
;
;	This is a modification of PEASECOLR for color-blind friendly plotting.
;	The suggested color palette in
;	https://www.nature.com/articles/nmeth.1618.pdf
;	is used to make the following base color assignments appropriate for a
;	white background:
;	X   colorname     	 r   g   b	hue    sat    val
;	1 = Blue:       	  0,114,178	191.46,1.0,  0.698
;	2 = Vermillion: 	213, 94,  0	 26.48,1.0,  0.835
;	3 = Bluish Green:	  0,158,115	163.67,1.0,  0.620
;	4 = Brown:      	179,109, 50	 27.44,0.721,0.702
;	5 = Reddish Purple:	204,121,167	326.75,0.407,0.800
;	6 = Sky Blue:   	 86,180,233	201.63,0.631,0.914
;	7 = Yellow:     	240,228, 66	 55.86,0.725,0.941
;	8 = Orange:     	230,159,  0	 41.48,1.0,  0.902
;	9 = Black:       	  0,  0,  0	  0.0, 0.0,  0.0
;
;	For a black background, each are rescaled such that X.val=1
;	also colors 4 and 7 are swapped
;
;	Each corresponding dex goes linearly from X0=X to
;	whiteX9.val=blackX.val and blackX9.val=whiteX.val
;
;	To obtain the complete set of color shades used, set the
;	keyword HELP in conjunction with VERBOSE>5
;	
;	color indices 100+ are not touched, but all else is
;	overwritten.
;
;subroutines
;	HSV2RGB()
;	HEXED()
;
;example
;	blindcolr,verbose=10
;	blindcolr,/help,charthick=3
;	!p.background=255 & blindcolr,verbose=10,/white,/help
;
;history
;	vinay kashyap (Jul21; based on peasecolr.pro)
;-

on_error,2

;	keywords
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
wbg=0 & if keyword_set(white) then wbg=1
if !d.NAME eq 'PS' then begin
  if vv gt 0 and (wbg eq 0 or n_elements(white) eq 0) then message,$
   'Assuming white background; set WHITE=0 explicitly for black background',/info
  if n_elements(white) ne 0 then wbg=1
endif
ih=0 & if keyword_set(help) then ih=1 & if vv gt 50 then ih=1
iret=2
if n_elements(retain) ne 0 then begin
  case fix(retain) of
    0: iret=0
    1: iret=1
    else: iret=2
  endcase
endif

;	set up the color table
hsvw=fltarr(3,9,10)	;[hue,sat,val],[X],[X0-X9]
chsvw=strarr(9,10)
X=1-1 & hsvw[0,X,*]=191.46 & hsvw[1,X,*]=1.0   & hsvw[2,X,*]=0.698 & chsvw[X,0]='BLUE'
X=2-1 & hsvw[0,X,*]= 26.48 & hsvw[1,X,*]=1.0   & hsvw[2,X,*]=0.835 & chsvw[X,0]='VERMILLION'
X=3-1 & hsvw[0,X,*]=163.67 & hsvw[1,X,*]=1.0   & hsvw[2,X,*]=0.620 & chsvw[X,0]='BLUISH GREEN'
X=4-1 & hsvw[0,X,*]= 27.44 & hsvw[1,X,*]=0.721 & hsvw[2,X,*]=0.702 & chsvw[X,0]='BROWN'
X=5-1 & hsvw[0,X,*]=326.75 & hsvw[1,X,*]=0.407 & hsvw[2,X,*]=0.800 & chsvw[X,0]='REDDISH PURPLE'
X=6-1 & hsvw[0,X,*]=201.63 & hsvw[1,X,*]=0.631 & hsvw[2,X,*]=0.914 & chsvw[X,0]='SKY BLUE'
X=7-1 & hsvw[0,X,*]= 55.86 & hsvw[1,X,*]=0.725 & hsvw[2,X,*]=0.941 & chsvw[X,0]='YELLOW'
X=8-1 & hsvw[0,X,*]= 41.48 & hsvw[1,X,*]=1.0   & hsvw[2,X,*]=0.902 & chsvw[X,0]='ORANGE'
X=9-1 & hsvw[0,X,*]=  0.00 & hsvw[1,X,*]=0.0   & hsvw[2,X,*]=0.000 & chsvw[X,0]='BLACK'
for k=1,9 do hsvw[2,*,k]=(1.-hsvw[2,*,0])*float(k)/9.+hsvw[2,*,0]
rXw=bytarr(9,10) & gXw=rXw & bXw=rXw
hsvb=hsvw & chsvb=chsvw & chsvb[8,0]='WHITE' & rXb=rXw & gXb=rXb & bXb=rXb
for X=0,8 do begin
  hsvb[2,X,*]=reverse(reform(hsvw[2,X,*]))
  for k=0,9 do begin
    zw=hsv2rgb(hsvw[0,X,k],hsvw[1,X,k],hsvw[2,X,k]) & rXw[X,k]=zw[0] & gXw[X,k]=zw[1] & bXw[X,k]=zw[2]
    if k gt 0 then chsvw[X,k]=hexed(zw)
    zb=hsv2rgb(hsvb[0,X,k],hsvb[1,X,k],hsvb[2,X,k]) & rXb[X,k]=zb[0] & gXb[X,k]=zb[1] & bXb[X,k]=zb[2]
    if k gt 0 then chsvb[X,k]=hexed(zb)
  endfor
endfor
tmp=reform(hsvb[*,3,*]) & hsvb[*,3,*]=hsvb[*,6,*] & hsvb[*,6,*]=tmp
tmp=reform(chsvb[3,*]) & chsvb[3,*]=chsvb[6,*] & chsvb[6,*]=tmp
tmp=reform(rXb[3,*]) & rXb[3,*]=rXb[6,*] & rXb[6,*]=tmp
tmp=reform(gXb[3,*]) & gXb[3,*]=gXb[6,*] & gXb[6,*]=tmp
tmp=reform(bXb[3,*]) & bXb[3,*]=bXb[6,*] & bXb[6,*]=tmp

;icol=[1,2,3,4,5,6,7,8,9,hci1,hci2,hci3,hci4,hci5,hci6,hci7,hci8,hci9]
icol=indgen(99)+1
if not keyword_set(wbg) then begin
  cols=[reform(chsvb[0:8,0]), reform(chsvb[0,*]), reform(chsvb[1,*]), reform(chsvb[2,*]), reform(chsvb[3,*]), reform(chsvb[4,*]), reform(chsvb[5,*]), reform(chsvb[6,*]), reform(chsvb[7,*]), reform(chsvb[8,*])]
  rr=[reform(rXb[0:8,0]), reform(rXb[0,*]), reform(rXb[1,*]), reform(rXb[2,*]), reform(rXb[3,*]), reform(rXb[4,*]), reform(rXb[5,*]), reform(rXb[6,*]), reform(rXb[7,*]), reform(rXb[8,*])]
  gg=[reform(gXb[0:8,0]), reform(gXb[0,*]), reform(gXb[1,*]), reform(gXb[2,*]), reform(gXb[3,*]), reform(gXb[4,*]), reform(gXb[5,*]), reform(gXb[6,*]), reform(gXb[7,*]), reform(gXb[8,*])]
  bb=[reform(bXb[0:8,0]), reform(bXb[0,*]), reform(bXb[1,*]), reform(bXb[2,*]), reform(bXb[3,*]), reform(bXb[4,*]), reform(bXb[5,*]), reform(bXb[6,*]), reform(bXb[7,*]), reform(bXb[8,*])]
endif else begin
  cols=[reform(chsvw[0:8,0]), reform(chsvw[0,*]), reform(chsvw[1,*]), reform(chsvw[2,*]), reform(chsvw[3,*]), reform(chsvw[4,*]), reform(chsvw[5,*]), reform(chsvw[6,*]), reform(chsvw[7,*]), reform(chsvw[8,*])]
  rr=[reform(rXw[0:8,0]), reform(rXw[0,*]), reform(rXw[1,*]), reform(rXw[2,*]), reform(rXw[3,*]), reform(rXw[4,*]), reform(rXw[5,*]), reform(rXw[6,*]), reform(rXw[7,*]), reform(rXw[8,*])]
  gg=[reform(gXw[0:8,0]), reform(gXw[0,*]), reform(gXw[1,*]), reform(gXw[2,*]), reform(gXw[3,*]), reform(gXw[4,*]), reform(gXw[5,*]), reform(gXw[6,*]), reform(gXw[7,*]), reform(gXw[8,*])]
  bb=[reform(bXw[0:8,0]), reform(bXw[0,*]), reform(bXw[1,*]), reform(bXw[2,*]), reform(bXw[3,*]), reform(bXw[4,*]), reform(bXw[5,*]), reform(bXw[6,*]), reform(bXw[7,*]), reform(bXw[8,*])]
endelse
ncol=n_elements(cols) & mcol=max(icol)

;	usage
if ih eq 1 then begin
  print,'Usage: blindcolr,/white,verbose=verbose,/help,$'
  print,'       oldr=oldr,oldg=oldg,oldb=oldb, XYOUTS keywords'
  print,'  set up a color table with some useful colors loaded in'
  if vv ge 5 then begin
    print,'INDEX','COLOR',form='(a5,a20)'
    bgcol='Black' & if keyword_set(wbg) then bgcol='White'
    print,0,bgcol+' (background)',form='(i3,2x,a20)'
    for i=0,ncol-1 do print,icol[i],cols[i],form='(i3,2x,a20)'
  endif else begin
    if keyword_set(wbg) then print,'White background' else print,'Black background'
    print,'BLUE:   1,10-19'
    print,'VERMILLION:    2,20-29'
    print,'BLUISH GREEN:  3,30-39'
    if keyword_set(wbg) then print,'BROWN:  4,40-49' else print,'YELLOW: 4,40-49'
    print,'REDDISH PURPLE:   5,50-59'
    print,'SKY BLUE:   6,60-69'
    if keyword_set(wbg) then print,'YELLOW: 7,70-79' else print,'BROWN:  7,70-79'
    print,'ORANGE: 8,80-89'
    print,'GREY: 9,90-99'
  endelse
endif

;	quit if colors cannot be loaded
if !D.NAME eq 'X' then begin
  if getenv('DISPLAY') eq '' then begin
    message,'DISPLAY environment variable is not set correctly; exiting',/info
    return
  endif
endif

;	pop up a window
;	(for 8-bit consoles, IDL does not assign the correct number of
;	!D.N_COLORS that are available until a window is actually created)
if vv gt 5 then begin
  if !D.NAME eq 'X' then begin
    if !D.WINDOW lt 0 then begin
      plot,[0] & wdelete
    endif
  endif
endif

;	warnings
dncol=!D.N_COLORS
if dncol gt 256 then begin
  if vv gt 2 then message,$
	'You have 24-color depth; set DECOMPOSED=0 in all calls to DEVICE',/informational
  device,decomposed=0,retain=iret
endif
if dncol lt mcol then message,'Maximum number of colors available: '+$
	strtrim(dncol,2),/info else dncol=mcol

;	get the existing color table
tvlct,r,g,b,/get
dnewcol=n_elements(r)
old_r=r & old_g=g & old_b=b	;and save

;	and overwrite
for i=0,ncol-1 do begin
  j=icol[i]
  if j lt dnewcol then begin
    r[j]=rr[i] & g[j]=gg[i] & b[j]=bb[i]
  endif
endfor

;	load new color table
tvlct,r,g,b

;	display
if vv ge 10 then begin
  ii=indgen(ncol) & xx=ii/20 & yy=ii mod 20
  cc=strtrim(icol,2)+':'+cols
  plot,[0,5],[0,21],/nodata,xstyle=4,ystyle=4,$
	position=[0.,0.,1.0,1.0],ymargin=[0,0],xmargin=[0,0]
  xyouts,xx+0.5,21-yy,cc,col=icol,charthick=1.5,charsize=1, _extra=e
  if vv gt 1000 then stop,'HALTing.  type .CON to continue'
endif

return
end
