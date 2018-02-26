pro peasecolr,white=white,brewer=brewer,verbose=verbose,help=help,$
	retain=retain,oldr=oldr,oldg=oldg,oldb=oldb, _extra=e
;+
;procedure	peasecolr
;	set up a color table by loading in some nice colors
;
;syntax
;	peasecolr,/white,/brewer,verbose=verbose,/help,$
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
;	this is based on Deron Pease's scheme of assigning specific
;	colors to color table indices and more to the point, being
;	able to remember which color goes where.  Deron's idea is
;	to have 1=BLUE, 2=RED, 3=GREEN and then have 10's to be
;	different shades of blue, 20's to be shades of red, and 30's
;	to be shades of greens, etc.  Further, within each dex, the
;	colors run as: X0=base, X1=dark, X9=light.
;	
;	This basic scheme has been enhanced as follows:
;	-- added 4=YELLOW, 5=PINK, 6=CYAN, 7=BROWN, 8=SAFFRON, 9=GREY
;	-- the shades are optimised for different backgrounds
;	-- if WHITE is set, then
;	   X1=light and X9=dark, and
;	   (4) <-> (7)
;	To obtain the complete set of color shades used, set the
;	keyword HELP in conjunction with VERBOSE>5
;	
;	color indices 100+ are not touched, but all else is
;	overwritten.
;
;example
;	peasecolr,verbose=10
;	peasecolr,/help,charthick=3
;	peasecolr,verbose=10,/white,/help & setkolor,'white',0
;
;history
;	vinay kashyap (Jul02)
;	returns to caller on display error (VK; Sep02)
;	!D.N_COLORS is useless on VNC? (VK; May04)
;	added saffron to 8X (VK; May05)
;	added grey to 9X (VK; Aug12)
;	added keyword BREWER (VK; Sep15)
;	added keyword RETAIN (VK; Dec15)
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

;	The default color shades used (overridden if BREWER is set)
;BLUE:
;black: MediumBlue,RoyalBlue,DodgerBlue,CornflowerBlue,DeepSkyBlue,SkyBlue,LightSkyBlue,PowderBlue,LightBlue
;white: DeepSkyBlue,LightSlateBlue,MediumSlateBlue,RoyalBlue,MediumBlue,DarkSlateBlue,MidnightBlue,DarkBlue,NavyBlue
;RED:
;black: red4,red3,red2,OrangeRed2,OrangeRed1,tomato2,tomato1,coral2,salmon1
;white: salmon3,coral3,tomato2,tomato3,OrangeRed2,OrangeRed3,red2,red3,red4
;GREEN:
;black: DarkGreen,DarkOliveGreen,green4,SeaGreen,MediumSeaGreen,SpringGreen,green2,LawnGreen,chartreuse
;white: chartreuse,LawnGreen,green2,SpringGreen,MediumSeaGreen,SeaGreen,green4,DarkOliveGreen,DarkGreen
;YELLOW:
;black: orange,yellow2,gold,LightGoldenrod,NavajoWhite,PaleGoldenrod,LightYellow2,beige,ivory
;white: bisque3,tan,peru,IndianRed,DarkGoldenrod,chocolate,sienna,SaddleBrown,firebrick (from BROWN)
;PINK:
;black: maroon,DeepPink,VioletRed,HotPink,orchid,violet,plum,LightPink1,LavenderBlush1
;white: RosyBrown,LightPink2,plum,violet,orchid,HotPink,VioletRed,DeepPink,maroon
;CYAN: 
;black: BlueViolet,MediumPurple1,cyan2,CadetBlue1,turquoise,PaleTurquoise3,MediumAquamarine,aquamarine,PaleTurquoise1
;white: Turquoise,Turquoise2,aquamarine2,MediumAquamarine,PaleTurquoise3,LightSeaGreen,CadetBlue2,MediumPurple1,BlueViolet
;BROWN:
;black: firebrick,SaddleBrown,sienna,chocolate,DarkGoldenrod,IndianRed,peru,tan,bisque1
;white: DarkGoldenrod,LightYellow4,gold3,wheat3,LightGoldenrod3,gold1,yellow3,khaki2,grey70 (from YELLOW)
;SAFFRON:
;black: #ee5511, #dd7711, #e07511, #e27111, #e46711, #e66311, #e85911, #ea5511, #ec5511, #ee5511
;white: #ee5511, #dd7711, #e07511, #e27111, #e46711, #e66311, #e85911, #ea5511, #ec5511, #ee5511
;GREY:
;black: #111111, #333333, #555555, #777777, #999999, #aaaaaa, #cccccc, #eeeeee, #f7f7f7
;white: #f7f7f7, #eeeeee, #cccccc, #aaaaaa, #999999, #777777, #555555, #333333, #111111

;	these are the colors used (overridden if BREWER is set)
hcn1b=['BLUE','MediumBlue','RoyalBlue','DodgerBlue','CornflowerBlue',$
	'DeepSkyBlue','SkyBlue','LightSkyBlue','PowderBlue','LightBlue']
hcn1w=['BLUE','DeepSkyBlue','LightSlateBlue','MediumSlateBlue','RoyalBlue',$
	'MediumBlue','DarkSlateBlue','MidnightBlue','DarkBlue','NavyBlue']
hcn2b=['RED','red4','red3','red2','OrangeRed2','OrangeRed1','tomato2',$
	'tomato1','coral2','salmon1']
hcn2w=['RED','salmon3','coral3','tomato2','tomato3','OrangeRed2',$
	'OrangeRed3','red2','red3','red4']
hcn3b=['GREEN','DarkGreen','DarkOliveGreen','green4','SeaGreen','MediumSeaGreen',$
	'SpringGreen','green2','LawnGreen','chartreuse']
hcn3w=['GREEN','chartreuse','LawnGreen','green2','SpringGreen','MediumSeaGreen',$
	'SeaGreen','green4','DarkOliveGreen','DarkGreen']
hcn4b=['YELLOW','orange','yellow2','gold','LightGoldenrod','NavajoWhite',$
	'PaleGoldenrod','LightYellow2','beige','ivory']
hcn4w=['BROWN','bisque3','tan','peru','IndianRed','DarkGoldenrod','chocolate',$
	'sienna','SaddleBrown','firebrick']
hcn5b=['PINK','maroon','DeepPink','VioletRed','HotPink','orchid','violet',$
	'plum','LightPink1','LavenderBlush1']
hcn5w=['PINK','RosyBrown','LightPink2','plum','violet','orchid','HotPink',$
	'VioletRed','DeepPink','maroon']
hcn6b=['CYAN','BlueViolet','MediumPurple1','cyan2','CadetBlue1','turquoise',$
	'PaleTurquoise3','MediumAquamarine','aquamarine','PaleTurquoise1']
hcn6w=['CYAN','Turquoise','Turquoise2','aquamarine2','MediumAquamarine',$
	'PaleTurquoise3','LightSeaGreen','CadetBlue2','MediumPurple1','BlueViolet']
hcn7b=['BROWN','firebrick','SaddleBrown','sienna','chocolate','DarkGoldenrod',$
	'IndianRed','peru','tan','bisque1']
hcn7w=['YELLOW','DarkGoldenrod','LightYellow4','gold3','wheat3','LightGoldenrod3',$
	'gold1','yellow3','khaki2','grey70']
hcn8b=['SAFFRON','#dd7711','#e07511','#e27111','#e46711','#e66311','#e85911','#ea5511','#ec5511','#ee5511']
hcn8w=hcn8b
hcn9b=['gray', '#111111', '#333333', '#555555', '#777777', '#999999', '#aaaaaa', '#cccccc', '#eeeeee']
hcn9w=['gray', '#eeeeee', '#cccccc', '#aaaaaa', '#999999', '#777777', '#555555', '#333333', '#111111']

;	at these indices
hci1=indgen(10)+10
hci2=indgen(10)+20
hci3=indgen(10)+30
hci4=indgen(10)+40
hci5=indgen(10)+50
hci6=indgen(10)+60
hci7=indgen(10)+70
hci8=indgen(10)+80
hci9=indgen(10)+90

;	with these rgb values (overridden if BREWER is set)
r1b=byte([  0,  0, 65, 30,100,  0,135,135,176,173]) & r1w=byte([  0,  0,132,123, 65,  0, 72, 25,  0,  0])
g1b=byte([  0,  0,105,144,149,191,206,206,224,216]) & g1w=byte([  0,191,112,104,105,  0, 61, 25,  0,  0])
b1b=byte([255,205,225,255,237,255,235,250,230,230]) & b1w=byte([255,255,255,238,225,205,139,112,139,128])
r2b=byte([255,139,205,238,238,255,238,255,238,255]) & r2w=byte([255,205,205,238,205,238,205,238,205,139])
g2b=byte([  0,  0,  0,  0, 64, 69, 92, 99,106,140]) & g2w=byte([  0,112, 91, 92, 79, 64, 55,  0,  0,  0])
b2b=byte([  0,  0,  0,  0,  0,  0, 66, 71, 80,105]) & b2w=byte([  0, 84, 69, 66, 57,  0,  0,  0,  0,  0])
r3b=byte([  0,  0, 85,  0, 46, 60,  0,  0,124,127]) & r3w=byte([  0,127,124,  0,  0, 60, 46,  0, 85,  0])
g3b=byte([255,100,107,139,139,179,255,238,252,255]) & g3w=byte([255,255,252,238,255,179,139,139,107,100])
b3b=byte([  0,  0, 47,  0, 87,113,127,  0,  0,  0]) & b3w=byte([  0,  0,  0,  0,127,113, 87,  0, 47,  0])
r4b=byte([255,255,238,255,238,255,238,238,245,255]) & r4w=byte([165,205,210,205,205,184,210,160,139,178])
g4b=byte([255,165,238,215,221,222,232,238,245,255]) & g4w=byte([ 42,183,180,133, 92,134,105, 82, 69, 34])
b4b=byte([  0,  0,  0,  0,130,173,170,209,220,240]) & b4w=byte([ 42,158,140, 63, 92, 11, 30, 45, 19, 34])
r5b=byte([255,176,255,208,255,218,238,221,255,255]) & r5w=byte([255,188,238,221,238,218,255,208,255,176])
g5b=byte([192, 48, 20, 32,105,112,130,160,174,240]) & g5w=byte([192,143,162,160,130,112,105, 32, 20, 48])
b5b=byte([203, 96,147,144,180,214,238,221,185,245]) & b5w=byte([203,143,173,221,238,214,180,144,147, 96])
r6b=byte([  0,138,171,  0,152, 64,150,102,127,187]) & r6w=byte([  0, 64,  0,118,102,150, 32,142,171,138])
g6b=byte([255, 43,130,238,245,224,205,205,255,255]) & g6w=byte([255,224,229,238,205,205,178,229,130, 43])
b6b=byte([255,226,255,238,255,208,205,170,212,255]) & b6w=byte([255,208,238,198,170,205,170,238,255,226])
r7b=byte([165,178,139,160,210,184,205,205,210,255]) & r7w=byte([255,184,139,205,205,205,255,205,238,179])
g7b=byte([ 42, 34, 69, 82,105,134, 92,133,180,228]) & g7w=byte([255,134,139,173,186,190,215,205,230,179])
b7b=byte([ 42, 34, 19, 45, 30, 11, 92, 63,140,196]) & b7w=byte([  0, 11,122,  0,150,112,  0,  0,133,179])
r8b=byte([238,221,224,226,228,230,232,234,236,238]) & r8w=r8b
g8b=byte([ 85,119,117,113,103, 99, 94, 89, 85, 85]) & g8w=g8b
b8b=byte([ 17, 17, 17, 17, 17, 17, 17, 17, 17, 17]) & b8w=b8b
r9b=byte([ 128, 17, 51, 85,119,153,170,204,238,240]) & r9w=[128, reverse(r9b[1:*])]
g9b=byte([ 128, 17, 51, 85,119,153,170,204,238,240]) & g9w=[128, reverse(g9b[1:*])]
b9b=byte([ 128, 17, 51, 85,119,153,170,204,238,240]) & b9w=[128, reverse(b9b[1:*])]

;	alternate if BREWER is set
if keyword_set(brewer) then begin
  r1b=byte([  0,  8,  8, 33, 66,107,158,198,222,247]) & r1w=[r1b[0],reverse(r1b[1:*])]
  g1b=byte([  0, 48, 81,113,146,174,202,219,235,251]) & g1w=[g1b[0],reverse(g1b[1:*])]
  b1b=byte([255,107,156,181,198,214,225,239,247,255]) & b1w=[b1b[0],reverse(b1b[1:*])]
  hcn1b=['BLUE','#08306b','#08519c','#2171b5','#4292c6','#6baed6','#92cae1','#c6dbef','#deebf7','#f7fbff'] & hcn1w=[hcn1b[0],reverse(hcn1b[1:*])]

  r2b=byte([255,103,165,203,239,251,252,252,254,255]) & r2w=[r2b[0],reverse(r2b[1:*])]
  g2b=byte([  0,  0, 15, 24, 59,106,146,187,224,245]) & g2w=[g2b[0],reverse(g2b[1:*])]
  b2b=byte([  0, 13, 21, 29, 44, 74,114,161,210,240]) & b2w=[b2b[0],reverse(b2b[1:*])]
  hcn2b=['RED','#67000d','#a50f15','#cb181d','#ef3b2c','#fb6a4a','#fc9272','#fcbba1','#fee0d2','#fff5f0'] & hcn2w=[hcn2b[0],reverse(hcn2b[1:*])]

  r3b=byte([  0,  0,  0, 35, 65,116,161,199,229,247]) & r3w=[r3b[0],reverse(r3b[1:*])]
  g3b=byte([255, 68,109,139,171,196,217,233,245,252]) & g3w=[g3b[0],reverse(g3b[1:*])]
  b3b=byte([  0, 27, 44, 69, 93,118,155,192,224,245]) & b3w=[b3b[0],reverse(b3b[1:*])]
  hcn3b=['GREEN','#00441b','#006d2c','#238b45','#41ab5d','#74c476','#a1d99b','#c7e9c0','#e5f5e0','#f7fcf5'] & hcn3w=[hcn3b[0],reverse(hcn3b[1:*])]

  r5b=byte([255, 63, 84,106,128,158,188,218,239,252]) & r5w=[r5b[0],reverse(r5b[1:*])]
  g5b=byte([192,  0, 39, 81,125,154,189,218,237,251]) & g5w=[g5b[0],reverse(g5b[1:*])]
  b5b=byte([203,125,143,163,186,200,220,235,245,253]) & b5w=[b5b[0],reverse(b5b[1:*])]
  hcn5b=['PINK','#3f007d','#54278f','#6a51a3','#807dba','#9e9ac8','#bcbddc','#dadaeb','#efedf5','#fcfbfd'] & hcn5w=[hcn5b[0],reverse(hcn5b[1:*])]

  r8b=byte([238,127,166,217,241,253,253,253,254,255]) & r8w=[r8b[0],reverse(r8b[1:*])]
  g8b=byte([ 85, 39, 54, 72,105,141,174,208,230,245]) & g8w=[g8b[0],reverse(g8b[1:*])]
  b8b=byte([ 17,  4,  3,  1, 19, 60,107,162,206,235]) & b8w=[b8b[0],reverse(b8b[1:*])]
  hcn8b=['SAFFRON','#7f2704','#a63603','#d94801','#f16913','#fd8d3c','#fdae6b','#fdd0a2','#fee6ce','#fff5eb'] & hcn8w=[hcn8b[0],reverse(hcn8b[1:*])]
endif

icol=[1,2,3,4,5,6,7,8,9,hci1,hci2,hci3,hci4,hci5,hci6,hci7,hci8,hci9]
if not keyword_set(wbg) then begin
  cols=['BLUE','RED','GREEN','YELLOW','PINK','CYAN','BROWN','SAFFRON','GREY',$
  	hcn1b,hcn2b,hcn3b,hcn4b,hcn5b,hcn6b,hcn7b,hcn8b,hcn9b]
  rr=[r1b[0],r2b[0],r3b[0],r4b[0],r5b[0],r6b[0],r7b[0],r8b[0],128b,$
  	r1b,r2b,r3b,r4b,r5b,r6b,r7b,r8b,r9b]
  gg=[g1b[0],g2b[0],g3b[0],g4b[0],g5b[0],g6b[0],g7b[0],g8b[0],128b,$
  	g1b,g2b,g3b,g4b,g5b,g6b,g7b,g8b,g9b]
  bb=[b1b[0],b2b[0],b3b[0],b4b[0],b5b[0],b6b[0],b7b[0],b8b[0],128b,$
  	b1b,b2b,b3b,b4b,b5b,b6b,b7b,b8b,b9b]
endif else begin
  cols=['BLUE','RED','GREEN','BROWN','PINK','CYAN','YELLOW','SAFFRON','GREY',$
  	hcn1w,hcn2w,hcn3w,hcn4w,hcn5w,hcn6w,hcn7w,hcn8w,hcn9w]
  rr=[r1w[0],r2w[0],r3w[0],r4w[0],r5w[0],r6w[0],r7w[0],r8w[0],128b,$
  	r1w,r2w,r3w,r4w,r5w,r6w,r7w,r8w,r9w]
  gg=[g1w[0],g2w[0],g3w[0],g4w[0],g5w[0],g6w[0],g7w[0],g8w[0],128b,$
  	g1w,g2w,g3w,g4w,g5w,g6w,g7w,g8w,g9w]
  bb=[b1w[0],b2w[0],b3w[0],b4w[0],b5w[0],b6w[0],b7w[0],b8w[0],128b,$
  	b1w,b2w,b3w,b4w,b5w,b6w,b7w,b8w,b9w]
endelse
ncol=n_elements(cols) & mcol=max(icol)

;	usage
if ih eq 1 then begin
  print,'Usage: peasecolr,/white,/brewer,verbose=verbose,/help,$'
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
    print,'RED:    2,20-29'
    print,'GREEN:  3,30-39'
    if keyword_set(wbg) then print,'BROWN:  4,40-49' else print,'YELLOW: 4,40-49'
    print,'PINK:   5,50-59'
    print,'CYAN:   6,60-69'
    if keyword_set(wbg) then print,'YELLOW: 7,70-79' else print,'BROWN:  7,70-79'
    print,'SAFFRON: 8,80-89'
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
  if vv gt 100 then stop,'HALTing.  type .CON to continue'
endif

return
end
