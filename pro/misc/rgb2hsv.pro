function rgb2hsv,rr,gg,bb, _extra=e
;+
;function	rgb2hsv
;	convert a (R,G,B) 3-tuple in the 0:255 range to (H,S,V) in the (0:360,0:1,0:1) range
;
;syntax
;	hsvarr=rgb2hsv(R,G,B)
;
;parameters
;	R	[INPUT; required] red
;		* must be in the range [0,255]
;		* if outside this range, is truncated
;		* may be an array
;	G	[INPUT; required] green
;		* must be in the range [0,255]
;		* if outside this range, is truncated
;		* size must match that of G
;	B	[INPUT; required] blue
;		* must be in the range [0,255]
;		* if outside this range, is truncated
;		* size must match that of G
;
;keywords	NONE
;
;example
;	print,hexed(rgb2hsv([240,0,120,60],[1,1,1,1]))	;blue,red,green,yellow
;	#0000FF #FF0000 #00FF00 #FFFF00
;	also, 
;		.run rgb2hsv
;	will produce a circular swatch showing how hue and saturation vary
;
;notes
;	based on the conversion described in
;	https://en.wikipedia.org/wiki/HSL_and_HSV#From_RGB
;
;history
;	Vinay Kashyap (2019Nov; based on hsv2rgb.pro)
;-

;	usage
ok='ok' & np=n_params() & nr=n_elements(rr) & ng=n_elements(gg) & nb=n_elements(bb)
if np lt 3 then ok='Insufficient parameters' else $
 if nr eq 0 then ok='R is not given' else $
  if ng eq 0 then ok='G is not given' else $
   if nb eq 0 then ok='B is not given' else $
    if nr ne ng then ok='R and G array sizes do not match' else $
     if nr ne nb then ok='R and B array sizes do not match'
if ok ne 'ok' then begin
  print,'Usage: rgb2hsv(red,green,blue)'
  print,'  converts (R,G,B) to (H,S,V)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	define the output
if nb eq 1 then hsv=fltarr(3) else hsv=fltarr(3,nb)

;	for each entry, compute (H,S,V)
for i=0L,nr-1L do begin
  zR=float(byte(rr[i]))/255. & zG=float(byte(gg[i]))/255. & zB=float(byte(bb[i]))/255.
  zmax=max([zR,zG,zB],imx) & zmin=min([zR,zG,zB],imn)

  ;	hue
  if imx eq imn then hsv[0,i]=0. else begin
    case imx of
      0: hsv[0,i]=60.*(0.+(zG-zB)/(zmax-zmin))
      1: hsv[0,i]=60.*(2.+(zB-zR)/(zmax-zmin))
      2: hsv[0,i]=60.*(4.+(zR-zG)/(zmax-zmin))
    endcase
  endelse
  hsv[0,i]=(hsv[0,i]+360.) mod 360.	;in case of -ves

  ;	saturation
  if zmax eq 0 then hsv[1,i]=0. else hsv[1,i]=(zmax-zmin)/zmax

  ;	value
  if zmax eq 0 or zmin eq 1 then hsv[2,i]=0. else hsv[2,i]=zmax
endfor

;	transpose HSV from [3,N] to [N,3]
;	(this helps preserve the structure of the array even if N=1)
if nb gt 1 then hsv=transpose(hsv)

return,hsv
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example
;convert some arbitary (rgb) to (hsv) and back again
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	usage
print,''
print,'--------------------------------------------------------------------------------'
print,rgb2hsv()
print,'--------------------------------------------------------------------------------'

;	example: red
print,''
print,'RED=(255,0,0) should go to (0,1,1)'
hsv=rgb2hsv(255,0,0) & rgb=fix(hsv2rgb(hsv[0],hsv[1],hsv[2]))
print,'(hsv)='+strtrim(strjoin(hsv,','),2) & print,'(rrggbb)='+strtrim(strjoin(rgb,','),2)

;	example: green
print,''
print,'GREEN=(0,255,0) should go to (120,1,1)'
hsv=rgb2hsv(0,255,0) & rgb=fix(hsv2rgb(hsv[0],hsv[1],hsv[2]))
print,'(hsv)='+strtrim(strjoin(hsv,','),2) & print,'(rrggbb)='+strtrim(strjoin(rgb,','),2)

;	example: blue
print,''
print,'BLUE=(0,0,255) should go to (240,1,1)'
hsv=rgb2hsv(0,0,255) & rgb=fix(hsv2rgb(hsv[0],hsv[1],hsv[2]))
print,'(hsv)='+strtrim(strjoin(hsv,','),2) & print,'(rrggbb)='+strtrim(strjoin(rgb,','),2)

;	example vector colors
colors=['grey',     'skyblue',       'khaki',        'gold',      'maroon',        'sienna1',    'sienna2',      'sienna3',    'sienna4']
irgb=[[127,127,127], [135,206,235],   [240,230,140],  [255,215,0], [176,48,96],     [255,130,71], [238,121,66],   [205,104,57], [139,71,38]]
ihsv=[[0.,0.,0.5],    [203.,0.46,0.98],[54.,0.42,0.94],[51.,1.,1.], [338.,0.73,0.69],[19.,0.72,1], [19.,0.72,0.93],[19.,0.72,0.80],[19.,0.72,0.55]]
hsv=rgb2hsv(irgb[0,*],irgb[1,*],irgb[2,*]) & rgb=hsv2rgb(hsv[*,0],hsv[*,1],hsv[*,2])
print,''
for i=0,n_elements(colors)-1 do print,colors[i]+'	'+strjoin(string(irgb[*,i],'(i3)'),' ')+' == '+strjoin(strtrim(ihsv[*,i],2),',')+' --> '+strjoin(strtrim(reform(hsv[i,*]),2),',')+' --> '+strjoin(string(fix(reform(rgb[i,*])),'(i3)'),' ')

end
