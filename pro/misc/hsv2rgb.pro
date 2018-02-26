function hsv2rgb,hh,ss,vv, _extra=e
;+
;function	hsv2rgb
;	convert (H,S,V) in the (0:360,0:1,0:1) range to an (R,G,B) 3-tuple in the 0:255 range
;
;syntax
;	rgbarr=hsv2rgb(H,S,V)
;
;parameters
;	H	[INPUT; required] Hue
;		* must be in the range [0,360]
;		* if outside this range, is taken to be mod 360
;		* may be an array
;	S	[INPUT; required] Saturation
;		* must be in the range [0,1]
;		* if outside this range, is truncated
;		* size must match that of H
;	V	[INPUT] Value
;		* must be in the range [0,1]
;		* if outside this range, is rescaled by the maximum
;		* if not given, assumed to be 1.0
;		* if size is greater than H, remaining elements are ignored
;
;keywords	NONE
;
;example
;	print,hexed(hsv2rgb([240,0,120,60],[1,1,1,1]))	;blue,red,green,yellow
;	#0000FF #FF0000 #00FF00 #FFFF00
;	also, 
;		.run hsv2rgb
;	will produce a circular swatch showing how hue and saturation vary
;
;notes
;	based on the conversion described in
;	http://en.wikipedia.org/wiki/HSL_and_HSV#Coverting_to_RGB
;
;history
;	Vinay Kashyap (2015Oct)
;	allowed (H,S,V) to be arrays (VK; 2015Nov)
;	allow V to be >1, in which case it just gets rescaled (VK; 2015Dec)
;-

;	usage
ok='ok' & np=n_params() & nh=n_elements(hh) & ns=n_elements(ss) & nv=n_elements(vv)
if np lt 2 then ok='Insufficient parameters' else $
 if nh eq 0 then ok='Hue is not given' else $
  if ns eq 0 then ok='Saturation is not given' else $
   if nh ne ns then ok='Hue and Saturation array sizes do not match'
if ok ne 'ok' then begin
  print,'Usage: hsv2rgb(hue,saturation,value)'
  print,'  converts (H,S,V) to (R,G,B)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	define the output
if nh eq 1 then rgb=bytarr(3) else rgb=bytarr(3,nh)

;	parse inputs
zH = (1.0*hh) mod 360.
zS = ((1.0*ss) > 0.) < 1.
zV = 0.*zh+1.
if nv gt 0 then zv[0L:(nv<nh)-1L]=vv[0L:(nv<nh)-1L]
if nv gt 1 then if max(zv,/nan) gt 1 then zv=za/max(zv,/nan)

;	see http://en.wikipedia.org/wiki/HSL_and_HSV#Coverting_to_RGB
zC = zV * zS
Hp = zH/60.
zX = zC*(1.-abs((Hp mod 2) -1.))
for i=0L,nh-1L do begin
  mm = zV[i]-zC[i]
  if 0 le Hp[i] and Hp[i] lt 1 then begin & rgb[0,i]=byte((zC[i]+mm)*255.+0.5) & rgb[1,i]=byte((zX[i]+mm)*255.+0.5) & rgb[2,i]=byte((0.+mm)*255.+0.5) & endif
  if 1 le Hp[i] and Hp[i] lt 2 then begin & rgb[0,i]=byte((zX[i]+mm)*255.+0.5) & rgb[1,i]=byte((zC[i]+mm)*255.+0.5) & rgb[2,i]=byte((0.+mm)*255.+0.5) & endif
  if 2 le Hp[i] and Hp[i] lt 3 then begin & rgb[0,i]=byte((0.+mm)*255.+0.5) & rgb[1,i]=byte((zC[i]+mm)*255.+0.5) & rgb[2,i]=byte((zX[i]+mm)*255.+0.5) & endif
  if 3 le Hp[i] and Hp[i] lt 4 then begin & rgb[0,i]=byte((0.+mm)*255.+0.5) & rgb[1,i]=byte((zX[i]+mm)*255.+0.5) & rgb[2,i]=byte((zC[i]+mm)*255.+0.5) & endif
  if 4 le Hp[i] and Hp[i] lt 5 then begin & rgb[0,i]=byte((zX[i]+mm)*255.+0.5) & rgb[1,i]=byte((0.+mm)*255.+0.5) & rgb[2,i]=byte((zC[i]+mm)*255.+0.5) & endif
  if 5 le Hp[i] and Hp[i] lt 6 then begin & rgb[0,i]=byte((zC[i]+mm)*255.+0.5) & rgb[1,i]=byte((0.+mm)*255.+0.5) & rgb[2,i]=byte((zX[i]+mm)*255.+0.5) & endif
endfor

;	transpose RGB from [3,N] to [N,3]
;	(this helps preserve the structure of the array even if N=1)
if nh gt 1 then rgb=transpose(rgb)

return,rgb
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example
;make circular swatch with varying hue (by degree) and saturation (by offaxis)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if not keyword_set(nhue) then nhue=360
if not keyword_set(nsat) then nsat=32
dhue=360./nhue & hues=findgen(nhue)*dhue
dsat=1./(nsat-1.) & sats=findgen(nsat)*dsat
if not keyword_set(val) then val=1.0

peasecolr & loadct,0 & peasecolr
if !d.name eq 'X' then window,0,title='Hue and Saturation colors',xsize=512,ysize=512,retain=2

plot,[0],/nodata,xstyle=5,ystyle=5,xr=[-1,1],yr=[-1,1],position=[0.1,0.1,0.9,0.9]
  for j=0L,nsat-1L do begin
    r1=sats[j] & r2=sats[j]+dsat
for i=0L,nhue-1L do begin
  th1=hues[i] & th2=th1+dhue
  c1=cos(th1*!pi/180.) & c2=cos(th2*!pi/180.)
  s1=sin(th1*!pi/180.) & s2=sin(th2*!pi/180.)
    xx=[r1*c1,r2*c1,r2*c2,r1*c2,r1*c1]
    yy=[r1*s1,r2*s1,r2*s2,r1*s2,r1*s1]
    setkolor,hexed(hsv2rgb(th1,r1,abs(val[0]))),99,/quiet
    polyfill,xx,yy,col=99
  endfor
endfor

end
