;+
;script	Zfip
;	plots the First Ionization Potential vs atomic number Z
;
;usage
;	white=0	;or 1
;	Zmax=30	;or 109
;	.run Zfip
;
;requires
;	INICON
;	PEASECOLR
;
;history
;	vinay kashyap (Jun2007)
;	slight modifications (VK; Jul2012)
;-

;	initialize colors
if not keyword_set(white) then white=0
peasecolr,white=white & loadct,3 & peasecolr,white=white

;	extract FIP
inicon,fip=fip,atom=atom

if not keyword_set(Zmax) then Zmax=30	;only look at these many elements
if Zmax gt n_elements(fip) then Zmax=n_elements(fip)

fip=fip[0:Zmax-1] & atom=atom[0:Zmax-1]
ymin=min(fip,max=ymax)

if !d.name eq 'PS' then th=5 else th=2
if !d.name eq 'X' then window,xsize=900,ysize=768
plot,findgen(Zmax)+1,fip,psym=3,xrange=[0,Zmax+1],/xstyle,yrange=[ymin-1,ymax+1],/ystyle,$
	xtitle='Z',ytitle='FIP [eV]',thick=th,xthick=th,ythick=th,charsize=1.3
for i=1,Zmax do $
  xyouts,i,fip[i-1],strtrim(atom[i-1],2),charsize=2,$
  charthick=th,align=0.5,color=1+(i mod 8)

end
