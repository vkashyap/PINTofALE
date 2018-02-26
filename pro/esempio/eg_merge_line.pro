;+
;EG_MERGE_LINE.PRO
;
;example program to call MERGE_LINE
;
;
;vinay kashyap (1998dec6)
;added call to getpoadef (VK; 2015aug5)
;-

dir1='$SPEX'		;SPEX directory
dir2='$CHIANTI'		;CHIANTI directory
;dir3=filepath('tmp',root_dir=getpoadef())	;output database -- comment out if
						;nothing should be written to disk
atom='Fe'		;stick to Fe lines
prec=0.01		;precision of data in Ang
doubt=3			;fuzziness of boundary
tdex=0.5		;match peaks in Temperature response
frat=4			;maximum ratio of matched summed emissivities
batch=0.0001		;ask no questions, wait but a jiffy

;	merge lines in DIR2 into lines in DIR1
merge_line,dir1,dir2,dir3,atom=atom,omatch=om,mish=mish,mash=mash,$
	prec=prec,doubt=doubt,tdex=tdex,frat=frat,batch=batch
;
;OM is a structure that contains the output matches
;MISH & MASH ensure that the next time MERGE_LINE is called, it does not
;	ask a whole lot of questions

merge_line,'../SPEX/code/emisspex','emisschianti',atom='Fe',omatch=om,$
	mish=mish,mash=mash,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001

w1=[om.(0)] & w2=[om.(2)] & f1=om.(1) & f2=om.(3)	;unpack OM
dw=abs(w1)-abs(w2) & nw=n_elements(w1)			;delta_Ang
rat=fltarr(nw)						;[1]/[2]

for i=0,nw-1 do begin			;{for each matched wavelength
  ;	compute ratios of summed emissivities over the region where
  ;	both lines have non-zero response
  ymax=max([f1(*,i),f2(*,i)]) & ymin=1e-10*ymax & yr=[ymin,ymax]
  hh=where(f1(*,i) gt ymin and f2(*,i) gt ymin,nhh)
  if nhh gt 0 then begin
    rat(i)=total(f1(hh,i))/total(f2(hh,i))
  endif else rat(i)=1e8
  ;
  ;	quality of match
  if om.q(i) eq 1 then c1='?' else c1=''
  print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
	' '+string(dw(i),'(f9.3)')+' '+strtrim(rat(i),2)+' '+c1
  ;print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
  ;	' '+string(dw(i),'(f9.3)')+' '+strtrim(tf1,2)+' '+$
  ;	strtrim(tf2,2)+' '+strtrim(rat(i),2)+' '+c1
endfor					;i=0,nw-1}

;	scatter plot
plot,dw,rat,psym=1,/ylog,title='wvl diff v/s intensity ratios'

;	histogram of ratios
rt_hrt=histogram(alog10(rat)+4,min=0,max=12,binsize=0.1)
rt_hdw=findgen(n_elements(rt_hrt))*0.1-4
plot,rt_hdw,rt_hrt,psym=10,xtitle='log!d10!n(Intensity Ratio)',ytitle='#',$
	title='SPEX v/s CHIANTI for Fe'

;	histogram of wavelength differences
dw_hdw=histogram(dw,min=-0.1,max=0.1,binsize=0.001)
dw_hrt=findgen(n_elements(dw_hdw))*0.001-0.1
plot,dw_hrt,dw_hdw,psym=10,xtitle='wavelength difference',ytitle='#',$
	title='SPEX v/s CHIANTI for Fe'

end
