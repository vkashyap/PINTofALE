!p.multi=[0,1,2]
printcap,'xinti',3
;-----------------------------------------------------------------------
merge_line,'../SPEX/code/emisspex','emisschianti',atom='Fe',omatch=om,$
	mish=mish,mash=mash,prec=0.01,doubt=3,frat=4,batch=0.001
w1=om.(0) & w2=om.(2) & dw=abs(w1)-abs(w2) & nw=n_elements(w1)
rat=fltarr(nw)
for i=0,nw-1 do begin
  f1=om.(1) & f2=om.(3)
  ymax=max([f1(*,i),f2(*,i)])>1. & ymin=1e-10*ymax
  oo=where(f1(*,i) gt ymin and f2(*,i) gt ymin,moo)
  if moo gt 0 then begin
    tf1=total(f1(oo,i)) & tf2=total(f2(oo,i))
    rat(i)=tf1/tf2
  endif else rat(i)=1e8
  if om.q(i) eq 1 then c1='?' else c1=''
  print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
	' '+string(dw(i),'(f9.3)')+' '+strtrim(tf1,2)+' '+$
	strtrim(tf2,2)+' '+strtrim(rat(i),2)+' '+c1
endfor
;
;plot,dw,rat,psym=1,/ylog,title='wvl diff v/s intensity ratios'
rt_hrt=histogram(alog10(rat)+4,min=0,max=12,binsize=0.05)
rt_hdw=findgen(n_elements(rt_hrt))*0.05-4
plot,rt_hdw,rt_hrt,psym=10,xtitle='log!d10!n(Intensity Ratio)',ytitle='#',$
	title='SPEX v/s CHIANTI for Fe',xr=[-1,1]
dw_hdw=histogram(dw,min=-0.05,max=0.05,binsize=0.001)
dw_hrt=findgen(n_elements(dw_hdw))*0.001-0.05
plot,dw_hrt,dw_hdw,psym=10,xtitle='wavelength difference',ytitle='#',$
	title='SPEX v/s CHIANTI for Fe'

;-----------------------------------------------------------------------
;this gets all wavelength matches with 2e-3 A, but essentially no temperature
;selection
merge_line,'../SPEX/code/emisspex','emisschianti',atom='Fe',omatch=om,$
	prec=1e-3,doubt=1e-3,tdex=1,frat=1e8,batch=0.001
;
w_w1=om.(0) & w_w2=om.(2) & w_dw=abs(w_w1)-abs(w_w2) & w_nw=n_elements(w_w1)
w_rat=fltarr(w_nw)
for i=0,w_nw-1 do begin
  f1=om.(1) & f2=om.(3)
  ymax=max([f1(*,i),f2(*,i)])>1. & ymin=1e-10*ymax
  oo=where(f1(*,i) gt ymin and f2(*,i) gt ymin,moo)
  if moo gt 0 then begin
    tf1=total(f1(oo,i)) & tf2=total(f2(oo,i))
    w_rat(i)=tf1/tf2
  endif else w_rat(i)=1e8
  if om.q(i) eq 1 then c1='?' else c1=''
  print,strtrim(i,2)+' '+string(w_w1(i),'(f8.3)')+' '+string(w_w2(i),'(f8.3)')+$
	' '+string(w_dw(i),'(f9.3)')+' '+strtrim(tf1,2)+' '+$
	strtrim(tf2,2)+' '+strtrim(w_rat(i),2)+' '+c1
endfor
;
;plot,w_dw,w_rat,psym=1,/ylog,title='wvl diff v/s intensity ratios'
w_rth=histogram(alog10(w_rat)+4,min=0,max=12,binsize=0.1)
w_dwh=findgen(n_elements(w_rth))*0.1-4
plot,w_dwh,w_rth,psym=10,xtitle='log!d10!n(Intensity Ratio)',ytitle='#',$
	title='perfect *wavelength* matching'
w_rth2=histogram(alog10(w_rat)+4,min=0,max=12,binsize=0.05)
w_dwh2=findgen(n_elements(w_rth2))*0.05-4
plot,w_dwh2,w_rth2,psym=10,xtitle='log!d10!n(Intensity Ratio)',ytitle='#',$
	title='perfect *wavelength* matching',xr=[-0.5,0.5]

;-----------------------------------------------------------------------
;this gets all temperature and profile matches exactly, but essentially
;no wavelength selection (anything within 10A is considered an OK match!)
merge_line,'../SPEX/code/emisspex','emisschianti',atom='Fe',omatch=om,$
	prec=5.,doubt=1,tdex=0.1,frat=1.1,batch=0.001
;
t_w1=om.(0) & t_w2=om.(2) & t_dw=abs(t_w1)-abs(t_w2) & t_nw=n_elements(t_w1)
t_rat=fltarr(t_nw)
for i=0,t_nw-1 do begin
  f1=om.(1) & f2=om.(3)
  ymax=max([f1(*,i),f2(*,i)])>1. & ymin=1e-10*ymax
  oo=where(f1(*,i) gt ymin and f2(*,i) gt ymin,moo)
  if moo gt 0 then begin
    tf1=total(f1(oo,i)) & tf2=total(f2(oo,i))
    t_rat(i)=tf1/tf2
  endif else t_rat(i)=1e8
  if om.q(i) eq 1 then c1='?' else c1=''
  print,strtrim(i,2)+' '+string(t_w1(i),'(f8.3)')+' '+string(t_w2(i),'(f8.3)')+$
	' '+string(t_dw(i),'(f9.3)')+' '+strtrim(tf1,2)+' '+$
	strtrim(tf2,2)+' '+strtrim(t_rat(i),2)+' '+c1
endfor
;
;plot,t_dw,t_rat,psym=1,/ylog,title='wvl diff v/s intensity ratios'
t_dwh=histogram(t_dw+10,min=0,max=20,binsize=1)
t_rth=findgen(n_elements(t_dwh))*1-10.
plot,t_rth,t_dwh,psym=10,xtitle='Wavelength Difference',ytitle='#',$
	title='(near) perfect temperature profile matching'
t_dwh2=histogram(t_dw+10,min=0,max=20,binsize=0.001)
t_rth2=findgen(n_elements(t_dwh2))*0.001-10.
plot,t_rth2,t_dwh2,psym=10,xtitle='Wavelength Difference',ytitle='#',$
	title='(near) perfect temperature profile matching',xr=[-0.1,0.1]

;-----------------------------------------------------------------------
save,file='xinti_mrg.sav'
device,/cl

end
