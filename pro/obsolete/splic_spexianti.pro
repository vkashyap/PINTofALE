;+
;SPLIC_SPEXIANTI
;	program to splice the CHIANTI spectral line database with that of SPEX
;
;history
;	based on mrg_spexianti (VK; Jan98)
;-

message,'OBSOLETE!',/informational

;	initialize
n_pres=8 & d_pres=1 & min_pr=13.		;pressure grid
n_den=14 & d_eden=2. & min_den=8.		;electron density grid
chianti='/data/drake7/vinay/SCAR/emissivity/chianti'
spex='/data/drake7/vinay/SCAR/SPEX/code/spex'
temp='/data/drake7/vinay/SCAR/emissivity/emisstmp'
xinti='/data/drake7/vinay/SCAR/emissivity/emissplice'	;comment out to
							;avoid output
elem=['H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al','Si','P',$
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;		first pass -- get SPEX < WCUT \AA
wcut=50.
wrange=[-1.,wcut]

;	emissivities at constant pressure
for iz=0,n_elements(elem)-1 do begin
  mish=0 & mash=0

  ;	first get the matches for a "good" logP
  merge_line,spex,'none',atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logP=16,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001,wrange=wrange
  w1=om.(0) & w2=om.(2) & dw=abs(abs(w1)-abs(w2)) & nw=n_elements(w1)
  print,'	',elem(iz)
  if w1(0) eq 0. and w2(0) eq 0. then nw=0L
  for i=0L,nw-1L do begin
    f1=om.(1) & f2=om.(3) & tf1=total(f1(*,i)) & tf2=total(f2(*,i))
    if om.q(i) eq 1 then c1='?' else c1=''
    print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
    ' '+string(dw(i),'(f8.3)')+' '+strtrim(tf1,2)+' '+$
    strtrim(tf2,2)+' '+strtrim(tf1/(tf2>1),2)+' '+c1
  endfor
  
  ;	now for all the pressures
  for ip=0,n_pres-1 do begin
    logp=min_pr+ip*d_pres		;[cm^-3 K]
    merge_line,spex,'none',xinti,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001,wrange=wrange
  endfor

  ;	and for all the densities
  for id=0,n_den-1 do begin
    n_e = 10.^(min_den+id/2)/d_eden		;[cm^-3]
    if d_eden eq 2 then d_eden=1. else d_eden=2.
    merge_line,spex,'none',xinti,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,n_e=n_e,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001,$
	wrange=wrange
  endfor

endfor

;		second pass -- get CHIANTI > WCUT \AA
wrange=[wcut,1e10]

;	emissivities at constant pressure
for iz=0,n_elements(elem)-1 do begin
  mish=0 & mash=0

  ;	first get the matches for a "good" logP
  merge_line,chianti,'none',atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logP=16,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001,wrange=wrange
  w1=om.(0) & w2=om.(2) & dw=abs(abs(w1)-abs(w2)) & nw=n_elements(w1)
  print,'	',elem(iz)
  if w1(0) eq 0. and w2(0) eq 0. then nw=0L
  for i=0L,nw-1L do begin
    f1=om.(1) & f2=om.(3) & tf1=total(f1(*,i)) & tf2=total(f2(*,i))
    if om.q(i) eq 1 then c1='?' else c1=''
    print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
    ' '+string(dw(i),'(f8.3)')+' '+strtrim(tf1,2)+' '+$
    strtrim(tf2,2)+' '+strtrim(tf1/(tf2>1),2)+' '+c1
  endfor
  
  ;	now for all the pressures
  for ip=0,n_pres-1 do begin
    logp=min_pr+ip*d_pres		;[cm^-3 K]
    merge_line,chianti,'none',temp,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001,wrange=wrange
  endfor

  ;	and for all the densities
  for id=0,n_den-1 do begin
    n_e = 10.^(min_den+id/2)/d_eden		;[cm^-3]
    if d_eden eq 2 then d_eden=1. else d_eden=2.
    merge_line,chianti,'none',temp,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,n_e=n_e,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001,$
	wrange=wrange
  endfor

endfor

;		third pass -- append CHIANTI > WCUT \AA

;	emissivities at constant pressure
for iz=0,n_elements(elem)-1 do begin
  mish=0 & mash=0

  ;	first get the matches for a "good" logP
  merge_line,xinti,temp,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logP=16,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001
  w1=om.(0) & w2=om.(2) & dw=abs(abs(w1)-abs(w2)) & nw=n_elements(w1)
  print,'	',elem(iz)
  if w1(0) eq 0. and w2(0) eq 0. then nw=0L
  for i=0L,nw-1L do begin
    f1=om.(1) & f2=om.(3) & tf1=total(f1(*,i)) & tf2=total(f2(*,i))
    if om.q(i) eq 1 then c1='?' else c1=''
    print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
    ' '+string(dw(i),'(f8.3)')+' '+strtrim(tf1,2)+' '+$
    strtrim(tf2,2)+' '+strtrim(tf1/(tf2>1),2)+' '+c1
  endfor
  
  ;	now for all the pressures
  for ip=0,n_pres-1 do begin
    logp=min_pr+ip*d_pres		;[cm^-3 K]
    merge_line,xinti,temp,xinti,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001
  endfor

  ;	and for all the densities
  for id=0,n_den-1 do begin
    n_e = 10.^(min_den+id/2)/d_eden		;[cm^-3]
    if d_eden eq 2 then d_eden=1. else d_eden=2.
    merge_line,xinti,temp,xinti,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,n_e=n_e,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001
  endfor

endfor

end
