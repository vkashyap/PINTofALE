;+
;MRG_SPEXIANTI
;	program to merge the CHIANTI spectral line database with that of SPEX
;
;history
;	vinay kashyap (Jan96)
;	slightly modified in no important way (VK; Nov98)
;-

;	initialize
pmin=13 & pmax=20 & dp=1		;pressure grid
dmin=8 & dmax=13 & dd=1			;e-density grid

;	input/output directories
chianti='/data/fubar/SCAR/emissivity/chianti'
spex='/data/fubar/SCAR/emissivity/spex'
xinti='/data/fubar/SCAR/emissivity/xinti'	;comment out to avoid output
elem=['H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al','Si','P',$
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

for iz=0,n_elements(elem)-1 do begin
  ;iz=0 & elem=['Fe13']		;to get only one set, for e.g.
  mish=0 & mash=0
  for ip=pmin,pmax,dp do begin
    logp=float(ip)
    merge_line,chianti,spex,xinti,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	logp=logp,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001
    w1=om.(0) & w2=om.(2) & dw=abs(abs(w1)-abs(w2)) & nw=n_elements(w1)
    print,'	',elem(iz),logp
    if w1(0) eq 0. and w2(0) eq 0. then nw=0L
    for i=0L,nw-1L do begin
      f1=om.(1) & f2=om.(3) & tf1=total(f1(*,i)) & tf2=total(f2(*,i))
      if om.q(i) eq 1 then c1='?' else c1=''
      print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
	' '+string(dw(i),'(f8.3)')+' '+strtrim(tf1,2)+' '+$
	strtrim(tf2,2)+' '+strtrim(tf1/(tf2>1),2)+' '+c1
    endfor
  endfor
endfor

for iz=0,n_elements(elem)-1 do begin
  mish=0 & mash=0
  for id=dmin,dmax,dd do begin
    eden=float(id)
    merge_line,chianti,spex,xinti,atom=elem(iz),omatch=om,mish=mish,mash=mash,$
	n_e=eden,prec=0.01,doubt=3,tdex=0.5,frat=4,batch=0.001
    w1=om.(0) & w2=om.(2) & dw=abs(abs(w1)-abs(w2)) & nw=n_elements(w1)
    print,'	',elem(iz),eden
    if w1(0) eq 0. and w2(0) eq 0. then nw=0L
    for i=0L,nw-1L do begin
      f1=om.(1) & f2=om.(3) & tf1=total(f1(*,i)) & tf2=total(f2(*,i))
      if om.q(i) eq 1 then c1='?' else c1=''
      print,strtrim(i,2)+' '+string(w1(i),'(f8.3)')+' '+string(w2(i),'(f8.3)')+$
	' '+string(dw(i),'(f8.3)')+' '+strtrim(tf1,2)+' '+$
	strtrim(tf2,2)+' '+strtrim(tf1/(tf2>1),2)+' '+c1
    endfor
  endfor
endfor

end
