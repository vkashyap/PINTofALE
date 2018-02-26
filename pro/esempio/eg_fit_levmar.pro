;+
;EG_FIT_LEVMAR.PRO
;	program to road test FIT_LEVMAR
;
;usage:
;	.run adjustie lmcoeff levmarq fit_levmar eg_fit_levmar
;
;vinay k
;-

;	first, decide how many counts in line
if not keyword_set(ct) then $
	ct=long(10.^(findgen(10)*0.4+1.5)+0.5)
	;ct=long(10.^(findgen(20)*0.2+1.5)+0.5)
	;ct=long(10.^(findgen(10)*0.3+1)+0.5)
nct=n_elements(ct)

;	what's the initial sigma
if not keyword_set(sig) then sig=1.

;	how many simulations/set?
if not keyword_set(nsim) then nsim=100L

;	is there a background?
if not keyword_set(bkg) then bkg=0.

;	initialize
simsig=fltarr(nct,nsim) & simsige=simsig
simpos=fltarr(nct,nsim) & simpose=simpos
avsig=fltarr(nct) & sigavsig=avsig & avsige=avsig
avpos=fltarr(nct) & sigavpos=avpos & avpose=avpos
if not keyword_set(hmin) then hmin=-10.
if not keyword_set(hmax) then hmax=10.
if not keyword_set(hbin) then hbin=0.1

;	for each specified counts in line, simulate line, fit gaussian,
;	store sigmas and errors
for i=0L,nct-1L do begin			;{CT
  ict=ct(i)
  kilroy,dot=strtrim(ict,2)
  for j=0L,nsim-1L do begin			;{NSIM
    kilroy
    sct=randomn(seed,ict)		;generate line counts
    kct = ict/2 > 10
    delp=4.
    sct2=randomn(seed,kct)+delp	;add an extra line
    hct=histogram([sct,sct2],min=hmin,max=hmax,binsize=hbin)	;make histogram
    mhct=max(hct) < 10.
    ;mhct=max(hct)
    nh=n_elements(hct) & oo=where(hct gt 0)
    x=findgen(nh)*hbin+hmin & xx=x(oo)
    wei=1./(sqrt(float(hct)+0.75)+1.)^2
    ;
    ;DEBUG
    if j eq 0 then plot,x,hct,psym=10 else oplot,x,hct,col=100
    ;
    ;initialize variables for FIT_LEVMAR
    a=[0.,1,mhct,1.,1.,mhct/10.]		;model parameters
    ties=['a2 = a2 > 0.1','a3=a0+'+strtrim(delp,2),'a4=a1','a5 = a5 > 0.01']
    ties=['a3=a0+'+strtrim(delp,2),'a4=a1'] & freeze=[3,4]
    ;ties=['a3=a0+'+strtrim(delp,2)] & freeze=[3]
    sig=sqrt(float(hct)+0.75)+1.
    funcs='x3model'
    ;type='lorentz' & beta=2.
    ;type=['gauss','lorentz']
    type='gauss'
    ;
    ;call FIT_LEVMAR
    fit_levmar,x,hct,a,freeze=freeze,erra=erra,chisq=chisq,$
	ties=ties,funcs=funcs,sig=sig,type=type,beta=beta,/dumb
    ;fit_levmar,x(oo),hct(oo),a,freeze=freeze,erra=erra,chisq=chisq,$
    ;	ties=ties,funcs=funcs,sig=sig(oo),type=type,beta=beta
    ;
    simpos(i,j)=a(0)
    simpose(i,j)=erra(0)
    simsig(i,j)=a(1)
    simsige(i,j)=erra(1)
    ;
    np=n_elements(a)/3
    jp=lindgen(np) & pp=a(3*jp) & ww=a(3*jp+1) & hh=a(3*jp+2)
    ymod=mk_3model(x,pp,ww,hh,_extra=e)
    ;
    ;DEBUG
    oplot,x,ymod,col=150
    ;
    ;DEBUG if a(2) gt 100 then stop,i,j,ict,simsig(i,j),simsige(i,j)
  endfor					;J=0,NSIM-1}
  wait,1
endfor						;I=0,NCT-1}

;	output
for i=0,nct-1 do avsig(i)=total(simsig(i,*))/nsim
for i=0,nct-1 do sigavsig(i)=sqrt(total(simsig(i,*)^2)/nsim-avsig(i)^2)
for i=0,nct-1 do avsige(i)=total(simsige(i,*))/nsim
for i=0,nct-1 do avpos(i)=total(simpos(i,*))/nsim
for i=0,nct-1 do sigavpos(i)=sqrt(total(simpos(i,*)^2)/nsim-avpos(i)^2)
for i=0,nct-1 do avpose(i)=total(simpose(i,*))/nsim

;	plot
;plot,ct,sigavpos,/xl,/yl & oplot,ct,avpose,line=1
;plot,ct,avpose,/xl,/yl & oplot,ct,sigavpos,line=1
;oplot,ct,float(ct(0))/ct,col=125
;oplot,ct,1./sqrt(ct),col=150
;oplot,ct,(float(ct(0))^2)/float(ct)^2,col=100
plot,ct,avsige,/xl,/yl,xtitle='ct',ytitle='avg(err(wdt))'

end
