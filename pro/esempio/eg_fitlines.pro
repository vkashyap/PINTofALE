;+
;EG_FITLINES
;	example program to exercize FITLINES.PRO
;
;usage:
;	.run fitlines_event fitlines eg_fitlines
;
;vinay kashyap
;-

peasecolr
ct=0 & pos=0 & wdt=0

;	first, decide how many counts in each line
nct=3L
if not keyword_set(ct) then $
	ct=long(randomu(seed,nct)*1e4)
nct=n_elements(ct)

;	what're the sigma
nsig=n_elements(wdt)
if nsig ne nct then wdt=abs(randomn(seed,nct))*3 > 0.5
nsig=n_elements(wdt)

;	what're the locations
npos=n_elements(pos)
if npos ne nct then pos=randomu(seed,nct)*20.-10.
npos=n_elements(pos)

;	initialize
if not keyword_set(hmin) then hmin=-10.
if not keyword_set(hmax) then hmax=10.
if not keyword_set(hbin) then hbin=0.1

;	make spectrum
nbin=long((hmax-hmin)/hbin+1)
x=findgen(nbin)*hbin+hmin & y=0.*x
for i=0,nct-1 do begin
  sct=randomn(seed,ct(i))*wdt(i)+pos(i)
  hct=histogram(sct,min=hmin,max=hmax,binsize=hbin)
  y=y+hct
endfor

sig=sqrt(abs(y)+0.75)+1.
;flx=ct*hbin
flx=ct/wdt/sqrt(2.*!pi)
thaw=intarr(3*nct)+1

ostr=fitlines(x,y,ysig=sig,pos=pos,wdt=wdt,flx=flx,ferrp=sqrt(flx),/dumb)

end
