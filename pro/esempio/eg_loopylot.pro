;+
;eg_loopylot
;	make plots of T and EM distributions using the output
;	of MCMC_DEM() for when the keyword LOOPY is set
;
;	the standard routine used to munge the output of MCMC_DEM(),
;	MCMC_PLOT, is not of much use when the LOOPY option is used.
;	this example routine shows which variables are of interest
;	and how they can be manipulated to get useful numbers out
;	of the output
;
;usage
;	mcmcsav='mcmc.save'
;	;set PSROOT to make hardcopy plots
;	;set NSMIN to avoid checking SIMPRB
;	.run eg_loopylot
;
;	NOTE:  SIMPRB are the probabilities at each point of the chain.
;	If a trend of any sort is discernible at the beginning parts,
;	do not include that part of the chain.  If there is a trend in
;	the later parts, do not use this chain: it hasn't converged.
;
;history
;	vinay k (jun2005)
;	now force contour axes arrays to match size of contour array
;	  (VK; apr2006)
;-

;	set up to run
ok='ok' & ssz=size(mcmcsav) & nssz=n_elements(ssz)
if not keyword_set(mcmcsav) then ok='MCMCSAV undefined: require name of output save file from MCMC_DEM()' else $
 if n_elements(mcmcsav) gt 1 then ok='MCMCSAV not understood: must be filename' else $
  if ssz[nssz-2] ne 7 then ok='MCMCSAV must be character string'
if ok ne 'ok' then message,ok
loadct,3 & peasecolr

;	restore savefile
fil=findfile(mcmcsav[0],count=nfil)
if nfil eq 0 then message,mcmcsav[0]+': not found'
restore,mcmcsav[0],verbose=verbose
sloop=1.5 & if n_tags(e) ne 0 then sloop=e.SLOOP
message,'DEM slope is '+strtrim(sloop,2),/informational

;	stop if doesn't make sense
ok='ok'
if n_elements(simdem) eq 0 then ok=mcmcsav[0]+': does not seem to be the right save file' else $
 if n_elements(simprb) ne nsim+1L then ok=mcmcsav[0]+': not in right format' else $
  if not keyword_set(loopy) then ok='LOOPY was not set during MCMC run; cannot work with this savefile'
if ok ne 'ok' then message,ok

if n_elements(nsmin) eq 0 then begin
  ;	pop up SIMPRB to see how much of it to use
  print,''
  print,'********************************************************************'
  print,'PICKRANGE(): choose the simulation number below which to discard all'
  print,'left-button-click+drag-right from simulation of interest to select'
  print,''
  print,'right-button-click to exit'
  print,''
  print,'********************************************************************'
  print,''
  osim=pickrange(simprb) & ismin=0L & if osim[0] gt 0 then ismin=osim[0]
endif else ismin=long(nsmin[0])>0
if ismin ge long(0.95*nsim) or ismin ge (nsim-20L) then $
	message,'something wrong? too few simulations left'

;	plot the DEMs, just for fun
ymin=min(simdem[where(simdem gt 0)],max=ymax)
if not keyword_set(psroot) then begin
  plot,logt,simdem[*,nsim],xtitle='logT',ytitle='DEM',/ylog,yrange=[ymin,ymax]
  for i=ismin,nsim-1 do oplot,logt,simdem[*,i],col=2
  wait,1
endif

;	set up some variables
nloopy=nt & if loopy[0] gt 1 then nloopy=long(loopy[0])
ksim=nsim-ismin
tcomp=fltarr(nloopy,ksim+1)
emcomp=fltarr(nloopy,ksim+1)

;	extract the T and EM components
for i=ismin,nsim do begin
  tmpdem=reform(simdem[*,i]) & demax=max(tmpdem)
  obreak=reverse(where(tmpdem[1:*]-tmpdem lt 0 and tmpdem gt 1e-7*demax,mobreak))
  if not keyword_set(psroot) then plot,logt,tmpdem,title=i,/ylog,yrange=[ymin/100.,ymax]
  ok='ok' & j=0L
  while ok eq 'ok' do begin
    logTmax=logt[obreak[0]] & EM=tmpdem[obreak[0]]
    tmp=loopem(logTmax,EM,logT=logT,sloop=sloop,verbose=verbose)
    tmp=EM*tmp/max(tmp) & tmpdem=tmpdem-tmp > 0
    tcomp[j,i-ismin]=logTmax & emcomp[j,i-ismin]=EM
    if not keyword_set(psroot) then oplot,logt,tmp,col=((j+1) mod 255)
    if not keyword_set(psroot) then oplot,logt,tmpdem,col=((j+1) mod 255)
    obreak=reverse(where(tmpdem[1:*]-tmpdem lt 0 and tmpdem gt 1e-7*demax,mobreak))
    ;
    j=j+1L
    if j eq nloopy then ok='j eq nloopy' else $
     if total(tmpdem) eq 0 then ok='total(tmpdem) eq 0' else $
      if mobreak eq 0 then ok='mobreak eq 0'
  endwhile
  ;for j=0L,mobreak-1L do begin
  ;  logTmax=logt[obreak[j]] & EM=tmpdem[obreak[j]]
  ;  tmp=loopem(logTmax,EM,logT=logT,sloop=sloop,verbose=verbose)
  ;  tmp=EM*tmp/max(tmp) & tmpdem=tmpdem-tmp > 0
  ;  tcomp[j,i-ismin]=logTmax & emcomp[j,i-ismin]=EM
  ;  oplot,logt,tmpdem,col=((j+1) mod 255)
  ;endfor
  if not keyword_set(psroot) then wait,0.001
endfor
;
message,'Temperature and EM components are in the variables:',/informational
help,tcomp,emcomp

;	make some useful histograms to plot
tmin=min(logt,max=tmax) & tbin=2*(median(logt[1:*]-logt))
ok=where(emcomp gt 0,mok) & if mok eq 0 then message,'No non-zero values found!'
emmax=max(emcomp[ok],min=emmin) & emnbin=(long(nsim/10L) > 10) < 101L
embin=(alog10(emmax)-alog10(emmin))/emnbin
ht=histogram(tcomp[ok],min=tmin,max=tmax,binsize=tbin)
xht=findgen(n_elements(ht))*tbin+tmin
hem=histogram(alog10(emcomp[ok]),min=alog10(emmin),max=alog10(emmax),binsize=embin)
xhem=findgen(n_elements(hem))*embin+alog10(emmin)
h2tem=hist_2d(tcomp[ok],alog10(emcomp[ok]),bin1=tbin,bin2=embin,min1=tmin,max1=tmax,min2=alog10(emmin),max2=alog10(emmax))
h2tem=h2tem/total(h2tem) & h2temmax=max(h2tem)
sh2tem=size(h2tem)
if sh2tem[1] ne n_elements(ht) then xht=findgen(sh2tem[1])*tbin+tmin
if sh2tem[2] ne n_elements(hem) then xhem=findgen(sh2tem[2])*embin+alog10(emmin)

;	hardcopy
if keyword_set(psroot) then begin
  szp=size(psroot) & nszp=n_elements(szp)
  outps='eg_loopylot' & if szp[nszp-2] eq 7 then outps=psroot[0]+'.ps'
  oldplot=!D.NAME & set_plot,'ps' & device,file=outps,/color
endif

if not keyword_set(psroot) and !D.NAME eq 'X' then window,0
  plot,xht,ht/total(ht),psym=10,xtitle='logTmax',ytitle='p(logTmax)',subtitle=mcmcsav[0]
  print,"plot,xht,ht/total(ht),psym=10,xtitle='logTmax',ytitle='p(logTmax)',subtitle=mcmcsav[0]"

if not keyword_set(psroot) and !D.NAME eq 'X' then window,1
  plot,xhem,hem/total(hem),psym=10,xtitle='logEMmax',ytitle='p(logEMmax)',subtitle=mcmcsav[0]
  print,"plot,xhem,hem/total(hem),psym=10,xtitle='logEMmax',ytitle='p(logEMmax)',subtitle=mcmcsav[0]"

if not keyword_set(psroot) and !D.NAME eq 'X' then window,2
  clevels=h2temmax*[0.001,0.01,0.1,0.3,0.6,0.9] & ccolors=indgen(n_elements(clevels))+1
  plot,tcomp[ok],alog10(emcomp[ok]),psym=1,col=81,xr=[tmin,tmax],yr=alog10([emmin,emmax]),/xs,/ys,subtitle=mcmcsav[0]
  print,"plot,tcomp[ok],alog10(emcomp[ok]),psym=1,col=81,xr=[tmin,tmax],yr=alog10([emmin,emmax]),/xs,/ys,subtitle=mcmcsav[0]"
  ;contour,h2tem,xht+tbin/2.,xhem+embin/2.,xr=[tmin,tmax],yr=alog10([emmin,emmax]),/xs,/ys,$
  contour,h2tem,xht+tbin,xhem+embin,xr=[tmin,tmax],yr=alog10([emmin,emmax]),/xs,/ys,$
	levels=clevels,c_colors=ccolors,/noerase,$
	xtitle='logTmax',ytitle='logEMmax',title='p(logTmax,logEMmax)'
  print,"contour,h2tem,xht+tbin,xhem+embin,xr=[tmin,tmax],yr=alog10([emmin,emmax]),/xs,/ys,levels=clevels,c_colors=ccolors,/noerase"

if keyword_set(psroot) then begin
  device,/close & set_plot,oldplot
  message,'Hardcopy output in file: '+outps,/informational
endif

end
