pro mcmc_plot,logt,simdem,demerr,simprb,schme,storidx,slect=slect,col_rng=col_rng,$
ps_fil=ps_fil,col_tabl=col_tabl,col_shft=col_shft,clev=clev,demslope=demslope,slopsig=slopsig,$
slopedt=slopedt,goodt=goodt,sampct=sampct,tbuff=tbuff,monochrom=monochrom, $
spshft=spshft,shftsig=shftsig,slopenv=slopenv,slpsprd=slpsprd,$
slope1=slope1,shift1=shift1,_extra = e 

;+ 
;procedure   mcmc_plot
;
;       This procedure plots Monte Carlo Markov-Chain DEM
;       reconstruction results. One sigma uncertanties are 
;       shown as a grey scale (or otherwise specified) with a gradient
;       determined by how the N simulated DEMs (which  ordain the one sigma
;       confidence intervals) are distributed about the best-fit DEM. 
;       Alternatively, the gradient may be determined by the p(D|M) of 
;       each DEM.
; 
;       NOTE: ALL REQUIRED INPUTS ARE STANDARD MCMC_DEM() OUTPUTS
;
;parameters
;       logt    [INPUT; required] array which determines the mid-bin values 
;               of the temperature grid on which DEM is defined
;       simdem  [INPUT; required] two dimensional array [NT, NSIM+1]
;               simulated arrays, the last of which is the best fit. 
;       demerr  [INPUT; required] confidence bounds on MAP estimates of DEM
;		* DEMERR(*,0) are lower bounds, DEMERR(*,1) are upper bounds
;       simprb  [INPUT; required] the p(D|M) of each DEM            
;       schme   [INPUT; optional] string designating plotting strategy. Can be:
;
;                   'ERRBAR'      - Just plot error bars using DEMERR (default)
;                   'DIFFERENCE'  - Colors are determined by
;                                   difference of each DEM simulation
;                                   bin value from the best fit DEM
;                                   values. Colors are polyfilled
;                                   between simulated values within
;                                   each  logT bin. 
;                   'INDEX'       - Colors are determined by the order in which
;                                   they appear (starting from the best
;                                   DEM outwards). Colors are polyfilled
;                                   between simulated values within each
;                                   logT bin. 
;                   'SPLINE INDEX'- Colors are determined by the order in which
;                                   they appear (starting from the best
;                                   DEM outwards). Colors are polyfilled
;                                   between simulated values within each
;                                   logT bin. Spline interpolated
;                                   curves are plotted in lieu of
;                                   histogram style bins
;                   'PROB'        - Colors are determined by the p(D|M)
;                                   of each DEM. Values are oploted
;                                   with thick bars between simulated
;                                   values within each  logT bin.
;                   'SPLINE'      - Simulated DEMS are spline
;                                   interpolated and colors determined
;                                   by the p(D|M) associated with
;                                   each. Each spline is oploted onto
;                                   a best-fit dem plot. 
;                   'NICE'        - For each temperature bin, a median simulation
;                                   value is calculated. The median values will be
;                                   plotted (dotted line) together with the DEM
;                                   simulation (solid line) that best fits them.
;                                   Uncertainty will be displayed by plotting the MIN and
;                                   MAX simulation values for each
;                                   bin. (dashed lines) Use the SLECT and CLEV
;                                   keywords to limit which simulations to include.
;                   'NICE DIFFERENCE' - A combination of DIFFERENCE and NICE
;                                   Colors are determined by
;                                   difference of each DEM simulation
;                                   bin value from the best fit DEM
;                                   simulation to the median value as
;                                   described in the 'NICE'
;                                   definition.  Colors are polyfilled
;                                   between simulated values within
;                                   each  logT bin. 
;                   'NICE INDEX'    - A combination of DIFFERENCE and NICE 
;                                   Colors are determined by the order in 
;                                   which they appear (starting from the 
;                                   best fit simulation to the median value 
;                                   as described in the 'NICE' definition. 
;                                   Colors are plyfilled between simulated
;                                   values within each logT value.                                 
;                   'SUPPORT'       - Just plot upper and lower bounds of the 
;                                   support in each bin given the set of 
;                                   selected DEMs (see SLECT keyword) 
;
;       storidx [INPUT ; optional]  from MCMC_DEM(), contains the indices of 
;                                   of all the sampled parameters stored in
;                                   STORPAR. If a temperature bin has not been 
;                                   adequately sampled ( default threshold is
;                                   nsim*0.10 ), then the bin is not
;                                   plotted. Only bins at the ends of the Log
;                                   T range are thrown out. (see SAMPCT keyword)
; 
;keywords
;       slect    [INPUT]  choose which DEMS to include default=[2] 
;            
;                   1:ALL     -Include all simulations
;                   2:CHI^2   -Limit simulations to best 50% p(D|M) DEMS
;                              (default = 50, use CLEV keyword to toggle)
;                   3:Errors  -Limit simulations to values in each bin that
;                              lie  within the specified confidence bounds (demerr)
;       clev     [INPUT]  To be used in conjunction with slect=2,
;                         determines what percentage by which to limit
;                         sims in p(D|M) 
;       col_rng  [INPUT]  index range in the color table to cover. 
;                         gradient will taper
;                         off from 0 (darkest) to 
;                         col_rng (brightest). (maximum is 255)
;       col_tabl [INPUT]  color table on which to base the gradient       
;       col_shft [INPUT]  number of color indices by which to shift
;                         brighter(or darker if black background), col_rng will
;                         be adjusted accordingly
;       ps_fil   [I/O]    name of ps file to output. If not specified 
;                         plot will be sent to current device
;       demslope [I/O]    if set, linear fits are made to the simulated DEMs 
;                         in a temperature range specefied by keyword SLOPEDT 
;                         and DEMSLOPE will return the mean slope. The
;                         distribution of slopes will be plotted
;       slopesig [OUTPUT] if set, returns standard deviation slopes of
;                            linear fits 
;       slopedt  [INPUT]  temperature range specified in Log(T[K]) within which
;                         to make linear fits to DEM 
;                             * if 2-element array, then SLOPEDT =
;                               [min(SLOPEDT),max(SLOPEDT)] 
;                             * if 1-element array, then SLOPEDT = 
;                               [SLOPEDT, logt(max(bestdem)))] or 
;                               [SLOPEDT, logt(min(bestdem))]
;                               if logt(max(bestdem)) eq SLOPEDT 
;                             * if not set, then SLOPEDT = [min(logt),
;                               logt(max(bestdem))] or [SLOPEDT,
;                               logt(min(bestdem))] if logt(max(bestdem)) eq SLOPEDT
;       sampct   [INPUT]  when STORIDX is input, set this to threshold for T
;                         bin to be included. Expressed as percent of nsim. (default = 10) 
;       goodt    [OUPUT]  if STORIDX is input, GOODT will containt the log T
;                         range used to plot 
;       tbuff    [INPUT]  if STORIDX is set, TBUFF can be set to the number of
;                         bins to add on either side of the temperature range 
;                         determined with STORIDX
;       monochrom[INPUT] set this to a value between 0 and 255 to make 
;                        all the plots in this color. 
;restrictions
;       The SPLINE and SPLINE INDEX plotting schemes may not be used with slect = 3
;       NOT tested yet with positional parameters e.g. !P.POSITION
;
;subroutines 
;       mid2bound
;       peasecolr
;       stample
;history
;    
;    (LL 6/03) 
;    ADDED 'Nice' and 'Nice Difference' option for schme parameter (LL 2/04)
;    ADDED call to STAMPLE (LL 2/04)
;    ADDED call to PEASECOLR to restore settings after loadct (LL 2/04)
;    ADDED recognition of weather plot background is black or white
;          and and adjustment to col_shft accordingly (LL 2/04) 
;    FIXED SPLINE INDEX schme so that it actually polyfills and doesn't just
;          oplot i.e. no more 'speckled' look. Also added a NOERASE
;          plot statement to prevent polyfills from hiding tickmarks (LL 2/04) 
;    FIXED color schemes so that the background color recognition
;          scheme would work with both psuedo color and true color (LL 2/04)
;    REMOVED PEASECOLR because IDL dynamically updates the x windows
;          display whith the current color table (LL 2/04) 
;    BUGFIX variables CUPBND and CLOBND should be of size nT not nkpt (LL 3/04)
;          if this affected you the routine would have crashed. 
;    ADDED 'NICE INDEX' scheme. (LL 10/04) 
;    CHANGED 'PROB' color schme so best-fit is plotted brightest for
;        'x' and darkest for 'ps', order in which bars plotted matters! (LL 9/05)
;    BUGFIX now compatible with multiple plots per page via !p.multi (LL 10/05)
;    ADDED 'SUPPORT' scheme and keywords MONOCHROM and OUTMED (LL 11/05) 
;    CHANGED col_tabl behavior, only loads col_tabl when set, so default is 
;           now whatever is currently loaded and not B/W Linear (LL 11/05)
;    ADDED keywords DEMSLOPE, SLOPSIG, SLOPEDT and parameter STORIDX (LL 1/06)
;    CHANGED behavior of output to be not in color if COL_TABL is 0 (VK 8/07)
;-

ok='ok' & np = n_params() & sze= size(simdem) & nT = n_elements(logt)
nERR = size(demerr)
if keyword_set(col_tabl) then loadct, col_tabl else col_tabl=0
if keyword_set(ps_fil) then begin
my_device='x'          ;assuming prior output is 'X'
set_plot, 'PS' 
if keyword_set(col_tabl) then device, filename =ps_fil, /color else $
	device, filename=ps_fil, color=0
endif 
if keyword_set(col_rng)  then col_rng  = col_rng<255L  else col_rng  = 255L
if keyword_set(col_shft) then col_shft =fix(col_shft) else col_shft = 0
if !p.background eq 0 then alt_cols=reverse(indgen(256L)) 
if keyword_set(slect) then slect = slect else slect = 2
if n_elements(schme) eq 0 then schme = 'ERRBAR' 
if np lt 3 then ok = 'insufficient parameters' else $ 
 if (schme eq 'SPLINE') and (slect eq 3) then ok = 'SPLINE may not be used with slect = 3' else $
  if (schme eq 'NICE INDEX') and (slect ne 1) then ok='NICE INDEX may only be used with slect = 1' else $
   if (schme eq 'SPLINE INDEX') and (slect eq 3) then ok = 'SPLINE INDEX may not be used with slect = 3' else $
    if (slect gt 3) or (slect le 0) then ok = 'slect keyword selection not valid' else $
     if (col_tabl lt 0 ) then ok =  'col_tabl keyword selection not valid' 
     
if ok ne 'ok' then begin 
   print, 'Usage:  mcmc_plot,logt,simdem,demerr,simprb,schme,storidx,slect=slect,$'
   print, 'col_rng=col_rng,ps_fil=ps_fil,col_tabl=col_tabl,col_shft=col_shft,tbuff=tbuff,$'
   print, 'demslope=demslope,slopsig=slopsig,slopedt=slopedt,goodt=goodt,sampct=sampct,monochrom=monochrom'
   if np ne 0 then message, ok, /info 
   return
endif 

nsim = n_elements(simdem(0,*))-1
demerrl = demerr[*,0]
demerru = demerr[*,1]

;      check stroidx and establish useful temperature range 
if n_elements(storidx) gt 0 then begin 
    strhstx = nT & strhst = hastogram(storidx, strhstx)     
    if n_elements(sampct) eq 0 then cutfrac = 0.10 else cutfrac = sampct(0)/100d
    if n_elements(tbuff) eq 0 then tbuff = 0 
    oo = where(strhst gt nsim*cutfrac and strhstx gt 0) 
    ologt =logt & osimdem = simdem 
    logt = logt[(min(oo)-tbuff)>0:(max(oo)+tbuff)<(nT-1)] & goodt = logt
    simdem = simdem[(min(oo)-tbuff)>0:(max(oo)+tbuff)<(nT-1),*]     
    demerrl = demerrl[(min(oo)-tbuff)>0:(max(oo)+tbuff)<(nT-1)]
    demerru =  demerru[(min(oo)-tbuff)>0:(max(oo)+tbuff)<(nT-1)]
    nT = n_elements(logt)
endif
 TBB = mid2bound(logt)  ;LOGT grid bin boundaries

;initialize container arrays for SPLINE INDEX 
if schme eq 'SPLINE INDEX' then begin 
 splndxu=[simdem-simdem]
 splndxl=[simdem-simdem]
endif

;initialize container arrays for NICE 
if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' or schme eq 'NICE INDEX' then begin 
 meanofeachbin = fltarr(nT) 
 maxofeachbin  = fltarr(nT) 
 minofeachbin  = fltarr(nT) 
endif 

;if not opting for SPLINE interpolation or 'ERRBAR' then  do this
if (schme ne 'SPLINE') and (schme ne 'ERRBAR') then begin 

 ;initial DEM plot to lay basis for polyfill bars
 ; _extra keywords are caught here
 if schme ne 'NICE' and schme ne 'NICE DIFFERENCE' and schme ne 'NICE INDEX' $
     and schme ne 'SUPPORT' then begin
     plot,logt,simdem(*,nsim),psym=10, /ylog,/ynoz,xtitle='logT'$
       ,ytitle = 'DEM', title = 'Best DEM + error bands  ', thick = 5,$
       yrange = [0.8*min(simdem), 1.2*max(simdem)],$
       xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1, /nodata,_extra = e
 endif

;loop through each bin in logT  
for dd = 0, nT-1 do begin 
 ;initialize some variables for polyfill/bar loop according to selection mechanism
     if slect eq 1 then begin 
         upbnd     = max(simdem(dd,*))
         lobnd     = min(simdem(dd,*)) 
         sortbyem  = sort(simdem(dd,*))
         sortdem   = simdem(dd,sortbyem)
         sortprb   = simprb(sortbyem)
     endif
     if slect eq 3 then begin 
         upbnd     = demerru(dd)
         lobnd     = demerrl(dd)
         sortbyem  = sort(simdem(dd,*))
         sortdem   = simdem(dd,sortbyem)
         sortprb   = simprb(sortbyem)
     endif
     if slect eq 2 then begin 
         if keyword_set(clev) then clev = clev else clev = 50
         cutnumber = fix(((clev/100d)<1.0)*nsim)
         cutarray  = simprb(sort(simprb))
         cutvalu   = cutarray(cutnumber)         
         sortndx   = where(simprb le cutvalu)
         tmpo      = simdem(dd,sortndx)      
         tmpoprb   = simprb(sortndx)
         sortbyem  = sort(tmpo) 
         sortdem   = tmpo(sortbyem)
         sortprb   = tmpoprb(sortbyem)
         upbnd     = max(tmpo(sortbyem))   
         lobnd     = min(tmpo(sortbyem)) 
     endif
 ;store initialized arrays   
     ;on first iteration initialize container arrays 
     if dd eq 0 then begin 
         nkpt = n_elements(sortprb) & cupbnd= dblarr(nT) & clobnd= dblarr(nT) 
         csortdem = dblarr(nT,nkpt) & csortprb=csortdem 
     endif 
     ;populate the container arrays 
     cupbnd(dd) = upbnd & clobnd(dd) = lobnd 
     csortdem(dd,*) = sortdem & csortprb(dd,*) = sortprb 
     ;if NICE or NICE DIFFERENCE set 
     if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' or schme eq 'NICE INDEX' then begin 
         meanofeachbin(dd)=median(sortdem)
         maxofeachbin(dd)=upbnd 
         minofeachbin(dd)=lobnd  
     endif ;the NICE/NICE DIFFERENCE setm if 
endfor ; the loop through each bin in logT space loop

;do the nice and nice difference stuff difference stuff
if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' or schme eq 'NICE INDEX' then begin     
  if slect eq 2 then ll=sortndx else ll=findgen(nsim) & cutsamp = simdem(*,ll) 
  ;find best fit of SIMDEM to MEANOFEACHBIN.
  ;use a statistic that takes the squared residuals weighted by (UPBOUND-LOBND)
    mystat = fltarr(n_elements(ll))
  for j = 0, n_elements(ll)-1 do mystat(j)= $
          total(((meanofeachbin-simdem[*,ll[j]])^2)/(maxofeachbin-minofeachbin))
    jnk = min(mystat,kk) 

  if schme eq 'NICE INDEX' then begin
       plot, logt, simdem[*,kk],linestyle=0, /ylog,/ynoz,xtitle='logT',$
       ytitle = 'DEM', title = 'Best DEM + error bands  ', thick = 1,$
       yrange = [0.8*min(simdem), 1.2*max(simdem)],$
       xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1 , _extra = e, psym = 10
  endif
  if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then begin 
       plot, logt, simdem[*,kk],linestyle=0, /ylog,/ynoz,xtitle='logT',$
       ytitle = 'DEM', title = 'Best DEM + error bands  ', thick = 5,$
       yrange = [0.8*min(simdem), 1.2*max(simdem)],$
       xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1 , _extra = e;, psym = 10
    if schme eq 'NICE' then begin 
        oplot, logt, meanofeachbin, linestyle = 1, thick = 3 ;, psym = 10
        oplot, logt, maxofeachbin, linestyle = 2, thick = 3 ;, psym = 10
        oplot, logt, minofeachbin, linestyle = 2, thick = 3 ; , psym = 10
    endif
  endif
endif ; the NICE/NICE DIFFERENCE if

;loop through each bin in logT space (again) and do the pollyfill et.al)
for dd = 0, nT-1 do begin 
 upbnd = cupbnd(dd) & lobnd = clobnd(dd) 
 sortdem = transpose(csortdem(dd,*)) & sortprb = transpose(csortprb(dd,*)) 
 ;prepare the simulated DEM arrays to be passed through the polyfill mechanism
 ;if schme 1 or to and then loop through sim values to polyfill
  if schme eq 'DIFFERENCE' or schme eq 'INDEX' or schme eq 'NICE DIFFERENCE' or schme eq 'NICE INDEX' then begin
         bestndx = where(sortdem eq simdem[dd, nsim]) ;  
         if schme eq 'NICE DIFFERENCE' or schme eq 'NICE INDEX' then bestndx = where(sortdem eq simdem[dd,kk])
         jnk1=min(abs(1-upbnd/sortdem(*)), undx)      ;
         jnk2=min(abs(1-lobnd/sortdem(*)), lndx)      ; 
         erorintu =sortdem(undx)-sortdem(bestndx(0))  ;
         erorintl =sortdem(bestndx(0))-sortdem(lndx)  ;
 ;first loop for sim values below best value
     for qq = lndx, bestndx(0)-1 do begin 

     ;determine the color for each sim according to the schme selected
         if schme eq 'DIFFERENCE' or schme eq 'NICE DIFFERENCE' then clr = fix( (col_rng-col_shft)*$
                                             ((sortdem(bestndx(0))-sortdem(qq)))/erorintl ) + col_shft
         if schme eq 'INDEX' or schme eq 'NICE INDEX' then clr = fix( (col_rng-col_shft)*$
                                             ((bestndx(0) - qq))/(bestndx(0)-lndx) ) + col_shft
         if !p.background eq 0 then clr = alt_cols(clr) 
         polyfill, [TBB(dd),TBB(dd+1),TBB(dd+1),TBB(dd)] $
           , [sortdem(qq),sortdem(qq), sortdem(qq+1), sortdem(qq+1)] $
           , color = clr(0),noclip = 0
     endfor; end loop through sim values -below

   ;then loop for sim values above best value
     for qq = bestndx(0), undx-1 do begin 
         if schme eq 'DIFFERENCE' or schme eq 'NICE DIFFERENCE' then clr = fix( (col_rng-col_shft)*$
                                             ((sortdem(qq)-sortdem(bestndx)))/erorintu ) + col_shft 
         if schme eq 'INDEX' or schme eq 'NICE INDEX' then clr = fix( (col_rng-col_shft)*$
                                             ((qq-bestndx(0)))/(undx-bestndx(0)) ) + col_shft 
         if !p.background eq 0 then clr = alt_cols(clr) 
         TBB = mid2bound(logt)  ; LOGT grid bin boundaries
         polyfill, [TBB(dd),TBB(dd+1),TBB(dd+1),TBB(dd)] $
           , [sortdem(qq),sortdem(qq), sortdem(qq+1), sortdem(qq+1)] $
           , color = clr(0), noclip = 0
     endfor ; end loop though sim values -above 

  endif; the schme eq 'DIFFERENCE' or 'INDEX' if 

  if schme eq 'PROB' then begin
      TBB= mid2bound(logt) 
      ;must plot lower confidence (high prb number) values 
      ;first so we sort again in prb space
      jnk1=min(abs(1-upbnd/sortdem(*)), undx)      
      jnk2=min(abs(1-lobnd/sortdem(*)), lndx)       
      sortdem = sortdem[lndx: undx]
      sortprb = sortprb[lndx: undx]
      sortndx = sort(sortprb)
      sortprb = sortprb(sortndx)
      sortdem = sortdem(sortndx)
      nsort = n_elements(sortndx)
      if !p.background eq 0 then begin
        icols=0>(col_rng-(fix((col_rng-col_shft)*(sortprb-min(sortprb))/(max(sortprb)-min(sortprb)))>0)-col_shft)<255  
        os=reverse(sort(icols))
        for i=0L,nsort-1 do $
        oplot,[TBB(dd),TBB(dd+1)],[replicate(sortdem(os(nsort-1-i)),2)], color=icols(os(nsort-1-i)) , _extra=e
      endif else begin 
        icols=((fix((col_rng-col_shft)*(sortprb-min(sortprb))/(max(sortprb)-min(sortprb)))>0)+col_shft)<255
        os=reverse(sort(icols))
        for i=0L,nsort-1 do $
        oplot,[TBB(dd),TBB(dd+1)],[replicate(sortdem(os(i)),2)], color=icols(os(i)) , _extra=e
      endelse
endif; the schme 'PROB' if 
  if schme eq 'SPLINE INDEX' then begin 
      bestndx = where(sortdem eq simdem[dd, nsim]) 
      jnk1=min(abs(1-upbnd/sortdem(*)), undx)      
      jnk2=min(abs(1-lobnd/sortdem(*)), lndx)       
      sortdem = sortdem[lndx: undx]
      splndxu[dd,findgen(undx-bestndx(0)+1)] = sortdem[bestndx(0):undx]
      splndxl[dd,findgen(bestndx(0)-lndx+1)] = sortdem[lndx:bestndx(0)]
  endif ;the SPLINY INDEX if 

endfor ;endfor loop through LOGT BINS
;oplot the 'best' DEM again over the painted stuff
if schme eq 'NICE INDEX' then  oplot, logt, simdem[*,kk],linestyle=0,psym=10, thick = 1.7
if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then $
    oplot, logt, simdem[*,kk],linestyle=0, thick = 2.0 else begin 
    if schme ne 'SPLINE INDEX' and schme ne 'NICE INDEX' and schme ne 'SUPPORT' then begin 
        oplot,logt,simdem(*,nsim),psym=10,thick=5
        oplot,logt,simdem(*,nsim),psym=10, thick = 5, linestyle=1, color=255
    endif
endelse

endif ;endif the not spliny if
if schme eq 'SPLINE INDEX' then begin 

    TBB= mid2bound(logt) & Tn = [logt, TBB] & sT = sort(Tn) & Tn=Tn(sT)

    ;get the  array w/ largest number of non zero elements (in DEM space)
    tst_gt0u=fltarr(nt)
    tst_gt0l=fltarr(nt)
    for qq=0, nt-1 do begin ;loop through temperature bins
        tst_gt0u(qq) = n_elements(where(splndxu[qq,*] ne 0))
        tst_gt0l(qq) = n_elements(where(splndxl[qq,*] ne 0))
    endfor
    lrgstu=max(tst_gt0u) & lrgstl=max(tst_gt0l)

    ;create array of size [nt,nlargest] and interpolate 
    ;the smaller arrays to an nlargest component array   
    nusplndxu=fltarr(nt,lrgstu) & nusplndxl=fltarr(nt,lrgstl)
    for ww=0, nt-1 do begin 
        nusplndxu[ww,*]=interpol(splndxu[ww,indgen(tst_gt0u(ww))],lrgstu)
        nusplndxl[ww,*]=interpol(splndxl[ww,indgen(tst_gt0l(ww))],lrgstl)
    endfor

    ;SPLINE interpolate to grid NT and polyfill plot
    splndu = 10^( spline(logt, alog10(nusplndxu[*,0]),Tn) )
    for ee=1,lrgstu-1 do begin 
        colru = fix( (col_rng-col_shft)*(ee/lrgstu) ) + col_shft
        if !p.background eq 0 then colru = alt_cols(colru)
        splndudn1=splndu
        splndu=10^( spline(logt, alog10(nusplndxu[*,ee]),Tn) )
          ;need this because Spline interpolation fuzzy:
          tmp = [[splndudn1],[splndu]]
          for j = 0,n_elements(Tn)-1 do splndu(j)=max(tmp(j,*)) 
          for j = 0,n_elements(Tn)-1 do splndudn1(j)=min(tmp(j,*))            
          polyfill,[Tn, reverse(Tn)],[splndudn1,reverse(splndu)],  color = colru,noclip=0
    ; oplot, Tn, splndu, color = colru
    endfor
    splndl = 10^( spline(logt, alog10(nusplndxl[*,0]),Tn) )
    for rr=1,lrgstl-1 do begin
        colrl = fix( (col_rng-col_shft)*((lrgstl-rr)/lrgstl) ) + col_shft
        if !p.background eq 0 then colrl =  alt_cols(colrl)
        splndldn1=splndl
        splndl=10^( spline(logt, alog10(nusplndxl[*,rr]),Tn) )
          ;need this because Spline interpolation fuzzy:
          tmp = [[splndldn1],[splndl]]
          for j = 0,n_elements(Tn)-1 do splndl(j)=max(tmp(j,*)) 
          for j = 0,n_elements(Tn)-1 do splndldn1(j)=min(tmp(j,*))      
        polyfill,[Tn, reverse(Tn)],[splndldn1,reverse(splndl)], color = colrl,noclip=0
;    oplot, Tn, splndl, color = colrl
    endfor

endif; the SPLNIY INDEX if

if schme eq 'SPLINE' then begin 

    if slect eq 1 then sortndx = findgen(nsim) 
    if slect eq 2 then begin 
        if keyword_set(clev) then clev = clev else clev = 50
        cutnumber= fix(((clev/100d)<1.0)*nsim)-1
        cutarray = simprb(sort(simprb))
        cutvalu = cutarray(cutnumber)         
        sortndx = where(simprb le cutvalu)
    endif 

    TBB= mid2bound(logt) & Tn = [logt, TBB] & sT = sort(Tn) & Tn=Tn(sT)
    sortprb=simprb(sortndx)
    sortdem=simdem[*,sortndx]
    icols=fix((col_rng-col_shft)*(sortprb-min(sortprb))/(max(sortprb)-min(sortprb)))+col_shft
    if !p.background eq 0 then icols = alt_cols(icols)
    os=reverse(sort(icols))

    ;initial DEM plot to lay basis for bars
    ; _extra keywords are caught here
    plot,Tn,10^( SPLINE(logt,alog10(simdem(*,nsim)),Tn )), /ylog, /ynoz,xtitle='logT'$
      ,ytitle='DEM', title = 'Best DEM + error bands  ', thick = 5,$
      yrange = [0.8*min(simdem), 1.2*max(simdem(*))],$
      xrange = [min(logt), max(logt)], xstyle = 1, _extra = e

    for i = 0L, n_elements(sortndx)-1L do $
      oplot, Tn,10^( SPLINE(logt,alog10(sortdem(*,os(i))),Tn) ),col=icols(os(i)), thick = 2

endif ;the SPLINY if 

if schme eq 'ERRBAR' then begin 
    plot,logt,simdem(*,nsim),psym=10, /ylog,/ynoz,xtitle='logT'$
      ,ytitle='DEM', title = 'Best DEM + error bands  ', thick = 5,$
      yrange = [0.8*min(simdem), 1.2*max(simdem)],$
      xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1 , _extra = e
    for j = 0,n_elements(logt)-1 do oplot, [logt(j),logt(j)], [demerrl(j), demerru(j)] 

endif ; the ERRBAR if

if schme eq 'SUPPORT' then begin 

    plot,logt,simdem(*,nsim),linestyle=1, /ylog,/ynoz,xtitle='logT'$
      ,ytitle='DEM', title = 'Best DEM + error bands  ', thick = 5,$
      yrange = [0.8*min(simdem), 1.2*max(simdem)],$
      xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1 , _extra = e
   hbnd = logt & for b = 0, nT -1 do hbnd(b) = max(csortdem[b,*]) 
   lbnd = logt & for b = 0, nT -1 do lbnd(b) = min(csortdem[b,*]) 
   oplot, logt, hbnd,linestyle = 1,_extra = e
   oplot, logt, lbnd,linestyle = 1,_extra = e
   if keyword_set(monochrom) then begin 
      polyfill,[logt, reverse(logt)],[lbnd,reverse(hbnd)],  color = monochrom,noclip=0
   endif  

endif 

;do some /noerase /nodata plots to ensure tick marks arn't painted over
     !p.multi[0] = !p.multi[0]+1. > 0      
     plot,logt,simdem(*,nsim),/yl,/ynoz,$
        yrange = [0.8*min(simdem), 1.2*max(simdem)],$
        xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1, /nodata,/noerase,_extra=e
     !p.multi[0] = !p.multi[0]-1. > 0
stample,_extra=e;put PoA stamp on corner of plot
;peasecolr,_extra=e;restore default PoA colors.. no don't do
;this.. because the x window plot is dynamically updated to these
;colors and mahem ensues

 if arg_present(demslope) then begin 
   bestdem = simdem[*,nsim] 
   if n_elements(slopedt) gt 0 then begin 
     nrng = n_elements(slopedt) 
     if nrng gt 1 then rng = [min(slopedt),max(slopedt)] 
     if nrng eq 1 then begin 
       jnk = max(bestdem,mxdx)
       lgtx = where(logt eq slopedt(0))
       if mxdx ne lgtx then  rng = [min(slopedt),logt(mxdx)] else begin 
         jnk = min(bestdem,mndx) & rng = [min(slopedt),logt(mndx)]         
       endelse
     endif
   endif else begin 
       jnk = max(bestdem,mndx)
       if mndx ne 0 then  rng = [min(logt),logt(mndx)] else begin 
         jnk = min(bestdem,mndx) & rng = [min(logt),logt(mndx)]         
       endelse   
   endelse
   shift1=fltarr(nsim+1)
   slope1=fltarr(nsim+1)
   if slect eq 2 then begin 
     if keyword_set(clev) then clev = clev else clev = 50
         cutnumber = fix(((clev/100d)<1.0)*nsim)
         cutarray  = simprb(sort(simprb))
         cutvalu   = cutarray(cutnumber)         
         sortndx   = where(simprb le cutvalu)
         ssimdem = simdem[*,sortndx]
   endif else ssimdem= simdem
   snsim = n_elements(ssimdem[0,*]) 
   shift1=fltarr(snsim+1)
   slope1=fltarr(snsim+1)
   oo = where(logt ge rng(0) and logt le rng(1),noo) 
   if keyword_set(slopenv) or arg_present(slopenv) then slopenv=fltarr(noo,snsim) 
   for j = 0, snsim-1 do begin 
      result = linfit(logt(oo),alog10(ssimdem(oo,j))) 
     if keyword_set(slopenv) or arg_present(slopenv) then $
      slopenv(*,j)=(10^result(0))*((10^logt(oo))^result(1))
     if keyword_set(slpsprd) then $ 
      oplot,logt(oo),(10^result(0))*((10^logt(oo))^result(1)), thick =0.6,linestyle=1 
     shift1(j)=result(0)
     slope1(j)=result(1)
   endfor 
;   mslop = mean(slope1) & mstv = stddev(slope1) 
   x1 = [-1*reverse(findgen(101)/10.0),findgen(101)/10.0] & y = hastogram(slope1,x1)
;   if !d.name eq 'X' then window,/free 
;      plot,[min(x1),x1],[0,y,0],title='mean: '+$
;      strcompress(string(mean(slope1),format='(F10.2)'),/remove_all)$
;      +'+-'+strcompress(string(stddev(slope1),format='(F10.2)'),/remove_all), $
;      xtitle='Slope'+'[range: '+ strmid(string(rng[0],format='(F10.2)'),6)+$
;      '-'+strmid(string(rng[1],format='(F10.2)'),6)+']',ytitle = 'N',psym=10,$
;      xr=[min(slope1),max(slope1)]
     slopsig = stddev(slope1)      
   ;  demslope = mean(slope1) 
     demslope = median(slope1)
     odope = where(slope1 eq demslope)     
     spshft = shift1(odope)
     slopedt = rng

     meanline= (10^spshft(0))*((10^logt(oo))^demslope)
  if keyword_set(slopenv) or arg_present(slopenv) then begin 
      slopenvu = fltarr(noo) & slopenvl = fltarr(noo)
      for zz = 0,noo-1 do slopenvu(zz)=max(slopenv[zz,*]) 
      for zz = 0,noo-1 do slopenvl(zz)=min(slopenv[zz,*]) 
      loadct,4 & peasecolr
      oplot,logt(oo), slopenvu,thick=!p.thick-0.7,linestyle=1,color=61
      oplot,logt(oo), slopenvl,thick=!p.thick-0.7,linestyle=1,color=61
      oplot,logt(oo), meanline,linestyle=0,color=61
      loadct,col_tabl
   endif
      if !d.name eq 'X' then window,/free 
      plot,[min(x1),x1],[0,y,0],title='mean: '+$
      strcompress(string(mean(slope1),format='(F10.2)'),/remove_all)$
      +'+-'+strcompress(string(stddev(slope1),format='(F10.2)'),/remove_all), $
      xtitle='Slope'+'[range: '+ strmid(string(rng[0],format='(F10.2)'),6)+$
      '-'+strmid(string(rng[1],format='(F10.2)'),6)+']',ytitle = 'N',psym=10,$
      xr=[mean(slope1)-4*stddev(slope1),mean(slope1)+4*stddev(slope1)],thick=4,yr =[0,1.05*max(y)]

 endif
if keyword_set(ps_fil) then begin 
    device, /close_file
    set_plot, my_device
endif 
if n_elements(storidx) gt 0 then begin 
logt   = ologt
simdem = osimdem
endif

end 
