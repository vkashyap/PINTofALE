pro mc_eror, g, sims, bndu, bndl, confidence=confidence, alas=alas, _extra = e 

;+
;procedure mc_eror 
;
;          simple procedure that calculates confidence bounds 
;          about a given parameter value g given a certain number 
;          of monte carlo simulations contained in array sims. 
;          assumes the simulations are distributed normally about 
;          the 'candidate' parameters.
;
;          by default,set to find smallest possible bounds. can also
;          be set to find confidence bounds by counting  upward 
;          and downwards by required fraction.        
;
;syntax
;           mc_eror, g, sims, erru, errl, confidence=confidence 
;
;parameters
;
;         g     [input, required] parameter about which to estimate bounds
;         sims  [input, required] array of monte carlo simulated g's 
;                                 minimum 10 simulations [though this 
;                                 seems a ridiculous number]
;         bndu  [output]          upper confidence bound
;         bndl  [output]          lower confidence bound 
; 
;keywords 
; 
;         confidence [default = 0.68269] confidence limit at which bounds
;                                    are desired. 
;         alas                    find confidence bounds by counting  upward 
;                                 and downwards by required
;                                 fraction. may also be set 
;                                 to find smallest possible bounds.               
;
;history 
;         LL Jan03  
;-        LL Apr03 add smallest interval mechanism and turn old into
;           frac keyword
;
ok='ok' & ng=n_elements(g) & nmc=n_elements(sims)
if ng eq 0 then ok='No g defined' else $
 if nmc lt 1 then ok='Please, more simulations.' 
if ok ne 'ok' then begin
  print,'Usage:mc_eror, g, sims, mcup, mcdown, confidence=confidence'  
if ok ne 'ok' then message,ok,/info
  return
endif
             
if keyword_set(confidence) then begin
  if (confidence gt 1.0) then begin 
    ok= 'Cannot be more than 100 percent confident unless W says so. ugh.'
    message, ok, /info
    return
  endif
conf=confidence/2 
endif else begin 
 conf=0.68269/2
endelse 

szsim = size(sims)
if (szsim(0) gt 1) and (szsim(1) eq 1) then sims = transpose(sims)

if keyword_set(alas) then begin 
  auxar = [g, sims] &  argg  = auxar(sort(auxar)) & a0w   = where(argg eq g) 
  a0w = median(a0w) & intt  = round(nmc*conf)    & nmcc = nmc-1 ;IDL index corr. 
   if (a0w + intt le nmcc) and (a0w - intt ge 0) then begin 
      erru = argg(a0w+intt) & errl = argg(a0w-intt)
   endif
   if (a0w + intt ge nmcc) then begin
     erru = argg(nmcc) & errl = argg(a0w-intt-(intt-(nmcc-a0w)))        
   endif
   if (a0w - intt le 0) then begin 
     errl = argg(0) & erru = argg(2*intt)
   endif
endif else begin 
  auxar = [g, sims] & argg = auxar(sort(auxar)) & a0w   = where(argg eq g) 
  a0w = median(a0w) & intt = round(nmc*conf) & nmcc = nmc;-1(inx corr)+1(g added)
  btbd = a0w - 2*intt > 0     ;lowest possible lower bound in IDL index space
  tpbd = a0w + 2*intt < nmcc  ;highest possible high bound in IDL index space
  testint = fltarr((tpbd-btbd-2*intt)+1) ;tricky.... 
  for j = btbd, tpbd-2*intt do testint(j-btbd) = argg(j+2*intt)-argg(j);calculate width of each interval
    jj = where(testint eq min(testint)); identify smallest width
                                       ; and corresponding bounds
    bndu = argg(jj+2*intt+btbd) & bndl = argg(jj+btbd)  
endelse  

;for j = 0, a0w do testint(j) = (auxar(j+2*intt)-auxar(j)) ;
;
; if (a0w - intt lt 0) then begin 
;    testint = fltarr(a0W)
;    for j = 0, a0w do  testint(j) = (auxar(j+2*intt)-auxar(j))
;     jj = where(testint eq min(testint))
;     erru = auxar(jj+2*intt) & errl = auxar(jj) 
; endif 

; if (a0w + intt gt nmcc) then begin 
;    testint = fltarr(nmcc-a0w)
;    for j = a0w, nmcc do begin 
;     testint(j-a0w) = (auxar(j)- auxar(j-2*intt))
;    endfor
;     jj = where(testint eq min(testint))
;     erru = auxar(aow+jj) & errl = auxar(a0w+jj-2*intt) 
; endif

; if (a0w + intt le nmcc) and (a0w - intt ge 0) then begin 
;     testint = fltarr(2*intt)
;     for j = a0w-intt, a0w+intt do begin 
;     testint(j) = (auxar(i)- auxar(i-2*intt))
;    endfor
;     jj = where(testint eq min(testint))
;     ii = nmcc - jj
;     erru = auxar(ii) & errl = auxar(ii - 2*intt)  
; endif 
end 

