function cat_idtresp,str1, str2, abund=abund, onlyrat=onlyrat,_extra = e

; + 
;function  cat_idtresp
;
;    Takes an ID structure from e.g. LINEID() or MUNGE_LIST() and 
;    converts it to a temperature response structure e.g. RSPSTR from 
;    SOLAR_RESP(). Conversely, if a response structure is entered, $
;    it will be converted to an ID structure.  
;
;    A second ID or temperature response structure may be entered as 
;    STR2. This will be converted and/or concatenated with STR1. Conversion 
;    direction is dictated by STR1. 
;
;syntax 
;    tid = cat_idtresp(str1,str2,abund=abund,lnstr=lnstr) 
;
;parameters 
;    str1     [INPUT;required] an ID structure or RSPSTR 
;             e.g. see LINEID or SOLAR_TRESP to be converted. 
;    str2     [INPUT] an ID structure or RSPSTR to be converted
;             as well and concatenated with STR1. 
;
;keywords 
;    abund    [INPUT] abundances relative to H (abund(0)=1)
;	        * abund(Z-1) contains the abundance for element Z
;		* if array size is smaller than the largest present Z,
;		  the last element is replicated to fill the gap
;    onlyrat  [OUTPUT] 
;    _extra   [INPUT ONLY] use this to pass defined keywords to subroutines
;
;restrictions 
;    WILL NOT WORK ON ID STRUCTURES that are products of CAT_IDTRESP()          
;subroutines 
;    CAT_ID 
;    SOLAR_TRESP
;histroy 
;    liwei lin (Sep05)  
;-

;    usage 
ok='ok' & np=n_params() & nid1=n_tags(str1) & nid2=n_tags(str2)

;    check and parse input
if np eq 0 then ok='Insufficient parameters' else $
 if nid1 eq 0 then ok='str1 should be a structure' else begin $
;  if nid2 eq 0 then ok='str2 should be a structure' else begin
    ;     check to see if idstr or rspstr 
    oki='ok' & okr = 'ok' 
    oki2='ok' & okr2 = 'ok' 
    nw1=n_elements(str1.(0))
    if nid1 ne nw1+1L then oki='STR1 in unknown format1'
    str1n = tag_names(str1)  
    jnk1 = where(strmatch(str1n,'INSTR') gt 0, o1)
    jnk2 = where(strmatch(str1n,'BANDPASS') gt 0, o2)
    if o1 le 0 or o2 le 0 then okr='STR1 in unknown format2' 
    if np eq 2 then begin 
       nw2=n_elements(str2.(0))     
       if nid2 ne nw2+1L then oki2='STR2 in unkown format3' 
       str1n = tag_names(str1)  
       jnk1 = where(strmatch(str1n,'INSTR') gt 0, o1)
       jnk2 = where(strmatch(str1n,'BANDPASS') gt 0, o2)
       if o1(0) le 0 or o2(0) le 0 then okr2='STR2 in unknown format4' 
       if oki2 ne 'ok' and okr2 ne 'ok' then ok = oki2 
   endif else begin 
       oki2 = 'none' & okr2 = 'none' 
   endelse 
    if oki ne 'ok' and okr ne 'ok' then ok = oki 
endelse
if ok ne 'ok' then begin 
print, 'Usage: tid = cat_idtresp(str1,str2,abund=abund,lnstr=lnstr)'
 if np ne 0 then message,ok,/info
  if np eq 0 then return,0L
return, -1 
endif

;     abundnaces 
nabu=n_elements(abund)
defabu=getabund('feldman') & abu=defabu
if nabu eq 0 then abu=defabu else begin
  abu(*)=abund(nabu-1L)
  if nabu lt n_elements(defabu) then abu(0L:nabu-1L)=abund else $
   abu(*)=abund(0L:n_elements(defabu)-1L)
endelse
inicon, atom=atom, roman=roman

;     for idstr entered in str1
if oki eq 'ok' then begin 
  for s    = 1,n_elements(str1.wvl) do begin  
     f     = str1.(s).emis 
     logt  = str1.(s).logt  &  wvl  = str1.(s).wvl 
     z     = str1.(s).z & ion  = str1.(s).ion 
     desig = ['-','-'] & econf = desig 
     src   = ['-'] & jon   = ['-'] 
     jnk   = max(str1.(s).flux,nmx) 
     wrange= [min(wvl)-0.01d,max(wvl)+0.01d] 
     lnstr = create_struct('LINE_INT',f,'LOGT',logt,'WVL',wvl,'Z',z,'ION',ion,$
        	  'DESIG',desig,'CONFIG',econf,'SRC',src,'JON',jon)
     nl = strcompress(n_elements(z),/remove_all)  
     mxid  = atom[z[nmx]-1]+' '+roman[ion[nmx]-1]+' '+strcompress(wvl[nmx],/remove_all)+$ 
             ', # of lines: '+nl
     jnk   = solar_tresp('-',mxid,tgrid=logt,lstr=lnstr,$ 
             /nocont,/allah,rspstr=orspstr,wrange=wrange,abund=abu,_extra=e) 
     if s eq 1 then nstr1=orspstr else nstr1=[nstr1,orspstr]       
     if s eq 1 then begin  
           onlyrat1 = replicate('+N1',nl) 
     endif else begin 
           onlyrat1 = [onlyrat1,replicate('+N'+strcompress(s,/remove_all),nl)] 
     endelse
  endfor 
endif
;      for idstr entered in str2
if oki2 eq 'ok' then begin 
  for s   = 1,n_elements(str2.wvl) do begin  
     f    = str2.(s).emis 
     logt = str2.(s).logt  &  wvl  = str2.(s).wvl 
     z    = str2.(s).z & ion  = str2.(s).ion 
     desig  = ['-','-'] & econf = desig 
     src   = ['-'] & jon   = ['-'] 
     wrange = [min(wvl)-0.01d,max(wvl)+0.01d] 
     lnstr = create_struct('LINE_INT',f,'LOGT',logt,'WVL',wvl,'Z',z,'ION',ion,$
        	  'DESIG',desig,'CONFIG',econf,'SRC',src,'JON',jon)
     jnk  = solar_tresp('_',atom[z[0]-1]+' '+roman[ion[0]-1],tgrid=logt,lstr=lnstr,$ 
            /nocont,/allah,rspstr=orspstr,wrange=wrange,_extra=e) 
     if s eq 1 then nstr2=orspstr else nstr2=[nstr2,orspstr]       
     if s eq 1 then begin  
           onlyrat2 = replicate('+N1',nl) 
     endif else begin 
           if oki eq 'ok' then ss=s+n_elements(str1.wvl) else ss=s+n_elements(nstr1.wvl) 
           onlyrat2 = [onlyrat2,replicate('+N'+strcompress(ss,/remove_all),nl)] 
     endelse
  endfor 
endif

;     for rspstr entered in str1 
if okr eq 'ok' then begin 
  nrsp  = n_elements(str1) & wvl = [-1]  
  for j = 0, nrsp-1 do begin  
  ;     check for presensce of effective area and construct wvl array 
     stags = tag_names(str1[j]) & tst = where(strmatch(stags,'EFF*'))
     if tst[0] eq -1 then wvl = [wvl,-1] else begin $ 
       jnk = max(str1[j].effar, nmx)  
       wvl = [wvl, str1[j].wvlar[nmx]] 
     endelse
  ;    construct an ID for each temperature response
     tst = where(strmatch(stags,'ZRESP')) 
     if tst gt 0 then tst2 = total(str1[j].zresp) else tst2=0  
     if tst2 eq 0 then begin 
         nw = 1 
     ;   if zresp not present then use tresp          
         nstr1c = {wvl:-1.0,$
              z:0.,$
              ion:0.,$
              labl:[str1.ldbdir+' '+str1.abref+' '+str1.eqfile,$
                    str1.instr+' '+str1.bandpass],$
              flux:0.,$
              fluxerr:0.,$
              logt:str1.tgrid,$
              emis:str1[j].resp} 
     endif else begin          
     ;   if zresp present then a seperate componenet for each z present 
         nabn = n_elements(str1[j].zresp[0,*]) & tstem = fltarr(nabn) 
         for zz = 0, nabn-1 do tstem[zz]=total(str1[j].zresp[*,zz])  
         oo = where(tstem gt 0,nw) 
         nstr1c = {wvl:replicate(-1.0,nw),$ 
                   z:oo,$
                   ion:replicate(0,nw),$
                   labl:transpose([ [replicate(str1.ldbdir+' '+str1.abref+' '+str1.eqfile,nw)],$
                            [replicate(str1.instr+' '+str1.bandpass,nw)]]),$
                   flux:replicate(0,nw),$
                   fluxerr:replicate(0,nw),$
                   logt:str1.tgrid,$
                   emis:str1[j].resp} 
     endelse
     if j eq 0 then nstr1 = create_struct('ID0',nstr1c) else $ 
     nstr1=creat_struct(nstr1,'ID'+strim(j,2),nstr1c)
     if j eq 0 then begin  
           onlyrat1 = replicate('+N1',nw) 
     endif else begin 
           onlyrat1 = [onlyrat1,replicate('+N'+strcompress(j+1,/remove_all),nw)] 
     endelse
  endfor 
     wvl=wvl[1:*] & nstr1=create_struct('wvl',wvl,nstr1)
endif 
;     for rspstr entered in str2 
if okr2 eq 'ok' then begin 
  nrsp  = n_elements(str2) & wvl = [-1]  
  for j = 0, nrsp-1 do begin  
  ;     check for presensce of effective area and construct wvl array 
     stags = tag_names(str2[j]) & tst = where(strmatch(stags,'EFF*'))
     if tst[0] eq -1 then wvl = [wvl,-1] else begin $ 
       jnk = max(str2[j].effar, nmx)  
       wvl = [wvl, str2[j].wvlar[nmx]] 
     endelse
  ;     construct an ID for each temperature response
     tst = where(strmatch(stags,'ZRESP')) 
     if tst gt 0 then tst2 = total(str2[j].zresp) else tst2=0  
     if tst2 eq 0 then begin 
     ;     if zresp not present then use tresp          
         nstr2c = {wvl:-1.0,$
              z:0.,$
              ion:0.,$
              labl:[str2.ldbdir+' '+str2.abref+' '+str2.eqfile,$
                    str2.instr+' '+str2.bandpass],$
              flux:0.,$
              fluxerr:0.,$
              logt:str2.logt,$
              emis:str2[j].resp} 
     endif else begin          
     ;     if zresp present then a seperate componenet for each z present 
         nabn = n_elements(str2[j].zresp[0,*]) & tstem = fltarr(nabn) 
         for zz = 0, nabn-1 do tstem[zz]=total(str2[j].zresp[*,zz])  
         oo = where(tstem gt 0,nw) 
         nstr2c = {wvl:replicate(-1.0,nw),$ 
                   z:oo,$
                   ion:replicate(0,nw),$
                   labl:transpose([ [replicate(str2.ldbdir+' '+str2.abref+' '+str2.eqfile,nw)],$
                            [replicate(str2.instr+' '+str2.bandpass,nw)]]),$
                   flux:replicate(0,nw),$
                   fluxerr:replicate(0,nw),$
                   logt:str2.tgrid,$
                   emis:str2[j].resp} 
     endelse
     if j eq 0 then nstr2 = create_struct('ID0',nstr2c) else $ 
     nstr2=creat_struct(nstr2,'ID'+strim(j,2),nstr2c)
     if oki eq 'ok' then jj=j+n_elements(str1.wvl) else jj=j+n_elements(nstr1.wvl) 
     if j eq 0 then begin  
           onlyrat2 = replicate('+N'+strcompress(jj+1,/remove_all),nw)
     endif else begin 
           if oki eq 'ok' then jj=j+n_elements(str1.wvl) else jj=j+n_elements(nstr1.wvl) 
           onlyrat2 = [onlyrat2,replicate('+N'+strcompress(jj+1,/remove_all),nw)] 
     endelse     
  endfor 
     wvl=wvl[1:*] & nstr2=create_struct('wvl',wvl,nstr2)
endif 
if np eq 2 then begin 
   if okr eq 'ok' and okr2 eq 'ok' then  outstr = cat_id(nstr1,nstr2,ask='k');do ask dont tell 
   if okr eq 'ok' and oki2 eq 'ok' then  outstr = cat_id(nstr1,str2,ask='k') 
   if oki eq 'ok' and oki2 eq 'ok' then  outstr = [nstr1,nstr2] 
   if oki eq 'ok' and okr2 eq 'ok' then  outstr  = [nstr1,str2] 
endif else begin 
   outstr = nstr1 
endelse
;if np eq 2 then outstr=[nstr1,nstr2] else outstr = nstr1
if np eq 2 then onlyrat=[onlyrat1,onlyrat2] else onlyrat=onlyrat1 
return, outstr
end
