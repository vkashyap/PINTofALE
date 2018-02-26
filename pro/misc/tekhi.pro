function tekhi, instring,nodollar=nodollar, _extra=e
;+ 
;function chitex
;       Converts spectroscopic terms and configurations to LaTeX
;       form, e.g.: 
;       	'1s.2p 3P2.0'      to '$1s2p \; ^3P_{2}$'
;       	'1s.2p 3P2.5'      to '$1s2p \; ^3P_{5/2}$'
;       	'(2s2.2p5) 2P_3'   to '$(2s^22p^5) \; ^2P_{3}$'
;       	'(2s2.2p5) 2P_3/2' to '$(2s^22p^5) \; ^2P_{3/2}$'
;
;syntax
;	outstring = tekhi(instring) 
;
;parameters 
;       instring [INPUT; required] single string that describes the 
;                electron configurations and/or level designations to be 
;                converted to LaTex form. 
;
;keywords 
;       nodollar [INPUT] if set then output will not be wrapped with '$'s 
;
;	_extra	 [JUNK] here only to prevent crashing the program
;
;subroutines
;          STREGEX
;          STRJOIN
;          STRMID
;
;history
;       LiWei Lin/Jeremy Drake (Nov04)     
;	added keyword _EXTRA and cleaned up usage (VK; Jul05)
;-

;       usage 
ok='ok' & np=n_params() & ns=n_elements(instring)
insz = size(instring) & ninsz=n_elements(insz)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='INSTRING is not defined' else $
  if ns gt 1 then ok='INSTRING cannot be an array' else $
   if insz[ninsz-2] ne 7 then ok='INSTRING is not a string'
if ok ne 'ok' then begin 
 print,'Usage: outstr=tekhi(instr)' 
 print,'  converts spectroscopic terms and configurations to LaTeX form'
 if np ne 0 then message,ok,/informational
 return,instring
endif 

;       initialize parameters
istr = '' 
nc = strlen(instring[0])

;       split string into array of single characters 
for j = 0, nc-1 do istr=[istr,strmid(instring[0],j,1)] & istr= istr[1:*] 

;       check if string is already tekhified
oo = where(stregex(istr,'[\^\$\{\}]') ge 0, ns) 
oo0= stregex(instring[0],'\;') 
if ns gt 0 or oo0(0) ge 0 then begin 
 message, 'Presence of ^${};\ characters indicates that no conversion is neccessary.',/info
 return, instring[0]
endif 

;       numbers after orbital designations need superscript '^' 
oo = where(stregex(istr,'[spdf]') ge 0, ns)
if ns gt 0 then for j = 0, ns-1 do $ 
if stregex(istr(oo[j]+1),'[[:digit:]]') ge 0 then istr(oo[j])=strjoin([istr(oo[j]),'^'])
       
;       numbers preceeding term designations need superscript '^'
oo = where(stregex(istr,'[SPDF]') ge 0, ns) 
if ns gt 0 then for j = 0, ns-1 do istr(oo[j]-1) = strjoin(['^',istr(oo[j]-1)])

;       for upper case characters followed by _ , the number or fraction  
;       immediately following needs to be enclosed in {}
oo = where(stregex(istr,'_') ge 0, ns) 
if ns gt 0 then begin 
  for j = 0, ns-1 do begin  
     ndx = oo[j] 
     if stregex(istr(ndx-1),'[SPDF]') ge 0 then begin 
        istr(ndx)=strjoin([istr(ndx),'{'])
        nrmain = n_elements(istr[ndx:*]) 
        if nrmain ge 4 then begin 
           if stregex(istr(ndx+2),'/') ge 0 then $   
           istr(ndx+3)=strjoin([istr(ndx+3),'}']) else $ 
           istr(ndx+1)=strjoin([istr(ndx+1),'}']) 
        endif else begin 
           istr(ndx+1)=strjoin([istr(ndx+1),'}']) 
        endelse
     endif
  endfor 
endif

oo = where(stregex(istr,'\.') ge 0, ns) 
if ns gt 0 then begin 
  for j = 0, ns-1 do begin 
     ndx = oo[j] 
     if ndx gt 1 and ndx+1 le nc-1 then begin 
       if stregex(istr(ndx-2),'[SPDF]') ge 0 $ 
       and stregex(istr(ndx-1),'[[:digit:]]') ge 0 then begin 
         if istr(ndx+1) eq '5' then begin  
           if istr(ndx-1) eq '4' then istr(ndx-1:ndx+1)=['_{9','/','2}'] 
           if istr(ndx-1) eq '3' then istr(ndx-1:ndx+1)=['_{7','/','2}'] 
           if istr(ndx-1) eq '2' then istr(ndx-1:ndx+1)=['_{5','/','2}']
           if istr(ndx-1) eq '1' then istr(ndx-1:ndx+1)=['_{3','/','2}'] 
           if istr(ndx-1) eq '0' then istr(ndx-1:ndx+1)=['_{1','/','2}'] 
         endif
         if istr(ndx+1) eq '0' then begin 
           istr(ndx:ndx+1)=''  
           istr(ndx-1)=strjoin(['_{',istr(ndx-1),'}']) 
         endif  
       endif else begin 
         istr[ndx]=''     
       endelse 
     endif
  endfor
endif

;       replace all spaces with \;  thick space in mathmode (med thin->\: \,)
oo = where(stregex(istr,' ') ge 0, ns) 
if ns gt 0 then for j = 0, ns-1 do istr(oo[j]) = ' \; ' 

;       wrap with '$'s by default
if not keyword_set(nodollar) then istr=['$',istr,'$'] 
return, strjoin(istr)
end
