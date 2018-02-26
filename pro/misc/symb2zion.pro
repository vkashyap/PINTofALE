pro symb2zion,symbol,z,ion,verbose=verbose, _extra=e
;+
;procedure	symb2zion
;		convert spectroscopic symbol to atomic number and ionic state
;
;syntax
;	symb2zion,symbol,z,ion
;
;parameters
;	symbol	[INPUT; required] string containing element designation
;		(e.g., "he_2", "He II", "HeII", "he 2", "he", "He", whatever)
;	z	[OUTPUT] atomic number
;	ion	[OUTPUT] ionic state (0=undefined, 1=neutral, 2=singly
;		ionized, etc.)
;
;keywords
;	verbose	[INPUT] higher the number, greater the chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	* cannot handle the "+/-" type notation (e.g., O^+5)
;	* currently coded only up to Es
;
;requires
;	INICON
;	LAT2ARAB
;
;history
;	written by vinay kashyap (Nov96), based on CHIANTI's convertname.pro
;	modified screen dump to also show input (VK; Dec97)
;	added call to INITSTUFF, include elements up to Es (VK; Aug98)
;	corrected bug that was returning garbage for single-letter elements
;	  (VK; Oct98)
;	was returning 0 for ions for, e.g., "O   VII"; corrected (VK; May99)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	added keyword VERBOSE (VK; MarMM)
;	added various trans-Uranium symbols (VK; AugMM)
;	now allow tab separation (VK; Mar03)
;-

np=n_params(0)
if np eq 0 then begin
  print,'Usage: symb2zion,symbol,z,ion,verbose=v'
  print,'  return atomic number and ionic state of given spectroscopic symbol'
  return
endif

v=0 & if keyword_set(verbose) then v=fix(verbose(0)) > 1

;just in case SYMBOL is not defined...
if n_elements(symbol) eq 0 then symbol='X'

;what is the separator?
symb=str_sep(strlowcase(symbol),' ')		;e.g., "He II"
ns=n_elements(symb)
if ns eq 1 then begin
  symb=str_sep(strlowcase(symbol),'	')	;e.g., "He	II"
endif
ns=n_elements(symb)
if ns eq 1 then begin
  symb=str_sep(strlowcase(symbol),'_')		;e.g., "He_2"
endif
if ns gt 2 then begin
  cc=['']
  for i=0L,ns-1L do begin
    if strtrim(symb(i),2) ne '' then cc=[cc,strtrim(symb(i),2)]
  endfor
  if n_elements(cc) gt 1 then symb=cc(1:*)
endif

ns=n_elements(symb)
if ns eq 2 then begin
  name=symb(0) & istate=symb(1)			;e.g., "fe" and "xii"
endif else begin
  slen=strlen(symb(0)) & symarr=strarr(slen)
  for i=0,slen-1 do symarr(i)=strmid(symb(0),i,1)
  k=1 & name=symarr(0) & if slen gt 1 then c1=symarr(1) else c1=''
  case symarr(0) of
    'h': begin					;H
      if slen gt 1 then begin
	if c1 eq 'e' or c1 eq 'o' or c1 eq 'f' or c1 eq 'g' or $
	   c1 eq 's' then begin
	  k=2 & name=name+c1
	endif					;He, Ho, Hf, Hg, Hs
      endif
    end
    'l': begin
      if slen eq 1 then name='li' else begin	;force Li
        if c1 eq 'i' or c1 eq 'a' or c1 eq 'u' or c1 eq 'r' then begin
	  k=2 & name=name+c1
	endif					;Li, La, Lu, Lr
      endelse
    end
    'b': begin					;B
      if slen gt 1 then begin
	if c1 eq 'i' then begin			;Bi
	  k=2 & name=name+c1
	  if v gt 0 then message,'assuming Bismuth, not Boron',/info
	endif
	if c1 eq 'e' or c1 eq 'r' or c1 eq 'a' or $
	   c1 eq 'k' or c1 eq 'h' then begin
	  k=2 & name=name+c1
	endif
						;Be,Br,Ba,Bk,Bh
      endif
    end
    'c': begin					;C
      if slen gt 1 then begin
	if c1 eq 'a' or c1 eq 'l' or c1 eq 'r' or c1 eq 'o' or $
	   c1 eq 'u' or c1 eq 'd' or c1 eq 's' or c1 eq 'e' or $
	   c1 eq 'm' or c1 eq 'f' then begin
	     k=2 & name=name+c1
	   endif				;Ca,Cl,Cr,Co,Cu,Cd,Cs,Ce,Cm,Cf
      endif
    end
    'n': begin					;N
      if slen gt 1 then begin
	if c1 eq 'i' then begin
	  k=2 & name=name+c1			;Ni
	  if v gt 0 then message,'assuming Nickel, not Nitrogen',/info
	endif
	if c1 eq 'e' or c1 eq 'a' or c1 eq 'b' or c1 eq 'd' or $
	   c1 eq 'p' or c1 eq 'o' then begin
	  k=2 & name=name+c1
	endif					;Ne,Na,Nb,Nd,Np,No
      endif
    end
    'o': begin					;O
      if slen gt 1 then begin
        if c1 eq 's' then begin
	  k=2 & name=name+c1
	endif					;Os
      endif
    end
    'f': begin					;F
      if slen gt 1 then begin
        if c1 eq 'e' or c1 eq 'r' or c1 eq 'm' then begin
	  k=2 & name=name+c1
	endif					;Fe, Fr, Fm
      endif
    end
    'm': begin
      if slen eq 1 then name='mg' else begin	;force Mg
        if c1 eq 'g' or c1 eq 'n' or c1 eq 'o' or c1 eq 'd' or $
	   c1 eq 't' then begin
	  k=2 & name=name+c1			;Mg, Mn, Mo, Md, Mt
	endif
      endelse
    end
    'a': begin
      if slen eq 1 then name='al' else begin	;force Al
        if c1 eq 'l' or c1 eq 'r' or c1 eq 's' or c1 eq 'g' or $
	 c1 eq 'u' or c1 eq 't' or c1 eq 'c' or c1 eq 'm' then begin
	  k=2 & name=name+c1			;Al,Ar,As,Ag,Au,At,Ac,Am
	endif
      endelse
    end
    's': begin					;S
      if slen gt 1 then begin
	if c1 eq 'i' then begin
	  k=2 & name=name+c1			;Si
	  if v gt 0 then message,'assuming Silicon, not Sulfur',/info
	endif
	if c1 eq 'c' or c1 eq 'e' or c1 eq 'r' or c1 eq 'n' or $
	   c1 eq 'b' or c1 eq 'm' or c1 eq 'g' then begin
	  k=2 & name=name+c1		
	endif					;Sc,Se,Sr,Sn,Sb,Sm,Sg
      endif
    end
    'p': begin					;P
      if slen gt 1 then begin
	k=2
	if c1 eq 'd' or c1 eq 'r' or c1 eq 'm' or c1 eq 't' or $
	   c1 eq 'b' or c1 eq 'o' or c1 eq 'a' or c1 eq 'u' then begin
	  k=2 & name=name+c1
	endif					;Pd,Pr,Pm,Pt,Pb,Po,Pa,Pu
      endif
    end
    'k': begin 					;K
      if slen gt 1 then begin
	if c1 eq 'r' then begin
	  k=2 & name=name+c1
	endif					;Kr
      endif
    end
    't': begin
      if slen eq 1 then name='ti' else begin	;force Ti
        if c1 eq 'i' or c1 eq 'c' or c1 eq 'e' or c1 eq 'b' or c1 eq 'm' or $
	   c1 eq 'a' or c1 eq 'l' or c1 eq 'h' then begin
	  k=2 & name=name+c1
        endif					;Ti,Tc,Te,Tb,Tm,Ta,Tl,Th
      endelse
    end
    'v':					;V
    'z': begin
      if slen eq 1 then name='zn' else begin	;force Zn
        if c1 eq 'n' or c1 eq 'r' then begin
	  k=2 & name=name+c1		
        endif					;Zn,Zr
      endelse
    end
    'g': begin
      if slen eq 1 then name='ge' else begin	;force Ge
        if c1 eq 'a' or c1 eq 'e' or c1 eq 'd' then begin
	  k=2 & name=name+c1
        endif					;Ga,Ge,Gd
      endelse
    end
    'r': begin
      if slen eq 1 then name='ra' else begin	;force Ra
        if c1 eq 'b' or c1 eq 'u' or c1 eq 'h' or c1 eq 'e' or $
	   c1 eq 'n' or c1 eq 'a' or c1 eq 'f' then begin
	  k=2 & name=name+c1
        endif					;Rb,Ru,Rh,Re,Rn,Ra,Rf
      endelse
    end
    'y': begin					;Y
      if slen gt 1 then begin
        if c1 eq 'b' then begin
	  k=2 & name=name+c1
	endif					;Yb
      endif
    end
    'i': begin					;I
      if slen gt 1 then begin			;In,Ir
        if c1 eq 'n' or c1 eq 'r' then begin
	  k=2 & name=name+c1
	endif
      endif
    end
    'x': begin					;Xe
      if slen gt 1 then begin
	if c1 eq 'e' then begin
	  k=2 & name='xe'
	endif
      endif
    end
    'e': begin
      if slen eq 1 then name='eu' else begin	;force Eu
        if c1 eq 'u' or c1 eq 'r' or c1 eq 's' then begin
	  k=2 & name=name+c1			;Eu,Er,Es
        endif
      endelse
    end
    'd': begin
      if slen eq 1 then begin
	message,'D for Deuterium',/info & name='H'	;force Deuterium
      endif else begin
	if c1 eq 'y' or c1 eq 'b' then begin
	  k=2 & name=name+c1
	endif					;Dy, Db
      endelse
      k=2 & name='dy'
    end
    'w': name='w'				;W
    'u': name='u'				;U
    else: begin
      print,symbol,' -- say again?'
      name=symbol & istate='0' & k=slen
    end
  endcase
  istate=''
  if slen ge k then for i=k,slen-1 do istate=istate+symarr(i)
endelse

;assign atomic numbers
atom=1 & inicon,atom=atom & atom=strlowcase(atom)
oz=where(name eq atom,moz) & z=oz(0)+1

;assign ionic state
c1=strmid(istate,0,1)
if c1 ne 'i' and c1 ne 'v' and c1 ne 'x' and c1 ne 'l' and c1 ne 'c' then begin
  ion=fix(istate) > 0			;already a number
endif else begin
  ion=lat2arab(istate)			;convert from latin
endelse

;trivial science check
if ion gt z+1 then begin
  message,'er..element with atomic number '+strtrim(z,2)+$
  ' cannot be ionized '+strtrim(ion-1,2)+' times',/info
  ion = 0
endif

case np of
  1: message,symbol+':= Z:'+strtrim(z,2)+' ION:'+strtrim(ion,2),/info
  2: message,'('+symbol+') ION:'+strtrim(ion,2),/info
  else:
endcase

return
end
