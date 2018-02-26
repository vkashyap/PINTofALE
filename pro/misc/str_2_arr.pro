function str_2_arr,str,delim,sep=sep,i1=i1,i4=i4,r4=r4,r8=r8,$
	squish=squish,array=array,delimit=delimit,nomult=nomult,$
	_extra=e
;+
;function	str_2_arr
;	converts a series of numbers written out as a string into an array.
;	the output is an array of numbers of specified type
;
;syntax
;	arr=str_2_arr(str,delim,sep=sep,/i1,/i4,/r4,/r8,/squish,$
;	array=array,delimit=delimit,/nomult)
;
;parameters
;	str	[INPUT; required] input string
;		* if itself an array, is concatenated into one huge string
;	delim	[INPUT] same as SEP, put here for consistency with SSWIDL
;		version of program with same name.
;		* if given, is used instead of SEP
;		* if given, the output is returned as a string array
;		  ignoring I1,I4,R4,R8
;		* if not given and DELIMIT is not set, and I1,I4,R4,R8
;		  are also not set, the output is an integer array
;
;keywords
;	sep	[INPUT] separator between entries in the string
;		* if not given as a string, first checks "," and
;		  if that doesn't seem to make any difference, then
;		  checks "<sp>"
;		* if explicitly set, checks only for SEP
;		* if array, uses only the first element
;	i1	[INPUT] if set, output is a BYTARR, default is I2
;	i4	[INPUT] if set, output is a LONARR
;	r4	[INPUT] if set, output is a FLTARR
;	r8	[INPUT] if set, output is a DBLARR
;		* precedence is R8 over R4 over I4 over I1 over I2
;	squish	[INPUT] if set, "squishes" adjoining SEPs, e.g.,
;		"1      2" is converted to "1 2"
;	array	[OUTPUT] string array version of the output
;	delimit	[INPUT] delimiter, same as DELIM, and ignored if either
;		DELIM or SEP are present
;		* but if set, the output is a string array, and I1,I4,R4,R8
;		  are ignored
;	nomult	[INPUT] same as SQUISH, just for compatibility with SSWIDL
;		version
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	SYZE
;
;a word
;	yes, I know about READS.  its limitations are that you still have
;	to define the input array first, cannot use other values for the
;	separator (e.g., ":"), and can't handle missing data (e.g., "1,2,")
;
;history
;	vinay kashyap (Nov98)
;	added parameter DELIM and keywords ARRAY, DELIMIT, and NOMULT
;	  for the sake of consistency with SSWIDL version; also added
;	  call to SYZE() (VK; Jun02)
;	changed name from STR2ARR to STR_2_ARR to avoid conflict with
;	  SSW routine because of compilation issues with SSW (VK; Apr05)
;-

;	usage
ns=n_elements(str)
if ns eq 0 then begin
  print,'Usage: arr=str_2_arr(str,delim,sep=sep,/I1,/I4,/R4,/R8,/squish,$'
  print,'       array=array,delimit=delimit,/nomult)'
  print,'  convert string to array'
  return,0
endif

;	SEP, DELIM, or DELIMIT?
if keyword_set(delimit) then sp=delimit(0)
if keyword_set(delim) then sp=delim(0) else $
 if keyword_set(sep) then sp=sep(0)

;	check input
ss='' & for i=0,ns-1 do ss=ss+string(str(i))

;	separate the output
szsep=size(sp) & nszsep=n_elements(szsep)
bss=byte(ss) & nbs=n_elements(bss)
if szsep(nszsep-2) eq 7 then begin		;(separator is a string
  bsep=(byte(sp(0)))(0)
  seplen=strlen(sp)
endif else begin				;SP is string)(not string
  bsep=(byte(','))(0)
  os=where(bss eq bsep,mos)
  if mos eq 0 then bsep=(byte(' '))(0)		;"," didn't work, try "<sp>"
  seplen=1
endelse						;SP is ignored)

if seplen gt 1 then begin		;(slow version, goes thru STR_SEP
  array=str_sep(ss,sp(0))
endif else begin			;)(faster?
  os=where(bss eq bsep,mos)
  if mos eq 0 then os=nbs	;all one string
  array=strarr(mos+1L)
  if os(0) gt 0 then array(0)=string(bss(0:os(0)-1L))
  for i=1L,mos-1L do $
    if os(i-1)+1L le os(i)-1L then array(i)=string(bss(os(i-1)+1L:os(i)-1L))
  if mos gt 0 then begin
    if os(mos-1)+1L le nbs-1L then array(mos)=string(bss(os(mos-1)+1L:*))
  endif
endelse					;SEPLEN)

;	how many elements in output?
narr=n_elements(array)

;	what type is output?
arr=strarr(narr) & I2=0
if keyword_set(R8) then arr=dblarr(narr) else $
 if keyword_set(R4) then arr=fltarr(narr) else $
  if keyword_set(I4) then arr=lonarr(narr) else $
   if keyword_set(I1) then arr=bytarr(narr) else $
    if not keyword_set(delim) and not keyword_set(delimit) then begin
      I2=1 & szarr=syze(array) & nszarr=n_elements(szarr)
      if szarr(nszarr-2L) gt 5 or szarr(nszarr-2L) lt 0 then I2=0
      if keyword_set(I2) then arr=intarr(narr)
    endif

;	"squish" if asked
if keyword_set(nomult) or keyword_set(squish) then begin
  oo=where(array ne '',moo)
  if moo eq 0 then return,arr(0)
  array=array(oo) & arr=arr(oo) & narr=moo
endif else begin
  oo=where(array eq '',moo)
  if moo gt 0 then array(oo)='0'
endelse

;	extract output array
for i=0L,narr-1L do $
  if keyword_set(R8) then arr(i)=double(array(i)) else $
   if keyword_set(R4) then arr(i)=float(array(i)) else $
    if keyword_set(I4) then arr(i)=long(array(i)) else $
     if keyword_set(I1) then arr(i)=(byte(fix(array(i))))(0) else $
      if keyword_set(I2) then arr(i)=fix(array(i)) else $
       arr(i)=array(i)

return,arr
end
