function num2str,num,len=len,nmax=nmax, _extra=e
;+
;FUNCTION	num2str
;	given an array of integers, prepends the appropriate number
;	of zeros to make all the strings of the same length, equal
;	to the maximum length needed.
;
;SYNTAX
;	sss=num2str(nnn,len=len,nmax=nmax)
;
;PARAMETERS
;	num	[INPUT; required] array of integers
;
;KEYWORDS
;	len	[INPUT] if set, final string length = LEN.
;		* if LEN < maximum length, defaults to maximum
;		* if LEN < 0, overrides above setting
;	nmax	[INPUT] the maximum possible value of NUM
;		* if set, computes LEN on the fly based on NUM
;		* if both LEN and NMAX are set, LEN takes precedence
;	_extra	[JUNK] here only to prevent crashing the program
;
;EXAMPLE
;	print, num2str(findgen(10),len=-4)
;
;HISTORY
;	written by vinay kashyap (1996?)
;	added keyword NMAX (VK; Jan2000)
;	now takes advantage of the In.n format (VK; Jul2005)
;-

;	usage
if n_params(0) eq 0 then begin
  print, 'Usage: str = num2str(num,len=string_length,nmax=maxnum)'
  print, '  converts integers to equal-length strings'
  return,''
endif

;	input
str = strtrim(num,2)

;	keywords
ll = strlen(str) & lmax=max(ll) & lfnl = lmax
if keyword_set(nmax) then begin
  if nmax(0) gt max([num]) then begin
    lnmax=strlen(strtrim(nmax(0),2)) > lfnl
    lfnl = lnmax
  endif else message,'Input number > alleged maximum! Ignoring.',/info
endif
;
if keyword_set(len) then begin
  lfnl = len
  if len gt 0 and len lt lmax then lfnl = lmax
  if len lt 0 then lfnl = -len
endif

;	make the strings
str = string(num,'(i'+strtrim(lfnl,2)+'.'+strtrim(lfnl,2)+')')
;
;for i=1,lfnl-1 do begin
;  h1 = where(ll eq i)
;  if h1(0) ne -1 then begin
;    prpnd = '' & n1 = 0
;    while n1+i lt lfnl do begin
;      prpnd = prpnd + '0' & n1 = strlen(prpnd)
;    endwhile
;    str(h1) = prpnd + str(h1)
;  endif
;endfor

if lfnl lt lmax then str = strmid(str,0,lfnl)

return,str
end
