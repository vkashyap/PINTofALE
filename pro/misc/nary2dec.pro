pro nary2dec,num,dec,binary=binary,octal=octal,hex=hex, _extra=e
;+
;procedure	nary2dec
;	converts num from n-ary base to decimal
;
;syntax
;	nary2dec,num,dec,/binary,/octal,/hex
;
;parameters
;	num	[INPUT] n-ary number or string
;	dec	[OUTPUT] corresponding decimal number
;
;keywords
;	binary	[INPUT] if set, NUM taken to be in binary notation
;	octal	[INPUT] if set, NUM taken to be in octal notation
;	hex	[INPUT] if set, NUM taken to be in hexadecimal notation
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (1996?)
;-

hinum = -1 & dec = 0L
if keyword_set(binary) then begin
  hinum = '1111111111' & base = 2. & head = 'BINARY' & hidbl = 14
endif
if keyword_set(octal) then begin
  hinum = '11111' & base = 8. & head = 'OCTAL' & hidbl = 5
endif
if keyword_set(hex) then begin
  hinum = '9' & base = 16. & head = 'HEXADECIMAL' & hidbl = 3
endif
if hinum eq -1 then begin
  decimal = 1 & hinum = '1111111111' & base = 10. & head = 'DECIMAL'
endif

n1 = n_params(0)
if n1 eq 0 then begin
  print, 'Usage: nary2dec,num,dec,/binary,/octal,/hex'
  print, '  converts number in n-ary notation to decimal'
  print, '  input num > ',hinum,' as character strings'
  return
endif

if keyword_set(decimal) then begin
  dec = num & return
endif

nary = strtrim(num,2) & lb = strlen(nary) & nn = strarr(lb)
for i=0,lb-1 do nn(i) = strmid(nary,i,1)
if keyword_set(hex) then begin
  h1 = where(nn eq 'a') & if h1(0) ne -1 then nn(h1) = '10'
  h1 = where(nn eq 'b') & if h1(0) ne -1 then nn(h1) = '11'
  h1 = where(nn eq 'c') & if h1(0) ne -1 then nn(h1) = '12'
  h1 = where(nn eq 'd') & if h1(0) ne -1 then nn(h1) = '13'
  h1 = where(nn eq 'e') & if h1(0) ne -1 then nn(h1) = '14'
  h1 = where(nn eq 'f') & if h1(0) ne -1 then nn(h1) = '15'
endif
nn = fix(nn)

if max(nn) gt base-1 then begin
  print, 'what?  ',nary,' is ',head,'?' & return
endif

if lb le hidbl then base = fix(base) & if lb gt hidbl then dec = 0.D
for i=lb-1,0,-1 do dec = dec + nn(i)*base^(lb-1-i)

if n1 eq 1 then print, dec

return
end
