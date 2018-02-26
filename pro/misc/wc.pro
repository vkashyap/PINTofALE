function wc,files,word=word,char=char,comm=comm,mat=mat,h=h, _extra=e
;+
;function	wc
;	obtains the number of lines, words or characters in the
;	specified file, using the UNIX wc command
;
;syntax
;	lines = wc(files,/word,/char,comm=comment,mat=mat,/h)
;
;parameters
;	files	[INPUT; required] name of file
;
;keywords
;	word	[INPUT] if set, returns number of words
;	char	[INPUT] if set, returns number of characters
;	comm	[INPUT] excludes lines where the first character is
;			the comment character, as specified
;	mat	[INPUT] operate only on text containing matching string
;	_extra	[JUNK] here only to prevent crashing the program
;
;notes
;	if input is an array, the output is also an array, but
;	with one extra element at the end, containing the total
;	number of lines, words, or characters found.
;
;restrictions
;	IDL wrapper to the UNIX wc command, so obviously works only for UNIX.
;
;history
;	vinay kashyap (3/25/93)
;-

if n_params(0) eq 0 then h = 1

if keyword_set(h) then begin
  print, 'Usage: lines=wc(file(s),/word,/char,comm=comment,mat=match,/h)'
  print, '  returns number of lines in file.  if file contains UNIX wildcards,'
  print, '  expands them and returns the output as an array, with the last'
  print, '  element containing the sum of all'
  print, '  /word: returns number of words'
  print, '  /char: returns number of characters'
  print, '  comm="...": excludes lines whose first character is a comment character'
  print, '  mat="...": operate only on lines with matching string'
  return,-1L
endif

sz = size(files) & c1 = files(0)
;if sz(0) eq 0 then c1 = files else c1 = ''
;if sz(0) gt 0 then for i=0,sz(1)-1 do c1 = c1 + files(i) + ' '
if sz(0) gt 0 then for i=1,sz(1)-1 do c1 = c1 + ' ' + files(i)

flag = ' -l '
if keyword_set(word) then flag = ' -w '
if keyword_set(char) then flag = ' -c '
cmnd = 'wc' + flag + c1 + "| awk '{print $1}'"

if keyword_set(comm) then cmnd = 'cat '+c1+' | grep -v "^'+comm+'" | wc'+$
	flag + "| awk '{print $1}'"
if keyword_set(mat) then cmnd = 'grep "'+mat+'" | wc' + flag +$
	"| awk '{print $1}'"

spawn,cmnd,wc & val = long(wc)

sz = size(val) & if sz(1) eq 1 then val = val(0)

return,val
end
