function is_keyword_set,key
;+
;function	is_keyword_set
;	a wrapper to KEYWORD_SET(), in order to produce the
;	same behavior with vector [0] as with IDL 5.6 and prior
;
;parameters
;	key	[INPUT] the keyword to test for
;
;keywords
;	none	isn't this ironic?
;
;description
;	prior to IDL 5.6, keyword_set(0) returned false and
;	keyword_set([0]) returned true.  this behavior changed
;	in IDL 5.6 to match the documentation of keyword_set,
;	so that both 0 and [0] now return false.  it is far too
;	difficult to change the code logic of all the programs
;	written with old IDLs, so better to have a wrapper
;	that duplicates the same behavior instead.
;
;history
;	Vinay Kashyap (Mar2006)
;-

if n_elements(key) eq 0 then return,0		;this is easy
if (size(key))[0] gt 0 then return,1	;all arrays exist

return,keyword_set(key)				;standard keyword check
end
