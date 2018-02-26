pro whee,spoke=spoke,moveit=moveit
;+
;procedure	whee
;	emulates a wheel rolling in place by printing the characters
;	"|", "/", "-", "\", "|", "/", "-", "\" in succession
;
;syntax
;	whee,spoke=spoke,/moveit
;
;parameters	NONE
;
;keywords
;	spoke	[I/O] a +ve integer that on input says at what state the
;		spoke in the wheel is in.  if not given, set to 1.
;		contains this value on output (see MOVEIT)
;		* if > 8, wraps back to 1
;	moveit	[INPUT] if set, increments SPOKE by 1 prior to output
;
;usage
;	IDL> for i=0,3000 do whee,spoke=spoke,/moveit
;
;history
;	vinay kashyap (Feb97)
;-

;	initialize
spk=['|','/','-','\','|','/','-','\'] & mspk=n_elements(spk)
if not keyword_set(spoke) then ist=1L else ist=long(spoke(0))
if ist lt 0 then ist=-ist
if ist gt mspk then begin				;get modulo
  fst=float(ist)/mspk & ist=(ist-long(fst)*mspk)
endif
if ist eq 0 then ist=mspk

;	wheeeeeeeeeeeeeeeeel
print,form="($,2a)",spk(ist-1),string(8b)

;	increment if necessary
if keyword_set(moveit) then spoke=ist+1

return
end
