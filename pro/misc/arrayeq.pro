function arrayeq,A,B,tol=tol, _extra=e
;+
;function	arrayeq
;	checks whether two input arrays are equal, and returns
;	1: if they are identical
;	0: if at least one element differs past tolerance
;	-1: if sizes differ, but one array is a simple subset of the other
;
;syntax
;	AeqB=arrayeq(A,B,tol=tol)
;
;parameters
;	A	[INPUT; required] first array
;	B	[INPUT; required] second array
;		* arrays must be of some >numerical< type
;
;keywords
;	tol	[INPUT; default: 1e-4] tolerance within which two numbers
;		are said to be identical
;	_extra	[JUNK] here only to avoid crashing the program
;
;history
;	vinay kashyap (Nov98)
;	converted to IDL5 notation (VK; OctMM)
;-

;	usage
np=n_params()
if np lt 2 then begin
  print,'Usage: AeqB=arrayeq(A,B,tol=tol)'
  print,'  checks whether input arrays are identical or not'
  return,0L
endif

;	check inputs
nA=n_elements(A) & nB=n_elements(B)
if nA eq 0 or nB eq 0 then return,0L		;missing array(s)?
if n_elements(tol) eq 0 then tol=1e-4

;	array type check
szA=size(A) & szB=size(B) & nszA=n_elements(szA) & nszB=n_elements(szB)
case szA(nszA-2) of
  1: ;byte
  2: ;int
  3: ;long
  4: ;float
  5: ;double
  else: return,0L
endcase
case szB(nszB-2) of
  1: ;byte
  2: ;int
  3: ;long
  4: ;float
  5: ;double
  else: return,0L
endcase

isgn=1L & if nA ne nB then isgn=-1L		;to check if one's a subset

iA=0 & iB=0
nn = nA < nB
if isgn lt 0 then begin			;(check if one's a subset of other
  if nA lt nB then begin		;(first is subset of second
    tmp=min(abs(A[0]-B),iB)
    iB = iB < (nB-nA)
  endif else begin			;)(second is subset of first
    tmp=min(abs(B[0]-A),iA)
    iA = iA < (nA-nB)
  endelse				;NA v/s NB)
endif					;ISGN < 0)

;	find difference between the two arrays
del=total(abs(A[iA:iA+nn-1]-B[iB:iB+nn-1]))

aeq=0L & if del lt tol[0] then aeq=isgn

return,aeq
end
