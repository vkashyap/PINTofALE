function contuti,y,mass=mass,vmink=vmink, _extra=e
;+
;function        contuti
;
;        returns the convective turnover time in days given a stellar
;        mass or a V-K colour.  Uses the Wright et al 2018,
;        MNRAS,479,2351 relations between V-K and mass. 
;        Input mass can be a vector or a scalar. 
;        MIST theoretical turnover times will be added as an option in 
;        the future, as will the option to give a spectral type.
;
;syntax
;        ctt=contuti(y,/mass,/vmink)
;
;parameters
;        y        [Input; required] vector or scalar of parameter(s)
;                 for which ch convective turnover times are to be calculated   
;                 * If /mass is set, assumed to be stellar mass [Default]
;                 * If /vmink is set, assumed to be V-K colours
;                 Note that mass keyword is currently redundant, but
;                 is expected to be used in future updates and
;                 expansion of options.
;
;keywords
;                
;        mass      [INPUT] specify that y is in units of solar masses
;        vmink     [INPUT] specify that y is V-K colour
;       _extra     [INPUT] junk -- here only to prevent program from crashing
;
;restrictions
;
;        The convective turnover times are only valid within the range of
;        masses for which data are available, and for unevolved
;        stars. These ranges are 0.08 < M/Msun < 1.36 and 1.1 < V-K < 7.0.
;        T Tauri stars are also not yet covered, but should be in the future
;        when MIST data are included.  Routine exits with error
;        message if input data are out of range.
;
;history
;        Jeremy Drake (5/4/19)
;-

;----------------------------------------------------------

;       usage
  
ny=n_elements(y)
if ny eq 0 then begin
  print,'Usage: convtime=contuti(y,mass=mass,vmink=vmink)
  print,'  return convective turnover times in days for each y[i]'
  return,[-1L]
endif

; Check if mass or V-K input and exit if input is outside validity range.
; Valid ranges are  0.08 < M/Msun < 1.36 and 1.1 < V-K < 7.0.

; If mass input specified (currently redundant), or no keywords set
; and mass assumed:

if keyword_set(mass) or ((keyword_set(mass) eq 0) and $
    (keyword_set(vmink) eq 0)) then begin
  if (min(y) lt 0.08 or max(y) gt 1.36) then begin
    print,'Input masses outside of valid range 0.08 < M/Msun < 1.36; returning'
    return,[-1L]
  endif else begin
    tauc=2.33-1.5*y + 0.31*y^2.  ; this is actually log10 tau
    tauc=10.^tauc
  endelse

; If V-K input specified:

endif else if keyword_set(vmink) then begin
  if (min(y) lt 1.1 or max(y) gt 7.0) then begin
    print,'Input V-K outside of valid range 1.1 < V-K < 7.0; returning'
    return,[-1L]
  endif else begin
    tauc=0.64+0.25*y    ; this is actually log10 tau
    tauc=10.^tauc
  endelse
endif

print,'mass keyword',keyword_set(mass)

return,tauc
end
