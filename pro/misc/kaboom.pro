pro kaboom,boom,bang=bang,beep=beep,flash=flash,help=help, _extra=e
;+
;procedure	kaboom
;	produces BOOMs and BANGs
;	(to be used only to indicate catastrophic failure, OK?  don't overdo it)
;
;syntax
;	kaboom,boom,bang=bang,/beep,/flash,/help
;
;parameters
;	boom	[INPUT; default='KABOOM!!'] could be 'BIFF'.  or 'POW'.  could
;		even be an array
;
;keywords
;	bang	[INPUT; default=42] number of "KABOOM"s to produce
;	beep	[INPUT] if set, also produces beeps
;	flash	[INPUT] if set, flashes the window
;	help	[INPUT] prints usage and quits
;	_extra	[JUNK] here only to prevent crashing the program
;
;usage examples
;	kaboom,boom,bang=bang,/beep,/flash,/help
;	kaboom,'?',bang=100,/beep,/flash
;
;side-effects
;	makes an almighty mess of the plot window; may also make noises
;
;history
;	vinay kashyap (rewritten from memory of 1993 program)
;	added keywords FLASH, HELP, and _EXTRA (VK; Nov98)
;	improved color-scale setting for 24-bit consoles (VK; FebMMI)
;	added a brief pause (VK; Apr04)
;-

;	usage
if keyword_set(help) then begin
  print,'Usage: kaboom,boom,bang=bang,/beep,/flash,/help'
  print,'  produces BANGs, BOOMs, and BEEPs'
  return
endif

;	check inputs
ka='KABOOM!!' & nbang=42L
nka=n_elements(boom) & if nka ge 1 then ka=string(boom) & nka=n_elements(ka)
if n_elements(bang) gt 0 then nbang=long(bang(0)) > 1L

;	initialize
dncolors=256. > !D.N_COLORS	;24-bit color screen temporary fix
ika=fix((randomu(seed,nbang)*nka)) < (nka-1L)		;BOOM?
icol=fix((randomu(seed,nbang)*(dncolors-1L)+1L)*!d.n_colors/dncolors)	;color?
iwrap=fix((randomu(seed,nbang)*4+1))			;color table wraps
ithk=randomu(seed,nbang)*5 > (0.5)			;thickness
iori=randomu(seed,nbang)*360				;orientation
isiz=randomu(seed,nbang)*5 > (0.5)			;size

;	set
if !d.window eq -1 then window,0
dx=!x.crange & dy=!y.crange
if dx(0) eq dx(1) or dy(0) eq dy(1) then begin
  plot,[0],/nodata,xs=4,ys=4
  dx=!x.crange & dy=!y.crange
endif
x=randomu(seed,nbang) & y=randomu(seed,nbang)
x=dx(0)+x*(dx(1)-dx(0)) & y=dy(0)+y*(dy(1)-dy(0))

for i=0L,nbang-1L do begin
  if keyword_set(flash) then multi,iwrap(i)
  xyouts,x(i),y(i),ka(ika(i)),color=icol(i),orientation=iori(i),$
	charsize=isiz(i),charthick=ithk(i)
  if keyword_set(beep) then print,format='($,a)',string("7b)
  if i eq 10L*long(i/10L) then wait,0.1	;pause a bit after every 10
endfor
multi,1

return
end
