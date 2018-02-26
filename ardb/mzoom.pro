function mzoom, img,xfac = xfac, yfac = yfac, factor=factor, _extra 
;+  
;procedure mzoom
;   rebins input image onto new grid and displays it in x window
;   
;5/22/04 
;7/24/04 add option of rebinning x and y differently and 
;        change name so it doesn't conflict with CHIANTI routine
;-

ox = n_elements(img[*,0,0]) 
oy = n_elements(img[0,*,0])

if keyword_set(factor) then xfac = 1.0d 
if keyword_set(factor) then yfac = 1.0d 
if not keyword_set(factor) then factor = 1d
newx = ox*factor*round(xfac)
newy = oy*factor*round(yfac)

newimg = rebin(img,newx,newy) 
;/sample keyword makes it use nearest neigbor sampling 
;rather than default and better but more time consuming 
;bi-linear interpolation
;window, 0
;tv, newimg

return, newimg 

end                

