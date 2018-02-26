;+
;EG_RD_LIST.PRO
;	example program to exercise RD_LIST.PRO and CAT_LN.PRO
;
;vinay k (1998dec6)
;changed pathnames and added call to getpoadef (2015aug5)
;-

;	initialize SCAR keywords
dbdir='$SCAR'				;line DB directory
n_e=1e10
chifil='ioneq/arnaud_raymond.ioneq'
if not keyword_set(lnlst) then lnlst=0
if not keyword_set(comm) then comm=0

;	list of lines to read in emissivities for:
lnlst=[	'Fe 23	132.831	'+getpoadef()+'/emissivity/emisschianti',$
	'Fe 18	110,112',$
	'C IV	(220,250)	'+getpoadef()+'emissivity/emisspex',$
	getpoadef()+'pro/essempio/eg_rd_list.lst']

;	call RD_LIST
f=rd_list(lnlst,/incieq,dbdir=dbdir,n_e=n_e,chifil=chifil)

;	view results
help,cat_ln(f,comm=comm,lnlst=lnlst)

end
