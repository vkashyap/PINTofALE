# USAGE: make arch=system where system is one of 
#          (osx=darwin, linux, solaris, sun, alpha, dec)

SHELL = /bin/sh
CC = gcc

SITEOBJ=	pimms.o \
		pms_index.o \
		pms_intlz.o \
		pms_rarea.o \
		pms_slmdl.o \
		pms_specl.o \
		pms_spec2.o \
		sis_obstm.o

OBJ=	pimms.o \
	pms_intlz.o \
	pms_slmdl.o \
	pms_selct.o \
	pms_docrt.o \
	pms_outpt.o \
	pms_rarea.o \
	pms_specl.o \
	pms_spec2.o \
	pms_write.o \
	pms_index.o \
	pms_afdir.o \
	pms_colmn.o \
	acis_pileup.o \
	bat_limit.o \
	epic_pileup.o \
	int_limit.o \
	sis_limit.o \
	sis_obstm.o \
	gis_limit.o \
	hxd_limit.o \
	hxt_limit.o \
	pca_limit.o \
	ha4_levin.o \
	nustar_limit.o \
	sax_limit.o \
	sxs_grades.o \
	sgd_limit.o \
	xrs_pileup.o \
	dcd_rrang.o \
	rd_real.o \
	simp32.o \
	spec.o \
	tranmm.o \
	tranuvo.o \
	intcom.o \
	getcom.o \
	set_param.o \
	get_param.o \
	force_comm.o \
	ark_help.o \
	comm_lib.o \
	sys.o \
	lentrim.o \
	upcase.o \
	getlun.o

all : linklib
	make ahmake xsing pimms

linklib :
	@if [ '$(arch)' = 'sun' ]; \
	then \
		rm -f sys.f; ln sun.sys sys.f; $(FC) -c sys.f; \
	elif [ '$(arch)' = 'solaris' ]; \
	then \
		rm -f sys.f; ln solaris.sys sys.f; $(FC) -c sys.f; \
	elif [ '$(arch)' = 'dec' ]; \
	then \
		rm -f sys.f; ln ultrix.sys sys.f; $(FC) -c sys.f; \
	elif [ '$(arch)' = 'alpha' ]; \
	then \
		rm -f sys.f; ln osf.sys sys.f; $(FC) -c sys.f; \
	elif [ '$(arch)' = 'linux' ]; \
	then \
		rm -f sys.f; ln linux.sys sys.f; $(FC) -c sys.f; \
	elif [ '$(arch)' = 'osx' ]; \
	then \
		rm -f sys.f; ln macosx.sys sys.f; $(FC) -c sys.f; \
	else \
		echo "Architecture not supported: arch=$(arch)";  exit 1 ; \
	fi

ahmake :	ahmake.o \
		sys.o \
		lentrim.o \
		upcase.o \
		getlun.o
		$(FC) -o ahmake ahmake.o sys.o lentrim.o upcase.o getlun.o
		./ahmake
		mv pimms.ahl ../data

xsing :	xsing.o \
	sys.o \
	lentrim.o \
	upcase.o \
	getlun.o
	$(FC) -o ../models/xsing xsing.o sys.o lentrim.o upcase.o getlun.o

pimms :	$(OBJ)
	$(FC) -o ../pimms $(OBJ)

clean :
	rm -f *.o sys.f ahmake

$(SITEOBJ):	./sitespec.inc
