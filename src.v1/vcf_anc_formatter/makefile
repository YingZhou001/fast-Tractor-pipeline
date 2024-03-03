
#CC=gcc -W -Wall -Werror -g #-std=c99
CC=gcc -W -g
HTSLIBDIR=/homes8/yingzhou/Tools/pub/htslib/htslib-1.17.install/lib
HTSINCDIR=/homes8/yingzhou/Tools/pub/htslib/htslib-1.17.install/include
LINKPATH= -I${HTSINCDIR} -L${HTSLIBDIR} -lhts -DHTSLIB

all: vcf_anc_formatter

vcf_anc_formatter: extractAncestry.v2.c
	${CC} -o $@ $(LINKPATH) $^
