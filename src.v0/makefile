
#CC=gcc -W -Wall -Werror -g #-std=c99
CC=gcc -W -g

htstag=/home/yzhou3/fast/tools/pub/htslib/build/lib/libhts.so

all: extractAncestry

extractAncestry: extractAncestry.v2.c ${htstag}
	${CC} -o $@ $^
