 CC        = /usr/local/Cellar/gcc48/4.8.2/bin/gcc-4.8 #/usr/bin/gcc
 CC_FLAGS  = -Wall -m64 -march=core2 -O3 -std=gnu99
 CC_PATHS  = 
 CC_LIBS   =

# uncomment to use the older, default GMP installation
# CC_PATHS +=
# CC_LIBS  +=              -lgmp

# uncomment to use the newer, bespoke GMP installation
 CC_PATHS += -I/usr/local/Cellar/gmp/5.1.3/include#/usr/local/gmp505/include/ 
 CC_PATHS += -L/usr/local/Cellar/gmp/5.1.3/lib#/usr/local/gmp505/lib/
 CC_LIBS  += #-Wl,-Bstatic -lgmp -Wl,-Bdynamic

all    : modmul

modmul : $(wildcard *.[ch])
	@${CC} ${CC_PATHS} ${CC_FLAGS} -o ${@} $(filter %.c, ${^}) ${CC_LIBS}

clean  : 
	@rm -f core modmul
