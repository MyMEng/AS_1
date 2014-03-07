#ifndef __MODMUL_H
#define __MODMUL_H

#include  <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <string.h>

#include   <math.h>

// #include  <fcntl.h>
#include <openssl/rand.h>
 #include "rdrand.h"

#include    <gmp.h>

// Reading from stdin constants
#define INPUT_YES 0
#define INPUT_NO 1
#define INPUT_LONG -1

// length of single line
#define IN_BUFF_SIZE 257

// system of input | 16<->HEX
#define INPUT_FORMAT 16

// sliding window exponentiation
#define WINDOW_SIZE 4

// yes / no structure
#define OK 1
#define NO 0

#endif
