#include "modmul.h"

// define lookup table for hex->bin
const char hexBin[16][4] = { "0000", "0001", "0010", "0011", "0100", "0101",
  "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110",
  "1111" } ;

// read in hex representation as binary into an array
void readHexToBin ( char *hex, char *bin ) {
  const char *temp;
  int end = 0;
  int i = 0;
  // each character in hex becomes 4 characters in bin
  while( !( hex[i] == '\0' || hex[i] =='\n' ) ) {
  // for ( int i = 0; hex[i+1] != '\0'; ++i ) {
    if (hex[i] >= '0' && hex[i] <= '9') {
      temp = hexBin[ hex[i] - '0' ];
    } else if ( hex[i] >= 'A' && hex[i] <= 'F' ) {
      temp = hexBin[ 10 + hex[i] - 'A' ];
    } else if ( hex[i] >= 'a' && hex[i] <= 'f' ) {
      temp = hexBin[ 10 + hex[i] - 'a' ];
    } else {
      fprintf( stderr,
        "Could not recognize number in hex during conversion: %c\n", hex[i] );
      ++i;
      continue;
    }

    for (int j = 0; j < 4; ++j) {

      // only for testing
      if (4*(IN_BUFF_SIZE-1) + 1 <= 4*i + j)
        fprintf(stderr, "Index out of range in conversion.\n");

      bin[4*i + j] = temp[j];
    }
    end = 4*i+3;
    ++i;
  }

  // append \0 at the end
  bin[end+1] = '\0';
}

// read in line of input with restriction to 256 characters
static int readLine (char *buffer, size_t buffSize) {

  // read buffSize characters; if unsuccessful report to caller
  if ( fgets (buffer, buffSize, stdin) == NULL ) {
        return INPUT_NO;
  }

  // get index of last character
  int lastCh = strlen(buffer) - 1;

  // check for overflow; if last character is not newline overflow occurred
  if ( buffer[lastCh] != '\n' ) {
    char singleton = '0';
    int overflow = 0;

    // clear input stream
    singleton = getchar();
    while( singleton != '\n' && singleton != EOF ) {
      overflow = 1;
      singleton = getchar();
    }

    if ( overflow )
      return INPUT_LONG;
    else
      return INPUT_YES;
  }

  // remove newline from string and return to caller
  buffer[lastCh] = '\0';
  return INPUT_YES;
}

// read input string to mpz_t structure
void readIn ( mpz_t rop, const char *inputStr, int base ) {

  int check = mpz_set_str (rop, inputStr, base);

  // check validity of input
  if ( check == -1 )
  {
    // couldn't read number in
    fprintf( stderr, "Could not import hex into mpz_t.\n" );
  } // else OK
}

// read n-tuple | together with binary representations of exponents
int readTuple ( const int n, mpz_t *reader, int ex1, char *exp1, int ex2, char *exp2 ) {
  // set the buffer for input for 256 characters + new line character
  // remember n-tuple
  char readBuffer[n][IN_BUFF_SIZE];
  // feedback from reader
  int feedback;

  for ( int i = 0; i < n; ++i ) {

    feedback = readLine ( readBuffer[i], sizeof ( readBuffer[i] ) );

    if ( feedback == INPUT_NO ) {
      // Finish reading
      break; // -> DONE
    } else if ( feedback == INPUT_LONG ) {
      // Ignore the rest of input and whole n-tuple
      fprintf( stderr, "Input line too long.\n" );
      break;
    } // else Memorize line
  }

  if ( feedback == INPUT_YES ) {
    // do_operation on numbers with GMP
    for (int i = 0; i < n; ++i) {
      readIn( reader[i], readBuffer[i], INPUT_FORMAT );

      // additionally read in exponent as a binary string
      if ( i == ex1 ) {
        readHexToBin( readBuffer[i], exp1 );
      } else if ( i == ex2 ) {
        readHexToBin( readBuffer[i], exp2 );
      }
    }
  }

  return feedback;
}


// MONTGOMERY
void zn_mont_rho_sq ( mpz_t output, mpz_t modulus ) {
  // initialize to 1
  mpz_set_ui ( output, 1 );
  // initialize variables
  long int l_N = mpz_size ( modulus );

  for ( long int i = 1; i <= 2 * l_N * mp_bits_per_limb; ++i ) {
    mpz_add( output, output, output );
    mpz_mod( output, output, modulus );
  }
}

void zn_mont_omega ( mpz_t output, mpz_t modulus, mpz_t base ) {
  // initialize to 1
  mpz_set_ui ( output, 1 );
  // Global Constant: const int mp_bits_per_limb = w
  for ( int i = 1; i <= mp_bits_per_limb-1; ++i ) {
    mpz_mul ( output, output, output );
    mpz_mod ( output, output, base );
    mpz_mul ( output, output, modulus );
    mpz_mod ( output, output, base );
  }
  mpz_neg( output, output );
  mpz_mod( output, output, base );
}

void zn_mont_mul ( mpz_t output, mpz_t x, mpz_t y, mpz_t N, mpz_t base, mpz_t omega ) {
  // initialize to 0
  mpz_set_ui ( output, 0 );
  long int l_N = mpz_size ( N );

  // initialize temporary variable
  mpz_t u;
  mpz_init( u );

  mp_limb_t *output_0 = output->_mp_d;
  mp_limb_t *x_0 = x->_mp_d;
  mp_limb_t *y_i = y->_mp_d;
  mp_limb_t temp0;
  mp_limb_t temp1;

  for ( long int i = 0; i <= l_N-1; ++i ) {
    output_0 = output->_mp_d;
    x_0 = x->_mp_d;
    y_i = y->_mp_d;
    // assign u
    // y_i * x_0
    // u = y_i[i] * x_0[0]
    temp0 = y_i[i];
    temp1 = x_0[0];
    mpz_set_ui( u, temp0 );
    mpz_mul_ui( u, u, temp1 );
    temp1 = output_0[0];
    mpz_add_ui ( u, u, temp1 );
    mpz_mul ( u, u, omega );
    mpz_mod ( u , u, base );

    // assign output
    mpz_mul( u, u, N );
    mpz_add( output, output, u );
    // define u as y_i * x
    //   put y_i into u
    mpz_set_ui ( u, y_i[i] );
    mpz_mul( u, u, x );
    mpz_add( output, output, u );
    // bitwise right shift with base = 2^ mp_bits_per_limb
    // mpz_cdiv_q ( output, output, base );
    mpz_tdiv_q_2exp ( output, output, mp_bits_per_limb );
  }
  if ( mpz_cmp( output, N ) >= 0 ) {
    mpz_sub( output, output, N );
  }

  mpz_clear( u );
}

// do not use --- testing purpose --- development version
void zn_mont_red ( mpz_t output, mpz_t t, mpz_t N, mpz_t base, mpz_t omega ) {
  // helper variable
  mpz_t u;
  mpz_init2( u, 1024 );
  mpz_t v;
  mpz_init2( v, 1024 );

  mp_limb_t *output_i = output->_mp_d;

  mpz_set( output, t );
  long int l_N = mpz_size ( N );
  for ( long int i = 0; i < l_N-1; ++i ) {
    output_i = output->_mp_d;
    //set u to output_i
    mpz_set_ui ( u, output_i[i] );
    mpz_mul( u, u, omega );
    mpz_mod( u, u, base );

    // set result
    mpz_mul( u, u, N );
    // base^i
    mpz_pow_ui( v, base, i );
    mpz_mul( u, u, v );
    mpz_add( output, output, u );
  }

  // r <- r/b^lN
  // set u to b^Nl
  mpz_pow_ui ( u, base, l_N);
  mpz_cdiv_q( output, output, u );


  if ( mpz_cmp( output, N ) >= 0 ) {
    mpz_sub( output, output, N );
  }
  mpz_clear( u );
  mpz_clear( v );
}

// do not use --- testing purpose --- development version
void zn_mont_exp ( mpz_t result, mpz_t xi, char *exp, mpz_t N, mpz_t base, mpz_t rho_sq, mpz_t omega ) {
  mpz_t t, x, one, temp;

  mpz_init( t ); mpz_init( x ); mpz_init( one );
  mpz_init( temp ); mpz_set_ui( one, 1 );

  zn_mont_mul( t, one, rho_sq, N, base, omega );
  zn_mont_mul( x, xi, rho_sq, N, base, omega );

  int y = strlen( exp );
  for (int i = 0; i < y; ++i) {
    zn_mont_mul( temp, t, t, N, base, omega );
    mpz_set( t, temp );
    if ( exp[i] == '1' ) {
      zn_mont_mul( temp, t, x, N, base, omega );
      mpz_set( t, temp );
    }
  }
  zn_mont_mul( result, t, one, N, base, omega );
  mpz_clear( t ); mpz_clear( x ); mpz_clear( one ); mpz_clear( temp );
}


// implementation of 2k-ary-slide-1exp
void slidingWindow ( mpz_t result, mpz_t base, char *exp, mpz_t mod,
  mpz_t omega, mpz_t b, mpz_t one, mpz_t rho ) {

  // string buffer for bin to int
  char buffer[WINDOW_SIZE];

  int ind = pow(2, WINDOW_SIZE);
  mpz_t lookup[ind/2];

  // temporary variable
  mpz_t temp; mpz_init( temp );
  mpz_t temp1; mpz_init( temp1 );

  // lookup table
  //   base multiplier
  // mpz_powm_ui( b_sq, base, 2, mod );
  mpz_set(temp1, base);
  zn_mont_mul ( temp, base, temp1, mod, b, omega );  

  // do first step outside the loop
  mpz_init( lookup[0] );
  mpz_set( lookup[0], base );

  for ( int i = 1; i < ind/2; ++i ) {
    mpz_init( lookup[i] );
    // precomputed lookup table
    // mpz_powm_ui ( lookup[i], base, (2*i+1), mod );
    // mpz_mul( lookup[i], lookup[i-1], b_sq );
    // mpz_mod( lookup[i], lookup[i], mod );
    zn_mont_mul ( lookup[i], lookup[i-1], temp, mod, b, omega );
  }

  // initialize 'result' to identity element a.k.a. 1
  //   also must be in Montgomery
  // mpz_set_ui ( result, 1 );
  zn_mont_mul( result, one, rho, mod, b, omega );

  // main loop
  int l = 0;
  int u = 0;
  int i = strlen( exp ) - 1;
  int beginning = i;
  while ( i >= 0 ) {
    if ( exp[beginning- i] == '0' ) {
      l = i;
      u = 0;
    } else {
      // l_+
      l = i - WINDOW_SIZE + 1;
      l = ( l < 0 ) ? 0 : l ;
      while ( exp[beginning -l] == '0' ) {
        ++l;
      }
      // must change subset of l from string to integer
      memcpy( buffer, &exp[beginning-i], i-l+1 );
      buffer[(i-l)+1] = '\0';
      // binary number in string to int
      u = strtol( buffer, NULL, 2 );
    }

    // Montgomery's Ladder to defend against side-channel attacks ??????????????

    // mpz_powm_ui ( result, result, (int)pow(2, (i-l+1)), mod );
    for ( int n = 0; n < i-l+1; ++n ) {
      // mpz_powm_ui ( result, result, 2, mod );
      //   MONTGOMERY --- square it
      mpz_set(temp1, result);
      zn_mont_mul ( temp, result, temp1, mod, b, omega );
      mpz_set(result, temp);
    }

    if ( u != 0 ) {
      // mpz_mul ( result, result, lookup[ (int)floor((u-1)/2) ] );
      // mpz_mod ( result, result, mod );
      // MONTGOMERY --- result mod n
      zn_mont_mul ( temp, result, lookup[ (int)floor((u-1)/2) ], mod, b, omega );
      mpz_set( result, temp );
    }
    i = l - 1;
  }

  for ( int i = 0; i < ind/2; ++i ) {
    mpz_clear( lookup[i] );
  } mpz_clear( temp ); mpz_clear(temp1);
}


// implementation of 2k-ary-slide-1exp --- before tuning in for Montgomery
void slidingWindowNOMO ( mpz_t result, mpz_t base, char *exp, mpz_t mod ) {
  // string buffer for bin to int
  char buffer[WINDOW_SIZE];

  int ind = pow(2, WINDOW_SIZE);
  mpz_t lookup[ind/2];

  // Montgomery's Ladder to defend against side-channel attacks

  // lookup table
  //   base multiplier
  mpz_t b_sq; mpz_init( b_sq );
  mpz_powm_ui( b_sq, base, 2, mod );
  // do first step outside the loop
  mpz_init( lookup[0] );
  mpz_set( lookup[0], base );
  for ( int i = 1; i < ind/2; ++i ) {
    mpz_init( lookup[i] );
    // precomputed lookup table
    // mpz_powm_ui ( lookup[i], base, (2*i+1), mod );
    mpz_mul( lookup[i], lookup[i-1], b_sq );
    mpz_mod( lookup[i], lookup[i], mod );
  }
  mpz_clear( b_sq );
  // initialize 'result' to identity element a.k.a. 1
  mpz_set_d ( result, 1 );

  // main loop
  int l = 0;
  int u = 0;
  int i = strlen( exp ) - 1;
  int beginning = i;
  while ( i >= 0 ) {
    if ( exp[beginning- i] == '0' ) {
      l = i;
      u = 0;
    } else {
      // l_+
      l = i - WINDOW_SIZE + 1;
      l = ( l < 0 ) ? 0 : l ;
      while ( exp[beginning -l] == '0' ) {
        ++l;
      }
      // must change sub set of l from string to integer
      memcpy( buffer, &exp[beginning-i], i-l+1 );
      buffer[(i-l)+1] = '\0';
      // binary number in string to int
      u = strtol( buffer, NULL, 2 );
    }

    // mpz_powm_ui ( result, result, (int)pow(2, (i-l+1)), mod );
    for ( int n = 0; n < i-l+1; ++n ) {
      mpz_powm_ui ( result, result, 2, mod );
    }
    if ( u != 0 ) {
      mpz_mul ( result, result, lookup[ (int)floor((u-1)/2) ] );
      // result mod n
      mpz_mod ( result, result, mod );
    }
    i = l - 1;
  }

  for ( int i = 0; i < ind/2; ++i ) {
    mpz_clear( lookup[i] );
  }
}


// generate random seed
int getSeed ( unsigned int bytes, unsigned char *randomBinStr ) {

  // unsigned int bytes = 128;
  // unsigned char seed[bytes];
  int feedback = rdrand_get_bytes( bytes, randomBinStr );

  return feedback;

  unsigned char *random_bytes = random_bytes = malloc((size_t)5 + 1);
  RAND_bytes(random_bytes, 5);
  *(random_bytes + 5) = '\0';
  fprintf(stderr, "%s\n", random_bytes);
  return 0;
}


/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c,
- then write the ciphertext c to stdout.
*/

void stage1() {

  // fill in this function with solution

  // n-tuple to read in
  int n = 3;
  // and output stream
  char *hexOut = NULL;

  // binary representation of exponent
  char binExp[4*(IN_BUFF_SIZE-1) + 1];
  // char *binExp = NULL;

  mpz_t rop[n];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  }
  
  mpz_t omega, rho, base, x_m, one, output, temp;
  mpz_init( omega ); mpz_init( rho ); mpz_init( base ); mpz_init( x_m );
  mpz_init( temp ); mpz_init( output );
  mpz_init2( one, 1024 );

  mpz_set_ui( one, 1 );
  mpz_ui_pow_ui( base, 2, mp_bits_per_limb );

  // for 3-tuple --- N, e, m
  //   last 4 parameters define which one need to be change into binary string
  int inputAvailable = readTuple( n, rop, 1, binExp, -1, NULL );
  while ( inputAvailable == INPUT_YES ) {

    // precompute omega and rho
    zn_mont_rho_sq( rho, rop[0] );
    zn_mont_omega( omega, rop[0], base );
    // get x in Montgomery
    zn_mont_mul( x_m, rop[2], rho, rop[0], base, omega );

    // raise to the power
    // mpz_powm ( output, rop[2], rop[1], rop[0] );
    // binExp = mpz_get_str (binExp, 2, rop[1]);
    slidingWindow( output, x_m, binExp, rop[0], omega, base, one, rho );

    // get back from Montgomery representation
    zn_mont_mul( temp, output, one, rop[0], base, omega );

    // convert to hex back again
    hexOut = mpz_get_str (hexOut, INPUT_FORMAT, temp);
    fprintf( stdout, "%s\n", hexOut );

    // check for another input
    inputAvailable = readTuple( n, rop, 1, binExp, -1, NULL );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  }
  mpz_clear( omega ); mpz_clear( rho ); mpz_clear( base ); mpz_clear( output );
  mpz_clear( one ); mpz_clear( temp ); mpz_clear( x_m );
  free( hexOut ); // free ( binExp );
}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m,
- then write the plaintext m to stdout.
*/

void stage2() {

  // fill in this function with solution

  // n-tuple to read in
  int n = 9;
  // and output stream
  char *hexOut = NULL;
  // string comparison
  int comparison;

  // binary representation of exponent
  char binExp0[4*(IN_BUFF_SIZE-1) + 1];
  char binExp1[4*(IN_BUFF_SIZE-1) + 1];
  // char *binExp = NULL;

  // GMP representation of numbers
  mpz_t rop[n];
  mpz_t output[3];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  } for (int i = 0; i < 3; ++i) {
    mpz_init( output[i] );
  }

  mpz_t omega0, rho0, omega1, rho1, base, one, temp, x_m0, x_m1;
  mpz_t a1, a2, b1, b2;
  mpz_init( a1 ); mpz_init( a2 ); mpz_init( b1 ); mpz_init( b2 );
  mpz_init( omega0 ); mpz_init( rho0 ); mpz_init( omega1 ); mpz_init( rho1 );
  mpz_init( base ); mpz_init( temp ); mpz_init( x_m0 ); mpz_init( x_m1 );
  mpz_init2( one, 1024 );

  mpz_set_ui( one, 1 );
  mpz_ui_pow_ui( base, 2, mp_bits_per_limb );

  // for 9-tuple --- N, d, p, q, d_p, d_q, i_p, i_q and c
  int inputAvailable = readTuple( n, rop, 4, binExp0, 5, binExp1 );
  while ( inputAvailable == INPUT_YES ) {

    // precompute omega and rho
    zn_mont_rho_sq( rho0, rop[2] );
    zn_mont_omega( omega0, rop[2], base );
    zn_mont_rho_sq( rho1, rop[3] );
    zn_mont_omega( omega1, rop[3], base );

    // get x in Montgomery
    zn_mont_mul( x_m0, rop[8], rho0, rop[2], base, omega0 );
    zn_mont_mul( x_m1, rop[8], rho1, rop[3], base, omega1 );

    // use CRT to decrypt message
    //   calculate first part
    // mpz_powm ( output[0], rop[8], rop[4], rop[2] );
    // binExp = mpz_get_str (binExp, 2, rop[4]);
      slidingWindowNOMO ( output[0], rop[8], binExp0, rop[2] );
    // slidingWindow ( temp, x_m0, binExp0, rop[2], omega0, base, one, rho0 );
    // zn_mont_mul( output[0], temp, one, rop[2], base, omega0 );

    //   calculate second part
    // mpz_powm ( output[1], rop[8], rop[5], rop[3] );
    // binExp = mpz_get_str (binExp, 2, rop[5]);
      slidingWindowNOMO ( output[1], rop[8], binExp1, rop[3] );
    // slidingWindow ( temp, x_m1, binExp1, rop[3], omega1, base, one, rho1 );
    // zn_mont_mul( output[1], temp, one, rop[3], base, omega1 );

    //   check which part is bigger
    comparison = mpz_cmp ( output[0], output[1] );
    if ( comparison > 0 ) { //  op0 > op1
      // op0 - op1
      mpz_sub ( output[2], output[0], output[1] );
      // multi
      mpz_mul ( output[2], output[2], rop[7]);
      // mod
      mpz_mod ( output[2], output[2], rop[2]);
      // reconstruct message
      //   multi
      mpz_mul ( output[2], output[2], rop[3] );
      //   add
      mpz_add( output[2], output[2], output[1] );

    } else if ( comparison == 0 ) { // op0 = op1
      // op0 - op1
      mpz_sub (output[2], output[0], output[1]);
      // multi
      mpz_mul ( output[2], output[2], rop[7]);
      // mod
      mpz_mod ( output[2], output[2], rop[2]);
      // reconstruct message
      //   multi
      mpz_mul ( output[2], output[2], rop[3] );
      //   add
      mpz_add( output[2], output[2], output[1] );

    } else { // comparison < 0 | op0 < op1
      // op1 - op0
      mpz_sub (output[2], output[1], output[0]);
      // multi
      mpz_mul ( output[2], output[2], rop[6]);
      // mod
      mpz_mod ( output[2], output[2], rop[3]);
      // reconstruct message
      //   multi
      mpz_mul ( output[2], output[2], rop[2] );
      //   add
      mpz_add( output[2], output[2], output[0] );
    }

    // convert to hex back again
    hexOut = mpz_get_str (hexOut, 16, output[2]);
    fprintf( stdout, "%s\n", hexOut );

    // check for another input
    inputAvailable = readTuple( n, rop, 4, binExp0, 5, binExp1 );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } for (int i = 0; i < 3; ++i) {
    mpz_clear( output[i] );
  }
  mpz_clear( omega0 ); mpz_clear( rho0 ); mpz_clear( omega1 ); mpz_clear( rho1 );
  mpz_clear( base ); mpz_clear( temp ); mpz_clear( one ); mpz_clear( x_m0 );
  mpz_clear( x_m1 );
  free( hexOut ); // free( binExp );
}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2),
- then write the ciphertext c to stdout.
*/

void stage3() {

  // fill in this function with solution

  // n-tuple to read in
  int n = 5;
  // and output stream
  char *hexOut = NULL;

  // binary representation of exponent
  // char binExp[4*(IN_BUFF_SIZE-1) + 1];
  char *binExp = NULL;

  mpz_t rop[n];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  }
  mpz_t output[2];
  for (int i = 0; i < 2; ++i) {
    mpz_init( output[i] );
  }

  mpz_t omega, rho, base, x_m, y_m, z_m, one, temp;
  mpz_init( omega ); mpz_init( rho ); mpz_init( base ); mpz_init( x_m );
  mpz_init( temp ); mpz_init( y_m ); mpz_init( z_m );
  mpz_init2( one, 1024 );

  mpz_set_ui( one, 1 );
  mpz_ui_pow_ui( base, 2, mp_bits_per_limb );

  // char randomState[(IN_BUFF_SIZE-1)*4 + 1];
  mpz_t y; mpz_init( y ); gmp_randstate_t randomState;

  // for 5-tuple
  int inputAvailable = readTuple( n, rop, -1, NULL, -2, NULL );
  while ( inputAvailable == INPUT_YES ) {

    // precompute omega and rho
    zn_mont_rho_sq( rho, rop[0] );
    zn_mont_omega( omega, rop[0], base );
    // get x in Montgomery
    zn_mont_mul( x_m, rop[2], rho, rop[0], base, omega );
    zn_mont_mul( y_m, rop[3], rho, rop[0], base, omega );
    zn_mont_mul( z_m, rop[4], rho, rop[0], base, omega );

    // generate ephemeral key 'y' in range 1 --- q-1  equivalent to(rop[1]-1)
    // mpz_set_ui ( y, 1 ); // testing purpose
    gmp_randinit_default ( randomState );
    // if random to big generate new one
    mpz_urandomm ( y, randomState, rop[1] );

    // converts his secret message m, into an element m, of G

    // compute first part of cipher
    // mpz_powm ( output[0], rop[2], y, rop[0] );
    binExp = mpz_get_str (binExp, 2, y);
    // slidingWindow ( output[0], rop[2], binExp, rop[0] );
    slidingWindow ( temp, x_m, binExp, rop[0], omega, base, one, rho );

    // convert back --- won't be needed any more
    zn_mont_mul( output[0], temp, one, rop[0], base, omega );

    // prepare base for second component
    // mpz_powm ( output[1], rop[3], y, rop[0] );
    // slidingWindow ( output[1], rop[3], binExp, rop[0] );
    slidingWindow ( output[1], y_m, binExp, rop[0], omega, base, one, rho );

    // calculate encryption
    // mpz_mul ( output[1], rop[4], output[1] );
    // mpz_mod( output[1], output[1], rop[0] );
    zn_mont_mul( temp, output[1], z_m, rop[0], base, omega );

    // convert back --- won't be needed any more
    zn_mont_mul( output[1], temp, one, rop[0], base, omega );

    // convert to hex back again
    for (int i = 0; i < 2; ++i) {
      hexOut = mpz_get_str (hexOut, INPUT_FORMAT, output[i]);
      fprintf( stdout, "%s\n", hexOut );
    }

    // check for another input
    inputAvailable = readTuple( n, rop, -1, NULL, -2, NULL );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } for (int i = 0; i < 2; ++i) {
    mpz_clear( output[i] );
  }
  mpz_clear( omega ); mpz_clear( rho ); mpz_clear( base ); mpz_clear( x_m );
  mpz_clear( temp ); mpz_clear( y_m ); mpz_clear( z_m ); mpz_clear( one );
  mpz_clear( y );
  gmp_randclear( randomState );
  free( hexOut ); free( binExp ); 
}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m,
- then write the plaintext m to stdout.
*/

void stage4() {

  // fill in this function with solution

  // n-tuple to read in
  int n = 6;
  // and output stream
  char *hexOut = NULL;

  // binary representation of exponent
  // char binExp[4*(IN_BUFF_SIZE-1) + 1];
  char *binExp = NULL;

  mpz_t rop[n];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  }

  mpz_t omega, rho, base, x_m, y_m, one, temp, output;
  mpz_init( omega ); mpz_init( rho ); mpz_init( base ); mpz_init( x_m );
  mpz_init( temp ); mpz_init( output ); mpz_init( y_m );
  mpz_init2( one, 1024 );

  mpz_set_ui( one, 1 );
  mpz_ui_pow_ui( base, 2, mp_bits_per_limb );

  // for 5-tuple
  int inputAvailable = readTuple( n, rop, -1, NULL, -2, NULL );
  while ( inputAvailable == INPUT_YES ) {

    // precompute omega and rho
    zn_mont_rho_sq( rho, rop[0] );
    zn_mont_omega( omega, rop[0], base );
    // get x in Montgomery
    zn_mont_mul( x_m, rop[4], rho, rop[0], base, omega );
    zn_mont_mul( y_m, rop[5], rho, rop[0], base, omega );

    // calculate shared secret --- avoid division | calculate (p-1)-x as exp
    mpz_sub_ui( output, rop[0], 1 );
    mpz_sub ( output, output, rop[3]);
    // mpz_powm ( output, rop[4], output, rop[0] );
    binExp = mpz_get_str ( binExp, 2, output );
    // slidingWindow( output, rop[4], binExp, rop[0] );
    slidingWindow( temp, x_m, binExp, rop[0], omega, base, one, rho );

    // decrypt
    // mpz_mul ( output, output, rop[5] );
    // mpz_mod ( output, output, rop[0] );
    zn_mont_mul( output, temp, y_m, rop[0], base, omega );

    // get back from Montgomery form
    zn_mont_mul( temp, output, one, rop[0], base, omega );

    // convert to hex back again
    hexOut = mpz_get_str (hexOut, INPUT_FORMAT, temp);
    fprintf( stdout, "%s\n", hexOut );

    // check for another input
    inputAvailable = readTuple( n, rop, -1, NULL, -2, NULL );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  }
  mpz_clear( output ); mpz_clear( omega ); mpz_clear( rho ); mpz_clear( base );
  mpz_clear( x_m ); mpz_clear( y_m ); mpz_clear( one ); mpz_clear( temp );
  free( hexOut ); free( binExp );
}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3();
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }

  // testing stage
  else if ( !strcmp( argv[1], "randomness" ) ) {
    unsigned int size = 128;
    unsigned char seed[size];
    int feedback = getSeed( size, seed );
    for( int i = 0; i < size; ++i ) {
      printf( "%x", seed[i] );
    } printf("\nFeedback: %d\n", feedback);

  char binary[1025];
  binary[1024] = '\0';
  // now bianry representation
  for ( int i = 0; i < size; ++i ) {
    for ( int j = 7; j >= 0; --j ) {
      printf("%x", (seed[i]>>j)&1  );
      binary[i*8 + 7-j] =  (char) ( ( (int)'0' ) + ((seed[i]>>j)&1));
    }
  }
  printf("\n%s\n", binary);
  

  }

  else if( !strcmp( argv[ 1 ], "test" ) ) {

    // x, N
    int n = 2;

    mpz_t rop[n];
    for (int i = 0; i < n; ++i) {
      mpz_init( rop[i] );
    }
    mpz_t rho_sq; mpz_init( rho_sq );
    mpz_t omega; mpz_init( omega );
    mpz_t output; mpz_init( output );
    mpz_t output1; mpz_init( output1 );
    mpz_t one; mpz_init2( one, 1024 );
    mpz_set_ui( one, 1 );
    mpz_t base; mpz_init( base );
    // mpz_set_ui( base, 2^mp_bits_per_limb );
    mpz_ui_pow_ui( base, 2, mp_bits_per_limb );
    int inputAvailable = readTuple( n, rop, -1, NULL, -2, NULL );
    while ( inputAvailable == INPUT_YES ) {

      zn_mont_rho_sq( rho_sq, rop[1] );
      zn_mont_omega( omega, rop[1], base );
      gmp_printf("Rho:\n%Zd\n", rho_sq);
      gmp_printf("Omg:\n%Zd\n", omega);
      zn_mont_mul ( output1, rop[0], rho_sq, rop[1], base, omega );
      zn_mont_mul ( output, one, output1, rop[1], base, omega );

      printf("Comparing...\n");
      if ( mpz_cmp(rop[0], output ) == 0 ) {
        fprintf ( stdout, "Montgomery multiplication OK!\n" );
      }
      if ( mpz_cmp(rop[0], output ) != 0 ) {
        fprintf ( stderr, "Montgomery multiplication not working!\n" );
      }

      inputAvailable = readTuple( n, rop, -1, NULL, -2, NULL );
    }
    for ( int i = 0; i < n; ++i ) {
      mpz_clear( rop[i] );
    } mpz_clear( rho_sq ); mpz_clear( omega ); mpz_clear( output );
    mpz_clear( output1 ); mpz_clear( one );




    // char readBuffer[IN_BUFF_SIZE];
    // int feedback = readLine ( readBuffer, sizeof ( readBuffer ) );
    // char binBuffer[4*(IN_BUFF_SIZE-1) + 1];
    // readHexToBin(readBuffer, binBuffer);

    // fprintf(stdout, "%s\n", binBuffer);


    // readBuffer[0] = 'A';
    // readBuffer[1] = 'b';
    // readBuffer[2] = 'c';
    // readBuffer[3] = 'd';
    // readBuffer[4] = 'E';
    // readBuffer[5] = '\0';

    // char buffer[WINDOW_SIZE];
    // memcpy( buffer, &readBuffer[1], (3-1)+1 );
    // buffer[(3-1)+1] = '\0';

    // int len = strlen(readBuffer);
    // fprintf(stdout, "%s, %d\n", buffer, len);

    // int alla = getSeed ( readBuffer );
    // fprintf(stderr, "%d\n", alla);
  }

  return 0;
}
