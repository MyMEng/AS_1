#include "modmul.h"

// define lookup table for hex->bin
const char hexBin[16][4] = { "0000", "0001", "0010", "0011", "0100", "0101", "0110",
  "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111" } ;

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
  // clean the rest of string ??????????????????????????????????????????????????
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

// how to quickly do exponentiation mod n
void readIn ( mpz_t rop, const char *inputStr, int base ) {

  int check = mpz_set_str (rop, inputStr, base);

  // check validity of input
  if ( check == -1 )
  {
    // couldn't read number in
    fprintf( stderr, "Could not import hex into mpz_t.\n" );
  } // else OK
}

// read n-tuple
// int readTuple ( const int n, mpz_t *reader, const int ex, char *exp ) {
int readTuple ( const int n, mpz_t *reader ) {
  // set the buffer for input for 256 characters + new line character
  // remember n-tuple
  char readBuffer[n][IN_BUFF_SIZE];
  // feedback from reader
  int feedback;

  for ( int i = 0; i < n; ++i ) {

    feedback = readLine ( readBuffer[i], sizeof ( readBuffer[i] ) ); // added [i] in sizeof

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
      // gmp_printf( "%Zd \n", reader[i] );

      // additionally read in exponent as a binary string
      // if ( i == ex ) {
        // readHexToBin(readBuffer[i], exp);
      // }
    }
  }

  return feedback;
}

// implementation of 2k-ary-slide-1exp
void slidingWindow ( mpz_t result, mpz_t base, char *exp, mpz_t mod ) {
  // string buffer for bin to int
  char buffer[WINDOW_SIZE];

  int ind = pow(2, WINDOW_SIZE);
  mpz_t lookup[ind/2];
  // Montgomery's Ladder to defend against side-channel attacks
  for ( int i = 0; i < ind/2; ++i ) {
    mpz_init( lookup[i] );
    // precomputed lookup table
    mpz_powm_ui ( lookup[i], base, (2*i+1), mod );
  }
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

    mpz_powm_ui ( result, result, (int)pow(2, (i-l+1)), mod );
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
  // char binExp[4*(IN_BUFF_SIZE-1) + 1];
  char *binExp = NULL;

  mpz_t rop[n];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  }
  mpz_t output;
  mpz_init( output );

  // for 3-tuple --- N, e, m
  int inputAvailable = readTuple( n, rop );
  // int inputAvailable = readTuple( n, rop, 1, binExp );
  // 1 is parameter which defines which one is exponent
  while ( inputAvailable == INPUT_YES ) {

    // raise to power
    // mpz_powm ( output, rop[2], rop[1], rop[0] );
    binExp = mpz_get_str (binExp, 2, rop[1]);
    slidingWindow( output, rop[2], binExp, rop[0] );

    // convert to hex back again | NOT SAFE ????????????????????????M+ NULL ??
    hexOut = mpz_get_str (hexOut, INPUT_FORMAT, output);
    fprintf( stdout, "%s\n", hexOut );

    // check for another input
    inputAvailable = readTuple( n, rop );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } mpz_clear( output ); free( hexOut ); free ( binExp );
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
  // char binExp0[4*(IN_BUFF_SIZE-1) + 1];
  // char binExp1[4*(IN_BUFF_SIZE-1) + 1];
  char *binExp = NULL;

  // GMP representation of numbers
  mpz_t rop[n];
  mpz_t output[3];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  } for (int i = 0; i < 3; ++i) {
    mpz_init( output[i] );
  }

  // for 9-tuple --- N, d, p, q, d_p, d_q, i_p, i_q and c
  int inputAvailable = readTuple( n, rop );
  while ( inputAvailable == INPUT_YES ) {

    // use CRT to decrypt message
    //   calculate first part
    // mpz_powm ( output[0], rop[8], rop[4], rop[2] );
    binExp = mpz_get_str (binExp, 2, rop[4]);
    slidingWindow ( output[0], rop[8], binExp, rop[2] );
    //   calculate second part
    // mpz_powm ( output[1], rop[8], rop[5], rop[3] );
    binExp = mpz_get_str (binExp, 2, rop[5]);
    slidingWindow ( output[1], rop[8], binExp, rop[3] );

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
      //  IS IT REALLY NEEDED ??????????????????????????????????????????????????
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
    // convert to hex back again | NOT SAFE ????????????????????????M+ NULL ????
    hexOut = mpz_get_str (hexOut, 16, output[2]);
    // gmp_printf( "%Zd \n", output );
    fprintf( stdout, "%s\n", hexOut );

    // check for another input
    inputAvailable = readTuple( n, rop );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } for (int i = 0; i < 3; ++i) {
    mpz_clear( output[i] );
  } free( hexOut ); free( binExp );
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
  mpz_t y; mpz_init( y );
  gmp_randstate_t randomState;

  // for 5-tuple
  int inputAvailable = readTuple( n, rop );
  while ( inputAvailable == INPUT_YES ) {

    // generate ephemeral key 'y' in range 1 --- q-1  equivalent to(rop[1]-1)
    // mpz_set_ui ( y, 1 ); // testing purpose
    gmp_randinit_default ( randomState );
    mpz_urandomm ( y, randomState, rop[1] );

    // converts his secret message m, into an element m, of G

    // compute first part of cipher
    // mpz_powm ( output[0], rop[2], y, rop[0] );
    binExp = mpz_get_str (binExp, 2, y);
    slidingWindow ( output[0], rop[2], binExp, rop[0] );

    // prepare base for second component
    // mpz_powm ( output[1], rop[3], y, rop[0] );
    slidingWindow ( output[1], rop[3], binExp, rop[0] );

    // calculate encryption
    mpz_mul ( output[1], rop[4], output[1] );
    mpz_mod( output[1], output[1], rop[0] );

    // gmp_printf( "%Zd \n", output );
    // convert to hex back again | NOT SAFE ????????????????????????M+ NULL ??
    for (int i = 0; i < 2; ++i) {
      hexOut = mpz_get_str (hexOut, INPUT_FORMAT, output[i]);
      fprintf( stdout, "%s\n", hexOut );
    }

    // check for another input
    inputAvailable = readTuple( n, rop );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } for (int i = 0; i < 2; ++i) {
    mpz_clear( output[i] );
  }  mpz_clear( y ); free( hexOut ); free( binExp );
  gmp_randclear( randomState );
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
  mpz_t output; mpz_init( output );

  // for 5-tuple
  int inputAvailable = readTuple( n, rop );
  while ( inputAvailable == INPUT_YES ) {

    // calculate shared secret --- avoid division | calculate (p-1)-x as exp
    mpz_sub_ui( output, rop[0], 1 );
    mpz_sub ( output, output, rop[3]);
    // mpz_powm ( output, rop[4], output, rop[0] );
    binExp = mpz_get_str (binExp, 2, output);
    slidingWindow( output, rop[4], binExp, rop[0] );

    // decrypt
    mpz_mul ( output, output, rop[5] );
    mpz_mod ( output, output, rop[0] );

    // gmp_printf( "%Zd \n", output );
    // convert to hex back again | NOT SAFE ????????????????????????M+ NULL ??
    hexOut = mpz_get_str (hexOut, INPUT_FORMAT, output);
    fprintf( stdout, "%s\n", hexOut );

    // check for another input
    inputAvailable = readTuple( n, rop );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } mpz_clear( output ); free( hexOut ); free( binExp );
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
  else if( !strcmp( argv[ 1 ], "test" ) ) {
    char readBuffer[IN_BUFF_SIZE];
    // int feedback = readLine ( readBuffer, sizeof ( readBuffer ) );
    // char binBuffer[4*(IN_BUFF_SIZE-1) + 1];
    // readHexToBin(readBuffer, binBuffer);

    // fprintf(stdout, "%s\n", binBuffer);


    readBuffer[0] = 'A';
    readBuffer[1] = 'b';
    readBuffer[2] = 'c';
    readBuffer[3] = 'd';
    readBuffer[4] = 'E';
    readBuffer[5] = '\0';

    char buffer[WINDOW_SIZE];
    memcpy( buffer, &readBuffer[1], (3-1)+1 );
    buffer[(3-1)+1] = '\0';

    int len = strlen(readBuffer);
    fprintf(stdout, "%s, %d\n", buffer, len);
  }

  return 0;
}
