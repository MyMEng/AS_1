#include "modmul.h"

/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c,
- then write the ciphertext c to stdout.
*/

// typedef struct {
//   char N[IN_BUFF_SIZE];
//   char e[IN_BUFF_SIZE];
//   char m[IN_BUFF_SIZE];
// } mpz_t;

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
    }
  }

  return feedback;
}





void stage1() {

  // fill in this function with solution

  // n-tuple to read in
  int n = 3;
  // and output stream
  char *hexOut = NULL;

  mpz_t rop[n];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  }
  mpz_t output;
  mpz_init( output );

  // for 3-tuple --- N, e, m
  int inputAvailable = readTuple( n, rop );
  while ( inputAvailable == INPUT_YES ) {

      // raise to power
      mpz_powm ( output, rop[2], rop[1], rop[0] );
      // gmp_printf( "%Zd \n", output );

      // convert to hex back again | NOT SAFE ????????????????????????M+ NULL ??
      hexOut = mpz_get_str (hexOut, INPUT_FORMAT, output);
      fprintf( stdout, "%s\n", hexOut );

      // check for another input
      inputAvailable = readTuple( n, rop );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } mpz_clear( output ); free( hexOut );
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
    mpz_powm ( output[0], rop[8], rop[4], rop[2] );
    //   calculate second part
    mpz_powm ( output[1], rop[8], rop[5], rop[3] );

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
  } free( hexOut );
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

  mpz_t rop[n];
  for (int i = 0; i < n; ++i) {
    mpz_init( rop[i] );
  }
  mpz_t output[2];
  for (int i = 0; i < 2; ++i) {
    mpz_init( output );
  }

  // for 5-tuple
  int inputAvailable = readTuple( n, rop );
  while ( inputAvailable == INPUT_YES ) {

      // raise to power
      mpz_powm ( output, rop[2], rop[1], rop[0] );
      // gmp_printf( "%Zd \n", output );
      // convert to hex back again | NOT SAFE ????????????????????????M+ NULL ??
      hexOut = mpz_get_str (hexOut, INPUT_FORMAT, output);
      fprintf( stdout, "%s\n", hexOut );
      // check for another input
      inputAvailable = readTuple( n, rop );
  }

  for ( int i = 0; i < n; ++i ) {
    mpz_clear( rop[i] );
  } mpz_clear( output ); free( hexOut );

}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m,
- then write the plaintext m to stdout.
*/

void stage4() {

  // fill in this function with solution

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

  return 0;
}
