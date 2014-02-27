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


void stage1() {

  // fill in this function with solution

  // set the buffer for input for 256 characters + new line character
  // remember 3-tuple
  char readBuffer[3][IN_BUFF_SIZE];
  int feedback;
  int inputAvailable = 1;
  int threeAgrees = 0;

  mpz_t rop[3];
  for (int i = 0; i < 3; ++i)
  {
    mpz_init( rop[i] );
  }

  while ( inputAvailable ) {
    // for 3-tuple --- N, e, m
    threeAgrees = 0;
    for (int i = 0; i < 3; ++i)
    {
      feedback = readLine ( readBuffer[i], sizeof ( readBuffer ) );
      if ( feedback == INPUT_NO )
      {
        // fprintf( stderr, "Couldn't read in a line.\n" );
        // Exit the program ???
        // Or just finish reading
        threeAgrees = 0;
        inputAvailable = 0;
        break;
      }
      else if ( feedback == INPUT_LONG )
      {
        // fprintf( stderr, "Input line too long.\n" );
        // Exit the program ???
        // Or just ignore the rest of input
        // ignore whole 3-tuple
        threeAgrees = 0;
      }
      else if ( feedback == INPUT_YES )
      {
        // Memorize line ??? -> DONE
        threeAgrees = 1;
      }
    }

    if ( threeAgrees )
    {
        // fprintf( stdout, "%s\n", readBuffer[0] );
        // fprintf( stdout, "%s\n", readBuffer[1] );
        // fprintf( stdout, "%s\n", readBuffer[2] );

        // do_operation on numbers with GMP
        for (int i = 0; i < 3; ++i)
        {
          readIn( rop[i], readBuffer[i], 16 );
          gmp_printf( "%Zd \n", rop[i] );
        }

    }

  }

  // sanity check
  // if you want check what is inside

  for ( int i = 0; i < 3; ++i )
  {
    mpz_clear( rop[i] );
  }

}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m,
- then write the plaintext m to stdout.
*/

void stage2() {

  // fill in this function with solution

}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2),
- then write the ciphertext c to stdout.
*/

void stage3() {

  // fill in this function with solution

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
