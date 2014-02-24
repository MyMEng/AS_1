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

void stage1() {

  // fill in this function with solution

  // set the buffer for input for 256 characters + new line character
  // remember 3-tuple
  char readBuffer[3][IN_BUFF_SIZE];
  int feedback;
  int inputAvailable = 1;

  while ( inputAvailable ) {
    // for 3-tuple --- N, e, m
    for (int i = 0; i < 3; ++i)
    {
      feedback = readLine ( readBuffer[i], sizeof ( readBuffer ) );
      if ( feedback == INPUT_NO )
      {
        // fprintf( stderr, "Couldn't read in a line.\n" );
        // Exit the program ???
        // Or just finish reading
        inputAvailable = 0;
        break;
      }
      else if ( feedback == INPUT_LONG )
      {
        fprintf( stderr, "Input line too long.\n" );
        // Exit the program ???
        // Or just ignore the rest of input
      }
      else if ( feedback == INPUT_YES )
      {
        // fprintf( stdout, "%s\n", readBuffer );
        // Memorize line ??? -> DONE
        // do_operation on numbers with GMP
      }
    }
    // fprintf( stdout, "%s\n", readBuffer[0] );
    // fprintf( stdout, "%s\n", readBuffer[1] );
    // fprintf( stdout, "%s\n", readBuffer[2] );
    // fprintf( stdout, "DONE\n");
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
