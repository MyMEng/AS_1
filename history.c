// implementation of 2k-ary-slide-1exp
// before tuning in for Montgomery
void slidingWindow ( mpz_t result, mpz_t base, char *exp, mpz_t mod ) {
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
