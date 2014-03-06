 
  // FILE *urandom = fopen ( "/dev/urandom", "r" );
  // FILE *random = fopen ( "/dev/random", "r" );
  // if ( urandom == NULL || random == NULL ) {
    // fprintf ( stderr, "Could not open random source in /dev\n" );
    // exit ( 1 );
    // return 1;
  // } 

  // read one byte at a time (8-bits) while not done
  int quantity = 5; //((IN_BUFF_SIZE-1)*4)/8;
  // char randomBytes[quantity];
  void *randomBytes = malloc(quantity);
  int isRandom = OK;
  int readSize;
  // do 10 reading trials for each source

  printf("Starting loop.\n");

  for ( int i = 0; i < 20; ++i ) {
    printf("Iteration: %d\n", i);
    if ( isRandom ) {
      // first try to read from random
      readSize = fread (&randomBytes, 1, quantity, random);
    } else {
      // if cannot read all bits from random fill with urandom
      readSize = fread (&randomBytes, 1, quantity, urandom);
    }
    if ( readSize == quantity ) {
      break;
    }

    if ( i == 9 ) {
      fprintf( stderr, "Not enough bytes in random; reading from urandom\n" );
      isRandom = NO;
    }
  }
  fclose(random);
  fclose(urandom);

  fprintf(stdout, "String:%u\n", randomBytes);
  for (int i = 0; i < quantity; ++i) {
    printf("%u\n", randomBytes[i]);
  }

  free( randomBytes );

  // check for correct reading
  if ( readSize == quantity ) {
    return 1;
  }

  // move data from byte string to bits string
  // blah
  // return 0;


  int random = open( "/dev/random", O_RDONLY );
  int urandom = open( "/dev/urandom", O_RDONLY );
  if ( random == -1 || urandom == -1 ) {
    fprintf( stderr, "Could not open random source in /dev\n" );
    return 1;
  }
  // int myRandomInteger;
  // size_t randomDataLen = 0;
  // while (randomDataLen < sizeof myRandomInteger) {
  //   ssize_t result = read(randomData,
  //     ((char*)&myRandomInteger) + randomDataLen,
  //     (sizeof myRandomInteger) - randomDataLen);

  //   if (result < 0){
  //     // error, unable to read /dev/random 
  //   }
  //   randomDataLen += result;
  // }
  close(random);
  close(urandom);

  return 0;