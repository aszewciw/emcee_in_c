/* MT stuff */


/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* seed by converting the input string to unsigned long */
int init_genrand_str(const char *str);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_uint32(void);

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void);

/* These real versions are due to Isaku Wada, 2002/01/09 added */
/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);

/*
   Seed the random number generator using /dev/random or /dev/urandom

   /dev/random is better but can block

   /dev/urandom will not block but can run out of entropy
*/
void randn_seed_devrand(void);
void randn_seed_devurand(void);


/*
   Generate random unsigned 32 but integers in [0,maxval)
*/

unsigned int genrand_uint32_max(unsigned int maxval);

/*
   Generate uniform random doubles between [0,1]
   currently calls genrand_res53
*/
double randu();

/*
   random numbers in the range [-1,1]
*/
double srandu();

/*
  Generate gaussian random deviates
*/

double randn();

/*
   Generate a poisson deviate.

   This is apparently from Knuth
*/
long poisson(double lambda);



