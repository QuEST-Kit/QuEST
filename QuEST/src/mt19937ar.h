#ifndef MT_RAND_H
#define MT_RAND_H

#ifdef __cplusplus
extern "C" {
#endif

void init_by_array(unsigned long init_key[], int key_length);

void init_genrand(unsigned long s);

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);

#ifdef __cplusplus
}
#endif

#endif // MT_RAND_H

