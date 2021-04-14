#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* INITIALISATION MATSUMOTO */

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
#define nbExp 30 /*Paramètre défini pour le nombre de cases du tableau*/ 

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* FIN INITIALISATION MATSUMOTO */

/* ------------------------- Partie TP ------------------------- */

/* ------ Question 1 ------ */

double simuPiDisk(long nbOfPoints)
{

	double xr = 0.; /* x et y les coordonnées de chaque point */
	double yr = 0.;
	double pi = 0.;
	long i;
	long intDisk = 0; /* Variable qui compte le nombre de points à l'intérieur du disque */
	
	/* On fait un boucle d'un nombre de tirages de points pour avoir un répartition à l'intérieur et à l'extérieur du disque */
	
	for(i = 0 ; i < nbOfPoints ; i++)
	{
		xr = genrand_real1();
		yr = genrand_real1(); /* Tirage aléatoire des 2 coordonnées */
		
		if( (xr * xr + yr * yr) <= 1)
		{
			intDisk++; /* On incrémente pour chaque point à l'intérieur du disque */
		}
	}
	
	/* On pense à faire un cast forcé pour utiliser les entiers en double dans le calcul */
	pi = ((double) intDisk / (double) nbOfPoints ) * 4.;
	
	return pi;
}
		
/* ------ Question 2 ------ */

/* Exercice où l'on réalise une simulation de pi pour un nombre d'expérience donné afin d'en faire ensuite une moyenne de pi */

double simuAvgPi(double tabPi[], long nbOfPoints)
{
	double addPi = 0.; /* Variable intermédiaire au calcul de la moyenne de pi */
	long i;
	
	for(i = 0 ; i < nbExp ; i++)
	{
		tabPi[i] = simuPiDisk(nbOfPoints); /* Chaque case du tableau tabPi prend une valeur de pi générée par la fontion simuPiDisk */
		addPi += tabPi[i]; /* On ajoute à chaque fois la valeur de pi pour ensuite faire la moyenne */
	}
	
	return addPi / (double) nbExp; /* On retourne la moyenne des valeurs de pi du tableau */
}
	

/* ------ Question 3 ------ */

/* On prend n=30 (correspond à la variable nbExp qui est en define, on aura donc une valeur de la Student Law égale à 2,042 */

double piComputing(double avgPi, double xPi[nbExp])
{
	long i;
	double inter = 0.;
	double sum = 0.;
	for(i = 0 ; i < nbExp ; i++) /* Dans cette boucle on calcule le numérateur de S²(n) */
	{
		inter = (xPi[i] - avgPi) * (xPi[i] - avgPi);
		sum = sum + inter;
	}
	
	double squaredS = (sum / (((double) nbExp) - 1)); /* Calcul de S²(n) */
	printf("S²(n) = %.10f \n", squaredS);
	double radius = 2.042 * sqrt((squaredS) / (double) nbExp); /* Calcul de la marge d'erreur */
	return radius;
}


int main(void)
{

/* ------ Initialisation Matsumoto ------ */

   //int i;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
    //printf("1000 outputs of genrand_int32()\n");
    /*for (i=0; i<1000; i++) {
      printf("%10lu ", genrand_int32());
      if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand_real2()\n");
    for (i=0; i<1000; i++) {
      printf("%10.8f ", genrand_real2());
      if (i%5==4) printf("\n");*/
      

/* ---------------------- Partie TP personnel ---------------------- */

	/*Initialisations pour la Question 1*/
		
		long nbOfPoints = 0; /* Le nombre de points utils pour la simulation de pi */

	/*Initialisations pour la Question 2*/
	
		double tabPi[nbExp] = { 0. }; /*On initialise ici le tableau du nombre d'expérience que l'on veut faire afin que ses valeurs initialisées dans les fonctions puissent être stockées et réutilisées */
		double avgPi = 0.; /* De même pour la moyenne de pi */

	/*Initialisations pour la Question 3*/
	
		double radius = 0.;

/* ----------------------------- MENU ------------------------------ */
	
	int choix; /* On demande le choix de l'exercice que l'on veut exécuter */
	printf("Quel exercice voulez-vous voir (1, 2 ou 3)?\n");
	scanf("%d", &choix);	

	/* Switch qui exécutera l'exercice voulu */
	
	switch(choix)
	{
	
		/* ------ Question 1 ------ */
		
		case 1 :
			
			printf("Combien voulez-vous de points pour la simulation ?\n");
			scanf("%ld", &nbOfPoints); /*Le nombre de points est demandé afin de donner la possibilité d'avoir une valeur de pi plus précise et plus proche de la réalité.*/
			
			printf("Réponse :\n");
			
			/* SimuPiDisk fait une simulation de pi en prenant un nombre de points voulus */
			printf("Nous avons pi = %.06f \n", simuPiDisk(nbOfPoints));
			break;
			
		/* ------ Question 2 ------ */
		
		case 2 :
			
			nbOfPoints = 1000000000; /* C'est le nombre de points voulus pour la simulation de Pi, stockés ici pour pouvoir changer la valeur si besoin */
			avgPi = simuAvgPi(tabPi, nbOfPoints); /*On stock la valeur de la moyenne de Pi réalisée dans la fonction simuAvgPi pour 30 expériences ici afin de ne pas réaliser d'autres expériences en amont.*/
			printf("Réponse :\n");
			printf("La moyenne de Pi pour %d expériences et %ld points est de : %.10f \n", nbExp, nbOfPoints, avgPi);
			break;
			
		/* ------ Question 3 ------ */
		
		case 3 :
			 
			nbOfPoints = 1000; 
			avgPi = simuAvgPi(tabPi, nbOfPoints); /* Ici on stock une nouvelle fois la variable pour ne pas fausser le résultat avec une moyenne qui pourrait être recalculée si on ne le faisait pas à cet endroit là.*/
			printf("Réponse :\n");
			radius = piComputing(avgPi, tabPi); /* Le radius correspond à la marge d'erreur du calcul de l'intervale */
			printf("R = %.10f \n", radius);
			printf("Nous obtenons un intervalle de Pi sur : [ %.10f , %.10f ]\n", avgPi-radius, avgPi+radius);
			break;
		
		/* Si la valeur du choix n'est pas correcte */
		
		default : 
		
			printf("Les exercices sont numérotés de 1 à 3.\n"); 
			break;
	}

    return 0;
}

