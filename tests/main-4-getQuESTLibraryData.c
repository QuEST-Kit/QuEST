#include <string.h>
#include <stdio.h>
#include "QuEST.h"

int main(int, char*[]) {
  QuESTEnv quenv = createQuESTEnv();

  {
       QuESTLibraryData libdat;
       printf("Testing `getQuESTLibraryData()`...\n"); fflush(stdout);
       getQuESTLibraryData(quenv, &libdat);

       printf("CUDA=%d OpenMP=%d MPI=%d threads=%d ranks=%d QuEST_PREC=%d \n",
               (int)libdat.CUDA, (int)libdat.OpenMP, (int)libdat.MPI,
               (int)libdat.threads, (int)libdat.ranks, (int)libdat.QuEST_Prec);

       printf("Test correctness of the precision info...\n"); fflush(stdout);
       if ( libdat.QuEST_Prec == (int)sizeof(qreal)/4 ) {
            printf("... fine.\n");
       } else {
            printf("`getQuESTLibraryData()` gave %d but truth is %d.\n",
                   (int)libdat.QuEST_Prec, (int)sizeof(qreal)/4);
       }

       printf("Test correctness of the rest via `getEnvironmentString()`...\n"); fflush(stdout);
       char str  [256] = "";
       char costr[256] = "";
       getEnvironmentString(quenv, str);
       sprintf(costr, "CUDA=%d OpenMP=%d MPI=%d threads=%d ranks=%d",
               (int)libdat.CUDA, (int)libdat.OpenMP, (int)libdat.MPI,
               (int)libdat.threads, (int)libdat.ranks);

       if ( strncmp(str,costr,256) ){
            str[255] = costr[255] = '\0';
            printf("`getQuESTLibraryData()` gave\n\t%s\nbut `getEnvironmentString()` gave\n\t%s\n",
                   costr, str);
       } else {
            printf("... fine.\n");
       }
  }

  destroyQuESTEnv(quenv);
  return 0;
}
