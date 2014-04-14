#include "classes.h"

int main(int argc, char const *argv[]) {
    // const char *method;
    int arg;
    if(argc>1) arg=atoi(&*argv[1]);
    else if(argc>1 && argv[1]>1000) arg=1000;
    else arg=5;

    ScoreMatrix A(arg);
    A.Print();
    
    return EXIT_SUCCESS;
}