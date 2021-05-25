#include "GoLCA.h"
#include <stdio.h>

void init()//initializing
{
    for(int i=0;i<N;i++)fa[i][0]=0,in[i]=0,head[i]=0,dep[i]=0,arr[i]=NULL;
    tot=1;
    n=1;
}
