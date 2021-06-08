#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int v=1;

void randtree(FILE *fp,int d)
{
    if(d>=10)return;
    int i,u;

    int n=rand()%4+1;
    for(i=1;i<=n;i++){
        u=v;
        v++;
        fprintf(fp,"%d\n",v);
        fprintf(fp,"%d\n",u);
        randtree(fp,d+1);
    }

    return;

}

int main()
{
    FILE *fp;
    srand(time(NULL));

    fp=fopen("random_tree.txt","w");
    randtree(fp,0);

    printf("successful!\n");

    return 0;
}
