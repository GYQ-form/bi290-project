#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GoLCA.h"

int n;
char *arr[N];
int fa[N][7];
int dep[N];
int head[N];
int nx[N];
int to[N];
int tot;
int in[N];


int main(int argc,char *argv[])
{
	if(argc!=3||strlen(argv[1])!=7||strlen(argv[2])!=7){
		printf("usage:<filename> GO_id1 GO_id2\n");
		exit(1);
	}
        int lca;
        init();
        save_plot();

        int l,r,i,flag=0;

         for(i=1;i<n;i++)
            if(!in[i]){
                dfs(i,0);//必须从入度为0的点开始dfs
                father0(i);
                father1(i);
            }

        for(i=1;i<=n;i++){
            if(strcmp(arr[i],argv[1])==0){l=i;flag++;}
            if(strcmp(arr[i],argv[2])==0){r=i;flag++;}
            if(flag==2)break;
        }

	if(flag!=2){
		printf("the GO id is not found!\n");
		exit(1);
	}

        if((lca=LCA(l,r))==0) printf("there is no LCA for the two GO terms!\n");
        else printf("the LCA between two GO terms is:%s\n",arr[lca]);

        for(i=1;i<n;i++)free(arr[i]);

        return 0;

}
