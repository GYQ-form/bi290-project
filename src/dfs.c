#include "GoLCA.h"
#include <stdio.h>

void dfs(int u,int d)//dfs处理深度,处理2^j的祖先关系
{
    int i;
    if(dep[u]<=d){
        dep[u]=d;

        for(i=head[u];i;i=nx[i]){
            int v=to[i];
            dfs(v,d+1);
        }
    }
    return;
}
