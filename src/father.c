#include "GoLCA.h"

void father0(int u)
{
    int i;
    for(i=head[u];i;i=nx[i]){
            int v=to[i];
            if(dep[u]>=dep[fa[v][0]])fa[v][0]=u;
            father0(v);
        }
}

void father1(int u)
{
    int i;
    for(i=1;i<M;i++)fa[u][i]=fa[fa[u][i-1]][i-1];
    for(i=head[u];i;i=nx[i]){
            int v=to[i];
            father1(v);
        }
}
