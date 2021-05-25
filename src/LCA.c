#include "GoLCA.h"

int LCA(int u,int v)
{
    int i;
    if(dep[u]<dep[v]){
        int tmp=u;
        u=v;
        v=tmp;
    }//保证u是深度大的那个节点
    int d=dep[u]-dep[v];//深度之差
    for(i=0;(1<<i)<=d;i++){//（1<<i）<=d是为了让u在保证不会跳过v的情况下进行跳跃
        if((1<<i)&d){   //（1<<i）&d其实就是转化二进制问那几个点可以跳跃
            u=fa[u][i];//例如差为5（101）只要在i==0和i==2的情况下跳跃(如果还是不懂就模拟一下）
        }
    }
    if(u==v)return u;//如果两个节点直接相等那么就如操作4所说，该点就是LCA，否则继续进行相同长度的跳跃
    for(i=M-1;i>=0;i--){
        if(fa[v][i]!=fa[u][i]){
            u=fa[u][i];
            v=fa[v][i];
        }
    }
    return fa[u][0];
}
