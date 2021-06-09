#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 442921
int fa[N][6];
int dep[N];
int head[N];
int nx[N];
int to[N];
int tot=1;
int vis[N];
int in[N];
int L,R;

void add()
{
    FILE *fp;
    if((fp=fopen("random_tree3.txt","r"))==NULL){
    printf("打开文件错误！\n");
    exit(1);
    }

    int u,v;
    char code[8];

    while(!feof(fp)){
        fgets(code,8,fp);
        v=atoi(code);
        fgets(code,8,fp);
        u=atoi(code);

        to[tot]=v;
        nx[tot]=head[u];
        head[u]=tot++;
        fa[v][0]=u;
        in[v]++;//处理入度

        }
    fclose(fp);
    }

void init(int n){//初始化
    for(int i=0;i<=n;i++)fa[i][0]=0,in[i]=0,vis[i]=0,head[i]=0,dep[i]=0;
    tot=1;
}

void dfs(int u,int d){//dfs处理深度,处理2^j的祖先关系
    vis[u]=1;
    dep[u]=d;
    for(int i=1;i<6;i++)fa[u][i]=fa[fa[u][i-1]][i-1];
    for(int i=head[u];i;i=nx[i]){
        int v=to[i];
        if(!vis[v])dfs(v,d+1);
    }
}

int LCA(int u,int v){
    if(dep[u]<dep[v])if(dep[u]<dep[v]){
        int tmp=u;
        u=v;
        v=tmp;
    }//保证u是深度大的那个节点

    int d=dep[u]-dep[v];//深度之差
    for(int i=0;(1<<i)<=d;i++){//（1<<i）<=d是为了让u在保证不会跳过v的情况下进行跳跃
        if((1<<i)&d){   //（1<<i）&d其实就是转化二进制问那几个点可以跳跃
            u=fa[u][i];//例如差为5（101）只要在i==0和i==2的情况下跳跃(如果还是不懂就模拟一下）
        }
    }
    if(u==v)return u;//如果两个节点直接相等那么就如操作4所说，该点就是LCA，否则继续进行相同长度的跳跃
    for(int i=5;i>=0;i--){
        if(fa[v][i]!=fa[u][i]){
            u=fa[u][i];
            v=fa[v][i];
        }
    }
    return fa[u][0];
}
int main()
{
        clock_t start,finish;
        start=clock();

        init(N);
        add();

        printf("save plot ... done\n");

        srand(time(NULL));

        L=rand()*10%442920+1;
        R=rand()*10%442920+1;
        for(int i=1;i<=tot;i++)
            if(!in[i]){
                dfs(i,0);//必须从入度为0的点开始dfs
                break;
            }

        printf("the LCA between %d and %d is %d\n",L,R,LCA(L,R));

        finish=clock();
        double duration;
        duration = (double)(finish - start) / CLOCKS_PER_SEC;

        printf("the excution time is:%f s\n",duration);

        return 0;
}
