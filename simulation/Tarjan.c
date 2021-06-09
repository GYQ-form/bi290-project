#include <stdio.h>
#include <stdlib.h>
#define N 5000000
#include <time.h>

int fa[N];
int head[N];
int nx[N];
int to[N];
int tot=1;
int vis[N];
int in[N];
int L,R;
void add(int u,int v){
    to[tot]=v;
    nx[tot]=head[u];
    head[u]=tot++;
}
int find(int x){
    if(fa[x]==x) return x;
    else {
        fa[x]=find(fa[x]);
        return fa[x];
    }

}
void tarjan(int u){
    vis[u]=1;
    for(int i=head[u];i;i=nx[i]){
        int v=to[i];
        if(!vis[v]){
            tarjan(v);//不停的向下走直到再往下就没有子节点
            fa[v]=u;//更新每个节点的父亲
        }
    }
    if(u==L){//如果u是需要查询两节点的任意一个节点，就要判断另一个节点是否已经被更新过，如果被更新过，那么另一个节点的祖先一定是两个节点的LCA
        if(fa[R]!=R)
            printf("LCA between %d and %d is %d\n",L,R,find(R));
    }
    else if(u==R){
        if(fa[L]!=L){
            printf("LCA between %d and %d is %d\n",L,R,find(L));
        }
    }
}
void init(int n){
    for(int i=0;i<=n;i++)fa[i]=i,in[i]=0,vis[i]=0,head[i]=0;
    tot=1;
}
int main(){
    clock_t start, finish;
    double  duration;
    start = clock();
    FILE *fp;
    if((fp=fopen("random_tree1.txt","r"))==NULL){
        printf("open file failed!\n");
        exit(1);
    }

    tot = 1;
        int max = 1;
        char code[8];
        int temp;

        while (!feof(fp)) {
            fgets(code, 8, fp);
            temp = atoi(code);
            max = temp > max ? temp : max;

        }
    printf("max=%d\n",max);
        init(max);
    rewind(fp);

        for (int i = 1; i < 22377; i++) {
            int l, r;
            fgets(code, 8, fp);
            r = atoi(code);
            fgets(code, 8, fp);
            l = atoi(code);

            add(l, r);
            in[r]++;//依旧还是从入度为0的点开始

        }

        srand(time(NULL));
        L = rand() % max + 1;
        R = rand() % max + 1;


        for (int i = 1; i <= 22377; i++)
            if (!in[i]) {
                tarjan(i);
                break;
            }
    tarjan(4);
    rewind(fp);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "execution time is %f seconds\n", duration );
    system("pause");
        return 0;

}
