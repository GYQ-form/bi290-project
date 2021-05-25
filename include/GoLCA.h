#define N 300
#define M 7
#ifndef GO_LCA
#define GO_LCA
void save_plot();
void init();
void dfs(int u,int d);
int LCA(int u,int v);
void father0(int u);
void father1(int u);
#endif // GO_LCA

extern int n;
extern char *arr[N];
extern int fa[N][M];
extern int dep[N];
extern int head[N];
extern int nx[N];
extern int to[N];
extern int tot;
extern int in[N];
