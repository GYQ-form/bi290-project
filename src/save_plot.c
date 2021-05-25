#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GoLCA.h"

void save_plot()//链式前向星存图
{

    FILE *fp;
    if((fp=fopen("slim_extract","r"))==NULL){
    printf("打开文件错误！\n");
    exit(1);
    }

    int i,j=0,k,flag=0,u=1,v=1;
    char line[100],code[8];

    while(!feof(fp)){
        fgets(line,100,fp);
        line[7]='\0';

        for(k=1;k<n;k++){
            if(strcmp(arr[k],line)==0){v=k;break;}
        }

        if(k==n){
            v=n;
            arr[n]=(char*)malloc(sizeof(char)*8);
            strcpy(arr[n++],line);
        }

        fgets(line,100,fp);

        for(i=0;line[i]!='\n';i++){
            if(line[i]!='\t') code[j++]=line[i];
            else{
                code[j]='\0';
                j=0;

                for(k=1;k<n;k++){
                    if(strcmp(arr[k],code)==0){u=k;break;}
                }

                if(k==n){
                    u=n;
                    arr[n]=(char*)malloc(sizeof(char)*8);
                    strcpy(arr[n++],code);
                }

                flag=1;
            }

            if(flag==1){
                to[tot]=v;
                nx[tot]=head[u];
                head[u]=tot++;
                in[v]++;//处理入度
                flag=0;
            }
        }
    }


    fclose(fp);
}
