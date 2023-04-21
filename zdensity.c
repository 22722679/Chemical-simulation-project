#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"

#include <math.h>
#define UTUBE 10.224  //cnt底  GO-Right-x
#define BTUBE 10.086      //cnt顶   GO-Left-x
#define RTUBE 0.401          //cnt半径（碳纳米管半径）
#define NOSOL 3936         //溶质数
#define TOTAL 22518      //总数
#define XBOX 6.39        //盒子的x轴长度    
#define YBOX 6.60			//盒子的y轴长度 
#define ZBOX  3.00			//盒子的z轴长度 
#define steplength 300		//z轴方向上共分为多少步 
#define accuracy 0.01		//z轴方向上每一步的精度 
#define Cnum 1620			//C原子的数量 
#define H2Onum 3461	//水分子的数量 
#define ionnum 89			//因为Na离子和Cl离子的数量相等，所以统称为离子数量 


int natoms, step,natom,read_return;
float time,prec;
matrix box;
rvec *x,*ox;
XDRFILE *xtc;
int u2b,b2u,note[TOTAL+100];
float r,or;
FILE *ff;

float pbcr(int a,int b)
{
	float X,Y,Z,R;
 	X=x[a][0]-x[b][0];
  	if(X>XBOX/2)X=X-XBOX;
  	if(X<-XBOX/2)X=X+XBOX;    //周期性边界条件
  	Y=x[a][1]-x[b][1];
  	if(Y>YBOX/2)Y=Y-YBOX;
  	if(Y<-YBOX/2)Y=Y+YBOX;
  	Z=x[a][2]-x[b][2];
  	if(Z>ZBOX/2)Z=Z-ZBOX;
  	if(Z<-ZBOX/2)Z=Z+ZBOX;
  	R=sqrt(X*X+Y*Y+Z*Z);
  	return R;                  //返回距离
}

main ()
{

	ff=fopen("md_regra2ion.txt","w");  	          //写f.dat文件
	xtc=xdrfile_open ("md_regra2ion.xtc","r");
   	read_xtc_natoms ("md_regra2ion.xtc",&natoms);
 	x=calloc(natoms, sizeof (x[0]));
	ox = calloc(natoms, sizeof (x[0]));
	int H2Ore=Cnum+H2Onum*3;
	int H2Opmax=0,Napmax=0,Clpmax=0;
	int Nare=Cnum+H2Onum*3+ionnum;
	int Clre=Nare+ionnum;
	int N=0, i,j,count1[1000]={0},count2[1000]={0},count3[1000]={0};
	read_return=read_xtc(xtc,natoms,&step,&time,box,x,&prec);
	
	printf ("%d\n",step);
	                                       //从trr文件中读取一个step time 的natoms个原子的坐标，先给ox，ox 就表示前一个原子
	while (1)
	{
	 	read_return=read_xtc (xtc,natoms,&step,&time,box,x,&prec); //这一次读取的时候就表示下一个原子，给X
	  	if (read_return!=0) break;
	  	if(step%1000000==0)printf ("%d %f\n",step,time);

	  	//fprintf(f1,"%f\n",time);
    	if(time>-1)
   		{
    		N++;
   			for(i=Cnum;i<H2Ore;i+=3)
 	   		{
	    		for(j=0;j<steplength;j++)
				{
					if(x[i][2]>j*accuracy && x[i][2]<=(j+1)*accuracy)
					{
						count1[j]++;
					}
				}
  	   		}

  			for(i=H2Ore;i<Nare;i++)
  			{
	  			for(j=0;j<steplength;j++)
	  			{
	  	  			if(x[i][2]>j*accuracy && x[i][2]<=(j+1)*accuracy)
					{
						count2[j]++;
					}
	  			}
	  		}
   			for(i=Nare;i<Clre;i++)
			{
	  			for(j=0;j<steplength;j++)
	  			{
	  	  			if(x[i][2]>j*accuracy && x[i][2]<=(j+1)*accuracy)
					{
						count3[j]++;
					}
	  			}
	  		}
 		}
  	}
  	//double ions=1/ionnum;
  	xdrfile_close(xtc);
	float res2=ionnum*23*0.01/(6.02214076*XBOX*YBOX*0.01),res3=ionnum*35.5*0.01/(6.02214076*XBOX*YBOX*0.01);
	fprintf(ff,"Z	  zH2Onum 	 H2O p/pbuik	Na+p/pbulk	  Cl-p/pbulk\nnm \n\n");
  	for(j=0;j<steplength;j++)
  	{
  		if(count2[steplength-10]==NULL && count3[steplength-10]==NULL){
		  	double res1=count1[j]/N;
  			fprintf(ff,"%4.3f	%10.5lf	 %10.10lf        \n",(j+1)*0.01,res1,count1[j]*18*0.01/(N*6.02214076*XBOX*YBOX*0.01*0.996782));
		  }
	  	else{
  			double res1=count1[j]/N;
  			fprintf(ff,"%4.3f	%10.5lf	 %10.6lf     %10.5lf   	%10.5lf   \n",(j+1)*0.01,res1,count1[j]*18*0.01/(N*6.02214076*XBOX*YBOX*0.01*0.996782),count2[j]*23*0.01/(N*6.02214076*XBOX*YBOX*0.01*res2),count3[j]*35.5*0.01/(N*6.02214076*XBOX*YBOX*0.01*res3));
  		}	
 	}
	fclose(ff);
//count2[j]*23*10/(N*6.02214076*XBOX*YBOX*0.01),count2[j]*ions,count3[j]*35.5*10/(N*6.02214076*XBOX*YBOX*0.01),count3[j]*ions
}
