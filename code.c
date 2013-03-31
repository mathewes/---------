/**
	Author: mathewes
	Mail  : mathewes@hotmail.com
**/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>

/**输入文件**/
#define fileName "config.ini" 
/**输出文件**/
#define outfileName "out.ini"
/**精度**/
#define delta 5

#define maxPoint 10000

#define MID 0
#define OUTMID 1
#define OUTOUT 2
#define PI 3.1415926

struct Point{
	double x;double y;
};

struct Circle{
	struct Point center;
	double radius;
	int startIndex;
	int endIndex;
	double mxError;
};

typedef struct Point Point_t;
typedef struct Circle Circle_t;

Point_t points[maxPoint];
Circle_t circles[maxPoint];

int pointNum;int cirNum;

void readFile(char*fileName1)
{
	FILE*rFile;char tmp1[10];char tmp2[10];
	Point_t tmp;
	
	pointNum=0;cirNum=0;
	rFile=fopen(fileName1,"r");
	if(!rFile)
	{
		printf("file Not exist");
		exit(-1);
	}
	while(fscanf(rFile,"%s %s",tmp1,tmp2)!=EOF){
		tmp.x=atof(tmp1);tmp.y=atof(tmp2);
		points[pointNum++]=tmp;
	}
}

/**参考 VB求三角形外界圆心与半径 一文**/
void GetCircle(Point_t p1,Point_t p2,Point_t p3,Circle_t*circle)
{
    double mat1,mat2,mat3;
	double a11=-2*p1.x+2*p2.x;
	double a12=-2*p1.y+2*p2.y;
	double b1=-p1.x*p1.x+p2.x*p2.x-p1.y*p1.y+p2.y*p2.y;
	double a21=-2*p1.x+2*p3.x;
	double a22=-2*p1.y+2*p3.y;
	double b2=-p1.x*p1.x+p3.x*p3.x-p1.y*p1.y+p3.y*p3.y;
	double D=a11*a22-a12*a21;
	double D1=b1*a22-a12*b2;
	double D2=a11*b2-b1*a21;
    if(fabs(mat3)<0.01)
    {
    	circle->center.x=(p1.x+p2.x)/2.0;
    	circle->center.y=(p1.y+p2.y)/2.0;    	
    }
    else{
    	circle->center.x = D1/D;
    	circle->center.y = D2/D;
    }
    circle->radius = sqrt((double)((p1.x-circle->center.x)*(p1.x-circle->center.x) + (p1.y-circle->center.y)*(p1.y-circle->center.y)));
}

/**通过点坐标获得角度**/
double getAngle(Point_t*p1)
{
	double atan1;
	if(p1->x==0&&p1->y==0)
		return 0;
	if(p1->x==0)
		return (p1->y)>0?PI/2:PI*3.0/2.0;
	else if(p1->y==0)
		return (p1->x)>0?0:PI;
	atan1=atan(p1->y/p1->x);
	if(atan1>=0)
	{
		if(p1->x>=0)
			return atan1;
		return PI+atan1;
	}
	else
	{
		if(p1->x>=0)
			return PI*1.5+atan1;
		return PI*0.5+atan1;
	}
}

void swap(double*d1,double*d2)
{
	double tmpAng;
	tmpAng=*d1;
	*d1=*d2;
	*d2=tmpAng;
}

/**通过平移的方式获得p1相对center的坐标**/
void GetNormalPoint(Point_t*p1,Point_t*center)
{
	p1->x=(p1->x)-(center->x);
	p1->y=(p1->y)-(center->y);
}

/**判断p3点位于弧的那一个位置（共三种情况)**/
int flagForPointToArc(Circle_t*circle,int p1Index)
{
	double ang1;double ang2;double ang3;double tmpAng;
	Point_t p1;Point_t p2;Point_t p3;Point_t center;
	center=circle->center;
	p1=points[circle->startIndex];
	p2=points[circle->endIndex];
	p3=points[p1Index];
	GetNormalPoint(&p1,&center);
	GetNormalPoint(&p1,&center);
	GetNormalPoint(&p1,&center);
	ang1=getAngle(&p1);ang2=getAngle(&p2);ang3=getAngle(&p3);
	if(ang1>ang2)
		swap(&ang1,&ang2);
	if(ang2-ang1<=PI)
	{
		ang2-=ang1;ang3-=ang1;ang1=0;
		if(ang3<0)
			ang3=2*PI-ang3;
	}
	else
	{
		double tmpRemain=2*PI-ang2;
		ang2=ang1+tmpRemain;ang1=0;ang3+=tmpRemain;
		if(ang3>=2*PI)
			ang3-=2*PI;
	}
	if(ang3<=ang2)
		return MID;
	else if(ang3>=PI&&ang3<ang2+PI)
		return OUTMID;
	else
		return OUTOUT;
}

double getDis(Point_t*p1,Point_t*p2)
{
	return sqrt(((p1->x)-(p2->x))*((p1->x)-(p2->x))+((p1->y)-(p2->y))*((p1->y)-(p2->y)));
}

/**获得p1Index对应的点相对圆弧的误差**/
double getErrorByCircle(Circle_t*circle,int p1Index)
{
	Point_t p1;double dis;
	p1=points[p1Index];
	dis=(p1.x-circle->center.x)*(p1.x-circle->center.x)+(p1.y-circle->center.y)*(p1.y-circle->center.y);
	dis=sqrt(dis);
	int rV=flagForPointToArc(circle,p1Index);
	if(rV==MID)
		return fabs(dis-(circle->radius));
	else
	{
		Point_t p2=points[circle->startIndex];
		Point_t p3=points[circle->endIndex];
		double d1=getDis(&p1,&p2);
		double d2=getDis(&p1,&p3);
		return d1>d2?d2:d1;
	}
}

/**核心函数**/
/**
这一个函数是判断step步， p1Index为起始点下标的时候能否达到要求
Return value: 
false 达不到要求
true 能达到要求，同时修改传入的circle参数，从而获得circle数据
**/
bool whetherOk(int p1Index,int step,Circle_t*circle)
{
	int midIndex;int midIndex2;Circle_t tmpCircle;

	double tmp11=-1;
	
	circle->startIndex=p1Index;
	circle->endIndex=p1Index+step-1;
	tmpCircle.startIndex=p1Index;
	tmpCircle.endIndex=p1Index+step-1;
	for(midIndex=p1Index+1;midIndex<p1Index+step-1;midIndex++)
	{
		GetCircle(points[p1Index],points[midIndex],points[p1Index+step-1],&tmpCircle);
		for(midIndex2=p1Index+1;midIndex2<p1Index+step-1;midIndex2++)
		{
			if(midIndex==midIndex2)
				continue;
			double vv=getErrorByCircle(&tmpCircle,midIndex2);
			if(vv>=delta){
				return false;
			}
			if(vv>tmp11)
			{
				circle->center=tmpCircle.center;
				circle->radius=tmpCircle.radius;
				circle->mxError=vv;
				tmp11=vv;
			}	
		}
	}
	return true;
}
/**
  在step为4的时候，如果也不满足则以首尾为直径。使用getSpecailCircle函数获得
**/
void getSpecialCircle(int p1Index,Circle_t*circle)
{
	double vv1;double vv2;
	circle->center.x=(points[p1Index].x+points[p1Index+3].x)/2.0;
	circle->center.y=(points[p1Index].y+points[p1Index+3].y)/2.0;
	circle->startIndex=p1Index;
	circle->endIndex=p1Index+3;
    circle->radius = sqrt((double)((points[p1Index].x-circle->center.x)*(points[p1Index].x-circle->center.x) + (points[p1Index].y-circle->center.y)*(points[p1Index].y-circle->center.y)));
    
    vv1=getErrorByCircle(circle,p1Index+1);
    vv2=getErrorByCircle(circle,p1Index+2);

    vv1=vv1>vv2?vv1:vv2;
    circle->mxError=vv1;
}
/**
	主要的函数，作用是以p1Index为起始点的坐标，计算最大可能的step,并生成圆弧。
**/
int function(int p1Index,Circle_t*circle)
{
	int pEndIndex;int step=4;Circle_t tmpCircle;
	if(p1Index+3>=pointNum)
	{
		return pointNum;
	}
	while(1)
	{
		if(p1Index+step>pointNum)
			break;
		if(whetherOk(p1Index,step,&tmpCircle)==false)
			break;

		*circle=tmpCircle;
		step=step+1;
	}
	if(step==4)
	{
		getSpecialCircle(p1Index,circle);
	}
	if(step!=4)
		step=step-1;
	return p1Index+step;
}

void writeFile(char*fileName1)
{
	FILE*wFile;int numTmp;Circle_t circle1;Point_t point;

	wFile=fopen(fileName1,"w");
	if(!wFile)
	{
		printf("can't write to file");
		exit(-1);
	}
	fprintf(wFile,"#每一段弧线起点下标 每一段弧线终点下标 圆心坐标x分量 圆心坐标y分量 半径 误差\n");
	for(numTmp=0;numTmp<cirNum;numTmp++)
	{
		circle1=circles[numTmp];
		point=circle1.center;
		fprintf(wFile,"%d %d %f %f %f %f\n",circle1.startIndex,circle1.endIndex,point.x,point.y,circle1.radius,circle1.mxError);
	}
}

int main(int argc, char const *argv[])
{
	int p1Index;Circle_t tmpCircle;
	
	p1Index=0;
	readFile(fileName);
	
	while(p1Index+3<pointNum)
	{
		p1Index=function(p1Index,&tmpCircle);
		circles[cirNum++]=tmpCircle;
	}
	writeFile(outfileName);
	return 0;
}