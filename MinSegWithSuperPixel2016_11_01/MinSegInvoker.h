#pragma once

#include <opencv2\opencv.hpp>
using namespace cv;

#define STEP 50//������Ŀ
#define GRADE 60//����
#define ADDPOINT 3
#define MAXDATA 999
#define MAXMAXDATA 999999
#define Dis 5
#define PARAM 0.95
typedef unsigned char UCHAR;

class MinSegInvoker
{
	int Initial();
	int RemovePoint();
	int CalcCost();
	int LocalNucCost();
	int BackgroundCost();
	int ModifyCost(Mat &_scanLineRet);
	void MinPath();
	int DrawEdge(Mat &_img);

	const Mat m_Img;
	const Mat m_SuperPixelRet;

	double m_dCost[GRADE][STEP];	//���۾���
	int m_nGray[GRADE][STEP];		//Rͨ����Gͨ��
	Point m_Point[GRADE][STEP];		//�����������
	int m_nMinroad[GRADE];			//��С����·��
	double m_dFinalcost;			//��С����

public:
	MinSegInvoker(const Mat &_img, const Mat &_superPixelRet);
	int Seg(Mat &_costRet, Mat &_superPixelCostRet, Mat &_scanLineRet, Mat &_img);
	virtual ~MinSegInvoker();
};

