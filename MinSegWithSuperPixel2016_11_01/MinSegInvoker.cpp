#include "MinSegInvoker.h"

MinSegInvoker::MinSegInvoker(const Mat &_img, const Mat &_superPixelRet): m_Img(_img), m_SuperPixelRet(_superPixelRet)
{
	Initial();
	RemovePoint();
	CalcCost();
	LocalNucCost();
	BackgroundCost();
}

MinSegInvoker::~MinSegInvoker()
{
}

int MinSegInvoker::Initial()
{
	int height = m_Img.rows;
	int width = m_Img.cols;
	/*提取灰度图像*/
	Mat resImg;
	medianBlur(m_Img, resImg, 5);
	//GaussianBlur(resImg, resImg, Size(5, 5), 0);

	/*计算起始坐标*/
	int t_nLength = std::min(width, height) * 14 / 30;	//较短的一条边,用于计算采样点坐标

	//图像中心点
	int CenterPointX = width / 2 + 1;
	int CenterPointY = height / 2 + 1;

	/*初始化数组*/
	for (int j = 0; j < GRADE; ++j)
	{
		for (int i = 0; i < STEP; ++i)
		{
			/*搜索点的坐标*/
			m_Point[j][i].x = (int)(t_nLength * cos(j * 2 * 3.14 / GRADE) * i / STEP) + CenterPointX;
			m_Point[j][i].y = (int)(t_nLength * sin(j * 2 * 3.14 / GRADE) * i / STEP) + CenterPointY;

			m_Point[j][i].x = std::min(m_Point[j][i].x, width - 1);
			m_Point[j][i].x = std::max(m_Point[j][i].x, 0);
			m_Point[j][i].y = std::min(m_Point[j][i].y, height - 1);
			m_Point[j][i].y = std::max(m_Point[j][i].y, 0);
			/*搜索点的灰度值*/

			auto pixel = resImg.at<Vec3b>(m_Point[j][i].y, m_Point[j][i].x);
			m_nGray[j][i] = static_cast<int>(pixel[2]) - pixel[1];	//红-绿
		}
	}

	memset(m_dCost, 0, sizeof(double) * GRADE * STEP);

	return 0;
}

//把深色的细胞核区域插值成浅色
int MinSegInvoker::RemovePoint()
{
	int threshold = 70;
	for (int j = 0; j < GRADE; ++j)
	{
		for (int i = 5; i < STEP; ++i)
		{
			/*搜索黑点*/
			int p, q;
			if (m_nGray[j][i] > threshold)
			{
				p = i;
				q = STEP - 1;
				for (i = p + 1; i < STEP; ++i)
				{
					if ((m_nGray[j][i - 1] < threshold) && (m_nGray[j][i] >= threshold))
					{
						q = i - 1;
						break;
					}
				}
				int len, g_len;
				len = q - p + 1;
				g_len = abs(m_nGray[j][p - 1] - m_nGray[j][q + 1]);
				for (i = p; i <= q; ++i)
				{
					if (m_nGray[j][p] < m_nGray[j][q])
					{
						m_nGray[j][i] = m_nGray[j][p - 1] + (i - p + 1) / len * g_len;
					}
					else
					{
						m_nGray[j][i] = m_nGray[j][p - 1] - (i - p + 1) / len * g_len;
					}
				}
				i = q;
			}
		}
	}

	return 0;
}

double Calcd(double d)
{
	double calcd;
	//calcd = 1 / ( d + 12);
	double a = -0.071;
	calcd = exp(a * d);
	return calcd;
}

//以梯度为基础，计算代价矩阵
int MinSegInvoker::CalcCost()
{
	/*计算梯度数组*/
	for (int j = 0; j < GRADE; ++j)
	{
		double sum = 0;
		for (int i = ADDPOINT; i < STEP - ADDPOINT; ++i)
		{
			int t_nAddUpA = 0;
			int t_nAddUpB = 0;

			for (int k = 0; k < ADDPOINT; ++k)
			{
				t_nAddUpA += m_nGray[j][i - k];
				t_nAddUpB += m_nGray[j][i + k];
			}
			//double d = abs(t_nAddUpA - t_nAddUpB);
			//double d = t_nAddUpB - t_nAddUpA;
			double d = t_nAddUpA - t_nAddUpB;
			m_dCost[j][i] = Calcd(d);
			sum += m_dCost[j][i];
		}
		for (int i = ADDPOINT; i < STEP - ADDPOINT; ++i)
		{
			m_dCost[j][i] /= sum;
		}
		for (int i = 0; i < ADDPOINT; ++i)
		{
			m_dCost[j][i] = MAXDATA;
		}

		for (int i = STEP - ADDPOINT; i < STEP; ++i)
		{
			m_dCost[j][i] = MAXDATA;
		}
	}

	return 0;
}

//屏蔽细胞核中心区域
int MinSegInvoker::LocalNucCost()
{
	/*自核范围15 = 150 * 9 / 10 * i / STEP,i = 7*/
	for (int j = 0; j < GRADE; ++j)
	{
		for (int i = 0; i < 11; ++i)
		{
			m_dCost[j][i] = MAXDATA;
		}
	}

	return 0;
}

//屏蔽边缘部分
int MinSegInvoker::BackgroundCost()
{
	/*背景代价设为无穷大*/
	for (int j = 0; j < GRADE; ++j)
	{
		for (int i = 0; i < STEP - 2; ++i)
		{
			if (m_nGray[j][i] == 0 && m_nGray[j][i + 1] == 0 && m_nGray[j][i + 2] == 0)
			{
				m_dCost[j][i] = MAXDATA;
				for (int m = i; m < STEP; ++m)
				{
					m_dCost[j][m] = MAXDATA;
				}
			}
		}
	}

	return 0;
}

double WeightFunc(int val)
{
	double ret;
	switch (val)
	{
	case 0:
		ret = 1.0;
		break;
	case 1:
		ret = 1.0;
		break;
	case 2:
		ret = 1.001;
		break;
	case 3:
		ret = 1.002;
		break;
	case 4:
		ret = 1.003;
		break;
	default:
		ret = 9999.;
	}

	return ret;
}

void DrawNop(Vec3b &vec, int flag)
{

	if (flag)
	{
		vec.val[0] = (UCHAR)0;
		vec.val[1] = (UCHAR)255;
		vec.val[2] = (UCHAR)255;
	}
	else
	{
		vec.val[0] = (UCHAR)255;
		vec.val[1] = (UCHAR)255;
		vec.val[2] = (UCHAR)255;
	}
}

void DrawYesp(Vec3b &vec, int flag)
{

	if (flag)
	{
		vec.val[0] = (UCHAR)0;
		vec.val[1] = (UCHAR)255;
		vec.val[2] = (UCHAR)0;
	}
	else
	{
		vec.val[0] = (UCHAR)255;
		vec.val[1] = (UCHAR)0;
		vec.val[2] = (UCHAR)255;
	}
}

//修正代价矩阵，生成扫描线图像
int MinSegInvoker::ModifyCost(Mat &_scanLineRet)
{
	_scanLineRet = m_Img.clone();

	for (int j = 0; j < GRADE; ++j)
	{
		double min_Cost = MAXMAXDATA;
		int min_position = -1;
		int p = 0;
		bool flag = true;
		for (int i = 1; i < STEP - 2; ++i)
		{
			if (m_nGray[j][i + 2] == 0 && m_nGray[j][i] == 0 && m_nGray[j][i + 1] == 0)
			{
				p = i;
				//_Cost[j][p] = 0;
				flag = false;
				break;
			}

			if (min_Cost > m_dCost[j][i])
			{
				min_Cost = m_dCost[j][i];
				min_position = i;
			}
		}

		if (p != 0)
		{
			/*line(m_ScanLineImg, Point(m_ScanLineImg.cols / 2 + 1, m_ScanLineImg.rows / 2 + 1)
			, Point(m_Point[j][p].x, m_Point[j][p].y), flag ? Scalar(0, 255, 0) : Scalar(255, 0, 255));*/
			//有p没找到边界绿色：找到边界跳出没有找p紫色
			/*Vec3b &vec = g_img.at<Vec3b>(m_Point[j][p].y, m_Point[j][p].x);
			vec.val[0] = (UCHAR) 0;
			vec.val[1] = (UCHAR) 255;
			vec.val[2] = (UCHAR) 0;*/
			for (int k = 0; k < p; ++k)
			{
				Vec3b &vec0 = _scanLineRet.at<Vec3b>(Point(m_Point[j][k].x, m_Point[j][k].y));
				DrawYesp(vec0, flag);
			}
		}
		else
		{
			/*line(m_ScanLineImg, Point(m_ScanLineImg.cols / 2 + 1, m_ScanLineImg.rows / 2 + 1)
			, Point(m_Point[j][STEP - 1].x, m_Point[j][STEP - 1].y), flag ? Scalar(0, 255, 255) : Scalar(255, 255, 255));*/
			//p=0没找到边界黄色：p=0找到背景跳出红色
			for (int k = 0; k < STEP; ++k)
			{
				Vec3b &vec0 = _scanLineRet.at<Vec3b>(Point(m_Point[j][k].x, m_Point[j][k].y));
				DrawNop(vec0, flag);
			}
		}
		if (p != 0)
		{
			if (abs(p - min_position) < 6)
			{
				m_dCost[j][p] = std::min(1 * min_Cost, m_dCost[j][p]);
			}
			else
			{
				m_dCost[j][p] = std::min(1. / PARAM * min_Cost, m_dCost[j][p]);
			}

		}
	}

	return 0;
}

void MinSegInvoker::MinPath()
{
	m_dFinalcost = 1e100;

	double cost[GRADE][STEP];
	int road[GRADE][STEP];
	memset(road, 0, sizeof(int) * GRADE * STEP);
	for (int k = 0; k < STEP; k++)
	{
		for (int i = 0; i < GRADE; i++)
		{
			for (int j = 0; j < STEP; j++)
			{
				cost[i][j] = 1e100;
			}
		}
		int begin, end;
		/*第0方向代价*/
		cost[0][k] = m_dCost[0][k];
		if (k < Dis)
		{
			begin = 0;
		}
		else
		{
			begin = k - (Dis - 1);
		}
		if (k < STEP - Dis)
		{
			end = k + Dis;
		}
		else
		{
			end = STEP;
		}
		/*第1方向代价*/
		for (int n = 0; n < STEP; n++)
		{
			cost[1][n] = 1e100;
		}
		for (int n = begin; n < end; n++)
		{
			cost[1][n] = cost[0][k] + WeightFunc(abs(n - k)) * m_dCost[1][n];
			road[1][n] = k;
		}

		for (int i = 2; i < GRADE; i++)
		{
			for (int j = 0; j < STEP; j++)
			{
				//计算起点到点[i][j]的最小代价
				double curMin = MAXMAXDATA;

				if (j < Dis)
				{
					begin = 0;
				}
				else
				{
					begin = j - (Dis - 1);
				}
				if (j < STEP - Dis)
				{
					end = j + Dis;
				}
				else
				{
					end = STEP;
				}
				for (int m = begin; m < end; m++)
				{
					double curCost = cost[i - 1][m] + WeightFunc(abs(m - j)) * m_dCost[i][j];
					if (curCost < curMin)
					{
						curMin = curCost;
						cost[i][j] = curMin;
						road[i][j] = m;
					}
				}
			}
		}

		//回到起点最小代价
		double curCostRoad;
		double curMinRoad = MAXMAXDATA;

		if (k < Dis)
		{
			begin = 0;
		}
		else
		{
			begin = k - (Dis - 1);
		}
		if (k < STEP - Dis)
		{
			end = k + Dis;
		}
		else
		{
			end = STEP;
		}
		for (int j = begin; j < end; j++)
		{
			curCostRoad = cost[0][k] + cost[GRADE - 1][j];
			if (curCostRoad < curMinRoad)
			{
				curMinRoad = curCostRoad;
				road[0][k] = j;
			}
		}
		/*计算最短路径*/
		if (curMinRoad < m_dFinalcost)
		{
			m_dFinalcost = curMinRoad;
			int l_row, l_col;
			l_row = 0;
			l_col = k;
			for (int i = 0; i < GRADE; i++)
			{
				int row, col;
				row = GRADE - 1 - i;
				col = road[l_row][l_col];
				m_nMinroad[i] = col;
				l_row = row;
				l_col = col;
			}
		}
	}
}

void DrawP(Vec3b &vec)
{
	vec.val[0] = (UCHAR)255;
	vec.val[1] = (UCHAR)0;
	vec.val[2] = (UCHAR)255;
}

int MinSegInvoker::DrawEdge(Mat &_img)
{
	_img = m_Img.clone();
	//cvtColor(_img, _color, CV_GRAY2RGB);

	for (int i = 0; i < GRADE; ++i)
	{
		Vec3b &vec0 = _img.at<Vec3b>(m_Point[GRADE - i - 1][m_nMinroad[i]].y, m_Point[GRADE - i - 1][m_nMinroad[i]].x);
		Vec3b &vec1 = _img.at<Vec3b>(m_Point[GRADE - i - 1][m_nMinroad[i]].y, m_Point[GRADE - i - 1][m_nMinroad[i]].x - 1);
		Vec3b &vec2 = _img.at<Vec3b>(m_Point[GRADE - i - 1][m_nMinroad[i]].y - 1, m_Point[GRADE - i - 1][m_nMinroad[i]].x);
		Vec3b &vec3 = _img.at<Vec3b>(m_Point[GRADE - i - 1][m_nMinroad[i]].y, m_Point[GRADE - i - 1][m_nMinroad[i]].x + 1);
		Vec3b &vec4 = _img.at<Vec3b>(m_Point[GRADE - i - 1][m_nMinroad[i]].y + 1, m_Point[GRADE - i - 1][m_nMinroad[i]].x);
		DrawP(vec0);
		DrawP(vec1);
		DrawP(vec2);
		DrawP(vec3);
		DrawP(vec4);
	}

	return 0;
}

int MinSegInvoker::Seg(Mat &_costRet, Mat &_superPixelCostRet, Mat &_scanLineRet, Mat &_img)
{
	return DrawEdge(_img);
}