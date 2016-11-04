#include <iostream>
#include <ctime>
using namespace std;

#include "MinSegInvoker.h"

int main()
{
	const string srcPath("C:/los/Proj/SuperPixelSrc/");
	const string dstPath("ret");

	clock_t a, b, c = 0;
	int count = 0;
	/*输入图片总数*/
	cin >> count;

	for (int i = 0; i < count; ++i)
	{
		//读取图像
		Mat superPixelRet;
		{
			ostringstream oss;
			oss << i << ".bmp";
			superPixelRet = imread(srcPath + oss.str(), CV_LOAD_IMAGE_GRAYSCALE);
		}

		Mat colorImg;
		{
			ostringstream oss;//绘制结果到彩色图像
			oss << "c" << i << ".bmp";
			colorImg = imread(srcPath + oss.str(), CV_LOAD_IMAGE_COLOR);
		}

		a = clock();
		MinSegInvoker invoker(colorImg, superPixelRet);
		Mat costRet, superPixelCostRet, scanLineRet, ret;
		invoker.Seg(costRet, superPixelCostRet, scanLineRet, ret);
		b = clock();

		c += b - a;
		//输出结果
		{
			ostringstream oss;
			oss << "cost" << i << ".jpg";
			imwrite(dstPath + oss.str(), costRet);
		}
		{
			ostringstream oss;
			oss << "scost" << i << ".jpg";
			imwrite(dstPath + oss.str(), costRet);
		}
		{
			ostringstream oss;
			oss << "line" << i << ".jpg";
			imwrite(dstPath + oss.str(), costRet);
		}
		{
			ostringstream oss;
			oss << "ret" << i << ".jpg";
			imwrite(dstPath + oss.str(), costRet);
		}		
	}

	cout << double(c) / CLOCKS_PER_SEC << "s for " << count << "imgs" << endl;
	return 0;
}