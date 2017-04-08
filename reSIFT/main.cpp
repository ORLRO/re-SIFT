//TODO: Do it in class
#define cimg_use_jpeg
#include "CImg.h"
#include <vector>
#include <cmath>
#include <iostream>
using namespace cimg_library;
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////Prototypes//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isExtreme(const CImg<float> &neighborhood);

CImg<float> GaussianKernel(float sigma, unsigned int HalfWidth);

float dominantOrientation(const CImg<float>& MagNeighborhood, const CImg<float>& DirNeighborhood,
	unsigned int histBinSize);

void computeScaleSpace();

void compute_DoG_pyramid();

void computeGradientPyramid();

void computeHarrisPyramid(float Threshold);

void detectExtremePoints(float threshold);

void computeOrientaions(unsigned int histBinSize);

void constructDescriptors(unsigned int histBinSize, unsigned int blockSize, unsigned int subBlockSize);

void markFeatures();

class SIFTFeature
{
public:
	CImg<float> histograms;
	float scale;
	float _x, _y;
};
SIFTFeature constructDescriptor(CImg<float> MagNeighborhood, CImg<float> DirNeighborhood,
	unsigned int histBinSize, unsigned int subBlockSize);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////Global Variables//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CImg<float> I, disp;
unsigned int quadrans = 4, nsigmas = 5; //number of quadrants
float sigmas[] = { sqrtf(2) / 2, 1.0f , sqrtf(2), 2.0f , 2 * sqrtf(2) };
//vector<CImgList<float> > scale_space;
CImgList<float> gaussian_pyramid;
CImgList<float> scale_space(quadrans);
CImgList<float> DoG_pyramid(quadrans);
CImgList<float> IxPyramid(quadrans), IyPyramid(quadrans), magPyramid(quadrans), dirPyramid(quadrans);
CImgList<float> Sx2(quadrans), Sy2(quadrans), Sxy(quadrans);
CImgList<bool> HarrisPyramid(quadrans);
CImgList<bool> intrestPyramid(quadrans);
CImgList<float> anglesPyramid(quadrans);
vector<SIFTFeature> descriptorList;
const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////Main//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
	I.assign("C:\\Users\\DRLR\\OneDrive\\_Master\\Computer Architecture\\Project\\sift_cimg_visualstudio\\cameraman.jpg");
	disp = I;
	I.channel(0);
	I.resize_halfXY();

	//scale space
	cout << "scale space...";
	computeScaleSpace();
	cout << "done\n";

	//difference of Gaussian pyramid
	cout << "Dog Pyramid...";
	compute_DoG_pyramid();
	cout << "done\n";

	//Derivative pyramid
	cout << "Gradient...";
	computeGradientPyramid();
	cout << "done\n";

	//harris edge detector
	cout << "Harris...";
	float Threshold = 5.0f;
	computeHarrisPyramid(Threshold);
	cout << "done\n";

	//free Ix and Iy pyramids
	IxPyramid.assign();
	IyPyramid.assign();


	//extema detedction & low contrast rejection & non corner rejection 
	cout << "extema detedction...";
	float threshold = 8.0f; //3% of 255
	detectExtremePoints(threshold);
	cout << "done\n";

	//free harris pyramid and DoG pyramid
	HarrisPyramid.assign();
	DoG_pyramid.assign();

	//orientation assignment
	cout << "orientation assignment...";
	unsigned int histBinSize = 36;
	computeOrientaions(histBinSize);
	cout << "done\n";

	//descriptor construction 
	histBinSize = 8; //from SIFT paper section 6.1
	unsigned int blockSize = 16;  // SIFT discripteor width
	unsigned int subBlockSize = 4;
	cout << "constructDescriptors...";
	constructDescriptors(blockSize, subBlockSize, histBinSize);
	cout << "done\n";

	cout << "#features: " << descriptorList.size() << endl;
	//display	
	markFeatures();
	CImgDisplay main_disp(disp, "Display");
	while (!main_disp.is_closed())
	{
		main_disp.wait();
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////Functions//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeScaleSpace()
{
	//constructing gaussian pyramid
	gaussian_pyramid.insert(I);
	for (unsigned int i = 1; i < quadrans; i++)
	{
		gaussian_pyramid.insert(gaussian_pyramid[i - 1].get_blur(1.0, true, true).get_resize(I.width() / 2, I.height() / 2, 1, 1, 1));
	}
	// multiple scales in each quadrant
	cimglist_for(gaussian_pyramid, level)
	{
		scale_space[level].assign(gaussian_pyramid[level]);
		for (unsigned int i = 1; i <= nsigmas; i++)
		{
			scale_space[level].append(scale_space[level].get_slice(i - 1).get_blur(sigmas[i - 1], true, true), 'z');
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_DoG_pyramid()
{
	for (unsigned int i = 0; i <quadrans; i++)
	{
		DoG_pyramid[i].assign(scale_space[i].get_slice(1) - scale_space[i].get_slice(0));
		for (unsigned int ii = 1; ii < nsigmas - 1; ii++)
		{
			DoG_pyramid[i].append(scale_space[i].get_slice(ii + 1) - scale_space[i].get_slice(ii), 'z');
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeGradientPyramid()
{
	//TODO: can be done without loop
	for (unsigned int i = 0; i <quadrans; i++)
	{
		IxPyramid[i].assign(scale_space[i].get_slice(0).get_gradient("x", 0)[0]);
		IyPyramid[i].assign(scale_space[i].get_slice(0).get_gradient("y", 0)[0]);

		magPyramid[i].assign((IxPyramid[i].get_slice(0).sqr() + IyPyramid[i].get_slice(0).sqr()).sqrt());
		dirPyramid[i].assign(IyPyramid[i].get_slice(0).get_atan2(IxPyramid[i].get_slice(0)));
		for (unsigned int ii = 1; ii < nsigmas - 1; ii++)
		{
			IxPyramid[i].append(scale_space[i].get_slice(ii).get_gradient("x", 0)[0], 'z');
			IyPyramid[i].append(scale_space[i].get_slice(ii).get_gradient("y", 0)[0], 'z');

			magPyramid[i].append((IxPyramid[i].get_slice(ii).sqr() + IyPyramid[i].get_slice(ii).sqr()).sqrt(), 'z');
			dirPyramid[i].append(IyPyramid[i].get_slice(ii).get_atan2(IxPyramid[i].get_slice(ii)), 'z');
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeHarrisPyramid(float Threshold)
{//TODO:try without assign then append
	for (unsigned int i = 0; i <quadrans; i++)
	{
		Sx2[i].assign(IxPyramid[i].get_slice(0).get_blur(sigmas[0], true, true).sqr().normalize(0.0f, 255.0f));
		Sy2[i].assign(IyPyramid[i].get_slice(0).get_blur(sigmas[0], true, true).sqr().normalize(0.0f, 255.0f));
		Sxy[i].assign(IxPyramid[i].get_slice(0).get_mul(IxPyramid[i].get_slice(0)).blur(sigmas[0], true, true).normalize(0.0f, 255.0f));
		CImg<bool> img = ((Sx2[i].get_slice(0).get_mul(Sy2[i].get_slice(0)) - Sxy[i].get_slice(0).get_sqr()).
			get_div(Sx2[i].get_slice(0) + Sy2[i].get_slice(0))).threshold(Threshold);
		HarrisPyramid[i].assign(img);
		for (unsigned int ii = 1; ii < nsigmas - 1; ii++)
		{
			Sx2[i].append(IxPyramid[i].get_slice(ii).get_blur(sigmas[ii], true, true).sqr().normalize(0.0f, 255.0f), 'z');
			Sy2[i].append(IyPyramid[i].get_slice(ii).get_blur(sigmas[ii], true, true).sqr().normalize(0.0f, 255.0f), 'z');
			Sxy[i].append(IxPyramid[i].get_slice(ii).get_mul(IxPyramid[i].get_slice(ii)).blur(sigmas[ii], true, true).normalize(0.0f, 255.0f), 'z');

			//cim = (Ix2.*Iy2 - Ixy. ^ 2). / (Ix2 + Iy2 + eps); % Harris corner measure
			img = ((Sx2[i].get_slice(ii).get_mul(Sy2[i].get_slice(ii)) - Sxy[i].get_slice(ii).get_sqr()).
				get_div(Sx2[i].get_slice(ii) + Sy2[i].get_slice(ii))).get_threshold(Threshold);
			HarrisPyramid[i].append(img, 'z');
		}
		Sx2[i].assign();
		Sy2[i].assign();
		Sxy[i].assign();
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void detectExtremePoints(float threshold)
{
	for (unsigned int i = 0; i < quadrans; i++)
	{
		intrestPyramid[i].assign(CImg<bool>(DoG_pyramid[i], "xyzc", false));
		CImg<float> neighborhood(3, 3, 3, 1);
		cimg_for3x3x3(DoG_pyramid[i], x, y, z, 0, neighborhood, float)
		{
			intrestPyramid[i](x, y, z) = isExtreme(neighborhood)
				&& DoG_pyramid[i](x, y, z) > threshold
				&& HarrisPyramid[i](x, y, z);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeOrientaions(unsigned int histBinSize)
{
	for (unsigned int i = 0; i < quadrans; i++)
	{
		anglesPyramid[i].assign(CImg<float>(intrestPyramid[i], "xyzc", 0.0f));
		for (unsigned int ii = 0; ii < intrestPyramid[i]._depth; ii++)
		{
			unsigned int HalfWidth = static_cast<unsigned int>(ceil(3 * sigmas[ii + 1])); //neighborhood half width
			CImg<float> GaussianWeight = GaussianKernel(sigmas[ii + 1], HalfWidth);
			cimg_for_insideXY(intrestPyramid[i], x, y, HalfWidth + 1)
			{
				if (intrestPyramid[i](x, y, ii))
				{
					CImg<float> MagNeighborhood = GaussianWeight.mul(magPyramid[i].slice(ii).get_crop(x - HalfWidth, y - HalfWidth, x + HalfWidth, y + HalfWidth));
					CImg<float> DirNeighborhood = dirPyramid[i].slice(ii).get_crop(x - HalfWidth, y - HalfWidth, x + HalfWidth, y + HalfWidth);
					anglesPyramid[i](x, y, ii) = dominantOrientation(MagNeighborhood, DirNeighborhood, histBinSize);
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void constructDescriptors(unsigned int histBinSize, unsigned int blockSize, unsigned int subBlockSize)
{
	//16 x 16 Neighborhood + margin to allow for safe rotation
	//(this is worst case -> can be reduced if dominant orientation is not multiple of 45deg)
	unsigned int HalfWidth = (unsigned int)ceil((blockSize + 1) / 2 * sqrt(2)); //neighborhood half width
	for (unsigned int i = 0; i < quadrans; i++)
	{
		//TODO: remove this loop
		for (unsigned int ii = 0; ii < intrestPyramid[i]._depth; ii++)
		{
			cimg_for_insideXY(intrestPyramid[i], x, y, HalfWidth + 1)
			{
				if (intrestPyramid[i](x, y, ii))
				{
					//16 x 16 Neighborhood + margin to allow for safe rotation
					CImg<float> Neighborhood = scale_space[i].slice(ii)
						.get_crop(x - HalfWidth, y - HalfWidth, x + HalfWidth, y + HalfWidth)
						//rotate Neighborhood
						.rotate(-1 * anglesPyramid[i](x, y, ii), 1, 0)
						//take center 16 x 16 region
						.crop(Neighborhood.width() / 2 - blockSize / 2 + 1,
							Neighborhood.width() / 2 + blockSize / 2,
							Neighborhood.height() / 2 - blockSize / 2 + 1,
							Neighborhood.height() / 2 + blockSize / 2);
					//gradiaent direction and magnitude
					CImgList<float> DxDy = Neighborhood.get_gradient("xy");
					CImg<float> MagNeighborhood = (DxDy[0].sqr() + DxDy[1].sqr()).sqrt();
					CImg<float> DirNeighborhood = DxDy[1].get_atan2(DxDy[0]);

					//construct the histogram array
					SIFTFeature descriptor = constructDescriptor(MagNeighborhood, DirNeighborhood, histBinSize, subBlockSize);
					descriptor._x = x*pow(2.0f, (float)i);
					descriptor._y = y*pow(2.0f, (float)i);
					descriptor.scale = pow(2.0f, sigmas[ii]);
					descriptorList.push_back(descriptor);
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void markFeatures()
{
	//TODO: add circles to the image
	for (int i = 0; i < descriptorList.size(); i++)
	{
		disp.draw_circle(2 * descriptorList[i]._x, 2 * descriptorList[i]._y,
			2 * 3 * descriptorList[i].scale, red, 1.0f, 2);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SIFTFeature constructDescriptor(CImg<float> MagNeighborhood, CImg<float> DirNeighborhood,
	unsigned int histBinSize, unsigned int subBlockSize)
{//returns the descriptor of sn intrest point SIFT paper section 6.1
	float minv = static_cast<float>(-cimg::PI);
	float maxv = static_cast<float>(+cimg::PI);
	CImg<float> histograms(MagNeighborhood.width() / subBlockSize, MagNeighborhood.width() / subBlockSize, histBinSize);
	for (unsigned int i = 0; i < MagNeighborhood.width() / subBlockSize; i++)
	{
		for (unsigned int ii = 0; ii < MagNeighborhood.width() / subBlockSize; ii++)
		{
			CImg<float> MagSubBlock = MagNeighborhood.get_crop(i, i, i + subBlockSize - 1, i + subBlockSize - 1);
			CImg<float> DirSubBlock = DirNeighborhood.get_crop(i, i, i + subBlockSize - 1, i + subBlockSize - 1);
			cimg_forXY(DirSubBlock, x, y)
			{
				int bin = int((DirSubBlock(x, y) - minv) / (maxv - minv) * histBinSize) % histBinSize;
				histograms(i, ii, bin) += MagSubBlock(x, y);
			}
		}
	}
	SIFTFeature feature;
	feature.histograms = histograms.normalize();
	return feature;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isExtreme(const CImg<float> &neighborhood)
{
	CImg<float> cube(neighborhood);
	float center = cube(2, 2, 2);
	cube(2, 2, 2) = static_cast<float>(cube.mean());
	return center > cube.max() || center < cube.min();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CImg<float> GaussianKernel(float sigma, unsigned int HalfWidth)
{
	// Compute gaussian kernel
	double color = 1;
	CImg<float> gaussian(2 * HalfWidth + 1, 2 * HalfWidth + 1);
	gaussian.draw_gaussian((float)HalfWidth, (float)HalfWidth, sigma, &color);
	return gaussian;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float dominantOrientation(const CImg<float>& MagNeighborhood, const CImg<float>& DirNeighborhood,
	unsigned int histBinSize)
{
	float minv = static_cast<float>(-cimg::PI);
	float maxv = static_cast<float>(+cimg::PI);
	vector<float> hist(histBinSize);
	cimg_forXY(DirNeighborhood, x, y)
	{
		int bin = int((DirNeighborhood(x, y) - minv) / (maxv - minv) * histBinSize) % histBinSize;
		hist[bin] += MagNeighborhood(x, y);
	}
	int index = distance(hist.begin(), max_element(hist.begin(), hist.end()));
	return ((index / histBinSize - 0.5f) - 1 / histBinSize) * 360.0f;
}