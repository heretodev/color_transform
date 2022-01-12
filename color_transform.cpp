#include <chrono>
#include <array>
#include <map>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "CImg.h"
using namespace cimg_library;

// TODO: great documentation for perceived luminance: https://stackoverflow.com/a/56678483
//  public domain function by Darel Rex Finley, 2006
//
//  This function expects the passed-in values to be on a scale
//  of 0 to 1, and uses that same scale for the return values.
//
//  See description/examples at alienryderflex.com/hsp.html
#define  Pr  .299
#define  Pg  .587
#define  Pb  .114
const std::array<double,3> rgb_to_hsp_pixel_map(const std::array<uint8_t,3>& in)
{
	double H,S,P, R = static_cast<double>(in[0])/255., G = static_cast<double>(in[1])/255., B = static_cast<double>(in[2])/255.;
  //  Calculate the Perceived brightness.
  P=sqrt(R*R*Pr+G*G*Pg+B*B*Pb);

  //  Calculate the Hue and Saturation.  (This part works
  //  the same way as in the HSV/B and HSL systems???.)
  if      (R==G && R==B) {
    H=0.; S=0.; return {H,S,P}; }
  if      (R>=G && R>=B) {   //  R is largest
    if    (B>=G) {
      H=6./6.-1./6.*(B-G)/(R-G); S=1.-G/R; }
    else         {
      H=0./6.+1./6.*(G-B)/(R-B); S=1.-B/R; }}
  else if (G>=R && G>=B) {   //  G is largest
    if    (R>=B) {
      H=2./6.-1./6.*(R-B)/(G-B); S=1.-B/G; }
    else         {
      H=2./6.+1./6.*(B-R)/(G-R); S=1.-R/G; }}
  else                   {   //  B is largest
    if    (G>=R) {
      H=4./6.-1./6.*(G-R)/(B-R); S=1.-R/B; }
    else         {
      H=4./6.+1./6.*(R-G)/(B-G); S=1.-G/B; }}
	return {H,S,P};
}

//  public domain function by Darel Rex Finley, 2006
//
//  This function expects the passed-in values to be on a scale
//  of 0 to 1, and uses that same scale for the return values.
//
//  Note that some combinations of HSP, even if in the scale
//  0-1, may return RGB values that exceed a value of 1.  For
//  example, if you pass in the HSP color 0,1,1, the result
//  will be the RGB color 2.037,0,0.
//
//  See description/examples at alienryderflex.com/hsp.html

const std::array<uint8_t,3> hsp_to_rgb_pixel_map(const std::array<double_t,3>& in)
{
	double R,G,B, H = in[0], S = in[1], P = in[2];
  double  part, minOverMax=1.-S ;

  if (minOverMax>0.) {
    if      ( H<1./6.) {   //  R>G>B
      H= 6.*( H-0./6.); part=1.+H*(1./minOverMax-1.);
      B=P/sqrt(Pr/minOverMax/minOverMax+Pg*part*part+Pb);
      R=(B)/minOverMax; G=(B)+H*((R)-(B)); }
    else if ( H<2./6.) {   //  G>R>B
      H= 6.*(-H+2./6.); part=1.+H*(1./minOverMax-1.);
      B=P/sqrt(Pg/minOverMax/minOverMax+Pr*part*part+Pb);
      G=(B)/minOverMax; R=(B)+H*((G)-(B)); }
    else if ( H<3./6.) {   //  G>B>R
      H= 6.*( H-2./6.); part=1.+H*(1./minOverMax-1.);
      R=P/sqrt(Pg/minOverMax/minOverMax+Pb*part*part+Pr);
      G=(R)/minOverMax; B=(R)+H*((G)-(R)); }
    else if ( H<4./6.) {   //  B>G>R
      H= 6.*(-H+4./6.); part=1.+H*(1./minOverMax-1.);
      R=P/sqrt(Pb/minOverMax/minOverMax+Pg*part*part+Pr);
      B=(R)/minOverMax; G=(R)+H*((B)-(R)); }
    else if ( H<5./6.) {   //  B>R>G
      H= 6.*( H-4./6.); part=1.+H*(1./minOverMax-1.);
      G=P/sqrt(Pb/minOverMax/minOverMax+Pr*part*part+Pg);
      B=(G)/minOverMax; R=(G)+H*((B)-(G)); }
    else               {   //  R>B>G
      H= 6.*(-H+6./6.); part=1.+H*(1./minOverMax-1.);
      G=P/sqrt(Pr/minOverMax/minOverMax+Pb*part*part+Pg);
      R=(G)/minOverMax; B=(G)+H*((R)-(G)); }}
  else {
    if      ( H<1./6.) {   //  R>G>B
      H= 6.*( H-0./6.); R=sqrt(P*P/(Pr+Pg*H*H)); G=(R)*H; B=0.; }
    else if ( H<2./6.) {   //  G>R>B
      H= 6.*(-H+2./6.); G=sqrt(P*P/(Pg+Pr*H*H)); R=(G)*H; B=0.; }
    else if ( H<3./6.) {   //  G>B>R
      H= 6.*( H-2./6.); G=sqrt(P*P/(Pg+Pb*H*H)); B=(G)*H; R=0.; }
    else if ( H<4./6.) {   //  B>G>R
      H= 6.*(-H+4./6.); B=sqrt(P*P/(Pb+Pg*H*H)); G=(B)*H; R=0.; }
    else if ( H<5./6.) {   //  B>R>G
      H= 6.*( H-4./6.); B=sqrt(P*P/(Pb+Pr*H*H)); R=(B)*H; G=0.; }
    else               {   //  R>B>G
      H= 6.*(-H+6./6.); R=sqrt(P*P/(Pr+Pb*H*H)); B=(R)*H; G=0.; }}
	return {static_cast<uint8_t>(round(R*255)),static_cast<uint8_t>(round(G*255)),static_cast<uint8_t>(round(B*255))};
}

const std::array<uint8_t,3> hsl_to_rgb_pixel_map(const std::array<float,3>& in)
{
	float h = in[0], s = in[1], l = in[2];
/* May not be necessary to validate hsl values:
	if(h < 0) h=0;
	if(s < 0) s=0;
	if(l < 0) l=0;
	if(h >= 360) h=359;
	if(s > 100) s=100;
	if(l > 100) l=100;
*/
	s /= 100.0f;
	l /= 100.0f;
	const float C = (1 - fabs(2.0f * l - 1)) * s;
	const float hh = h / 60;
	const float float_remainder_of_2 = hh - (floor(hh / 2) * 2); // it's % 2 but keep the fraction
	const float X = C * (1 - fabs(float_remainder_of_2 - 1));
	const float m = l - C / 2.0f;
	std::array<float,3> rgb;
	if(h < 60)
		rgb = {C,X,0};
	else if(h < 120)
		rgb = {X,C,0};
	else if(h < 180)
		rgb = {0,C,X};
	else if(h < 240)
		rgb = {0,X,C};
	else if(h < 300)
		rgb = {X,0,C};
	else
		rgb = {C,0,X};
	return {static_cast<uint8_t>(round((rgb[0] + m) * 255)), static_cast<uint8_t>(round((rgb[1] + m) * 255)), static_cast<uint8_t>(round((rgb[2] + m) * 255))};
}

// Returns float HSL.  Since H ranges from 0 to 360, it doesn't fit into uint8_t range.
const std::array<float,3> rgb_to_hsl_pixel_map(const std::array<uint8_t,3>& in)
{
	float r = static_cast<float>(in[0]);
	float g = static_cast<float>(in[1]);
	float b = static_cast<float>(in[2]);
/* may not be needed:
	if( r<0 ) r=0;
	if( g<0 ) g=0;
	if( b<0 ) b=0;
	if( r>255 ) r=255;
	if( g>255 ) g=255;
	if( b>255 ) b=255;
*/
	r/=255.0f;
	g/=255.0f;
	b/=255.0f;
	const std::array<float,3> pixel = {r,g,b};
	const std::pair<const float*,const float*> min_max = std::minmax_element(pixel.begin(), pixel.end());
	const float min = *min_max.first, max = *min_max.second;
	const float delta = max - min;
	std::array<float,3> hsl = {0,0,(max + min) / 2.0f};
	if(max == 0) // TODO: floating point equality!
		hsl[0] = 0;
	else if(max == r)
		hsl[0] = ((g - b) / delta) + (g < b ? 6 : 0); // justification: https://stackoverflow.com/a/39147465
	else if(max == g)
		hsl[0] = ((b - r) / delta) + 2;
	else if(max == b)
		hsl[0] = ((r - g) / delta) + 4;
	hsl[0] *= 60; // Note: 60 is in degrees
	if(delta > 0.000001) // is delta != 0
		hsl[1] = 100.0f * (delta / (1 - fabs(2 * hsl[2] - 1)));
	hsl[2] *= 100.0f;
	return hsl;
}

template<typename T1, typename T2, typename F>
void cimg_color_transform(const CImg<T1>& input, CImg<T2>& output, F f)
{
	cimg_forXY(output,x,y)
	{
		std::array<T1,3> in_pixel{input(x,y,0),input(x,y,1),input(x,y,2)};
		std::array<T2,3> out_pixel = f(in_pixel);
		output(x,y,0) = out_pixel[0];
		output(x,y,1) = out_pixel[1];
		output(x,y,2) = out_pixel[2];
	}
}

// Compare outputs to: https://www.rapidtables.com/convert/color/hsl-to-rgb.html
void verify_hsl_to_rgb(const std::array<float,3>& hsl_pixel)
{
	const std::array<uint8_t,3> rgb_pixel = hsl_to_rgb_pixel_map(hsl_pixel);
	std::cout << "RGB from HSL:" << static_cast<int>(rgb_pixel[0]) << "," << static_cast<int>(rgb_pixel[1]) << "," << static_cast<int>(rgb_pixel[2]) << std::endl;
}

// Compare outputs to: https://www.rapidtables.com/convert/color/rgb-to-hsl.html
void verify_rgb_to_hsl(const std::array<uint8_t,3>& rgb_pixel = {228,224,221}) // the a whiteboard hsl to rgb
{
	const std::array<float,3> hsl_pixel = rgb_to_hsl_pixel_map(rgb_pixel);
	std::cout << "HSL from RGB:" << hsl_pixel[0] << "," << hsl_pixel[1] << "," << hsl_pixel[2] << std::endl;
}

void verify_rgb_to_hsp(const std::array<uint8_t,3>& rgb_pixel = {228,224,221}) // the a whiteboard hsl to rgb
{
	const std::array<double,3> hsp_pixel = rgb_to_hsp_pixel_map(rgb_pixel);
	std::cout << "HSP from RGB:" << hsp_pixel[0] << "," << hsp_pixel[1] << "," << hsp_pixel[2] << std::endl;
}

int main(int argc, char * argv[])
{
    // image loading
	if(argc < 2)
	{
		std::cout << "Usage: specify an input image." << std::endl << "./color_transform <image_location>" << std::endl;
		return 1;
	}

	CImg<uint8_t> image(argv[1]);
	std::cerr << "Read the image" << ": width:" << image.width() << " x " << image.height() << " and channels:" << image.spectrum() << " and pixel_type:" << image.pixel_type() << std::endl;
	// Convert to HSL:
	if(image.spectrum() > 3)
	{
		image.channels(0,2);
	}
	CImg<double> hsp_image(image.width(),image.height(),1,image.spectrum());
	cimg_color_transform(image, hsp_image, rgb_to_hsp_pixel_map);

	double mean_p = hsp_image.get_channel(2).mean(); // or .median()

/*
	CImg<float> quantized_hsl_image(image.width(),image.height(),1,image.spectrum());
    auto start = std::chrono::high_resolution_clock::now();
	int niters = 1;
	for(int iter = 0; iter < niters; iter++)
	{
		cimg_color_transform(hsl_image, quantized_hsl_image, hsl_pixel_map);
	}
    // stats
    auto end = std::chrono::high_resolution_clock::now();
    float elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() / static_cast<float>(niters);
    std::cerr << "time " << elapsed << "ms" << std::endl;
*/
    // save
	CImg<uint8_t> rgb_image(image.width(),image.height(),1,image.spectrum());
	cimg_color_transform(hsp_image, rgb_image, hsp_to_rgb_pixel_map);
	rgb_image.save("out.jpg");

    return 0;
}
