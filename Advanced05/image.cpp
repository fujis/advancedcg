/*!
@file image.cpp

@brief 画像関係

@author Makoto Fujisawa
@date   2025-02
*/


//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------
#include "image.h"

// 画像読み書き
#define STB_IMAGE_WRITE_STATIC
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


/*!
* 画像データ読み込み
*/
int ReadImage(const std::string &fn, int &w, int &h, int &c, unsigned char* &data)
{
	int r = stbi_info(fn.c_str(), &w, &h, &c);
	data = stbi_load(fn.c_str(), &w, &h, &c, c);
	return r;
}

/*!
* 画像データをファイル保存
*/
int SaveImage(const std::string &fn, const std::string &ext, int w, int h, int c, const unsigned char* data, int quality)
{
	stbi_flip_vertically_on_write(1);

	int r = 0;
	if(ext == "bmp"){
		r = stbi_write_bmp(fn.c_str(), w, h, c, data);
	}
	else if(ext == "jpg" || ext == "jpeg"){
		r = stbi_write_jpg(fn.c_str(), w, h, c, data, quality);
	}
	else if(ext == "png"){
		r = stbi_write_png(fn.c_str(), w, h, c, data, w*c);
	}
	else if(ext == "tga"){
		r = stbi_write_tga(fn.c_str(), w, h, c, data);
	}
	return r;
}

/*!
* 画像データをファイル保存(HDR型)
*/
int SaveImageHDR(const std::string &fn, int w, int h, int c, const float* data)
{
	stbi_flip_vertically_on_write(1);
	return stbi_write_hdr(fn.c_str(), w, h, c, data);
}

