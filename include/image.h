/*!
@file image.h

@brief 画像関係

@author Makoto Fujisawa
@date 2025-02
*/

#ifndef _RX_IMAGE_H_
#define _RX_IMAGE_H_


//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------
#include <string>



//-----------------------------------------------------------------------------
// Function declarations
//-----------------------------------------------------------------------------
/*!
* 画像データ読み込み
*  - stb_image : https://github.com/nothings/stb/blob/master/stb_image.h
*    を使用
*  - BMP, PNG, JPEG, TGA, GIF, PIC, PNM, PSDに対応
* @param[in] fn 画像ファイル名
* @param[out] w,h 画像の横と縦のピクセル数
* @param[out] c 画像の色数(RGBで3,RGBAで4)
* @param[out] data 画像データ(unsigned char型)
*/
int ReadImage(const std::string &fn, int &w, int &h, int &c, unsigned char* &data);



/*!
* 画像データをファイル保存
*  - stb_image_write : https://github.com/nothings/stb/blob/master/stb_image_write.h
*    を使用
*  - BMP, PNG, JPEG, TGAに対応
* @param[in] fn,ext 保存ファイル名と拡張子(拡張子は保存ファイル名にも含まれている)
* @param[in] w,h 画像の横と縦のピクセル数
* @param[in] c 画像の色数(RGBで3,RGBAで4)
* @param[in] data 画像データ(unsigned char型)
* @param[in] quality JPEG保存時の画質([0,100]で値が大きいほど高品質,JPEG以外の時には使わない)
*/
int SaveImage(const std::string &fn, const std::string &ext, int w, int h, int c, const unsigned char* data, int quality = 95);

/*!
* 画像データをファイル保存(HDR型)
*  - stb_image_write : https://github.com/nothings/stb/blob/master/stb_image_write.h
*    を使用
* @param[in] fn 保存ファイル名
* @param[in] w,h 画像の横と縦のピクセル数
* @param[in] c 画像の色数(RGBで3,RGBAで4)
* @param[in] data 画像データ(float型)
*/
int SaveImageHDR(const std::string &fn, int w, int h, int c, const float* data);


#endif // #ifndef _RX_IMAGE_H_