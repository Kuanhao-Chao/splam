#ifndef _GDIMG_
#define _GDIMG_
#include "GBase.h"
#include "gd.h"


class GDImg {
 protected:
     gdImagePtr img;
     int imgW;
     int imgH;
     FILE* fout;
     int currentColor;
     int bgColor; //preallocated white by default
     void setFile(const char* fname);
     void setFile(FILE* f);
     static int defaultBg;
 public:

  void init(int w=64, int h=64, const char* fname=NULL, int bg_rgb=defaultBg);
  GDImg(int w=64, int h=64, const char* fname=NULL, int bg_rgb=defaultBg) {
	  init(w,h, fname, bg_rgb);
  }
  GDImg(int w,int h, int bg_rgb) {
     init(w,h, (const char*)NULL, bg_rgb);
     }
  ~GDImg();
  void write(const char* fname=NULL); //automatically write GIF
  void setTransparent(int cidx); // -1 means 'no transparency'
  void setTransparent(bool v=true) {
    setTransparent(v ? (int) bgColor : (int)-1);
    }
  int color(byte r, byte g, byte b);
  int color(int rgb) {
    return color( (byte)(rgb>>16) & 255,
                 (byte)(rgb>>8) & 255,
                 (byte)(rgb & 255)); }
  int colorAllocate(byte r, byte g, byte b) { return color(r,g,b); }
  int colorAllocate(int rgb) { return color(rgb); }
  void setColorIdx(int color) { currentColor=color; } //current color for drawing operations
  int setColor(int r, int g, int b) {
     currentColor=this->color(r,g,b);
     return currentColor;
     }
  int setColor(int rgb) {
    currentColor=this->color(rgb);
    return currentColor;
    }
  void setPixel(int x, int y, int color=-1) {
      if (color==-1) color=currentColor;
      gdImageSetPixel(img, x,y,color);
      }
  int getPixel(int x, int y) { return gdImageGetPixel(img, x, y); }
  void setBg(int rgb);
  void clear(int color=-1) {
    if (color==-1) color=bgColor;
    fillRectangle(0,0,imgW,imgH,color);
    }
  void line(int x1, int y1, int x2, int y2, int color=-1);
  void rectangle(int x1, int y1, int x2, int y2, int color=-1);
  void fillRectangle(int x1, int y1, int x2, int y2, int color=-1);
  void fillPolygon(gdPointPtr points, int ptotal, int color=-1);
};

#endif
