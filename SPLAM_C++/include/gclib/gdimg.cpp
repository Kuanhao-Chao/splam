#include "gdimg.h"

int GDImg::defaultBg=0xFFFFFF;

void GDImg::init(int w, int h, const char* fn, int bg_rgb) {
 fout=NULL;
 imgH=h;
 imgW=w;
 img = gdImageCreate(imgW, imgH);
 currentColor=-1;
 if (fn!=NULL) setFile(fn);
 //allocate white background color by default
 bgColor=color(bg_rgb);
 gdImageColorTransparent(img, -1); //default is non-transparent white background
}

void GDImg::setBg(int rgb) {
    //if (rgb!=defaultBg) bgColor=color(rgb);
    gdImageColorDeallocate(img, bgColor);
    bgColor=this->color(rgb);
    }

GDImg::~GDImg() {
 gdImageDestroy(img);
 if (fout!=NULL && fout!=stdout) fclose(fout);
}

int GDImg::color(byte r, byte g, byte b) {
  return gdImageColorAllocate(img,(int)r,(int)g,(int)b);
}

void GDImg::line(int x1, int y1, int x2, int y2, int color) {
  if (color==-1) color=currentColor;
  gdImageLine(img,x1,y1,x2,y2,color);
}

void GDImg::rectangle(int x1, int y1, int x2, int y2, int color) {
  if (color==-1) color=currentColor;
  gdImageRectangle(img,x1,y1,x2,y2,color);
}

void GDImg::fillRectangle(int x1, int y1, int x2, int y2, int color) {
  if (color==-1) color=currentColor;
  gdImageFilledRectangle(img,x1,y1,x2,y2,color);
}


void GDImg::fillPolygon(gdPointPtr points, int ptotal, int color) {
  if (color==-1) color=currentColor;
  gdImageFilledPolygon(img, points,ptotal,color);
}

void GDImg::setTransparent(int cidx) {
  //cidx must be the color index of a color allocated previously for img!
  gdImageColorTransparent(img,cidx);
}

void GDImg::setFile(const char* fname) {
  if (fout!=NULL && fout!=stdout) fclose(fout);
  if (fname[0]=='-' && fname[1]==0) {
    //special "-" file name means stdout
    fout=stdout;
    }
  else {
    fout=fopen(fname, "wb");
    if (fout==NULL) GError("Error: cannot open file %s for writing!\n",fname);
    }
}

void GDImg::setFile(FILE* f) {
  if (fout!=NULL && fout!=stdout) fclose(fout);
  fout=f;
}

void GDImg::write(const char* fname) {
 if (fname==NULL && fout==NULL)
     GError("Error at GDImg::writeGIF() - no destination file given!\n");
 if (fname!=NULL) setFile(fname);
    gdImageGif(img,fout);
}
