#ifndef __GCDB_H
#define __GCDB_H
#include "GBase.h"
#include <stddef.h>
#include <fcntl.h>

#ifdef _WIN32
  #define PROT_READ  1
  #define PROT_WRITE  2
  #define PROT_READWRITE  3
  #define MAP_SHARED  1
  #define MAP_PRIVATE  2
  #define F_OK 0
  #define R_OK 4
  #define W_OK 2
  #define RW_OK 6

  #ifndef MAP_FAILED
  #define MAP_FAILED      ((void *) -1)
  #endif
  void *mmap(char *,size_t,int,int,int,off_t);
  int   munmap(void *,size_t);
#else
  #include <sys/mman.h>
#endif

//=====================================================
//-------------     buffer stuff    -------------------
//=====================================================
#define GCDBUFFER_INSIZE 8192
#define GCDBUFFER_OUTSIZE 8192

typedef int (*opfunc)(int, char*, size_t);

//typedef unsigned long gcdb_seek_pos;
typedef off_t gcdb_seek_pos;
typedef unsigned int (*uint_conv_func)(void*); //uint conversion function pointer
typedef off_t (*offt_conv_func)(void*); //uint conversion function pointer
typedef int16_t (*int16_conv_func)(void*); //int16 conversion function pointer


//conversion function --> to platform independent uint
extern uint_conv_func gcvt_uint;
extern offt_conv_func gcvt_offt;
extern int16_conv_func gcvt_int16;
/*
unsigned int uint32_sun(void* x86int);
unsigned int uint32_x86(void* x86int);
//for file offsets: off_t runtime conversions:
off_t offt_sun(void* offt);
off_t offt_x86(void* offt);
int16_t int16_sun(void* i16);
int16_t int16_x86(void* i16);
*/

void gcvt_endian_setup();

class GCDBuffer {
 public:
  char *x;
  unsigned int p;
  unsigned int n;
  int fd;
  opfunc op;
//methods:
  GCDBuffer():x(NULL),p(0),n(0),fd(0),op(NULL) {
    }
  GCDBuffer(opfunc aop,int afd,char *buf,unsigned int len) {
    //check endianness
    gcvt_endian_setup();
    init(aop, afd, buf, len);
    }
  void init(opfunc aop,int afd,char *buf,unsigned int len) {
     x=buf;
     fd=afd;
     op=aop;
     p=0;
     n=len;
     }
  int  flush();
  int  write_all(char* buf, unsigned int pt);
  int  put(char* buf,unsigned int len);
  int  putalign(char* buf,unsigned int len);
  int  putflush(char* buf,unsigned int len);
  int  puts(char *buf);
  int  putsalign(char *buf);
  int  putsflush(char *buf);
  int  oneRead(char* buf, unsigned int len);
  int  getthis(char* buf,unsigned int len);
  int  get(char* buf,unsigned int len);
  int  bget(char* buf,unsigned int len);
  int  feed();
  char *peek();
  void seek(unsigned int len);
  int copy(GCDBuffer* bin);
};


//=====================================================
//-------------     cdb utils       -------------------
//=====================================================
#ifndef _WIN32
 extern int errno;
#endif

extern int error_intr;
extern int error_nomem;
extern int error_proto;

//additional data to be appended to the cdb file:
#define CDBMSK_OPT_MULTI    0x00000001
#define CDBMSK_OPT_C        0x00000002
#define CDBMSK_OPT_CADD     0x00000004
#define CDBMSK_OPT_COMPRESS 0x00000008
#define CDBMSK_OPT_GSEQ     0x00000010
//creates a compressed version of the database
//uses plenty of unions for ensuring compatibility with
// the old 'CIDX' info structure

//trying to prevent [64bit] machines to align this to 64bit -- sizeof() gets it wrong!
#pragma pack(4)
// eek, gcc 2.95.3 alpha-decosf version does not
// recognize this pragma directive

//32 bit limits for index file size
struct cdbInfo {
    uint32 num_keys;
    union {
     uint32 num_records;
     char oldtag[4]; // 'CIDX' for old tag style
     };
    // data file size -- used to be uint32, now it could be 64bit
    union {
     int64_t dbsize;
     uint32 oldnum[2]; //num_keys, num_records
     };
    union {
     uint32 idxflags;
     uint32 old_dbsize;
     };
    union {
     int dbnamelen;
     int old_idxflags;
     };
      // -- the actual db name precedes this fixed-size record
    union {
     char tag[4]; //'CDBX' for new files with LFS
     uint32 old_dbnamelen;
     };
   };

// for passing around index data:
struct CIdxData32 {
   uint32 fpos;
   uint32 reclen;
   };
/*
struct CIdxSeqData32 { //4+4+2+1 = 11 bytes
   uint32 fpos;
   uint32 reclen;
   uint16_t linelen; //line length for FASTA-formatted seq
   byte elen; //length of end-of-line delimiter: 1 (unix/mac) or 2 (Windows)
   };
*/
struct CIdxData {
   off_t fpos; //64bit value on Linux
   uint32 reclen;
 };
/*
struct CIdxSeqData { //8+4+2+1 = 15 bytes
   off_t fpos; //64bit value on Linux
   uint32 reclen;
   uint16_t linelen; //line length for FASTA-formatted seq
   byte elen; //length of end-of-line delimiter: 1 (unix/mac) or 2 (Windows)
  };
*/
#pragma pack()

extern int cdbInfoSIZE;
extern int IdxDataSIZE;
extern int IdxDataSIZE32;
/*
extern int IdxSeqDataSIZE;
extern int IdxSeqDataSIZE32;
*/

void uint32_pack(char *,uint32);
void uint32_pack_big(char *,uint32);
void uint32_unpack(char *,uint32 *);
void uint32_unpack_big(char *,uint32 *);

//=====================================================
//-------------     cdb index       -------------------
//=====================================================

#define CDB_HPLIST 1000

struct cdb_hp { uint32 h; uint32 p; } ;

struct cdb_hplist {
  struct cdb_hp hp[CDB_HPLIST];
  struct cdb_hplist *next;
  int num;
  };


//the index file should always be smaller than 4GB !

class GCdbWrite {
   GCDBuffer* cdbuf;
   char bspace[8192];
   char fname[1024];
   char final[2048];
   uint32 count[256];
   uint32 start[256];
   struct cdb_hplist *head;
   struct cdb_hp *split; /* includes space for hash */
   struct cdb_hp *hash;
   uint32 numentries;
   uint32 pos; //file position
   int posplus(uint32 len);
   int fd; //file descriptor
  public:
  //methods:
   GCdbWrite(int afd); //was: init
   GCdbWrite(char* fname);
   ~GCdbWrite();
   int addbegin(unsigned int keylen,unsigned int datalen);
   int addend(unsigned int keylen,unsigned int datalen,uint32 h);
   int addrec(const char *key,unsigned int keylen,char *data,unsigned int datalen);
   int add(const char *key, char *data, unsigned int datalen);
   int getNumEntries() { return numentries; }
   int finish();
   int close();
   int getfd() { return fd; }
   char* getfile() { return fname; }
};


//=====================================================
//-------------        cdb          -------------------
//=====================================================

#define CDB_HASHSTART 5381

uint32 cdb_hashadd(uint32,unsigned char);
uint32 cdb_hash(const char *,unsigned int);

#define MCDB_SLOT_BITS 8                  /* 2^8 = 256 */
#define MCDB_SLOTS (1u<<MCDB_SLOT_BITS)   /* must be power-of-2 */
#define MCDB_SLOT_MASK (MCDB_SLOTS-1)     /* bitmask */
#define MCDB_HEADER_SZ (MCDB_SLOTS<<4)    /* MCDB_SLOTS * 16  (256*16=4096) */
#define MCDB_MMAP_SZ (1u<<19)             /* 512KB; must be >  MCDB_HEADER_SZ */
#define MCDB_BLOCK_SZ (1u<<22)            /*   4MB; must be >= MCDB_MMAP_SZ */

class GCdbRead {
  //struct mcdb_mmap *map;
  char *map;         // ptr, mmap pointer
    uintptr_t size; // mmap size, initialized if map is nonzero
    uint32_t b;   // hash table stride bits: (data < 4GB) ? 3 : 4
    uint32_t n;   // num records in mcdb

  uint32 loop; // number of hash slots searched under this key
  uint32 hslots; // initialized if loop is nonzero

  uintptr_t kpos; // initialized if loop is nonzero
  uintptr_t hpos; // initialized if loop is nonzero
  uintptr_t dpos; // initialized if cdb_findnext() returns 1

  uint32 dlen; // initialized if cdb_findnext() returns 1
  uint32 klen; // initialized if cdb_findnext() returns 1

  uint32 khash; // initialized if loop is nonzero

  char fname[1024];
  //char *map; // 0 if no map is available
  int fd;
 public:
//methods:
  GCdbRead(int fd); //was cdb_init
  GCdbRead(char* afname); //was cdb_init
  ~GCdbRead(); //was cdb_free
  int read(char *,unsigned int,uint32);
  int match(const char *key, unsigned int len, uint32 pos);
  void findstart() { loop =0; }
  int findnext(const char *key,unsigned int len);
  int find(const char *key);
  int datapos() { return dpos; }
  int datalen() { return dlen; }
  int getfd() { return fd; }
  char* getfile() { return fname; }
};

class GReadBuf {
 protected:
  FILE* f;
  uchar* buf;
  int buflen;
  int bufused; //
  int bufpos;
  off_t fpos;
  bool eof;
  bool eob;

  int refill(bool repos=false) {
   //refill the buffer-----------
   if (repos && bufpos==0) return 0; //no need to repos
   if (eof) return 0;
   int fr=0;
   if (repos && bufpos<bufused) {
      int kept=bufused-bufpos;
      memmove((void*)buf, (void*)(buf+bufpos),kept);
      fr=(int)fread((void *)(buf+kept), 1, buflen-kept, f);
      if (fr<buflen-kept) eof=true;
      buf[kept+fr]='\0';
      bufused=kept+fr;
      }
     else {
      fr=(int)fread((void *)buf, 1, buflen, f);
      if (fr<buflen) eof=true;
      buf[fr]='\0'; //only for text record parsers
      bufused=fr;
      }
   if (feof(f)) eof=true;
   if (ferror(f)) {
     GMessage("GReadBuf::refill - error at fread!\n");
     eof=true;
     }
   bufpos=0;
   fpos+=fr; //bytes read from file so far
   return fr;
   }
 public:
  GReadBuf(FILE* fin, int bsize=4096) {
    f=fin;
    buflen=bsize;
    GMALLOC(buf,buflen+1);
    bufpos=0; //current pointer for get function
    bufused=0;
    fpos=0;
    eof=false;
    eob=false;
    refill();
    }
  ~GReadBuf() { GFREE(buf); }

  //reads len chars from stream into the outbuf
  //updates bufpos
  //->returns the number of bytes read
  int get(uchar *outbuf, int len) {
    if (eob) return 0;
    int rd=0; //bytes read
    while (!eob && rd<len) {
      int to_read=GMIN((bufused-bufpos),(len-rd));
      memcpy((void*)(outbuf+rd),(void*)(buf+bufpos), to_read);
      bufpos+=to_read;
      rd+=to_read;
      if (bufpos>=bufused) {
        if (eof) eob=true;
           else refill();
        }
      }//while
    return rd;
    }

  uchar* getStr(uchar *outbuf, int len) {
    int rd=get(outbuf,len);
    if (rd==0) return NULL;
      else {
       outbuf[rd]='\0';
       return outbuf;
       }
    }

  // getc equivalent
  int getch() {
    if (eob) return -1;
    int ch=(int)(uchar)buf[bufpos];
    bufpos++;
    if (bufpos>=bufused) {
        if (eof) eob=true;
           else refill();
        }
    return ch;
    }

  //---
  bool isEof() { return eob; }
  bool ended() { return eob; }
  off_t getPos() {
  //returns the virtual file position
  // = the actual file offset of the byte at bufpos
    return fpos-(bufused-bufpos);
    }
  //skip into the stream the specified number of bytes
  int skip(int skiplen) {
   if (eob) return 0;
   int r=0; //the actual number of bytes skipped
   while (skiplen && !eob) {
     int dif=GMIN(bufused-bufpos,skiplen);
     skiplen-=dif;
     bufpos+=dif;
     r+=dif;
     if (bufpos>=bufused) {
       if (eof) { eob=true; return r; }
       refill();
       }
     }
    return r;
   }
  //look ahead without updating the read pointer (bufpos)
  //Cannot peek more than buflen!
  int peek(uchar* outbuf, int len) {
    if (eob) return -1;
    //if (eob || len>buflen) return -1;
    if (len>bufused-bufpos) refill(true);
    int mlen=GMIN((bufused-bufpos),len);
    memcpy((void*)outbuf, (void*)(buf+bufpos), mlen);
    return mlen;
    }
  char peekChar() {
    if (eob) return -1;
    //if (eob || len>buflen) return -1;
    if (1>bufused-bufpos) refill(true);
    return *(buf+bufpos);
    }
  uchar* peekStr(uchar* outbuf, int len) {
   int rd=peek(outbuf,len);
   if (rd>0) { outbuf[rd]='\0'; return outbuf; }
        else return NULL;
   }
  //looks ahead to check if what follows matches
  int peekCmp(char* cmpstr, int cmplen=-1) {
    if (cmplen==0) return 0;
    if (eob) //GError("GReadBuf::peekcmp error: eob!\n");
         return -2;
    if (cmplen<0) cmplen=strlen(cmpstr);
    if (cmplen>bufused-bufpos) {
       refill(true);
       if (cmplen>bufused-bufpos) return -2;
       }
    //use memcmp
    return memcmp((void*)(buf+bufpos), cmpstr, cmplen);
    }

};

//circular line buffer, with read-ahead (peeking) capability
class GReadBufLine {
  protected:
    struct BufLine {
        off_t fpos;
        int len;
        char* chars;
        };
    int bufcap; //total number of lines in the buf array
    int bufidx; // the "current line" index in buf array
    bool isEOF;
    int lno;
    FILE* file;
    off_t filepos; //current file/stream offset for the first char of buf[bufidx]
    BufLine* buf; //array of bufferred lines
    char* readline(int idx);//read line from file into the buffer
    int fillbuf();
    bool isEOB;
  public:
    const char* line(); //gets current line and advances the "current line" pointer
                     //use putLine() to revert/undo this advancement
    off_t fpos(); //gets current line's byte offset in the file
                        // does NOT advance the "current line" pointer
    int   len();  //gets current line's length
                       // does NOT advance the "current line" pointer
    bool isEof() { return isEOB; }
    bool eof() { return isEOB; }
    off_t getfpos() { return fpos(); }
    const char* getline() { return line();  }
    const char* getLine() { return line();  }
    int getLen() { return len();  }
    int linenumber() { return lno; }
    int lineno() { return lno; }
    int getLineNo()  { return lno; }
    void putLine();
    GReadBufLine(FILE* stream, int bcap=20) {
      if (bcap<2) bcap=2; //at least 1 prev line is needed for putLine()
      bufcap=bcap;
      bufidx=-1;
      isEOB=false;
      isEOF=false;
      lno=0;
      GMALLOC(buf, bufcap * sizeof(BufLine));
      for (int i=0;i<bufcap;i++) {
          buf[i].chars=NULL;
          buf[i].fpos=-1;
          buf[i].len=0;
          }
      file=stream;
      fillbuf();
      }
    ~GReadBufLine() {
      for (int i=0;i<bufcap;i++) {
          GFREE(buf[i].chars);
          }
      GFREE(buf);
      }
};

#endif
