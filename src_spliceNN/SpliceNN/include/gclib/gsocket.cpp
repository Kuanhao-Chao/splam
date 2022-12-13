#include "gsocket.h"
#include <errno.h>             // For errno

#ifdef _WIN32
static bool initialized = false;
#endif

void GSocketErr(GStr message, bool inclSysMsg) {
  if (inclSysMsg) {
    message.append(": ");
    message.append(strerror(errno));
  }
  GError("%s\n",message.chars());
}

// Function to fill in address structure given an address and port
static void fillAddr(const GStr &address, unsigned short port, 
                     sockaddr_in &addr) {
  memset(&addr, 0, sizeof(addr));  // Zero out address structure
  addr.sin_family = AF_INET;       // Internet address

  hostent *host;  // Resolve name
  if ((host = gethostbyname(address.chars())) == NULL) {
    // strerror() will not work for gethostbyname() and hstrerror() 
    // is supposedly obsolete
    GSocketErr("Failed to resolve name (gethostbyname())");
  }
  addr.sin_addr.s_addr = *((unsigned long *) host->h_addr_list[0]);

  addr.sin_port = htons(port);     // Assign port in network byte order
}

// GSocket Code

GSocket::GSocket(int type, int protocol) {
  #ifdef _WIN32
    if (!initialized) {
      WORD wVersionRequested;
      WSADATA wsaData;

      wVersionRequested = MAKEWORD(2, 0);              // Request WinSock v2.0
      if (WSAStartup(wVersionRequested, &wsaData) != 0) {  // Load WinSock DLL
        GSocketErr("Unable to load WinSock DLL");
      }
      initialized = true;
    }
  #endif

  // Make a new socket
  if ((sockDesc = socket(PF_INET, type, protocol)) < 0) {
    GSocketErr("GSocket creation failed (socket())", true);
  }
}

GSocket::~GSocket() {
  #ifdef _WIN32
    ::closesocket(sockDesc);
  #else
    ::close(sockDesc);
  #endif
  sockDesc = -1;
}

GStr GSocket::getLocalAddress() {
  sockaddr_in addr;
  unsigned int addr_len = sizeof(addr);

  if (getsockname(sockDesc, (sockaddr *) &addr, (socklen_t *) &addr_len) < 0) {
    GSocketErr("Fetch of local address failed (getsockname())", true);
  }
  return inet_ntoa(addr.sin_addr);
}

unsigned short GSocket::getLocalPort() {
  sockaddr_in addr;
  unsigned int addr_len = sizeof(addr);

  if (getsockname(sockDesc, (sockaddr *) &addr, (socklen_t *) &addr_len) < 0) {
    GSocketErr("Fetch of local port failed (getsockname())", true);
  }
  return ntohs(addr.sin_port);
}

void GSocket::setLocalPort(unsigned short localPort) {
  // Bind the socket to its port
  sockaddr_in localAddr;
  memset(&localAddr, 0, sizeof(localAddr));
  localAddr.sin_family = AF_INET;
  localAddr.sin_addr.s_addr = htonl(INADDR_ANY);
  localAddr.sin_port = htons(localPort);

  if (bind(sockDesc, (sockaddr *) &localAddr, sizeof(sockaddr_in)) < 0) {
    GSocketErr("Set of local port failed (bind())", true);
  }
}

void GSocket::setLocalAddressAndPort(const GStr &localAddress,
    unsigned short localPort) {
  // Get the address of the requested host
  sockaddr_in localAddr;
  fillAddr(localAddress, localPort, localAddr);

  if (bind(sockDesc, (sockaddr *) &localAddr, sizeof(sockaddr_in)) < 0) {
    GSocketErr("Set of local address and port failed (bind())", true);
  }
}

void GSocket::cleanUp() {
  #ifdef _WIN32
    if (WSACleanup() != 0) {
      GSocketErr("WSACleanup() failed");
    }
  #endif
}

unsigned short GSocket::resolveService(const GStr &service,
                                      const GStr &protocol) {
  struct servent *serv;        /* Structure containing service information */

  if ((serv = getservbyname(service.chars(), protocol.chars())) == NULL)
    return atoi(service.chars());  /* Service is port number */
  else 
    return ntohs(serv->s_port);    /* Found port (network byte order) by name */
}

// GCommSocket Code
void GCommSocket::setTimeout(int microsecs) {
 #ifdef _WIN32
   DWORD timeout = microsecs;
   setsockopt(sockDesc, SOL_SOCKET, SO_RCVTIMEO, (const char*)&timeout, sizeof(timeout));
 #else
   struct timeval tv;
   if (microsecs>1000) {
     tv.tv_sec=microsecs / 1000;
     tv.tv_usec=microsecs % 1000;
   } else {
     tv.tv_sec=0;
     tv.tv_usec=microsecs; 
   }
   setsockopt(sockDesc, SOL_SOCKET, SO_RCVTIMEO, (const char*)&tv,sizeof(struct timeval));
 #endif
}

void GCommSocket::connect(const GStr &foreignAddress,
    unsigned short foreignPort) {
  // Get the address of the requested host
  sockaddr_in destAddr;
  fillAddr(foreignAddress, foreignPort, destAddr);

  // Try to connect to the given port
  if (::connect(sockDesc, (sockaddr *) &destAddr, sizeof(destAddr)) < 0) {
    GSocketErr("Connect failed (connect())", true);
  }
}

void GCommSocket::send(const void *buffer, int bufferLen)  {
  if (::send(sockDesc, (raw_type *) buffer, bufferLen, 0) < 0) {
    GSocketErr("Send failed (send())", true);
  }
}

int GCommSocket::recv(void *buffer, int bufferLen) {
  int rtn;
  if ((rtn = ::recv(sockDesc, (raw_type *) buffer, bufferLen, 0)) < 0) {
    GSocketErr("Received failed (recv())", true);
  }
  return rtn;
}

GStr GCommSocket::recvline() {
  GStr r;
  char buf[1024];
  char* p=NULL;
  while (p==NULL) {
    int rtn = ::recv(sockDesc, (raw_type *) buf, 1024, 0);
    if (rtn<0) GSocketErr("Received failed (recv())", true);
    if (rtn==0) return r;
    p=(char*)memchr((void*)buf, '\n', rtn);
    if (p) {
      r.appendmem(buf, p-buf);
      return r;
    }
    r.appendmem(buf, rtn);
  }
 return r;
}


GStr GCommSocket::getForeignAddress() {
  sockaddr_in addr;
  unsigned int addr_len = sizeof(addr);
  if (getpeername(sockDesc, (sockaddr *) &addr,(socklen_t *) &addr_len) < 0) {
    //GSocketErr("Fetch of foreign address failed (getpeername())", true);
    return "";
  }
  return inet_ntoa(addr.sin_addr);
}




unsigned short GCommSocket::getForeignPort() {
  sockaddr_in addr;
  unsigned int addr_len = sizeof(addr);
  if (getpeername(sockDesc, (sockaddr *) &addr, (socklen_t *) &addr_len) < 0) {
    return 0;
  }
  return ntohs(addr.sin_port);
}

// GTCPServerSocket Code

GTCPSocket *GTCPServerSocket::accept() {
  int newConnSD;
  if ((newConnSD = ::accept(sockDesc, NULL, 0)) < 0) {
    GSocketErr("Accept failed (accept())", true);
  }

  return new GTCPSocket(newConnSD);
}

void GTCPServerSocket::setListen(int queueLen) {
  if (listen(sockDesc, queueLen) < 0)
    GSocketErr("Set listening socket failed (listen())", true);
}

// GUDPSocket Code

void GUDPSocket::setBroadcast() {
  // If this fails, we'll hear about it when we try to send.  This will allow 
  // system that cannot broadcast to continue if they don't plan to broadcast
  int broadcastPermission = 1;
  setsockopt(sockDesc, SOL_SOCKET, SO_BROADCAST, 
             (raw_type *) &broadcastPermission, sizeof(broadcastPermission));
}

void GUDPSocket::disconnect() {
  sockaddr_in nullAddr;
  memset(&nullAddr, 0, sizeof(nullAddr));
  nullAddr.sin_family = AF_UNSPEC;

  // Try to disconnect
  if (::connect(sockDesc, (sockaddr *) &nullAddr, sizeof(nullAddr)) < 0) {
   #ifdef _WIN32
    if (errno != WSAEAFNOSUPPORT) {
   #else
    if (errno != EAFNOSUPPORT) {
   #endif
      GSocketErr("Disconnect failed (connect())", true);
    }
  }
}

void GUDPSocket::sendTo(const void *buffer, int bufferLen, 
    const GStr &foreignAddress, unsigned short foreignPort) {
  sockaddr_in destAddr;
  fillAddr(foreignAddress, foreignPort, destAddr);
  // Write out the whole buffer as a single message.
  if (sendto(sockDesc, (raw_type *) buffer, bufferLen, 0,
             (sockaddr *) &destAddr, sizeof(destAddr)) != bufferLen) {
    GSocketErr("Send failed (sendto())", true);
  }
}

int GUDPSocket::recvFrom(void *buffer, int bufferLen, GStr &sourceAddress,
    unsigned short &sourcePort) {
  sockaddr_in clntAddr;
  socklen_t addrLen = sizeof(clntAddr);
  int rtn;
  if ((rtn = recvfrom(sockDesc, (raw_type *) buffer, bufferLen, 0, 
                      (sockaddr *) &clntAddr, (socklen_t *) &addrLen)) < 0) {
    GSocketErr("Receive failed (recvfrom())", true);
  }
  sourceAddress = inet_ntoa(clntAddr.sin_addr);
  sourcePort = ntohs(clntAddr.sin_port);

  return rtn;
}

void GUDPSocket::setMulticastTTL(unsigned char multicastTTL) {
  if (setsockopt(sockDesc, IPPROTO_IP, IP_MULTICAST_TTL, 
                 (raw_type *) &multicastTTL, sizeof(multicastTTL)) < 0) {
    GSocketErr("Multicast TTL set failed (setsockopt())", true);
  }
}

void GUDPSocket::joinGroup(const GStr &multicastGroup) {
  struct ip_mreq multicastRequest;

  multicastRequest.imr_multiaddr.s_addr = inet_addr(multicastGroup.chars());
  multicastRequest.imr_interface.s_addr = htonl(INADDR_ANY);
  if (setsockopt(sockDesc, IPPROTO_IP, IP_ADD_MEMBERSHIP, 
                 (raw_type *) &multicastRequest, 
                 sizeof(multicastRequest)) < 0) {
    GSocketErr("Multicast group join failed (setsockopt())", true);
  }
}

void GUDPSocket::leaveGroup(const GStr &multicastGroup) {
  struct ip_mreq multicastRequest;

  multicastRequest.imr_multiaddr.s_addr = inet_addr(multicastGroup.chars());
  multicastRequest.imr_interface.s_addr = htonl(INADDR_ANY);
  if (setsockopt(sockDesc, IPPROTO_IP, IP_DROP_MEMBERSHIP, 
                 (raw_type *) &multicastRequest, 
                 sizeof(multicastRequest)) < 0) {
    GSocketErr("Multicast group leave failed (setsockopt())", true);
  }
}
