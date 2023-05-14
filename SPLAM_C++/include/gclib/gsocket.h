#ifndef GSOCKET_DEFINED
#define GSOCKET_DEFINED
#include "GBase.h"
#include "GStr.h"

#ifdef _WIN32
  #include <winsock.h>         // For socket(), connect(), send(), and recv()
  typedef int socklen_t;
  typedef char raw_type;       // Type used for raw data on this platform
#else
  #include <sys/types.h>       // For data types
  #include <sys/socket.h>      // For socket(), connect(), send(), and recv()
  #include <netdb.h>           // For gethostbyname()
  #include <arpa/inet.h>       // For inet_addr()
  #include <unistd.h>          // For close()
  #include <netinet/in.h>      // For sockaddr_in
  typedef void raw_type;       // Type used for raw data on this platform
#endif


/**
 *   Signals a problem with the execution of a socket call.
 *   @param message explanatory message
 *   @param incSysMsg true if system message (from strerror(errno))
 *   should be postfixed to the user provided message
 */
void GSocketErr(GStr message, bool inclSysMsg = false);

/**
 *   Base class representing basic communication endpoint
 */
class GSocket {
public:
  // Close and deallocate this socket
  ~GSocket();
  /**
  *   Get the local address
  *   @return local address of socket
  */
  GStr getLocalAddress();

  /**
   *   Get the local port
   *   @return local port of socket
  */
  unsigned short getLocalPort();

  /**
   *   Set the local port to the specified port and the local address
   *   to any interface
   *   @param localPort local port
   */
  void setLocalPort(unsigned short localPort);

  /**
   *   Set the local port to the specified port and the local address
   *   to the specified address.  If you omit the port, a random port 
   *   will be selected.
   *   @param localAddress local address
   *   @param localPort local port
   */
  void setLocalAddressAndPort(const GStr &localAddress, 
    unsigned short localPort = 0);

  /**
   *   If WinSock, unload the WinSock DLLs; otherwise do nothing.  We ignore
   *   this in our sample client code but include it in the library for
   *   completeness.  If you are running on Windows and you are concerned
   *   about DLL resource consumption, call this after you are done with all
   *   Socket instances.  If you execute this on Windows while some instance of
   *   Socket exists, you are toast.  For portability of client code, this is 
   *   an empty function on non-Windows platforms so you can always include it.
   *   @param buffer buffer to receive the data
   *   @param bufferLen maximum number of bytes to read into buffer
   *   @return number of bytes read, 0 for EOF, and -1 for error
   */
  static void cleanUp();

  /**
   *   Resolve the specified service for the specified protocol to the
   *   corresponding port number in host byte order
   *   @param service service to resolve (e.g., "http")
   *   @param protocol protocol of service to resolve.  Default is "tcp".
   */
  static unsigned short resolveService(const GStr &service,
                                       const GStr &protocol = "tcp");

private:
  // Prevent the user from trying to use value semantics on this object
  GSocket(const GSocket &sock);
  void operator=(const GSocket &sock);

protected:
  int sockDesc;              // Socket descriptor
  GSocket(int type, int protocol);
  GSocket(int sockDesc) { this->sockDesc = sockDesc; }
};

/**
 *   Socket which is able to connect, send, and receive
 */
class GCommSocket : public GSocket {
public:
  /**
   *   Establish a socket connection with the given foreign
   *   address and port
   *   @param foreignAddress foreign address (IP address or name)
   *   @param foreignPort foreign port
  */
  void setTimeout(int microsecs);
  void connect(const GStr &foreignAddress, unsigned short foreignPort);

  /**
   *   Write the given buffer to this socket.  Call connect() before
   *   calling send()
   *   @param buffer buffer to be written
   *   @param bufferLen number of bytes from buffer to be written
  */
  void send(const void *buffer, int bufferLen);
  void send(const GStr& str) { send(str.chars(), str.length()); }

  /**
   *   Read into the given buffer up to bufferLen bytes data from this
   *   socket.  Call connect() before calling recv()
   *   @param buffer buffer to receive the data
   *   @param bufferLen maximum number of bytes to read into buffer
   *   @return number of bytes read, 0 for EOF, and -1 for error
   */
  int recv(void *buffer, int bufferLen);
  GStr recvline();
  /**
   *   Get the foreign address.  Call connect() before calling recv()
   *   @return foreign address
   */
  GStr getForeignAddress();

  /**
   *   Get the foreign port.  Call connect() before calling recv()
   *   @return foreign port
   */
  unsigned short getForeignPort();

protected:
  GCommSocket(int type, int protocol) : GSocket(type, protocol) { }
  GCommSocket(int newConnSD) : GSocket(newConnSD) { }
};


//   TCP socket for communication with other TCP sockets
class GTCPSocket : public GCommSocket {
public:
  //   Construct a TCP socket with no connection
  GTCPSocket() : GCommSocket(SOCK_STREAM, IPPROTO_TCP) { }

  /**
   *   Construct a TCP socket with a connection to the given foreign address
   *   and port
   *   @param foreignAddress foreign address (IP address or name)
   *   @param foreignPort foreign port
   */
  GTCPSocket(const GStr &foreignAddress, unsigned short foreignPort) 
           : GCommSocket(SOCK_STREAM, IPPROTO_TCP) {
     connect(foreignAddress, foreignPort);
  }

private:
  // Access for TCPServerSocket::accept() connection creation
  friend class GTCPServerSocket;
  GTCPSocket(int newConnSD) : GCommSocket(newConnSD) { }
};

//  TCP socket class for servers
class GTCPServerSocket : public GSocket {
public:
  /**
   *   Construct a TCP socket for use with a server, accepting connections
   *   on the specified port on any interface
   *   @param localPort local port of server socket, a value of zero will
   *                   give a system-assigned unused port
   *   @param queueLen maximum queue length for outstanding 
   *                   connection requests (default 5)
   */
  GTCPServerSocket(unsigned short localPort, int queueLen = 5)
          : GSocket(SOCK_STREAM, IPPROTO_TCP) {
     setLocalPort(localPort);
     setListen(queueLen);
  }


  /**
   *   Construct a TCP socket for use with a server, accepting connections
   *   on the specified port on the interface specified by the given address
   *   @param localAddress local interface (address) of server socket
   *   @param localPort local port of server socket
   *   @param queueLen maximum queue length for outstanding 
   *                   connection requests (default 5)
   */
  GTCPServerSocket(const GStr &localAddress, unsigned short localPort,
      int queueLen = 5) : GSocket(SOCK_STREAM, IPPROTO_TCP) {
    setLocalAddressAndPort(localAddress, localPort);
    setListen(queueLen);
  }

  //  Blocks until a new connection is established on this socket or error
  //  @return new connection socket
  GTCPSocket *accept();

private:
  void setListen(int queueLen);
};

/**
  *   UDP socket class
  */
class GUDPSocket : public GCommSocket {
public:
  //   Construct a UDP socket
  GUDPSocket() : GCommSocket(SOCK_DGRAM, IPPROTO_UDP) { setBroadcast(); }

  //   Construct a UDP socket with the given local port
  //   @param localPort local port
  GUDPSocket(unsigned short localPort) : 
      GCommSocket(SOCK_DGRAM, IPPROTO_UDP) {
    setLocalPort(localPort);
    setBroadcast();
  }

   //   Construct a UDP socket with the given local port and address
   //   @param localAddress local address
   //   @param localPort local port
  GUDPSocket(const GStr &localAddress, unsigned short localPort) : 
      GCommSocket(SOCK_DGRAM, IPPROTO_UDP) {
    setLocalAddressAndPort(localAddress, localPort);
    setBroadcast();
  }

  //   Unset foreign address and port
  //   @return true if disassociation is successful
  void disconnect();
   /*   Send the given buffer as a UDP datagram to the
    *   specified address/port
    *   @param buffer buffer to be written
    *   @param bufferLen number of bytes to write
    *   @param foreignAddress address (IP address or name) to send to
    *   @param foreignPort port number to send to
    *   @return true if send is successful
    */
  void sendTo(const void *buffer, int bufferLen, const GStr &foreignAddress,
            unsigned short foreignPort);

  /*   Read read up to bufferLen bytes data from this socket.  The given buffer
   *   is where the data will be placed
   *   @param buffer buffer to receive data
   *   @param bufferLen maximum number of bytes to receive
   *   @param sourceAddress address of datagram source
   *   @param sourcePort port of data source
   *   @return number of bytes received and -1 for error
   */
  int recvFrom(void *buffer, int bufferLen, GStr &sourceAddress, 
               unsigned short &sourcePort);

   //   Set the multicast TTL
   //   @param multicastTTL multicast TTL
  void setMulticastTTL(unsigned char multicastTTL);

  //   Join the specified multicast group
  //   @param multicastGroup multicast group address to join
  void joinGroup(const GStr &multicastGroup);

  //   Leave the specified multicast group
  //   @param multicastGroup multicast group address to leave
  void leaveGroup(const GStr &multicastGroup);

private:
  void setBroadcast();
};


#endif /* GSOCKET_DEFINED */
