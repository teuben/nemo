/*  TCPpipe - creates a TCP connection and (monodirectionally) sends any kind of data through it
 *  Copyright (C) 2001  Grischa Weberstaedt
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *        2001   Original version                  Grischa Weberstaedt
 *  2-jul-2001   Adapted for NEMO                  Peter Teuben
 * 10-feb-2004   added report= and changed host= to better reflect its use PJT
 *  
 */ 

#include <nemo.h>

#include <unistd.h>
#include <errno.h>
#include <sys/socket.h>
#include <sys/times.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

string defv[] = {    /* NEMO's command line keyword, and their defaults */
  "host_sender=\n          hostname, only needed in sender mode",
  "port=2801\n             TCP/IP port to communicate with",
  "bufsize=1024\n          buffersize to use",
  "report=f\n              Report I/O speed info (to stderr)",
  "VERSION=0.4\n           10-feb-04 PJT",
  NULL,
};

string usage = "create TCP connection to send data from/to pipe";

#define LISTEN  1
#define CONNECT 2
#define SEND    4
#define RECEIVE 8

int do_listen(int port, char *rhost);
int do_connect(int port, char *rhost);
void do_transfer(int from, int to, int buffer_size, bool Qreport, int port);
int get_ip(char *string, struct in_addr *ip);

int nemo_main() {
    int sock, mode = 0, 
      direction = 0,
      port = getiparam("port"),
      buffer_size = getiparam("bufsize");
    string rhost = NULL;
    bool Qreport = getbparam("report");

    if (hasvalue("host_sender")) {
      mode = CONNECT;
      direction = RECEIVE;
      rhost = getparam("host_sender");
    } else {
      mode = LISTEN;
      direction = SEND;
    }

    switch (mode) {
	case LISTEN:	sock = do_listen(port, rhost);
			break;
	case CONNECT:	sock = do_connect(port, rhost);
			break;
        default:	error("bad mode");
                        break;
    }
    switch (direction) {
        case SEND:	do_transfer(STDIN_FILENO, sock, buffer_size, FALSE, port);
			break;
        case RECEIVE:	do_transfer(sock, STDOUT_FILENO, buffer_size, Qreport, port);
			break;
        default:	error("bad direction");
			break;
    }
    return 0;
}

int do_listen(int port, char *rhost) {
    struct sockaddr_in address;
    int address_size = sizeof(address);
    int listen_socket;
    int data_socket;

    address.sin_family=AF_INET;
    address.sin_addr.s_addr=INADDR_ANY;
    address.sin_port=htons(port);

    if ((listen_socket = socket(PF_INET, SOCK_STREAM, 0)) == -1) error("socket-1 port=%d",port);
    if (bind(listen_socket, &address, address_size) == -1) error("bind-1 port=%d",port);
    if (listen(listen_socket, 1) == -1) error("listen-1 port=%d",port);
    if ((data_socket=accept(listen_socket, &address, &address_size)) == -1) error("accept-1 port=%d",port);
    if (close(listen_socket) == -1) error("close-1 port=%d",port);
    return data_socket;
}

int do_connect(int port, char *rhost) {
    struct sockaddr_in address;
    int data_socket;

    address.sin_family=AF_INET;
    if (!get_ip(rhost, &address.sin_addr)) {
	warning("`%s': ", rhost);
	herror(NULL);
	exit(EXIT_FAILURE);
    }
    address.sin_port=htons(port);

    if ((data_socket=socket(PF_INET, SOCK_STREAM, 0)) == -1) 
      error("socket-2 port=%d",port);
    if (connect(data_socket, &address, sizeof(address)) == -1) 
      error("connect-2 port=%d",port);
    return data_socket;
}

void do_transfer(int from, int to, int buffer_size, bool Qreport, int port) {
    char *buffer;
    int used_length;
    clock_t ct0, ct1;
    int n = 0;
    struct tms buf;
    long clk_tck;

    if ((buffer=malloc(buffer_size)) == NULL) error("malloc-3");
    if (Qreport) ct0 = times(&buf);
    while ( (used_length=read(from, buffer, buffer_size)) ) {
	int written_length=0;
	if (used_length == -1) error("read-3");
	while (used_length - written_length > 0) {
	    int write_ret;
	    if ((write_ret=write(to, buffer+written_length, used_length-written_length)) == -1) 
	      error("write-3");
	    written_length += write_ret;
	}
	n += used_length;
    }
    if (Qreport) {
      double r1, r2;
      ct1 = times(&buf);
      clk_tck = sysconf(_SC_CLK_TCK);
      r1 = (double)n / (double)(ct1-ct0);
      r2 = r1/1024/1024*clk_tck;
      dprintf(0,"tcppipe(port %d) : %d bytes/ %ld ticks => %g bytes/clock_tick = %g MB/sec\n",
	      port,n,ct1-ct0,r1,r2);
    }
}

int get_ip(char *string, struct in_addr *ip) {
    struct hostent *host;
    if (inet_aton(string, ip)) return 1;
    else {
	if ((host=gethostbyname2(string, AF_INET)) == NULL) return 0;
	*ip = *(struct in_addr *)host->h_addr;
	return 1;
    }
}
