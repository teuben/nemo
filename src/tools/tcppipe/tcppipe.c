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
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#define VERSION "0.2"
#define DEFAULT_BUFFER_SIZE 1024

char *argv0;

int do_listen(int port, char *rhost);
int do_connect(int port, char *rhost);
void do_transfer(int from, int to, int buffer_size);
int get_ip(char *string, struct in_addr *ip);

void print_usage() {
    printf("Usage: %s [options] [host]\n", argv0);
    printf("Send/receive any kind of data via a TCP connection

  -l, --listen		wait for a connection
  -c, --connect		connect to host
			Either -l or -c must be used.
  -s, --send		read from stdin and send, default if -l
  -r, --receive		receive and write to stdout, default if -c
  -p, --port=P		listen on / connect to port P, necessary
  -b, --buffer=P	use buffer of size P, default: %d
  -h, --help		print this
  -V, --version		print version

Report bugs to <epibrator@gmx.net>.\n", DEFAULT_BUFFER_SIZE);
    exit(EXIT_SUCCESS);
}

void print_version() {
    printf("TCPpipe %s\n", VERSION);
    printf("Copyright (C) 2001 Grischa Weberstaedt
 TCPpipe comes with NO WARRANTY.
 You may redistribute copies of TCPpipe
 under the terms of the GNU General Public License.\n");
    exit(EXIT_SUCCESS);
}

void param_error(const char *string) {
    if (string==NULL) fprintf(stderr, "%s: There is a problem with the arguments.\n", argv0);
    else if (*string!=0) fprintf(stderr, "%s: %s\n", argv0, string);
    fprintf(stderr, "%s: Try `%s --help' for usage information.\n", argv0, argv0);
    exit(EXIT_FAILURE);
}

// incomplete !
void error(const char *text) {
    perror(argv0);
    exit(EXIT_FAILURE);
}

void oops() {
    fprintf(stderr, "Ooops!\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    struct option longopts[] = {
	{ "listen",	no_argument,		NULL,	'l' },
	{ "connect",	no_argument,		NULL,	'c' },
	{ "send",	no_argument,		NULL,	's' },
	{ "receive",	no_argument,		NULL,	'r' },
	{ "port",	required_argument,	NULL,	'p' },
	{ "buffer",	required_argument,	NULL,	'b' },
	{ "help",	no_argument,		NULL,	'h' },
	{ "version",	no_argument,		NULL,	'V' },
	{ NULL, 0, NULL, 0 }
    };
    int opt;
    enum { LISTEN=1, CONNECT } mode = 0;
    enum { SEND=1, RECEIVE } direction = 0;
    int port = 0;
    int buffer_size = 0;
    char *rhost = NULL;
    int sock;

    argv0=argv[0];

    while ((opt=getopt_long(argc, argv, "-lcsrp:b:hV", longopts, NULL)) != -1) {
	switch (opt) {
	    char *tailptr;
	    case 'l':	if (mode==0) mode=LISTEN;
			else param_error(NULL);
			break;
	    case 'c':	if (mode==0) mode=CONNECT;
			else param_error(NULL);
			break;
	    case 's':	if (direction==0) direction=SEND;
			else param_error(NULL);
			break;
	    case 'r':	if (direction==0) direction=RECEIVE;
			else param_error(NULL);
			break;
	    case 'p':	if (port==0) {
			    errno=0;
			    port = strtol(optarg, &tailptr, 0);
			    if (errno || *tailptr!=0) param_error(NULL);
			}
			else param_error(NULL);
			break;
	    case 'b':	if (buffer_size==0) {
			    errno=0;
			    buffer_size = strtol(optarg, &tailptr, 0);
			    if (errno || *tailptr!=0) param_error(NULL);
			}
			else param_error(NULL);
			break;
	    case 'h':	print_usage();
			break;
	    case 'V':	print_version();
			break;
	    case 1:	if (rhost==NULL) rhost=optarg;
			else param_error(NULL);
			break;
	    case '?':	param_error("");
			break;
	    default:	oops();
			break;
	}
    }
    if (mode==0) param_error(NULL);
    if (direction==0) direction = mode==LISTEN ? SEND : RECEIVE;
    if (port==0) param_error(NULL);
    if (port<0 || port>65535) param_error(NULL);
    if (buffer_size==0) buffer_size = DEFAULT_BUFFER_SIZE;
    if (buffer_size<0) param_error(NULL);
    if (mode==CONNECT && rhost==NULL) param_error(NULL);

    switch (mode) {
	case LISTEN:	sock = do_listen(port, rhost);
			break;
	case CONNECT:	sock = do_connect(port, rhost);
			break;
	default:	oops();
                        break;
    }
    switch (direction) {
	case SEND:	do_transfer(STDIN_FILENO, sock, buffer_size);
			break;
	case RECEIVE:	do_transfer(sock, STDOUT_FILENO, buffer_size);
			break;
	default:	oops();
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

    if ((listen_socket = socket(PF_INET, SOCK_STREAM, 0)) == -1) error("");
    if (bind(listen_socket, &address, address_size) == -1) error("");
    if (listen(listen_socket, 1) == -1) error("");
    if ((data_socket=accept(listen_socket, &address, &address_size)) == -1) error("");
    if (close(listen_socket) == -1) error("");
    return data_socket;
}

int do_connect(int port, char *rhost) {
    struct sockaddr_in address;
    int data_socket;

    address.sin_family=AF_INET;
    if (!get_ip(rhost, &address.sin_addr)) {
	fprintf(stderr, "%s: `%s': ", argv0, rhost);
	herror(NULL);
	exit(EXIT_FAILURE);
    }
    address.sin_port=htons(port);

    if ((data_socket=socket(PF_INET, SOCK_STREAM, 0)) == -1) error("");
    if (connect(data_socket, &address, sizeof(address)) == -1) error("");
    return data_socket;
}

void do_transfer(int from, int to, int buffer_size) {
    char *buffer;
    int used_length;

    if ((buffer=malloc(buffer_size)) == NULL) error("");
    while (used_length=read(from, buffer, buffer_size)) {
	int written_length=0;
	if (used_length == -1) error("");
	while (used_length - written_length > 0) {
	    int write_ret;
	    if ((write_ret=write(to, buffer+written_length, used_length-written_length)) == -1) error("");
	    written_length += write_ret;
	}
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

