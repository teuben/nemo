/*
 * bwfreq 66|60|50|48
 * 
 * Jim Gasbarro
 * Xerox PARC
 * 415-494-4803
 * 
 * This program changes the vertical refresh rate of the B&W monitor on a
 * SPARCStation1 using the ECL frame buffer board.  This is useful for
 * reducing flicker when videotaping or filming a screen image. 60 Hz mode
 * should be used in the US, 50 Hz in Europe.  48 Hz should be used for film.
 * 66 Hz is the normal rate.
 * 
 * To change the the refresh rate you must either be logged in as root, or
 * change the protection of /dev/sbus<n> to 666.  n = [1, 2, or 3] is the
 * number of the SBus slot where the frame buffer is installed.  Slot 1 is
 * next to the power switch.
 * 
 * WARNING: Slowing the vertical rate increases the current in the vertical
 * deflection yoke and drive circuits.  This may cause overheating and damage
 * to the monitor if continued for extended periods of time.  This feature is
 * not supported by Sun Microsystems Inc.  Use at your own risk.
 * 
 * Note: To use 48 or 50 Hz mode you must adjust the V.Freq pot (R306) inside
 * the monitor.  Take the outside plastic cover off and remove the top metal
 * EMI shield.  R306 is in the upper left quadrant of the rear board.  Ignore
 * the red glpt and go for it!  The extra trash at the bottom of the screen
 * in 48 and 50 Hz modes is alas unavoidable due to hardware limitations.
 */

#include <stdio.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>

#define VREG_OFFSET 0x400000

#define R18_66 0x5a03a824
#define R18_60 0x5a04128e
#define R18_50 0x5a04e0ff
#define R18_48 0x5a0515ff
#define R1c_66 0x0007ff01
#define R1c_60 0x3037ff01
#define R1c_50 0x3037ff01
#define R1c_48 0x3037ff01

main(argc, argv)
    char          **argv;
{
    char           *addr, *base_addr, *vreg_base;
    int             size = 0xc00000;
    int             result, speed;
    int             fd, slot;
    char            buf[80], buf1[80];
    FILE           *fp;

    if (argc < 2) {
        fprintf(stderr, "Usage: bwfreq 66|60|50|48\n");
        fprintf(stderr,
                "Changes vertical frame rate to allow videotaping.\n");
        fprintf(stderr,
                "Works only with SPARCStation1 low res B&W frame buffer.\n");
        exit(1);
    }
    speed = atoi(argv[1]);
    /* look in /usr/adm/messages to find out what sbus slot bwtwo0 is in */
    fp = popen("/etc/dmesg | grep bwtwo0 | tail -1 | awk \'{print $5}\'", "r");
    if (fp) {
        slot = -10;
        fscanf(fp, "%d", &slot);
        if (slot != -10)
            printf("using slot %d\n");
    } else {
        fprintf(stderr, "couldn't figure out which SBus slot to use\n");
        exit(1);
    }
    sprintf(buf, "/dev/sbus%d", slot);
    fd = open(buf, O_RDWR);
    if (fd < 0) {
        sprintf(buf1, "open: %s", buf);
        perror(buf1);
        exit(1);
    };

    base_addr = (char *) mmap(0, size, (PROT_READ | PROT_WRITE),
        (MAP_SHARED), fd, 0);
    if ((int) base_addr == -1) {
        perror("mmap failed");
        exit(1);
    };

    /* base addresses for video control registers and frame buffer */
    vreg_base = base_addr + VREG_OFFSET;

    addr = vreg_base + 0x18;
    switch (speed) {
    case 66:
        *(int *) addr = R18_66;
        break;
    case 60:
        *(int *) addr = R18_60;
        break;
    case 50:
        *(int *) addr = R18_50;
        break;
    case 48:
        *(int *) addr = R18_48;
        break;
    default:
        fprintf(stderr,
                "%d is not a valid frequency\n", speed);
        exit(1);
    }

    addr = vreg_base + 0x1c;
    switch (speed) {
    case 66:
        *(int *) addr = R1c_66;
        break;
    case 60:
        *(int *) addr = R1c_60;
        break;
    case 50:
        *(int *) addr = R1c_50;
        break;
    case 48:
        *(int *) addr = R1c_48;
        break;
    }
}
