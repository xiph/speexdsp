#!/bin/sh
gcc -Wall -I../include speexclient.c alsa_device.c -o speexclient -lspeex -lspeexvoip -lasound -lm
