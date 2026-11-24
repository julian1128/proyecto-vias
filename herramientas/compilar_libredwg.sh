#!/bin/sh
# Compila LibreDWG 0.13.3 (lector de DWG de código abierto) en /home/claude/dwg. Tarda ~6-8 min con un solo núcleo.
# Resultado: /home/claude/dwg/libredwg-0.13.3/programs/dwg2dxf
set -e
mkdir -p /home/claude/dwg && cd /home/claude/dwg
[ -f libredwg.tar.xz ] || curl -sL --max-time 180 -o libredwg.tar.xz https://github.com/LibreDWG/libredwg/releases/download/0.13.3/libredwg-0.13.3.tar.xz
[ -d libredwg-0.13.3 ] || tar xf libredwg.tar.xz
cd libredwg-0.13.3
[ -f Makefile ] || ./configure --disable-bindings --disable-shared --enable-static > /home/claude/dwg/configure.log 2>&1
make CFLAGS="-O0 -g0" > /home/claude/dwg/make.log 2>&1      # -O0: decode.c compila en minutos en lugar de >20 min
echo FIN >> /home/claude/dwg/make.log
