#!/bin/bash
set -e

sudo apt-get install libboost-system-dev libboost-iostreams-dev libboost-program-options-dev

if [ -d ~/clustgun ]; then
	cd ~/clustgun
	git pull
	make clean
else
	cd
	git clone git://github.com/MG-RAST/clustgun.git
fi

cd ~/clustgun

make

sudo ln -s /home/ubuntu/clustgun/clustgun /usr/bin/clustgun

if [ -d ~/apps/bin/ ]; then
	cd ~/apps/bin/
	ln -sf /home/ubuntu/clustgun/clustgun
fi