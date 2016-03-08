#!/bin/sh

fs=$1; #file string (the variable bit, at least)

numpops=$2; # number of populations (default = 2)
repeats=$3; # number of repeats (default = 1)
outdir="wtccc";
#program="~/install/structure_linux_v21/bin/structure"; # version 2.1
program="/home/gringer/install/structure2.2.3/bin/structure.exe"; # version 2.2.3


if [ -e $fs ]
    then echo "error: no file name given"
    echo "usage: $0 <file> [without structure_, .txt]";
    exit
fi

if [ ! -f runs/structure_${fs}.txt ]
		then echo "File not found: runs/structure_${fs}.txt";
		exit
fi

if [ -e $numpops ]
    then numpops=2;
fi

if [ -e $repeats ]
    then repeats=1;
fi

mkdir -p results/${outdir};

ni=$(( $(cat runs/structure_${fs}.txt | wc -l) - 1 )); # number of individuals
nm=$(head -n1 runs/structure_${fs}.txt | wc -w); # number of markers

## runs structure for random samplings

for markers in $(echo -e "$nm")
do
    markstr=$(printf '%04d' ${markers})
    for iter in $(seq 1 ${repeats})
    do
	[ ! -f results/${outdir}/out_${fs}_${markstr}_${iter}_K${numpops}_f ] &&
	${program} \
	    -L ${markers} -N ${ni} -K ${numpops} -i \
	    runs/structure_${fs}.txt \
	    -o results/${outdir}/out_${fs}_${markstr}_${iter}_K${numpops}
    done
done

