#!/bin/sh

inFile=$1;
Cstart=$(( $(grep -ni "probability intervals" ${inFile} | awk -F ':' '{print $1}') + 1 ));
Cend=$(( $(grep -ni "estimated allele" ${inFile} | awk -F ':' '{print $1}') - 3 ));

echo "No ID   Q     Qlow  Qhigh";
head -n ${Cend} ${inFile} | tail -n +${Cstart} | perl -pe \
    's/^ //; s/^ /0/; s/[0-9]\.[0-9]{3} +\((.*?),(.*?)\).*$/$1 $2/;' | \
    perl -pe 's/ +\(.*?\) +://;s/ +/ /g;';
