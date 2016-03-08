#!/bin/sh

from="$1"
to="$2"
file="$3"

usage() {
    echo -e "Usage: $0 FROM TO FILE\n\n"
    echo "FROM: marker to start extraction of genotypes"
    echo "TO  : marker to stop extraction of genotypes"
    echo "FILE: file to extract genotypes from (standard input will be"
    echo "      assumed if this is omitted)"
}


if [[ ( -z "$from" ) || ( -z "$to" ) ]]
    then
    echo -e "error: from/to markers have not been specified\n\n"
    usage
    exit 1
fi

# get line numbers of markers
if [[ -z "$file" ]]
    then
    echo -e "error: file has not been specified\n\n"
    usage
    exit 1
fi

# use zless instead of zcat, because it works with
# files that are not gzipped
zless $file | grep -n -e "$from" -e "$to" | sed 's/^\([0-9]\+\):.*$/\1/' > line_nums.tmp

top=$(head -n 1 line_nums.tmp)
bottom=$(tail -n 1 line_nums.tmp)

if [[ ( -z "$top" ) || ( "$top" == "$bottom" ) ]]
    then
    echo -e "error: input does not contain both specified markers, or"
    echo -e "       markers are the same\n\n"
    usage
    exit 1
fi

if [[ $top > $bottom ]]
    # if first marker occurrs after second in the file, switch round
    then
    tmp=$top
    top=$bottom
    bottom=$tmp
fi

difference="$(echo "$bottom - $top" | bc)"

#echo -e "$bottom and $top are $difference lines away from each other\n\n"

zless $file | head -n 1
#echo "---"
zless $file | head -n $bottom | tail -n $(( $difference + 1 ))
