#!/bin/sh

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <A> <B>"
    exit 1
fi

for i in $(seq $1 $2)
do
	mkdir -p result$i
	output/main.out result$i &
done
