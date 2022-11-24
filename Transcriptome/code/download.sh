#!/bin/bash
while read -r link
do
	echo $link
	wget $link
done
