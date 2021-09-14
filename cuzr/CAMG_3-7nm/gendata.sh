#!/bin/bash
##################################################################
##GENDATA.SH							##
##								##
##Dependencies	:						##
##Influences	:						##
##################################################################
## ver.	: 2019--, Syamal Praneeth Chilakalapudi, KIT, INT	##
##Author Email    :syamalpraneeth@gmail.com			##
##################################################################

for item in {'3nm','7nm'}
do
  cd $item
  pwd
  ./exec_post.py
  cd ../
done

[[ -d "collated" ]] || mkdir collated
./collate.py

