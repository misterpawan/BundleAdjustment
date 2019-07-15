#!/bin/bash

clear

rm -rf build

mkdir build && cd build && cmake .. && make

clear

cd Apps

./bundle_large ~/Dataset/problem-356-226730-pre.txt
