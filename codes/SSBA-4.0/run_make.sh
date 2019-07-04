#!/bin/bash

clear

rm -rf build

mkdir build && cd build && cmake .. && make

clear

cd Apps

./bundle_large ~/Dataset/problem-885-97473-pre.txt
