#!/bin/bash

clear

rm -rf build

mkdir build && cd build && cmake .. && make

clear

cd Apps

./bundle_large ~/Dataset/problem-49-7776-pre.txt
