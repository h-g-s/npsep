#!/bin/bash

g++ -O2 cutchecker.cpp `pkg-config --cflags cbc` `pkg-config --libs cbc` -o cutChecker
