#!/bin/bash

g++ cutchecker.cpp -lOsiCbc -lCbc `pkg-config --cflags cbc` `pkg-config --libs cbc` -o cutChecker
