#!/bin/bash

g++ cutChecker.cpp -lOsiCbc -lCbc `pkg-config --cflags cbc` `pkg-config --libs cbc` -o cutChecker
