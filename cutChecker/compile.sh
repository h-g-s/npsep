#!/bin/bash

g++ cutChecker.cpp `pkg-config --cflags cbc` `pkg-config --libs cbc` -o cutChecker
