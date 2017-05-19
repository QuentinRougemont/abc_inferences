#!/bin/bash

threshold=95
folder=11.isol.selection."$threshold"
arch=../00.architecture.simul #path to the full architecture downloaded from github
cp -r "$arch"/* "$folder" 
