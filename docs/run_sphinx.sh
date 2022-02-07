#!/usr/bin/bash

cd ~/code/depert/diatom/docs

#sphinx-apidoc -f -o ./source ../diatomic

make clean

make html
