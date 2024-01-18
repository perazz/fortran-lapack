#!/bin/sh

# run modularization script
python3 modularize_blas.py 

# prettify
fprettify -i 4 -l 132 -w 2 --disable-indent --strip-comments --c-relations --enable-replacements --enable-decl --recursive ../src/ 
