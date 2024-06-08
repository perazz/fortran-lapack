#!/bin/bash
#
#
# preprocess linear algebra sources with fypp
#
declare -a linalg_sources=("solve" "inverse" "least_squares" "determinant" "eye" "svd" "eigs" "qr" "cholesky")

fypp_path="../fypp"
declare -a operations=("src" "test")
declare -a oper_prefix=("stdlib" "test")

# Preprocess all sources
for str in "${linalg_sources[@]}"; do
   for i in "${!operations[@]}"; do 
      operation=${operations[i]}
      pref=${oper_prefix[i]}
      source="$fypp_path/$operation/${pref}_linalg_$str.fypp"
      dest="../$operation/${pref}_linalg_$str.f90"
      fypp -I ../include "$source" > "$dest"

      # prettify
      fprettify -i 4 -l 132 -w 2 --disable-indent --strip-comments --c-relations --enable-replacements --enable-decl --whitespace-comma 0 $dest

   done
done

