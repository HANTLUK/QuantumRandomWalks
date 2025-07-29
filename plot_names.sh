#!/bin/bash

for file in Figures/*.pdf; do
    mv $file ${file%%.*}_pdf.pdf
    export trimmedone="${file%%.*}"
    export filename="${trimmedone##*/}"
    export subst="s/${filename}/${filename}_pdf/g"
    sed -i $subst Figures/${filename}.tex
done
