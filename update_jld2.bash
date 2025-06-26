#!/bin/bash

echo "copying dynamics files for" $1 "from" $2
rm ~/research/ContactImplicitMPC.jl/src/dynamics/$1/dynamics/*
cp -r ~/.julia/packages/ContactImplicitMPC/$2/src/dynamics/$1/dynamics/* ~/research/ContactImplicitMPC.jl/src/dynamics/$1/dynamics/

echo "copying simulation files for" $1 "from" $2
rm -rf ~/research/ContactImplicitMPC.jl/src/simulation/$1/*
cp -r ~/.julia/packages/ContactImplicitMPC/$2/src/simulation/$1/* ~/research/ContactImplicitMPC.jl/src/simulation/$1/

echo "git add"
git add -A
git commit -am "added jld2 files"
git push

echo "done"
