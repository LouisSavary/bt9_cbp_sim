#!/bin/bash

declare -a predlist=("JB" "AS" )

TARGET="eval"
path=$(pwd)
ORIGIN="/home/lsavary/tages/ltage"
ORIGIN=$path
EXECDIR=$ORIGIN
#cd "/home/lsavary/tages/ltage"

for pred in ${predlist[@]}
do 
	TARGETDIR="compare_preds_"$pred"pred"
	$("cp "$pred"predictor.h sim/predictor.h")
	$("cp "$pred"predictor.cc sim/predictor.cc")
	$("cp "$pred"main.cc sim/main.cc")
	
	cd sim
	echo $(pwd)
#	tmp=$(make 2>&1) #mute
	make
#	cp predictor $EXECDIR/sim
	cd $EXECDIR/scripts
	tmp=""
	$("./runall.pl -s $ORIGIN/sim/predictor -w $TARGET -f 8 -d $TARGETDIR")
	$("./getdata.pl -w $TARGET -d $TARGETDIR > res_"$pred"pred")
	cd $ORIGIN
done
cd $path
