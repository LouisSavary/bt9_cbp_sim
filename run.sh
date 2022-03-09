#!/bin/bash

declare -a predlist=("predictor" "ASpredictor" )

TARGET="eval"
path=$(pwd)
ORIGIN="/home/lsavary/tages/ltage"
ORIGIN=$path
EXECDIR=$ORIGIN+"/cbp2016.final/trois"
cd "/home/lsavary/tages/ltage"

for pred in ${predlist[@]}
do 
	TARGETDIR="compare_preds_"+$pred
	$(cp $pred.h sim/predictor.h)
	$(cp $pred.cc sim/predictor.cc)
	
	cd sim
	echo $(pwd)
#	tmp=$(make 2>&1) #mute
	make
#	cp predictor $EXECDIR/sim
	cd $EXECDIR/scripts
	tmp=""
	tmp=$(./runall.pl -s $($ORIGIN/sim/predictor) -w $TARGET -f 8 -d $TARGETDIR)
	./getdata.pl -w $TARGET -d $TARGETDIR > "res_"+$pred
	cd $ORIGIN
done
cd $path
