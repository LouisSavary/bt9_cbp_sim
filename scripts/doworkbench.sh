#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************
#  This shell file gives examples of launching simulations for all 4 categories of workloads


#set this variable to NumCores in your cluster machine for faster sims
num_parallel_jobs=4


###########  HOW TO RUN JOBS?  ################

DUMP_FILE="/tmp/doworkbench_dump_file"
RESULT_DIR="../results/trace_workbench"
wsuite="workbench"
trace_dir="~/traces/"
gen_graph=0
while getopts w:d:f:t: flag
do
    case "${flag}" in
        w) wsuite=${OPTARG};;
        d) dest_dir=${OPTARG};;
        f) DUMP_FILE=${OPTARG};;
        t) trace_dir=${OPTARG};;
        g) gen_graph=1;;
        G) gen_graph=0;;
    esac
done

# The following line will launch sims for all workloads when you run ./doit.sh (comment it if you dont want it to) 

# ./runall.pl -s ../sim/predictor -w all -f  $num_parallel_jobs -d ../results/SEZNEC2014.08KB

./runall2.pl -s ../sim/predictor -w $wsuite -f  $num_parallel_jobs -d $RESULT_DIR 


###########  HOW TO GET STATS?  ################

# This scripts creates stats, after all the earlier jobs finish

./getdata.pl -w $wsuite -d $RESULT_DIR > $DUMP_FILE
if [ $gen_graph ]
then 
    python graph_gen_script.py $DUMP_FILE
else 
    cat $DUMP_FILE
fi

#./getdata.pl -w all -d ../results/SEZNEC2014.08KB

# To compare MPKI numbers against GSHARE for the provided benchmarks , uncomment this line 
# ./getdata.pl -w all -d ../results/MYRESULTS ../results/GSHARE.04KB  ../results/GSHARE.08KB ../results/GSHARE.16KB ../results/GSHARE.32KB



###########  PAST RUNS FOR 4KB, 8KB, and 16KB ###########

# The results are already in "../results" directory

#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.04KB"
#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.08KB"
#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.16KB"

################## GOOD LUCK! ##################
