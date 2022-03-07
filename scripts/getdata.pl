#!/usr/bin/perl -w
#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************

require ( "./bench_list.pl");

$stat      = "MISPRED_PER_1K_INST";
$statscnd  = "MISPREPRED_PER_1K_INST";
$statthrd  = "MISPRED_TRACE";
$wsuite    = "all";
$max_dirs  = 1024;
$amean     = 1;  # always print amean
$prepred   = 6;

#####################################
######### USAGE OPTIONS      ########
#####################################

$debug     = 0;
$noxxxx    = 0;


#####################################
######### USAGE OPTIONS      ########
#####################################

sub usage(){

$USAGE = "Usage:  '$0 <-options> -d <dir1> ... <dirN>'";
 
print(STDERR "$USAGE\n");
print(STDERR "\t-h                     : help -- print this menu. \n");
print(STDERR "\t-d <statdirs>          : directory for stats (more than one is ok, if -d is last knob) \n");
print(STDERR "\t-s <statname>          : name of the stat (wild card ok) \n");
print(STDERR "\t-w <workload/suite>    : name of the workload suite from bench_list\n");
print(STDERR "\t-noxxxx                : Print 0 for no data instead of xxxx. \n");
print(STDERR "\n");

exit(1);
}

######################################
########## PARSE COMMAND LINE ########
######################################
 
while (@ARGV) {
    $option = shift;
 
    if ($option eq "-h"){
	usage();
    }elsif ($option eq "-d") {
        while(@ARGV){
	    $mydir = shift;
	    $mydir .= "/";
	    push (@dirs, $mydir);
	}
    }elsif ($option eq "-s") {
        $stat = shift;
    }elsif ($option eq "-w") {
        $wsuite = shift;
    }elsif ($option eq "-debug") {
        $debug = 1;
    }elsif ($option eq "-noxxxx") {
        $noxxxx = 1;
    }else{
	usage();
        die "Incorrect option ... Quitting\n";
    }
}
             

##########################################################
# get the suite names, num_w, and num_p from bench_list

die "No benchmark set '$wsuite' defined in bench_list.pl\n"
        unless $SUITES{$wsuite};
    
my (@w) = split(/\s+/, $SUITES{$wsuite});
$num_w = scalar @w;

init_stats();
get_stats();

print_header();    
print_stats();    
print_amean()     if($amean);

# $stat = "MISPRED_TRACE";
# init_stats();
# get_stats();

# print_header();    
# print_stats();    
# print_amean()     if($amean);


########################## INIT STATS ####################

sub init_stats{
  for($ii=0; $ii< $max_dirs; $ii++){
    for($jj=0; $jj< $num_w; $jj++){
      for($kk=0; $kk<= $prepred+1; $kk++){
        $data[$ii][$jj][$kk]=0;

      }
    }
  }

}


########################## GET STATS ####################

sub get_stats{
    
    for($dirnum=0; $dirnum<@dirs; $dirnum++){
	$dir = $dirs[$dirnum];
	
	for($ii=0; $ii< $num_w; $ii++){
	    $wname   = $w[$ii];
	    $fname   = $dir . $wname . ".res";
	    $readable = 1;
	    
	    open(IN, $fname) or $readable=0;
	    
	    print "cannot open $fname for read\n" unless $readable;
	    
	    $found=0;
	    
	    while( ($readable) &&  ($_ = <IN>) ){
            @words = split/\s+/, $_;
            $num_words = scalar(@words);
            
            for($jj=0; $jj<$num_words-1; $jj++){
                if($words[$jj] && ($words[$jj] eq $stat) ){
                    $pos = $jj+2;
                    $val = $words[$pos];
                    $data[$dirnum][$ii] += $val;
                    printf "stat match for $stat found for $_" if($debug);
                    $found++;
                }
                if($words[$jj] && ($words[$jj] eq $statscnd) ){
                    $pos = $jj+3;
                    $pos2= $jj+1;
                    $val = $words[$pos];
                    $val2= $words[$pos2];
                    $data[$dirnum][$ii][$val2] += $val;
                    printf "stat match for $stat found for $_" if($debug);
                    #$found++;
                }
                if($words[$jj] && ($words[$jj] eq $statthrd) ){
                    $pos = $jj+2;
                    $val = $words[$pos];
                    $data[$dirnum][$ii][$prepred] += $val;
                    printf "stat match for $stat found for $_" if($debug);
                    $found++;
                }
            }
	    }
	    close(IN);
	}
    }
    
}


########################## INIT HEADER ####################

sub print_header{

  $header = "";

  for($dirnum=0; $dirnum< @dirs; $dirnum++){
    $dir = $dirs[$dirnum];
    chop $dir;
    $len = length($dir);
    $mystring = substr($dir, ($len-14), 14);
    $header .= $mystring."\t";
  }

  printf("\n%-30s\t", "ResultDirs ==>" );
  print "$header\n";
}


################# PRINT WORKLOAD STATS ######################

sub print_stats{


    for($ii=0; $ii< $num_w; $ii++){
        $wname   = $w[$ii];
        $wstring = substr($wname, 0, 20);
        printf("\n%-20s\t", $wstring);

        for($dirnum=0; $dirnum < @dirs ; $dirnum++){
            $val     = $data[$dirnum][$ii];
            print_val();
        }
    }

    print "\n";    
}

################# PRINT VAL ######################

sub print_val{

    if($val){ 
	printf("%12.3f \t", $val);  
    }
    else { 
	printf("xxxxxxxxxxx \t") unless($noxxxx); 
	printf("0           \t") if($noxxxx); 
    }
}


################# PRINT AMEAN STATS ######################

sub print_amean{


    for($dirnum=0; $dirnum < @dirs; $dirnum++){
	    for($statnum=0; $statnum <= $prepred+1; $statnum++ ){
            $dir_sums[$dirnum][$statnum]=0;
            for($ii=0; $ii< $num_w; $ii++){
                $dir_sums[$dirnum][$ii][$statnum] += $data[$dirnum][$ii][$statnum];
            }
        }
    }
    
    printf("\n%-20s\t", "AMEAN");
    
    for($dirnum=0; $dirnum < @dirs; $dirnum++){
        for($statnum=0; $statnum <= $prepred+1; $statnum++){
            $val = $dir_sums[$dirnum][$statnum]/$num_w;
            print_val();
            
            if ($statnum < $prepred){
                printf("\n%-19s %1d\t", "AMEAN", $statnum+1);
            }
            if ($statnum == $prepred){
                printf("\n%-20s\t", "AMEAN TRACE");
            }
        }
    }
    
    
    # printf("\n%-20s\t", "AMEAN");
    
    # for($dirnum=0; $dirnum < @dirs; $dirnum++){
    #     $val = $dir_sums[$dirnum]/$num_w;
    #     print_val();
    # }
    
    print "\n\n";
}

###########################################################
