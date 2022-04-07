#!/usr/bin/perl -w
#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************

require ( "./bench_list.pl");

$statsize  = 8;
my @stat   = ("MISPRED_PER_1K_INST", "MEAN_TRACE_PRECISION", "MEAN_TRACE_USAGE", "MEAN_TRACE_INSTRUCTION", "TRACE_NUMBER", "CONSTRUCTION_FAILS", "UNUSED_TRACE_PER", "OVERLAP_COEF");
# $statscnd  = "MEAN_TRACE_PRECISION";
# $statthrd  = "MEAN_TRACE_USAGE";
$wsuite    = "all";
$max_dirs  = 1024;
$amean     = 1;  # always print amean
$prepred   = 6;

#####################################
######### USAGE OPTIONS      ########
#####################################

$debug     = 0;
$noxxxx    = 1;


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
    # }elsif ($option eq "-s") {
    #     $stat = shift;
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
      for($kk=0; $kk<= 2*$prepred+1; $kk++){
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
                for($jjj=0; $jjj<$statsize; $jjj++){
                    if($words[$jj] && ($words[$jj] eq $stat[$jjj]) ){
                        $pos = $jj+2;
                        $val = $words[$pos];
                        $data[$dirnum][$ii][$jjj] += $val;
                        printf "stat match for $stat[$jjj] found for $_" if($debug);
                        $found++;
                    }
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
            for ($iii = 0; $iii < $statsize; $iii++){
                $val     = $data[$dirnum][$ii][$iii];
                print_val();
            }
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
	    for($statnum=0; $statnum < $statsize; $statnum++ ){
            $dir_sums[$dirnum][$statnum]=0;
            for($ii=0; $ii< $num_w; $ii++){
                $dir_sums[$dirnum][$statnum] += $data[$dirnum][$ii][$statnum] if !($data[$dirnum][$ii][$statnum] eq 'nan');
            }
        }
        
        printf("\n%-20s\t", "AMEAN MISPRED");
        $val = $dir_sums[$dirnum][0]/$num_w;
        print_val();
        
        printf("\n%-20s\t", "AMEAN PRECISION");
        $val = $dir_sums[$dirnum][1]/$num_w;
        print_val();
        
        printf("\n%-20s\t", "AMEAN USAGE");
        $val = $dir_sums[$dirnum][2]/$num_w;
        print_val();

        printf("\n%-20s\t", "AMEAN INSTRUCTION");
        $val = $dir_sums[$dirnum][3]/$num_w;
        print_val();
        printf("\n%-20s\t", "AMEAN TRACES NB");
        $val = $dir_sums[$dirnum][4]/$num_w;
        print_val();
        printf("\n%-20s\t", "AMEAN FAILS");
        $val = $dir_sums[$dirnum][5]/$num_w;
        print_val();
        printf("\n%-20s\t", "UNUSED TRACE PERCENT");
        $val = $dir_sums[$dirnum][6]/$num_w;
        print_val();
        printf("\n%-20s\t", "OVERLAP");
        $val = $dir_sums[$dirnum][7]/$num_w;
        print_val();

    }

    # for($dirnum=0; $dirnum < @dirs; $dirnum++){
    #         $val = $dir_sums[$dirnum][$statnum]/$num_w;
    #         print_val();
            
    #         if ($statnum < $preped){
    #             printf("\n%-18s %d\t", "AMEAN PRECISION", $statnum+1);
    #         }
    #         if ($statnum == $prepred){
    #             printf("\n%-20s\t", "AMEAN TRACE");
    #         }
        
    # }
    
    
    # printf("\n%-20s\t", "AMEAN");
    
    # for($dirnum=0; $dirnum < @dirs; $dirnum++){
    #     $val = $dir_sums[$dirnum]/$num_w;
    #     print_val();
    # }
    
    print "\n\n";
}

###########################################################
