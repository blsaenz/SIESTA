#!/usr/local/bin/perl

    # Time Commands
    # ------------------------------------------------------------------

    #@timeData = localtime(time);
    #print join(' ', @timeData);

    # You should see something similar to this, although the number could
    # be very different: 20 36 8 27 11 105 2 360 0
    # These elements of the current time are, in order:

    # Seconds past the minute
    # Minutes past the hour
    # Hours past midnight
    # Day of the month
    # Months past the start of the year
    # Number of years since 1900
    # Number of days since the start of the week (Sunday)
    # Number of days since the start of the year
    # Whether or not daylight savings is active

    # File/System Commands
    # ------------------------------------------------------------------

    #chdir "/u/rminer/tmp";
    #mkdir("test", 0755);
    #rename("test","gunk");
    #system("ps | grep rminer");

# check to see if this is restart job or not
$constants_file  = "constants.txt";
open(DB, $constants_file) || die "\nCan't open \"$constants_file\n";
read(DB,$c_data,10000) || die "Can't read file\n";
close(DB);

@c_data_lines = split('\n',$c_data);
@restart = split('\t',$c_data_lines[0]);

if ($restart[0] == 0) {

    # create date string to be used for output directory and filenames
    @timeData = localtime(time);

    $timeData[4] = $timeData[4] + 1;

    if ($timeData[4] == 1) { $mo = "Jan"; }
    elsif ($timeData[4] == 2) { $mo = "Feb"; }
    elsif ($timeData[4] == 3) { $mo = "Mar"; }
    elsif ($timeData[4] == 4) { $mo = "Apr"; }
    elsif ($timeData[4] == 5) { $mo = "May"; }
    elsif ($timeData[4] == 6) { $mo = "Jun"; }
    elsif ($timeData[4] == 7) { $mo = "Jul"; }
    elsif ($timeData[4] == 8) { $mo = "Aug"; }
    elsif ($timeData[4] == 9) { $mo = "Sep"; }
    elsif ($timeData[4] == 10) { $mo = "Oct"; }
    elsif ($timeData[4] == 11) { $mo = "Nov"; }
    elsif ($timeData[4] == 12) { $mo = "Dec"; }
    $yr = $timeData[5]+1900;


    # create output directory
    if ($timeData[3] < 10) {
        $out_dir = "data/output/run_0" . $timeData[3] . $mo . $yr;
    } else {
        $out_dir = "data/output/run_" . $timeData[3] . $mo . $yr;
    }
    if ($timeData[2] < 10) {
        $out_dir = $out_dir . "_0" . $timeData[2];
    } else {
        $out_dir = $out_dir . "_" . $timeData[2];
    }
    if ($timeData[1] < 10) {
        $out_dir = $out_dir . "0$timeData[1]";
    } else {
        $out_dir = $out_dir . $timeData[1];
    }
    mkdir($out_dir,0744);
    mkdir("$out_dir/code",0744);
    print "output directory: $out_dir\n";

    # copy options files, source code, & exececutable to output directory
    system("cp aph.txt $out_dir/code");
    system("cp aozone.txt $out_dir/code");
    system("cp aice.txt $out_dir/code");
    system("cp awater.txt $out_dir/code");
    system("cp awater_tc.txt $out_dir/code");
    system("cp awater_sc.txt $out_dir/code");
    system("cp awater_v.txt $out_dir/code");
    system("cp surfacelight.txt $out_dir/code");
    system("cp constants.txt $out_dir/code");
    system("cp stations.txt $out_dir/code");
    system("cp Makefile* $out_dir/code");
    system("cp comp.pl $out_dir/code");
    system("cp run.pl $out_dir/code");
    system("cp *.F* $out_dir/code");
    system("cp *.f90 $out_dir/code");
    system("cp *.cu $out_dir/code");
    system("cp sia2temp $out_dir/code");
    system("cp data/Boundary/boundary.nc $out_dir/code");

}

# call executable
if ($ARGV[0] ne "preonly") {
    if ($ARGV[0] eq "idb") {
        system("idb sia2temp -dbx");
    } elsif ($ARGV[0] eq "idbc") {
        system("idbc sia2temp -dbx");
    } elsif ($ARGV[0] eq "gdb") {
        system("gdb sia2temp");
    } elsif ($ARGV[0] eq "lldb") {
        system("lldb sia2temp");
    } elsif ($ARGV[0] eq "dbx") {
        system("dbx sia2temp");
    } else {
        system("./sia2temp");
    }
}