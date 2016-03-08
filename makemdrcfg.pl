#!/usr/bin/perl

# makemdrcfg.pl -- generates a configuration file for the MDR/pMDR
# programs

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use strict;
use warnings;
use English '-no_match_vars';
use File::Spec;

our ($dummy1, $dummy2, $pName) = File::Spec->splitpath($PROGRAM_NAME);

sub usage {
    print("usage: ./makemdrcfg.pl <genotype file> [options]");
}

my %options = ();

# write out default values
$options{"CALCTHRESHOLD"} = "SET";
$options{"CROSSVALINTERVAL"} = "1";
$options{"MODELSTOKEEP"} = "1";
$options{"NOTRECOGNIZEDRESPONSE"} = "1";
$options{"TIECELLVALUE"} = "1";
$options{"BESTPARTITION"} = "OFF";
$options{"MISSING"} = "-1";
$options{"POWERMODEL"} = "OFF";
$options{"PTEST-N"} = "0";

# parse command line arguments (so while(<>) works)
    while (@ARGV) {
        my $test=shift @ARGV;
        if(-f $test){
            if($options{"FILENAME"}){
                die("Error: no more than one file may be specified\n".usage());
            }
            $options{"FILENAME"} = $test;
        }
        else{
            $options{uc($test)} = shift @ARGV;
        }
}


print "// Configuration file for MDR (pMDR)\n";
print "// Generated using $pName [written by David Eccles (gringer)]\n\n";
print "// Comments are formatted in the following way:\n";
print "// OPTION        [default] Description\n";
print "// [Taken from pMDR readme File]\n\n";

# There's lots of repeated code here... a better way to do this would
# be to have (variable, default, description) triples (i.e.
# $default{"variable"}, $description{"variable"}) and iterate through
# the list.


print "// FILENAME           [none]  Specifies the location and name of the\n";
print "//                            data input file.\n";
if($options{"FILENAME"}){
    print("FILENAME ".$options{"FILENAME"}."\n");
}
print "\n";
print "// MAXLOCIVALUE       [none]  Maximum value for an input variable.  The\n";
print "//                            minimum is always 0.\n";
if($options{"MAXLOCIVALUE"}){
    print("MAXLOCIVALUE ".$options{"MAXLOCIVALUE"}."\n");
}
print "\n";
print "// CALCTHRESHOLD      [SET]   Threshold for assigning high/low risk\n";
print "//                            for a genotype during evaluation of a\n";
print "//                            model.  It is either calculated by the\n";
print "//                            overall ratio in the dataset (SET) or it\n";
print "//                            can be calculated for each model (MODEL).\n";
print "//                            It is faster to calculate on the SET, but\n";
print "//                            may be necessary to calculate per model\n";
print "//                            if there are missing genotypes in the dataset.\n";
if($options{"CALCTHRESHOLD"}){
    print("CALCTHRESHOLD ".$options{"CALCTHRESHOLD"}."\n");
}
print "\n";
print "// COMBOSTART         [none]  Minimum number of loci making up\n";
print "//                            a combination to be tested.  A '2' means\n";
print "//                            2-loci models will be tested initially\n";
print "//                            and larger models will be tested exhaustively\n";
print "//                            until COMBOEND is reached.\n";
if($options{"COMBOSTART"}){
    print("COMBOSTART ".$options{"COMBOSTART"}."\n");
}
print "\n";
print "// COMBOEND           [none]  Maximum number of loci making up\n";
print "//                            a combination to be tested. A '3' means\n";
print "//                            no models involving more than 3 loci\n";
print "//                            will be tested.\n";
if($options{"COMBOEND"}){
    print("COMBOEND ".$options{"COMBOEND"}."\n");
}
print "\n";
print "// MODELBUILDINTERVAL [none]  Specifies how many models will be evaluated\n";
print "//                            by each process before requesting more\n";
print "//                            models to evaluate.\n";
if($options{"MODELBUILDINTERVAL"}){
    print("MODELBUILDINTERVAL ".$options{"MODELBUILDINTERVAL"}."\n");
}
print "\n";
print "// CROSSVALINTERVAL   [1]     Number of cross-validation intervals to\n";
print "//                            perform.  Setting it to 1 means no\n";
print "//                            cross-validation will be done.\n";
if($options{"CROSSVALINTERVAL"}){
    print("CROSSVALINTERVAL ".$options{"CROSSVALINTERVAL"}."\n");
}
print "\n";
print "// MODELSTOKEEP       [1]     Number of models that will be kept per\n";
print "//                            loci combination and per cross-validation.\n";
print "//                            Set to 'ALL' to keep every model. WARNING:\n";
print "//                            If testing many combinations this can\n";
print "//                            greatly slow the analysis and create\n";
print "//                            a very large output file.\n";
if($options{"MODELSTOKEEP"}){
    print("MODELSTOKEEP ".$options{"MODELSTOKEEP"}."\n");
}
print "\n";
print "// RANDOMSHUFFLE     [none]   Values are 'ON' or 'OFF'.  When 'ON', the\n";
print "//                            individuals are shuffled before splitting\n";
print "//                            them to do cross-validation.\n";
if($options{"RANDOMSHUFFLE"}){
    print("RANDOMSHUFFLE ".$options{"RANDOMSHUFFLE"}."\n");
}
print "\n";
print "// RANDOMSEED        [none]   Seeds the random number generator for the\n";
print "//                            random shuffle of individuals.\n";
if($options{"RANDOMSEED"}){
    print("RANDOMSEED ".$options{"RANDOMSEED"}."\n");
}
print "\n";
print "// NOTRECOGNIZEDRESPONSE [1]  Sets how a genotype is evaluated that is\n";
print "//                            encountered during testing but not classified\n";
print "//                            during training.\n";
print "//                            Values are: -1 (Unknown), 0 (Unaffected),\n";
print "//                            1 (Affected)\n";
if($options{"NOTRECOGNIZEDRESPONSE"}){
    print("NOTRECOGNIZEDRESPONSE ".$options{"NOTRECOGNIZEDRESPONSE"}."\n");
}
print "\n";
print "// TIECELLVALUE         [1]   Sets a genotype risk when equal numbers of\n";
print "//                            affected and unaffected individuals during\n";
print "//                            training.\n";
print "//                            Values are: -1 (Unknown), 0 (Unaffected),\n";
print "//                            1 (Affected)\n";
if($options{"TIECELLVALUE"}){
    print("TIECELLVALUE ".$options{"TIECELLVALUE"}."\n");
}
print "\n";
print "// BESTPARTITION      [OFF]   Shows data partition for best model.\n";
print "//                            Values are: OFF, ON\n";
if($options{"BESTPARTITION"}){
    print("BESTPARTITION ".$options{"BESTPARTITION"}."\n");
}
print "\n";
print "// MISSING             [-1]   Sets value for missing data in dataset.\n";
if($options{"MISSING"}){
    print("MISSING ".$options{"MISSING"}."\n");
}
print "\n";
print "// POWERMODEL         [OFF]   Specifies whether to determine and display\n";
print "//                            the best model across all combination\n";
print "//                            sizes.\n";
print "//                            Values are:  OFF, ON\n";
if($options{"POWERMODEL"}){
    print("POWERMODEL ".$options{"POWERMODEL"}."\n");
}
print "\n";
print "// PTEST-N              [0]   Specifies number of permutation tests to perform\n";
print "//                            to assign significance to the models.  This will\n";
print "//                            increase the running time as each permutation test\n";
print "//                            takes as long as the initial run of the data.\n";
if($options{"PTEST-N"}){
    print("PTEST-N ".$options{"PTEST-N"}."\n");
}
print "\n";
