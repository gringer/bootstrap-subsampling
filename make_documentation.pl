#!/usr/bin/perl

# make_documentation.pl -- creates LaTeX documentation for script files.

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./make_documentation.pl <file name>\n");
  print("\n");
}

my @fileNames = ();

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        push(@fileNames, $argument);
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        else{
            printf(STDERR "Error: Unknown command-line option, '$argument'\n");
            usage();
            exit(1);
        }
    }
}

if(!@fileNames){
    printf(STDERR "Error: No files specified on command line\n");
    usage();
    exit(2);
}

while(@fileNames){
    my $currentFile = shift(@fileNames);
    $currentFile =~ s#^.*(\|/)##; # remove dir names from argument
    print(STDERR "Now checking '$currentFile'\n");
    open(INFILE, "< $currentFile");
    my %foundItems = ();
    my $codeLines = 0;
    my $oneLinerContinue = 0; # false
    my $usageSection = 0; # false
    my $usageContinue = 0; # false
    my $nameCredit = 0; # false
    if($currentFile =~ /\.pl$/){
        # Perl specific tests
        while(<INFILE>){
            if(/David Eccles \(gringer\), 20[0-9][0-9]/){
                $nameCredit = 1; # true
            }
            if(/sub usage {/){
                $foundItems{"usage"} = 1;
                $usageSection = 1; # true
            }
            if($usageSection){
                if($usageContinue){
                    if(/^\s+"(.*?)\\n"/){
                        $foundItems{"fileUsage"} .= $1;
                    }
                }
                if(/^\s+print.*?"usage: (.*?)\\n"/){
                    $foundItems{"fileUsage"} = $1;
                }
                if(/^\s+print.*?"usage: (.*?)"\./){
                    $foundItems{"fileUsage"} = $1;
                    $usageContinue = 1; # true
                }
                if(/^\s+print.*?"-(.*?)\s*:\s*(.*?)\\n"/){
                    my $optionName = $1;
                    my $optionDesc = $2;
                    if(!exists($foundItems{"commandLineOptions"})){
                        $foundItems{"commandLineOptions"} = "";
                    }
                    $foundItems{"commandLineOptions"} .=
                        "\\item[-$optionName] $optionDesc\n";
                }
                if(/^\s+}/){
                    $usageSection = 0; #false
                }
            }
            if($oneLinerContinue && (/^# (.*)$/)){
                $foundItems{"oneLiner"}.= " ".$1;
            } else {
                $oneLinerContinue = 0; # false
            }
            if(/^# $currentFile -- (.*)$/){
                $foundItems{"oneLiner"} = $1;
                $oneLinerContinue = 1; # true
            }
            if(/^use warnings;/){
                $foundItems{"warningUsed"} = 1;
            }
            if(/^use strict;/){
                $foundItems{"strictUsed"} = 1;
            }
            if(!(/^\s*$/) && !(/^#/)){
                $codeLines++;
            }
        }
        close(INFILE);
        # Perl specific warnings
        if(!exists($foundItems{"warningUsed"})){
            print(STDERR "Error: no 'use warnings;' line found in file\n");
            usage();
            exit(3);
        }
        if(!exists($foundItems{"strictUsed"})){
            print(STDERR "Error: no 'use strict;' line found in file\n");
            usage();
            exit(3);
        }
    }
    if($currentFile =~ /\.r$/){
        # R specific tests, uses 'cat' instead of print
        while(<INFILE>){
            if(/David Eccles \(gringer\), 20[0-9][0-9]/){
                $nameCredit = 1; # true
            }
            if(/usage <- function\(\){/){
                $foundItems{"usage"} = 1;
                $usageSection = 1; # true
            }
            if($usageSection){
                if($usageContinue){
                    if(/^\s+"(.*?)\\n"/){
                        $foundItems{"fileUsage"} .= $1;
                    }
                }
                if(/^\s+cat.*?"usage: (.*?)\\n"/){
                    $foundItems{"fileUsage"} = $1;
                }
                if(/^\s+cat.*?"usage: (.*?)",/){
                    $foundItems{"fileUsage"} = $1;
                    $usageContinue = 1; # true
                }
                if(/^\s+cat.*?"-(.*?)\s*:\s*(.*?)\\n"/){
                    my $optionName = $1;
                    my $optionDesc = $2;
                    if(!exists($foundItems{"commandLineOptions"})){
                        $foundItems{"commandLineOptions"} = "";
                    }
                    $foundItems{"commandLineOptions"} .=
                        "\\item[-$optionName] $optionDesc\n";
                }
                if(/^\s+}/){
                    $usageSection = 0; #false
                }
            }
            if($oneLinerContinue && (/^## (.*)$/)){
                $foundItems{"oneLiner"}.= " ".$1;
            } else {
                $oneLinerContinue = 0; # false
            }
            if(/^## $currentFile -- (.*)$/){
                $foundItems{"oneLiner"} = $1;
                $oneLinerContinue = 1; # true
            }
            if(!(/^\s*$/) && !(/^#/)){
                $codeLines++;
            }
        }
        close(INFILE);
        # R specific warnings
    }
    # generic warnings and final printing
    if(!exists($foundItems{"oneLiner"})){
        print(STDERR "Error: no one Liner description found in file\n");
        usage();
        exit(3);
    }
    if(!exists($foundItems{"usage"})){
        print(STDERR "Error: no usage section found in file\n");
        usage();
        exit(3);
    }
    if(!$nameCredit){
        print(STDERR "David Eccles (gringer), <date> is not acknowledged as author in file\n");
        usage();
        exit(3);
    }
    my $fileSize = qx{wc -c $currentFile};
    $fileSize =~ s/ .*$//;
    if($fileSize > 1024){
        $fileSize = sprintf("%d",$fileSize / 1024)." KiB";
    } else {
        $fileSize = sprintf("%d",$fileSize)." B";
    }
    printf("\\section{%s}\n", $currentFile);
    printf("\\label{sec:%s}\n\n", $currentFile);
    printf("\\subsection{Overview}\n");
    printf("\\label{sec:%s-overview}\n\n", $currentFile);
    printf("\\begin{description}\n");
    printf("\\item[Usage] %s\n", $foundItems{"fileUsage"});
    printf("\\item[Purpose] %s\n", $foundItems{"oneLiner"});
    printf("\\item[Lines of Code] %d\n", $codeLines);
    printf("\\item[File size] %s\n", $fileSize);
    printf("\\end{description}\n\n");
    printf("\\subsection{Command Line Options}\n");
    printf("\\label{sec:%s-command-line}\n\n", $currentFile);
    if(exists($foundItems{"commandLineOptions"})){
        printf("\\begin{description}\n");
        printf($foundItems{"commandLineOptions"});
        printf("\\end{description}\n\n");
    } else {
        printf("This script has no additional command line options.\n\n");
    }
    printf("\\emph{[Additional Comments]}\n\n");
    print(STDERR "Finished checking '$currentFile'\n");
}
