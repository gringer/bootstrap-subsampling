#!/usr/bin/perl

use warnings;
use strict;

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Author:  Zev Hochberg
#
# File Description:
# taxidToGeneLocs.pl -- Accepts tax_id as an argument,
# Gets name information from gene via e-utils
#

# Modified 2008-08-23 by David Eccles (gringer)
# Now outputs chromosome and location for genes as a CSV file
# Note: CSV lines will only output if ChrLoc data is present
# Modified 2009-04-22 by David Eccles (gringer)
# - help/usage converted to my own format for automated documentation

use strict;

# ---------------------------------------------------------------------------
#
use LWP::Simple;    # Define library for the 'get' function.

# ---------------------------------------------------------------------------
#
# Global variables

my $verbose = 0;
my $taxId;
my $dopt;
my $beginTime = time();

my $startTag = "<eSummaryResult>";
my $stopTag = "</eSummaryResult>";

# ---------------------------------------------------------------------------
#
# Main program
#

&processCommandLineParameters;

&verifyInput;

&printHeader if $dopt eq 'csv';

&process;

&end;


################################################################################
#
sub processCommandLineParameters {

    $taxId = "9606"; # default to Homo sapiens (human)
    $dopt = "csv"; # default to CSV output

    while(defined ($_= shift @ARGV)) {
        if ($_ eq '-h') {
            usage();
            exit(0);
        } elsif ($_ eq '-v') {
            $verbose = 1;
        } elsif ($_ eq '-t') {
            $taxId = shift @ARGV;
        } elsif ($_ eq '-o') {
            $dopt = shift @ARGV;
        } else {
            # Ignore extra input
        }
    }
    if (($taxId eq "") || ($dopt eq "")){
        usage();
        exit(1);
    }
}

################################################################################
#
sub verifyInput {
    unless (($taxId =~ m/^\d*$/)  and  ($taxId != 0)) {
        print STDERR "taxonomyId (-t) must be numeric and > 0\n";
        usage();
        exit(2);
    }
    unless (($dopt eq 'xml') or  ($dopt eq 'csv')) {
        print STDERR "output (-o) is required and must equal 'xml' or 'csv'\n";
        usage();
        exit(3);
  }
}

################################################################################
#
sub usage {
    print("usage: ./taxidToGeneLocs [option] -t taxonomyId -o xml|csv\n");
    print("Command-line options:\n\n");
    print("-h            : Display this usage help information\n");
    print("-v            : Verbose\n");
    print("-o <xml/csv>  : output options\n");
    print("   xml  - XML\n");
    print("   csv  - comma-delimited\n");
}


################################################################################
#
sub printHeader {
    print STDOUT "geneId,name,description,chromosome,cytoloc,chrstart,chrstop\n";
}


################################################################################
#
sub process {
    my $geneStart = 0;
    my $maxGenes = 250;
    my $totalGenes;
    my $haveTotal = 0;

    my $geneId = "";
    my $geneCount = 0;

    my $qs;
    my $esearch_result;
    my $esummary_result;

    my $first = 1;

    # Main loop: get up to maxGenes Gene ID's for requested taxId from Gene
    $totalGenes = $geneStart+1; # to get started

    GeneLoop:
    for (; $geneStart<$totalGenes; $geneStart+=$maxGenes) {
        if ($verbose) {print STDOUT "Processing genes $geneStart to ", $geneStart+$maxGenes-1, "\n";}

        #defining the query
        #this option looks for GeneIDs with the $taxId value of interest that are 
        #encoded on the mitochondrion. Note use of the field restriction [properties] 
        #and the boolean AND
        #$qs = "db=gene&retstart=$geneStart&retmax=$maxGenes&term=" . $taxId . "[taxid]+AND+source_mitochondrion[properties]";

        #this option looks for GeneIDs with the $taxId value of interest that are 
        #NOT encoded on the mitochondrion and do not have RefSeqs of the type model.
        #Note use of the boolean NOT and field restriction
        $qs = "db=gene&retstart=$geneStart&retmax=$maxGenes&term=" . $taxId . "[taxid]+NOT+source_mitochondrion[properties]+NOT+srcdb_refseq_model[properties]";
        
        $esearch_result = &Eutil ("esearch", $qs);
        if ($esearch_result =~ m/<Count>(\d+)<\/Count>/) {
            if ($haveTotal) {
                if ($totalGenes != $1) {
                    die "esearch reported new total genes: was $totalGenes; now $1\nFor request $qs\n";
                }
            } else {
                $totalGenes = $1; # extract total
                $haveTotal = 1;
            }
        } else {
            die "$esearch_result did not contain expected <Count>,\n for request $qs\n";
        }


        # Build querystring for GENE search
        $qs = "db=gene";
        while ($esearch_result =~ m/<Id>(\d+)<\/Id>/g) {
            $geneCount++;
            $geneId = $1;   # extract a geneId
            $qs .= "&id=$geneId";
        }

        # Get Gene summary
        $esummary_result = &Eutil ("esummary", $qs);
        
        # Extract and output information for all genes
        if ($dopt eq 'csv'){
            &extractAndOutput ($esummary_result);
        } else {
            &extractAndOutputXml ($esummary_result, $first);
        }

        $first = 0;
        
        if ($verbose) {
            print STDOUT "\ngeneCount: $geneCount\n";
        }
    }

    # Close xml output
    if ($dopt eq 'xml') {print STDOUT $stopTag}
}

################################################################################
#
sub extractAndOutput {
    my $xml = $_[0];

    my $match;

    while ($xml =~ m/<DocSum>(.*?)<\/DocSum>/gs) {
        $match = $1;

        my $geneId = "";
        my $name = "";
        my $description = "";
        my $chromosome = "";
        my $cytoloc = "";
        my $chrstart = "";
        my $chrstop = "";

        if ($match =~ /<Id>(\d+)<\/Id>/) {$geneId = $1}
        if ($match =~ m/<Item Name=.Name.*>(.+)<\/Item>/) {$name = $1}
        if ($match =~ m/<Item Name=.Chromosome.*>(.+)<\/Item>/) {$chromosome = $1}
        if ($match =~ m/<Item Name=.MapLocation.*>(.+)<\/Item>/) {$cytoloc = $1}
        if ($match =~ m/<Item Name=.Description.*>(.+)<\/Item>/) {$description = $1}
        if ($match =~ m/<Item Name=.ChrStart.*>(.+)<\/Item>/) {$chrstart = $1}
        if ($match =~ m/<Item Name=.ChrStop.*>(.+)<\/Item>/) {$chrstop = $1}
        if ($chrstart){
            print STDOUT "$geneId,\"$name\",\"$description\",$chromosome,$cytoloc,$chrstart,$chrstop\n";
        }
    }
}

################################################################################
#
sub extractAndOutputXml {
    my $xml = shift;
    my $first = shift;
    
    # Strip leading <!DOCTYPE, <?xml, and bracketing  <eSummaryResult> tags,
    # from all but first set

    # Strip closing <eSummaryResult> tag from all
    
    if ($first) {
        my $start = 0;
        my $stop = index ($xml, $stopTag);
        print STDOUT substr ($xml, $start, $stop-$start);
    } else {
        my $start = index ($xml, $startTag) + length($startTag);
        my $stop = index ($xml, $stopTag);
        print STDOUT substr ($xml, $start, $stop-$start);
    }
}

################################################################################
#
# Subroutine to handle all eutil calls
#
# Create a BEGIN block to provide effect of static data for sub
BEGIN {
    # static data
    my $lastEutilTime = 0; # init to avoid delay on first Eutil 
    
    # storing local constants here too.
    my $eutilBaseUrl = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils";
    my $minUtilPeriod = 3; # seconds
    
    sub Eutil {
        my ($eutil, $querystringData) = @_;
        
        my $elapsed = time() - $lastEutilTime;
        my $delay = $minUtilPeriod - $elapsed;
        if ($delay < 0) {$delay = 0}
        sleep ($delay);
        $lastEutilTime = time();  # save for next time

        my $eutilUrl = "$eutilBaseUrl/$eutil.fcgi?$querystringData";
        if ($verbose) {print STDOUT "\neutilUrl: $eutilUrl\n";}
        
        my $result = get($eutilUrl);    # get result of the eutil for return
        if ((not defined $result)  or  ($result eq ""))
        {
            $result = ""; # Simplify error testing on return
            print STDERR "$eutil failed for request: $querystringData\n\n"; 
        } 
                
        if ($verbose) {
            print STDOUT "\neutil result: $result\n";
            my $elapsedTime = time() - $lastEutilTime;
            print STDOUT "$eutil took $elapsedTime seconds\n";
        }
        
        $result; # for return
    }
}

################################################################################
#
sub end {
    if ($verbose) {
        my $elapsedTime = time() - $beginTime;
        print STDOUT "\nElapsed time: $elapsedTime seconds\n"; 
    }
    exit;
}
