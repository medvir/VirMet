#!/usr/bin/perl

#===============================================================================
#   Author: Robert SCHMIEDER, Computational Science Research Center @ SDSU, CA
#
#   File: deconseq
#   Date: 2013-05-05
#   Version: 0.4.3
#
#   Usage:
#      deconseq [options] -f <file> -dbs <list> -dbs_retain <list> ...
#
#      Try 'deconseq -h' for more information.
#
#    Purpose: DeconSeq will help you remove unwanted sequences (contaminations)
#             from your sequence data sets by using the BWA-SW algorithm.
#
#===============================================================================

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";

use DeconSeqConfig;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path); #requires version 2.07 or newer
use Cwd;

$| = 1; # Do not buffer output

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions(\%params, 'help|h', 'man', 'no_seq_out', 'keep_tmp_files', 'dbs=s', 'dbs_retain=s', 'f=s', 'out_dir=s', 'i=i', 'c=i', 'group=i', 'id=s', 'version' => sub { print VERSION_INFO."\n"; exit; }, 'show_dbs' => sub { print $_." - ".DBS->{$_}->{name}."\n" foreach(sort keys %{(DBS)}); exit; }, 'S=i', 'z=i', 'T=i') or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

DeconSeq - DECONtamination of SEQuence data

=head1 VERSION

0.4.3

=head1 SYNOPSIS

deconseq [options] -f <file> -dbs <list> -dbs_retain <list> ...

=head1 DESCRIPTION

DeconSeq will help you remove unwanted sequences (contaminations) from your sequence datasets by using the altered BWA-SW algorithm.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Prints the help message and exists.

=item B<-man>

Prints the full documentation.

=item B<-version>

Prints the version of the program.

=item B<-show_dbs>

Prints a list of available databases.

=item B<-f> <file>

Input file in FASTA or FASTQ format that contains the query sequences.

=item B<-dbs> <list>

Name of database(s) to use (default: hsref). Names are according to their definition in the config file. Separate multiple database names by comma without spaces.

Example: -dbs hs1,hs2,hsref

=item B<-dbs_retain> <list>

Name of database(s) to use for cross-check. Query sequences with hit against any dbs will be compared to these databases. Databases have to be different from names in dbs. Names are according to their definition in the config file. Separate multiple database names by comma without spaces.

Example: -dbs_retain bact,vir

=item B<-out_dir> <dir>

Directory where the results should be written (default: .). If the directory does not exist, it will be created.

=item B<-i> <integer>

Alignment identity threshold in percentage (integer from 1-100 without %) used to define matching sequences as similar. The identity is calculated for the part of the query sequence that is aligned to a reference sequence. For example, a query sequence of 100 bp that aligns to a reference sequence over the first 50 bp with 40 matching positions has an identity value of 80%.

=item B<-c> <integer>

Alignment coverage threshold in percentage (integer from 1-100 without %) used to define matching sequences as similar. The coverage is calculated for the part of the query sequence that is aligned to a reference sequence. For example, a query sequence of 100 bp that aligns to a reference sequence over the first 50 bp with 40 matching positions has an coverage value of 50%.

=item B<-group> <integer>

If dbs_retain is set, then this option can be used to group the sequences similar to dbs and dbs_retain databases with either the clean or the contamination output file. If group is not set and dbs_retain is set, then three separate files will be generated.

Use B<-group 1> for grouping "Clean + Both" and use B<-group 2> for grouping "Contamination + Both".

=item B<-no_seq_out>

Prevents the generation of the fasta/fastq output file for the given coverage and identity thresholds. This feature is e.g. useful for the web-version since -i and -c are set interactively and not yet defined at the data processing step.

=item B<-keep_tmp_files>

Prevents from unlinking the generated tmp files. These usually include the id file and the .tsv file(s). This feature is e.g. useful for the web-version since .tsv files are used to dynamically generate the output files.

=item B<-id> <string>

Optional parameter. If not set, ID will be automatically generated to prevent from overwriting previous results. This option is useful if integrated into other tools and the output filenames need to be known.

=item B<-S> <integer>

Chunk size of reads in bp as used by BWA-SW (default: 10000000).

=item B<-z> <integer>

Z-best value as used by BWA-SW (default: 1).

=item B<-T> <integer>

Alignment score threshold as used by BWA-SW (default: 30).

=back

=head1 AUTHOR

Robert SCHMIEDER, C<< <rschmieder_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <rschmieder_at_gmail_dot_com> >> so that I can make DeconSeq better.

=head1 COPYRIGHT

Copyright (C) 2010-2011  Robert SCHMIEDER

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

#
################################################################################
## DATA AND PARAMETER CHECKING
################################################################################
#

#Check if input file exists and estimate file format
my $format = '';
if(exists $params{f}) {
    if(-e $params{f}) {
        #check for file format
        $format = &checkFileFormat($params{f});
        unless($format eq 'fasta' || $format eq 'fastq') {
            &printError('input file for -f is in '.uc($format).' format not in FASTA or FASTQ format');
        }
    } else {
        &printError("could not find input file \"".$params{f}."\"");
    }
} else {
    &printError("you did not specify an input file containing the query sequences");
}

#Check if output dir is defined and if not, if it can be created
if(exists $params{out_dir}) {
    make_path($params{out_dir}, { verbose => 0, mode => 0711, error => \my $err });
    if(@$err) {
        &printError("problems while trying to create directory \"".$params{out_dir}."\"");
    }
} else {
    $params{out_dir} = cwd();
}
$params{out_dir} .= '/' unless($params{out_dir} =~ /\/$/);

#Check if databases are defined, if their names are valid, and if they exist
if(exists $params{dbs}) {
    my @tmp = split(/\,/,$params{dbs});
    foreach my $db (@tmp) {
        unless(exists DBS->{$db}) {
            &printError("database \"".$db."\" does not exist in config file");
        }
        unless(&checkForDbFiles($db)) {
            &printError("cannot find all database files for database \"".$db."\" in dir \"".DB_DIR."\"");
        }
    }
} else {
    $params{dbs} = DB_DEFAULT;
    unless(&checkForDbFiles(DB_DEFAULT)) {
        &printError("cannot find all database files for database \"".DB_DEFAULT."\" in dir \"".DB_DIR."\"");
    }
}

#Check if keep databases are defined, if their names are valid, and if they exist
if(exists $params{dbs_retain}) {
    my @tmp = split(/\,/,$params{dbs_retain});
    foreach my $db (@tmp) {
        unless(exists DBS->{$db}) {
            &printError("database \"".$db."\" does not exist in config file");
        }
        unless(&checkForDbFiles($db)) {
            &printError("cannot find all database files for database \"".$db."\" in dir \"".DB_DIR."\"");
        }
    }
}

#Check if databases and keep databases are not the same
if(exists $params{dbs_retain}) {
    my @tmp = split(/\,/,$params{dbs});
    my @tmp_keep = split(/\,/,$params{dbs_retain});
    foreach my $keep (@tmp_keep) {
        if(grep {/^$keep$/} @tmp) {
            &printError("cannot use the same database(s) for dbs AND dbs_retain");
        }
    }
}

#Check if thresholds are defined and valid
if(exists $params{i}) {
    unless($params{i} > 0 && $params{i} <= 100) {
        &printError("invalid value for identity threshold");
    }
} else {
    $params{i} = 1;
}
if(exists $params{c}) {
    unless($params{c} > 0 && $params{c} <= 100) {
        &printError("invalid value for coverage threshold");
    }
} else {
    $params{c} = 1;
}

#Check if group is set correctly
if(exists $params{group}) {
    unless($params{group} == 1 || $params{group} == 2) {
        &printError("invalid group option");
    }
    unless(exists $params{dbs_retain}) {
        &printError("group can only be used when dbs_retain is set");
    }
}

#
################################################################################
## DATA PROCESSING
################################################################################
#

## Get new id for data
unless(exists $params{id}) {
    $params{id} = &generateNewIdFile($params{out_dir});
}

## Run modified BWA-SW for input file on all selected databases
my @dbs = split(/\,/,$params{dbs});
my @dbs_retain;
@dbs_retain = split(/\,/,$params{dbs_retain}) if(exists $params{dbs_retain});
my %data;
my ($time,$seqnum1,$seqnum2,$seqnum3);
foreach my $db (@dbs) {
    my @refdbs = split(/\,/,DBS->{$db}->{db});
    foreach my $refdb (@refdbs) {
        my $tsvfile = $params{out_dir}.$params{id}.'_'.$db.'_'.$refdb.'.tsv';
        my $cmd = PROG_DIR.PROG_NAME.' bwasw -t 24 -A -f '.$tsvfile.(exists $params{S} ? ' -S '.$params{S} : '').(exists $params{z} ? ' -z '.$params{z} : '').(exists $params{T} ? ' -T '.$params{T} : '').' '.DB_DIR.$refdb.' '.$params{f};
        unless(-e $tsvfile) {
            open(TSV, ">$tsvfile") or &printError("Could not create tsv file $tsvfile: $!");
            close(TSV);
        }
        $time = &runCmd($cmd);
        print STDERR "[deconseq] CMD finished in $time seconds\n" if(DEBUG);
        $data{$db.$refdb}->{tsv} = $tsvfile;
    }
}


## If dbs_retain is defined
if(@dbs_retain) {
    # Parse TSV file(s) to get unique query ids
    my %ids;
    foreach my $db (@dbs) {
        my @refdbs = split(/\,/,DBS->{$db}->{db});
        foreach my $refdb (@refdbs) {
            &parseTsvFiles($data{$db.$refdb}->{tsv},\%ids);
        }
    }
    if(scalar(keys %ids)) {
        #generate FASTA file for query ids
        my $fasta = $params{out_dir}.$params{id}.'_tmp.fa';
        ($time,$seqnum1) = &generateFileFromIds(\%ids,$params{f},'fasta',$fasta);
        print STDERR "[deconseq] Generated file ".$fasta." with $seqnum1 sequences in $time seconds\n" if(DEBUG);

        # Run BWA-SW for input file on all selected databases
        foreach my $db (@dbs_retain) {
            my @refdbs = split(/\,/,DBS->{$db}->{db});
            foreach my $refdb (@refdbs) {
                my $tsvfile = $params{out_dir}.$params{id}.'_'.$db.'_'.$refdb.'.tsv';
                my $cmd = PROG_DIR.PROG_NAME.' bwasw -A -f '.$tsvfile.(exists $params{S} ? ' -S '.$params{S} : '').(exists $params{z} ? ' -z '.$params{z} : '').(exists $params{T} ? ' -T '.$params{T} : '').' '.DB_DIR.$refdb.' '.$params{f};
                 unless(-e $tsvfile) {
                    open(TSV, ">$tsvfile") or &printError("Could not create tsv file $tsvfile: $!");
                    close(TSV);
                }
                $time = &runCmd($cmd);
                print STDERR "[deconseq] CMD finished in $time seconds\n" if(DEBUG);
                $data{$db.$refdb}->{tsv} = $tsvfile;
            }
        }

        unlink($fasta);
    } else {
        print STDERR "[deconseq] No sequences for dbs_retain runs\n" if(DEBUG);
    }
}


## Generate output files
if(exists $params{no_seq_out}) {
    print STDERR "[deconseq] Did not generate sequence output files, since no_seq_out was set as parameter.\n" if(DEBUG);
} else {
    ## Parse TSV files and combine results
    my %ids;
    foreach my $db (@dbs) {
        my @refdbs = split(/\,/,DBS->{$db}->{db});
        foreach my $refdb (@refdbs) {
            &parseTsvFiles($data{$db.$refdb}->{tsv},\%ids);
        }
    }
    my %ids_keep;
    foreach my $db (@dbs_retain) {
        my @refdbs = split(/\,/,DBS->{$db}->{db});
        foreach my $refdb (@refdbs) {
            &parseTsvFiles($data{$db.$refdb}->{tsv},\%ids_keep) if(exists $data{$db.$refdb}->{tsv});
        }
    }
    my $filterids = &combineResults(\%ids,\%ids_keep,$params{c},$params{i},$params{id});


    ## Generate output files with query sequence subsets
    ($time,$seqnum1,$seqnum2,$seqnum3) = &generateFileFromIds($filterids,$params{f},$format,$params{out_dir}.$params{id}.'_cont.f'.($format eq 'fasta' ? 'a' : 'q'),$params{out_dir}.$params{id}.'_clean.f'.($format eq 'fasta' ? 'a' : 'q'),$params{out_dir}.$params{id}.'_both.f'.($format eq 'fasta' ? 'a' : 'q'));
    print STDERR "[deconseq] Generated files ".$params{out_dir}.$params{id}.'_cont.f'.($format eq 'fasta' ? 'a' : 'q').", ".$params{out_dir}.$params{id}.'_clean.f'.($format eq 'fasta' ? 'a' : 'q').", ".$params{out_dir}.$params{id}.'_both.f'.($format eq 'fasta' ? 'a' : 'q')." with $seqnum1, $seqnum2, and $seqnum3 sequences in $time seconds\n" if(DEBUG);
}


## Clean up temp files generated
#tsv files and id file
if(exists $params{keep_tmp_files}) {
    print STDERR "[deconseq] Did not unlink files, since keep_tmp_files was set as parameter.\n" if(DEBUG);
} else {
    print STDERR "[deconseq] Unlink ".$params{out_dir}.$params{id}."\n" if(DEBUG);
    unlink($params{out_dir}.$params{id});
    foreach my $db (keys %data) {
        print STDERR "[deconseq] Unlink ".$data{$db}->{tsv}."\n" if(DEBUG);
        unlink($data{$db}->{tsv});
    }
}

#TODOs
#add timing to all parts for later benchmarks and analytics

#
################################################################################
## MISC FUNCTIONS
################################################################################
#

sub printError {
    my $msg = shift;
    print STDERR "ERROR: ".$msg.".\n\nTry \'deconseq -h\' for more information.\nExit program.\n";
    exit(0);
}

sub checkForDbFiles {
    my $db = shift;
    my @exts = qw(amb ann bwt pac rbwt rpac rsa sa);
    my $status = 1;
    my @tmp = split(/\,/,DBS->{$db}->{db});
    foreach my $tmpdb (@tmp) {
        foreach my $ext (@exts) {
            unless(-e DB_DIR.$tmpdb.'.'.$ext) {
                $status = 0;
                last;
            }
        }
    }
    return $status;
}

sub generateNewIdFile {
    my $dir = shift;
    my $id;
    $id = time();
    while(-e $dir.$id) {
	$id = time();
    }
    `touch $dir$id`;
    return $id;
}

sub runCmd {
    my $cmd = shift;
    my $time = time();
    print STDERR "[deconseq] Run CMD: $cmd\n" if (DEBUG);
    system($cmd) == 0 or &printError("system call \"$cmd\" failed: $?");
    return (time()-$time);
}

sub parseTsvFiles {
    my ($file,$ids) = @_;
    my (@lineargs);
    open(IN,"<$file") or die "ERROR: could not open file $file: $! \n";
    while(<IN>) {
        next if(/^\#/);
        chomp();
        @lineargs = split(/\t/);
        #query id, dbid, db start, db aligned, coverage, identity
        next if($lineargs[4] < 1 || $lineargs[5] < 1); #problems with hits from bwasw (hit over more than 2 reference sequences);
        next if(exists $ids->{$lineargs[0]} && exists $ids->{$lineargs[0]}->{$lineargs[4]} && $ids->{$lineargs[0]}->{$lineargs[4]} >= $lineargs[5]);
        $ids->{$lineargs[0]}->{$lineargs[4]} = $lineargs[5];
    }
    close(IN);
    return $ids;
}

sub generateFileFromIds {
    my ($ids,$in,$format,$out1,$out2,$out3) = @_;
    my $time = time();
    my $group;
    if(exists $params{dbs_retain}) {
        if(exists $params{group}) {
            $group = $params{group};
        } else {
            $group = 0;
        }
    } else {
        $group = -1;
    }
    my ($id,$seqnum1,$seqnum2,$seqnum3,$type);
    $seqnum1 = $seqnum2 = $seqnum3 = $id = $type = 0;
    open(IN,"<$in") or die "ERROR: could not open file $in: $! \n";
    open(OUT1,">$out1") or die "ERROR: could not write to file $out1: $! \n"; #cont
    open(OUT2,">$out2") or die "ERROR: could not write to file $out2: $! \n" if(defined $out2); #clean
    open(OUT3,">$out3") or die "ERROR: could not write to file $out3: $! \n" if($group == 0 && defined $out3); #both

    if($format eq 'fastq') { #fastq input
        my $count = 0;
        while (<IN>) {
            if($count == 0 && /^\@(\S+)/) {
                $id = $1;
                if(exists $ids->{$id}) {
                    if($ids->{$id} == 1) { #cont
                        $seqnum1++;
                        $type = 1;
                    } elsif($group == 0) { #both
                        $seqnum3++;
                        $type = 3;
                    } elsif($group == 1) { #clean+both
                        $seqnum2++;
                        $type = 2;
                    } elsif($group == 2) { #cont+both
                        $seqnum1++;
                        $type = 1;
                    }
                } else { #clean
                    $type = 2;
                    $seqnum2++;
                }
            } elsif($count == 3) {
                $count = -1;
            }
            $count++;
            if($type == 1) {
                print OUT1 $_;
            } elsif($type == 2) {
                print OUT2 $_ if(defined $out2);
            } elsif($type == 3) {
                print OUT3 $_ if(defined $out3);
            }
        }
    } elsif($format eq 'fasta') {
        while(<IN>) {
            if(/^>(\S+)/) {
                $id = $1;
                if(exists $ids->{$id}) {
                    if($ids->{$id} == 1) { #cont
                        $seqnum1++;
                        $type = 1;
                    } elsif($group == 0) { #both
                        $seqnum3++;
                        $type = 3;
                    } elsif($group == 1) { #clean+both
                        $seqnum2++;
                        $type = 2;
                    } elsif($group == 2) { #cont+both
                        $seqnum1++;
                        $type = 1;
                    }
                } else { #clean
                    $type = 2;
                    $seqnum2++;
                }
            }
            if($type == 1) {
                print OUT1 $_;
            } elsif($type == 2) {
                print OUT2 $_ if(defined $out2);
            } elsif($type == 3) {
                print OUT3 $_ if(defined $out3);
            }
        }
    } else {
        die "ERROR: unknow file format: $format\n";
    }
    close(OUT1);
    close(OUT2) if(defined $out2);
    close(OUT3) if($group == 0 && defined $out3);
    close(IN);
    return (time()-$time,$seqnum1,$seqnum2,$seqnum3);
}

sub combineResults {
    my ($ids,$ids_keep,$c,$i) = @_;
    my (%filterids,$filter);
    foreach my $qid (keys %$ids) {
        $filter = 0;
        while(my ($cov,$ident) = each(%{$ids->{$qid}}) ) {
            if($cov >= $c && $ident >= $i) {
                $filter = 1;
                #check if in $ids_keep with same or better alignment values
                if(exists $ids_keep->{$qid}) {
                    while(my ($cov_keep,$ident_keep) = each(%{$ids_keep->{$qid}}) ) {
                        if($cov_keep >= $cov && $ident_keep >= $ident) {
                            $filter = 2;
                            last;
                        }
                    }
                }
                last;
            }
        }
        if($filter) {
            $filterids{$qid} = $filter;
        }
    }
    return \%filterids;
}

sub checkFileFormat {
    my $file = shift;

    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    while (<FILE>) {
        chomp();
        next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        }
    }
    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}
