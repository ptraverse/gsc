#!/usr/bin/perl

# barnacle.pl
#
# Created by Lucas Swanson
# Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.

# TODO
# write intermediate output between grouping and adding support/annotations

# add path to directory holding CPAN perl libraries required for Barnacle
use FindBin;
use lib "$FindBin::Bin/perl_libs";

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Set::IntSpan;
use Benchmark;
use Config::General;
use File::Basename;
use File::Path 2.07 qw(make_path remove_tree);
use File::Spec::Functions qw(catdir catfile file_name_is_absolute splitpath);
use Cwd qw(abs_path getcwd);
use POSIX qw(ceil);

# store the arguments passed in, to write to log file
my $barnacle_cmd = $0." ".join(" ", @ARGV);

my %OPT;
GetOptions(\%OPT,
     "lib=s",
     "lib_dir=s",
     "config=s",
     #"outdir=s",
     "identify_candidates",
     "add_support",
     "cid_memory=s",
     "use_smart_chooser",
     "use_quick_chooser",
     "max_match_fraction=f",
     "min_merge_overlap=f",
     "min_ctg_represented=f",
     "num_aligns=i",
     "min_identity=f",
     "mito_prefilter_on",
     "mito_prefilter_off",
     "flag_exon_boundary_junctions",
     "no_flag_exon_boundary_junctions",
     "prefer_exons",
     "no_exon_preference",
     "max_num_groups=i",
     "pair_to_genome",
     "no_pair_to_genome",
     "pair2gen_split=i",
     "p2g_memory=s",
     "frag_len=i",
     "frag_diff=f",
     "min_mapq=i",
     "read_to_contig",
     "no_read_to_contig",
     "r2c_memory=s",
     "r2c_min_overlap=i",
     "r2c_ctgs_per_job=i",
     "r2c_short_ctgs",
     "r2c_no_short_ctgs",
     "flag_RNA",
     "no_flag_RNA",
     "noMito",
     "flag_repeats",
     "no_flag_repeats",
     "repeat_search_size",
     "breakpoint_genes",
     "no_breakpoint_genes",
     "read_length=i",
     "no_split_check",
     "no_gap_check",
     "check_split_and_gap",
     "include_gap_events",
     "min_gap_size=i",
     "min_gap_identity=f",
     "min_gap_fraction=f",
     "gap_max_len=i",
     "max_sam_retries=i",
     "cluster=s",
     "hostname=s",
     "queue=s",
     "use_wall_time",
     "email=s",
     "assembly_ver=s",
     "assembler",
     "ver=s",
     "disable_profiling_timer",
     "debug",
     "debug_high",
     "help|h",
     "man|m",
     "ignore_cmds=s",
     "robust_r2c",
     );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{lib_dir} || !$OPT{config} ||
  (!$OPT{identify_candidates} && !$OPT{add_support}));

=pod

=head1 SYNOPSIS

barnacle.pl -lib LIB -lib_dir LIB_DIR -config FILE (-identify_candidates or -add_support) [ OPTIONS ]

NOTE:

=over

=item LIB: library name

=item LIB_DIR: the base directory for the library assembly

=back

=head1 OPTIONS

 -identify_candidates
                    Find putative trans-splicing events from
                    contigs-to-genome alignment files.
 -add_support
                    Group putative trans-splicing events and add support.
                    Use this option only after using the identify_candidates
                    option.
 -cid_memory N      Amount of memory to request when submitting candidate
                    identification jobs to the cluster [ default: 8G ]
 -use_smart_chooser
                    Use the smarter, but slower, algorithm for choosing the
                    best alignments for a contig [ default ]
 -use_quick_chooser
                    Use the more naive, but quicker, algorithm for choosing
                    the best alignments for a contig
 -max_match_fraction F
                    If any single alignment for a contig uses at least F%
                    of the contig length, do not look for any split
                    alignments for that contig. [ default: 0.999 ]
 -min_merge_overlap F
                    If at least F% of two alignments represent the same
                    region of the contig, consider only the higher scoring
                    of the two alignments. [ default: 0.80 ]
 -min_ctg_represented F
                    Only consider split alignments for which the two
                    alignments represent at least F% of the total contig.
                    [ default: 0.85 ]
 -num_aligns N      Extract the best N alignments for each contig from the
                    alignments file. [ default: 500 ]
 -min_identity PID  Extract from the alignments file only alignments with a
                    percent identity of at least PID%. [ default: 40.0 ]
 -mito_prefilter_on
                    Filter out events involving mitochondrial DNA when first
                    searching the alignments [ default ]
 -mito_prefilter_off
                    Do not filter out events involving mitochondrial DNA
                    when first searching the alignments
 -flag_exon_boundary_junctions
                    Mark events that have junctions matching up with exon
                    boundaries [ default ]
 -no_flag_exon_boundary_junctions
                    Do not check which events have junctions matching up
                    with exon boundaries
 -exon_bound_buffer N
                    Count the junction coordinate as matching an exon
                    coordinate if it is within Nbp to either side of it.
                    [ default:4 ]
 -prefer_exons      Prefer alignments that overlap exons when paring
                    alignment groups
 -no_exon_preference
                    Do not consider exon overlap when paring alignment
                    groups [ default ]
 -max_num_groups N  At least one contig showing the event must appear in no
                    more than N groups [ default: 0 (off) ]
 -pair_to_genome    Calculate the number of read-pairs aligning across the
                    genomic event region
 -no_pair_to_genome Do not calculate the number of read-pair to genome
                    alignments supporting the event [ default ]
 -pair2gen_split N  When splitting the read-pair to genome jobs, examine
                    N events in each job [ default: 5000 ]
 -p2g_memory N      Amount of memory to request when submitting pair-to-genome
                    support jobs to the cluster [ default: 5G ]
 -frag_len N        When looking for read-pairs spanning the genomic event
                    region, use an expected fragment length of N bp.
                    [ default: 200 ]
 -frag_diff F       When looking for read-pairs spanning the genomic event
                    region, only count a pair as being significantly
                    different from the expected length if the fractional
                    difference between the observed and expected fragment
                    length is more than F [ default: 0.10 ]
 -min_mapq N        When looking for read-pairs spanning the genomic event
                    region, only count pairs with a mapping quality greater
                    than N. [ default: 10 ]
 -read_to_contig    Filter events based on the number of reads aligning
                    to the event region of the contigs in the
                    event [ default ]
 -no_read_to_contig Do not calculate the number of read to contig alignments
                    supporting the event
 -r2c_memory N      Amount of memory to request when submitting read-to-contig
                    support jobs to the cluster [ default: 5G ]
 -r2c_min_overlap N When calculating read-to-contig alignment support,
                    only count reads that overlap the event region by at
                    least N bp [ default: 5 ]
 -r2c_ctgs_per_job N
                    When creating read-to-contig support jobs,
                    put N contigs in each job [ default: 1500 ]
 -r2c_short_ctgs    When creating read-to-conting support jobs, create
                    jobs for all contigs
 -r2c_no_short_ctgs When creating read-to-conting support jobs, do not
                    create jobs when all contigs are shorter than the read
                    length [ default ]
 -flag_RNA          Flag events overlapping smal structural RNAs [ default ]
 -no_flag_RNA       Do not flag events overlapping small structural RNAs
 -noMito            Fail events occuring purely in mitochondrial DNA
 -flag_repeats
                    Use the repeat coordinate file(s) in the config file
                    to check for and mark alignment blocks overlapping
                    repeat sequences within the genome [ default ]
 -no_flag_repeats   Do not check whether any of the alignment blocks overlap
                    repeat sequences
 -repeat_search_size N
                    Look for repeats in a region N bp wide around the event
                    coordinates. [default: 100]
 -breakpoint_genes  Mark each event with the gene features in which the
                    breakpoints occur. [ default ]
 -no_breakpoint_genes
                    Do not mark events with the gene features in which the
                    breakpoints occur.
 -read_length N     Each read is N bp long [ default: 75 ]
 -check_split_and_gap
                    Examine contig-to-genome alignments to identify contigs
                    matching either a split-alignment or a gapped-alignment
                    event signature [ default ]
 -no_split_check    Do not attempt to identify contigs matching a split-
                    alignment event signature
 -no_gap_check      Do not attempt to identify contigs matching a gapped-
                    alignment event signature
 -include_gap_events
                    Parse results of examining contig-to-genome alignments
                    for contigs matching gapped-alignment signatures
 -min_gap_size N    Only consider gaps of size N or greater [ default: 4 ]
 -min_gap_identity F
                    Require that gap sequences map to the rest of the
                    contig with at least F% identity [ default: 0.95 ]
 -min_gap_fraction F
                    Require that at least F% of the gap sequence bases map
                    to the rest of the contig [ default: 0.30 ]
 -gap_max_len N     Only run the gap-realigner on contigs shorter than
                    Nbp long. [ default: 50000 ]
 -max_sam_retries N
                    If there is an error running a samtools view command,
                    retry the command a maximum of N times. [ default: 5 ]
 -cluster CLUSTER_HEAD
                    Specify cluster head node if running jobs on cluster
 -hostname HOSTNAME The hostname(s) to submit to [default: read from config]
 -queue QUEUE       The queue(s) to submit to [default: read from config]
 -use_wall_time     Use the wall-time option when submitting jobs to the
                    cluster (not all versions of mqsub support this option)
 -email             When submitting jobs to the cluster, e-mail status
                    updates to the given email address
 -assembly_ver VER  The version of the assembly to use when constructing paths
 -assembler NAME    The name of the assembler used
 -ver VER           Use VER as the output version
 -disable_profiling_timer
                    Sometimes this script can hang when trying to spawn
                    child processes, due to the profiling timer used by
                    the kernel. Use this option to disable the profiling
                    timer if the script seems to be hanging
 -debug             Print status info while running
 -debug_high        Print even more status info
 -h, -help          Display a brief help message
 -man               Display the full documentation
 -ignore_cmds FILE  Ignore alignment files listed in FILE (for continuing
                    interrupted runs)

=head1 NAME

barnacle.pl -> Wrapper to identify candidate contigs from contig-to-genome
               alignments and add support and annotations

=head1 DESCRIPTION

December 2010

Script that runs all the BARNACLE code for a library and generates a single summary file with a PASS/FAIL status

=head1 AUTHOR

Lucas Swanson

=cut

# -outdir DIR        Output_directory [ default: INPUT_DIR/barnacle/ver_VER ]

my $run_start = new Benchmark;

# Candidate Identification Exit Statuses
my $CI_SUCCESS = 0;
my $CI_NO_ALIGNMENTS = 3;
# candidate contig output columns
my $CC_CTG_ID        = 0;
my $CC_TOPOLOGY      = 1;
my $CC_TARGET        = 2;
my $CC_CTG_COORDS    = 3;
my $CC_JUNCTION      = 4;
my $CC_ALIGN_METRICS = 5;
my $CC_BLOCKS        = 6;
my $CC_GENES         = 7;
my $CC_NEARBY        = 8;
my $CC_META          = 9;
my $CC_SEQ           = 10;
# Event group field indices
my $EGI_JUNCT_FROM    = 0;
my $EGI_JUNCT_TO      = 1;
my $EGI_CTG_TYPE      = 2;
my $EGI_ALIGNER       = 3;
my $EGI_CTG_ID        = 4;
my $EGI_TOPOLOGY          = 5;
my $EGI_TARGET        = 6;
my $EGI_CTG_COORDS    = 7;
my $EGI_JUNCTION      = 8;
my $EGI_ALIGN_METRICS = 9;
my $EGI_BLOCKS        = 10;
my $EGI_GENES         = 11;
my $EGI_NEARBY        = 12;
my $EGI_META          = 13;
my $EGI_SEQ           = 14;
# SamTools constants
my $SAM_BOTH       = 0;
my $SAM_UPSTREAM   = 1;
my $SAM_DOWNSTREAM = 2;
# Sam output constants
my $SAM_READ_ID = 0;
my $SAM_FLAG    = 1;
my $SAM_CTG_ID  = 2;
my $SAM_LEFT    = 3;

die ("missing library name: -lib LIB\n") unless $OPT{lib};
#die ("missing project name: -project PROJECT\n") unless $OPT{project};
die ("must use exactly one of -identify_candidates or -add_support\n") unless
  (  $OPT{identify_candidates} && ! $OPT{add_support}) ||
  (! $OPT{identify_candidates} &&   $OPT{add_support});

$OPT{barnacle_src_dir} = abs_path(dirname(__FILE__));
$OPT{barnacle_dir} = dirname($OPT{barnacle_src_dir});
print "BARNACLE script found in directory: \"".$OPT{barnacle_src_dir}."\"\n";
&CheckOptionsForConflicts();
&SetOptionDefaults();

# create hash for executables
my %EXEC;
$EXEC{summarize_script}=catfile($OPT{barnacle_src_dir}, "summarize.py");

# setup input and output directories
#print "Rel LIB DIR: \"".$OPT{lib_dir}."\".\n";
&error("Cannot find library directory: \"".$OPT{lib_dir}."\"") unless
  (-d $OPT{lib_dir});
my $raw_lib_dir = $OPT{lib_dir};
my $abs_lib_dir = abs_path($OPT{lib_dir});
#print "Abs LIB DIR: \"".$abs_lib_dir."\".\n";
chomp $abs_lib_dir;
$abs_lib_dir =~ s/\/$//; # remove final slash(/) if necessary
$OPT{lib_dir} = $abs_lib_dir;

# setup configuration options from config file
$OPT{config} = abs_path($OPT{config});
&CheckPath($OPT{config}, "config");
my %CFG;
&ReadConfigFile();

#my $input_dir = catdir($OPT{lib_dir}, "Assembly", $OPT{assembly_ver});
&error("Input dir does not exist: \"".$OPT{input_dir}."\"") unless (-d $OPT{input_dir});
my $out_base_dir = &SetupOutputBaseDirectory();
&GetOutputVersion($out_base_dir);
($OPT{results_dir}, $OPT{outdir}) = &SetupOutputDirectory($out_base_dir);
# Open the log file
my ($log_file_name, $logref) = &OpenLogFile();
&status("BARNACLE VERSION ".$OPT{code_version});
&status("OUTPUT VERSION ".$OPT{ver});
&status("COMMAND: ".$barnacle_cmd."\n");
&debug("Raw Lib dir: ".$raw_lib_dir);
&debug("Abs Lib dir: ".$OPT{lib_dir});
&debug("Input dir: ".$OPT{input_dir});

# need output directory to finish setting-up pair-to-genome info
unless ($OPT{no_pair_to_genome} || $OPT{identify_candidates}) {
  &FinishPairToGenomeSetup();
}

#my $merge_dir = catdir($input_dir, "merge");
#&debug("Merge dir: ".$merge_dir);
#unless (-e $merge_dir) {
#  &error("Can't find merge directory: \"$merge_dir\"");
#}

# if using a cluster
if ($OPT{cluster}) {
  # get the submit jobs script path and check that it is valid
  $EXEC{submit_script}=catfile($OPT{barnacle_src_dir}, "utils", "submit.py");
  &CheckPath($EXEC{submit_script}, "job submission script");
}

# if identifying candidates
if ($OPT{identify_candidates}) {
  # get the identify_candidates script path and check that it is valid
  $EXEC{identify_candidates}=catfile($OPT{barnacle_src_dir}, "alignment_processing",
    "identify_candidate_contigs.py");
  &CheckPath($EXEC{identify_candidates}, "candidate identifier script");
  # if using the gap realigner
  unless ($OPT{no_gap_check}) {
    # get the gap realigner script path and check that it is valid
    $EXEC{gap_realigner}=catfile($OPT{barnacle_src_dir}, "alignment_processing",
      "gap_realigner");
    if ($OPT{cluster}) {
      $EXEC{gap_realigner} .= "_cluster";
    }
    &CheckPath($EXEC{gap_realigner}, "gap realigner tool");
  }
}

# if adding support
if ($OPT{add_support}) {
  # get the check events script path and check that it is valid
  $EXEC{check_format}=catfile($OPT{barnacle_src_dir}, "utils",
    "check_format.py");
  &CheckPath($EXEC{check_format}, "candidate format checking script");
  # if using the results of jobs submitted to the cluster
  if ($OPT{cluster}) {
    # get the check cluster jobs status script path and check that it is valid
    $EXEC{check_status}=catfile($OPT{barnacle_src_dir}, "alignment_processing",
      "check_status.py");
    &CheckPath($EXEC{check_status},
      "check candidate identification status script");
  }
  # if getting pair2genome support
  unless ($OPT{no_pair_to_genome}) {
    # if submitting pair2genome support jobs to the cluster
    if ($OPT{cluster_pair2gen}) {
      # get the pair2genome support script path and check that it is valid
      $EXEC{p2g_script}=catfile($OPT{barnacle_src_dir}, "support", "pair_to_genome",
        "calculate.py");
      &CheckPath($EXEC{p2g_script}, "pair-to-genome support script");
    }
  }
  # If calculating read-to-contig support
  unless ($OPT{no_read_to_contig}) {
    # get the read-to-contig support path and check that it is valid
    $EXEC{r2c_submit}=catfile($OPT{barnacle_src_dir}, "support", "read_to_contig",
      "submit.py");
    &CheckPath($EXEC{r2c_submit}, "read-to-contig support script");
  }
  # If flagging exon boundary junctions
  unless ($OPT{no_flag_exon_boundary_junctions}) {
    # get the exon-boundary junctions script path and check that it is valid
    $EXEC{exon_bounds_script}=catfile($OPT{barnacle_src_dir}, "annotation",
      "exon_bounds.py");
    &CheckPath($EXEC{exon_bounds_script}, "exon-boundary junctions script");
  }
  # if checking for repeat sequences
  if ($OPT{flag_repeats}) {
    # get the repeats script path and check that it is valid
    $EXEC{repeats_script}=catfile($OPT{barnacle_src_dir}, "annotation", "repeats.py");
    &CheckPath($EXEC{repeats_script}, "repeat-overlaps script");
  }
  # if marking breakpoint genes
  if ($OPT{breakpoint_genes}) {
    # get the repeats script path and check that it is valid
    $EXEC{breakpoint_genes_script}=catfile($OPT{barnacle_src_dir}, "annotation",
      "breakpoint_genes.py");
    &CheckPath($EXEC{breakpoint_genes_script}, "breakpoint genes script");
  }
}

my $warning;
my $msg;

# determine which types of contigs to consider
my %contig_types = map {$_ => 1} split(/,/, $OPT{contig_types});
#for my $contig_type (sort keys %contig_types) {
#  print "checking $contig_type\n";
#}
# determine which aligners types to consider
my %aligner_types = ();
&debug_high("Setting alignment file extension using: ".$OPT{c2g_template});
($aligner_types{"split"}{align_ext}) = $OPT{c2g_template} =~ /\.([^.]+$)/;
&debug_high("  Extension set to: ".$aligner_types{"split"}{align_ext});
#unless ($OPT{only_exon_aligns}) {
#  $aligner_types{blat}{align_ext} = "psl";
#}
#if ($OPT{only_exon_aligns} || $OPT{blat_and_exon_aligns}) {
#  $aligner_types{exonerate}{align_ext} = "exon";
#}

# if processing alignments, output library info to a file
if ($OPT{identify_candidates}) {
  # attempt to get proper assembly version
  #&debug("Attempting to get assembly version from ".$outdir);
  #my ($assembly_ver) = $outdir =~ /Assembly[^\/]*\/([^\/]+)\//;
  #&debug("  Assembly version: ".$assembly_ver);
  #if ("current" eq $assembly_ver) {
  #  my ($current_path) = $outdir =~ /^(.*Assembly[^\/]*\/current)\//;
  #  &debug("  Current path: ".$current_path);
  #  my $current_target = `ls -ld $current_path`;
  #  &debug("  Current target: ".$current_target);
  #  my $current_ver = "";
  #  ($current_ver) = $current_target =~
  #    /([^ ]+-[0-9]+\.[0-9]+(\.[0-9]+)?)\/?$/;
  #  &debug("  Current version: ".$current_ver);
  #  if (($OPT{assembler} && $current_ver =~ /^$OPT{assembler}/) ||
  #      (! $OPT{assembler} && $current_ver =~ /\w/)) {
  #    $assembly_ver = $current_ver;
  #    # replace "current" in the path
  #    $outdir =~ s/current/$assembly_ver/;
  #  } else {
  #    &debug("Could not get assembly version pointed to by \"current\".");
  #  }
  #}
  my $lib_info_path = catfile($OPT{outdir}, "lib_info");
  if (open(LIB_INFO, ">", $lib_info_path)) {
    print LIB_INFO $OPT{lib}.",".$OPT{assembly_ver}.",".$OPT{ver}."\n";
    close LIB_INFO;
  } else {
    &warn("Can't create library info file $lib_info_path $!\n");
  }
}

# Run the code
&Main();

# Count the number of warnings
my $num_warnings = `grep -i -c warning $log_file_name`;
chomp $num_warnings;
if (0 < $num_warnings) {
  &status("$num_warnings warnings encountered, ".
          "see $log_file_name for details.");
}

my $run_end = new Benchmark;
my $td = timediff($run_end, $run_start);
my $ts = timestr($td);
&status("Total running time: $ts");
&status("Complete");
close $logref;

sub Main {
  # No parameters
  $| = 1; # REMINDER: turn on autoflush

  # Generate raw candidate groups (either initial cid jobs, or
  # grouping and structural RNA overlaps)
  if ("stop" eq &GenerateRawCandidateGroups()) {
    # stop here if creating candidate identification jobs
    return undef;
  }

  # Remove unparsable events
  &RemoveUnparsableEvents();
  my $stage = 2;
  my $next_input_dir = catdir($OPT{results_dir}, $stage."_parsable_candidates");

  # after creating external python p2g-submit script, call it here
  # If calculating pair-to-genome support
  unless ($OPT{no_pair_to_genome}) {
    unless ($OPT{cluster_pair2gen} || !$OPT{cluster}) {
      $stage += 1;
      $next_input_dir = catdir($OPT{results_dir}, $stage."_with_p2g");
    }
  }

  # If calculating read-to-contig support
  unless ($OPT{no_read_to_contig}) {
    &CalculateReadToCtgSupportExternal($next_input_dir);
    #&GetReadToCtgSupport($groups, $contigs_in_events, $next_input_dir);
  }
  #&PrintGroupCtgIDs($groups); # DEBUG

  # If flagging exon boundary junctions
  unless ($OPT{no_flag_exon_boundary_junctions}) {
    &FlagExonBoundaryJunctions($next_input_dir);
    $stage += 1;
    $next_input_dir = catdir($OPT{results_dir}, $stage."_with_exon_bounds");
  }

  # If flagging repeats
  if ($OPT{flag_repeats}) {
    &CheckForRepeatsExternal($next_input_dir);
    $stage += 1;
    $next_input_dir = catdir($OPT{results_dir}, $stage."_with_repeats");
  }

  # If marking breakpoint genes
  if ($OPT{breakpoint_genes}) {
    &MarkBreakpointGenes($next_input_dir);
    $stage += 1;
    $next_input_dir = catdir($OPT{results_dir}, $stage."_breakpoint_genes");
  }

  $| = 0; # REMINDER: turn off autoflush
  return undef;
}

sub GenerateRawCandidateGroups {
  # No parameters
  # If using cluster results, check that the cluster jobs are all complete
  if ($OPT{cluster} && $OPT{add_support}) {
    &CheckClusterJobStatus();
  }

  # Get the BARNACLE event groups
  my ($groups, $contigs_in_events, $fasta_data) = &GetEventGroups();
  #&PrintGroupCtgIDs($groups); # DEBUG

  if ($OPT{identify_candidates}) {
    &FinishCandidateIdentification();
    return "stop";
  }

  # if the flag_RNA option is set parse the RNA repeats file and
  # integrate into the data structure
  #unless ($OPT{no_flag_RNA}) {
  #  &Parse_RNA_Repeats($groups);
  #}
  #&PrintGroupCtgIDs($groups); # DEBUG

  # If locally calculating pair-to-genome support
  unless ($OPT{no_pair_to_genome} || $OPT{cluster_pair2gen}) {
    # in AddPairSupportToGroups function
    my $next_input_dir = "NOT_USED";
    &AddPairSupportToGroups($groups, $next_input_dir);
  }

  #Summarize the events
  &Summary($groups, $contigs_in_events, $fasta_data);

  #Generate the fasta file with candidate contigs
  &Fasta($fasta_data);

  # If calculating pair-to-genome support on a cluster
  unless ($OPT{no_pair_to_genome} || !$OPT{cluster_pair2gen}) {
    # for now, we do not use "next_input_dir/previous_output_dir"
    # in AddPairSupportToGroups function
    my $next_input_dir = "NOT_USED";
    &AddPairSupportToGroups($groups, $next_input_dir);
    #unless ($OPT{cluster_pair2gen}) {
    #  $stage += 1;
    #  $next_input_dir = catdir($results_dir, $stage."_with_p2g");
    #}
  }
  #&PrintGroupCtgIDs($groups); # DEBUG

  return "continue";
}

sub CheckClusterJobStatus {
  # No parameters
  my $complete = 1;
  my $cluster_dir = &GetClusterDir();
  if ($OPT{debug}) {
    system("echo \$BARNACLE_PATH");
  }
  &debug("CLUSTER: ".$cluster_dir);
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $check_status_cmd = join(" ",
    $python, $EXEC{check_status}, $cluster_dir);
  if ($OPT{disable_profiling_timer}) {
    $check_status_cmd .= " --disable-profiling-timer";
  }
  if ($OPT{debug_high}) {
    $check_status_cmd .= " -d";
  }
  &status($check_status_cmd);
  my $status = `$check_status_cmd`;
  &status($status);
  chomp $status;
  # use "m" at the end of the match to use multi-line matching
  my ($num_complete) = $status =~ /^complete:\s*(\d+)$/m;
  my ($num_jobs)     = $status =~ /^total:\s*(\d+)$/m;
  if (! $num_jobs or (! $num_complete and 0 != $num_complete)) {
    &error("could not get job status from \"".$status."\"");
  }
  &debug("JOBS: $num_jobs, COMPLETE: $num_complete");
  if ($num_jobs > $num_complete) {
    my $num_incomplete = $num_jobs - $num_complete;
    &status($num_incomplete." cluster job(s) not complete.\n");
    $complete = 0;
  }
  unless ($complete) {
    &error("Run the following command to check cluster status:\n".
      $check_status_cmd);
  }
}

sub GetEventGroups {
  # No parameters
  # Load the candidate contigs into memory
  my ($ungrouped_members, $contigs_in_events, $fasta_files) =
    &IdentifyOrLoadCandidates();
  if ($OPT{identify_candidates}) {
    return undef;
  }
  my $num_contigs = keys %{$contigs_in_events};
  if (0 == $num_contigs) {
    &error("No candidate contigs found");
  }
  #Group the events
  my $groups =
    &GroupEvents($ungrouped_members, $contigs_in_events);
  return ($groups, $contigs_in_events, $fasta_files);
}

sub IdentifyOrLoadCandidates {
  # No parameters
  # Setup the input directories
  my $input_paths = &SetupInputPaths();

  # Get the alignment files and populate the data structure
  &Find_Files($input_paths);

  # Run the candidate contig identification code
  my %num_gaps_found = ();
  for my $contig_type (keys %{$input_paths}) {
    $num_gaps_found{$contig_type} = 0;
  }
  if ($OPT{identify_candidates}) {
    my $cluster_jobs_file =
      &IdentifyCandidates($input_paths, \%num_gaps_found);
    if ($OPT{cluster}) {
      close $cluster_jobs_file->{handle};
      &SubmitJobsToCluster(
        $cluster_jobs_file->{dir},
        $cluster_jobs_file->{path},
        "candidate-identification",
        "cid", $OPT{cid_memory},
      );
    }
    return undef;
    #} elsif ($OPT{cluster}) {
  } else {
    &FindClusterResults($input_paths, \%num_gaps_found);
  }
  return &ParseAllCandidateFiles($input_paths, \%num_gaps_found);
}

sub ParseAllCandidateFiles {
  # get function parameters
  my ($input_paths, $num_gaps_found) = @_;
  # open a file for concatenating the psl alignment files
  my $psl_file = &OpenEventAlignmentsFile();
  # Create an empty dictionary to hold all the candidate contig output files
  # keyed first by chromosome, then genomic junction coordinate
  my %ungrouped_members = ();
  # ungrouped_members{$chr_from}{$chr_to}{$junction_from}{$junction_to}
  #                  {$abs_member_id} = $event_grouping_info
  # Create an empty dictionary for holding fasta file names
  my %fasta_file = ();
  # Create an empty dictionary for holding contigs to reads data
  my %contigs_in_events = ();
  # Every candidate contig/aligner type pair should have a unique id
  my $abs_member_id = 0;
  #Parse the candidate contig output files
  &status("Parsing candidate identification output files...");
  my $parse_candidate_start = new Benchmark;
  for my $contig_type (keys %{$input_paths}) {
    &debug("  Parsing candidate ".&Desc($contig_type)."...");
    #&print_full($logref, "Dir: ".$contig_type."\n"); # DEBUG
    #&print_full($logref,
    #  "Adding fasta file: $input_paths{$contig_type}{contig_seq_file}\n");
    $fasta_file{$contig_type}{seq_file} =
      $input_paths->{$contig_type}{contig_seq_file};
    my %aligner_paths = %{$input_paths->{$contig_type}{aligners}};
    for my $aligner (sort keys %aligner_paths) {
      &debug("    Parsing candidate contigs found from ".
             $aligner." alignments...");
      # check that candidate contigs were found in the current alignment group
      if (exists $aligner_paths{$aligner}{candidates}) {
        for my $candidate_path (@{$aligner_paths{$aligner}{candidates}}) {
          &debug("      Candidate $aligner file: \"$candidate_path\"");
          if (-e $candidate_path) {
            &ParseCandidateContigs
              ($candidate_path, $aligner, $contig_type, \$abs_member_id,
               \%ungrouped_members, \%contigs_in_events);
          } else {
            if ("gap" eq $aligner) {
              my $num_gaps = $num_gaps_found->{$contig_type};
              &warn("gap candidates file not found: $candidate_path ".
                    "($num_gaps gap candidate contigs reported).");
            } else {
              my $counts_path = $candidate_path;
              $counts_path =~ s/(gap|split)\.candidates/counts/;
              if (! -f $counts_path) {
                &error("Candidate identification counts file not found: ".
                  $counts_path);
              }
              my $split_line = `grep -m1 "Split:" $counts_path`;
              my ($num_split) = $split_line =~ /Split: ([0-9]+)$/;
              if (0 == $num_split) {
                &warn("No split candidates found for ".$candidate_path);
              } else {
                my $align_check_log = $log_file_name;
                $align_check_log =~ s/\.continue$//;
                system("grep \"Number of split candidates found\" ".
                       $align_check_log);
                &error($aligner." candidate contigs file not found: ".
                    "\"".$candidate_path."\". Perhaps -identify_candidates ".
                    "did not find any events? Check log file.");
              }
            }
          }
          my $small_psl_path = $candidate_path.".psl";
          if (-e $small_psl_path) {
            &AddAlignmentsToFile($small_psl_path, $psl_file);
          }
        }
      } elsif ($OPT{debug}) {
        &status($contig_type." ".$aligner.
          " candidate_paths not in input_paths hash");
      }
    }
  }
  if ($OPT{cluster} && $OPT{add_support}) {
    &status("");
  }
  close($psl_file);
  &status($abs_member_id." event members found");
  &TimeSpent("parsing identify_candidates output", $parse_candidate_start, 1);
  &debug(""); # a blank newline
  return (\%ungrouped_members, \%contigs_in_events, \%fasta_file);
}

sub OpenEventAlignmentsFile {
  # No parameters
  my $psl_file_path = &GetOutputPath($OPT{outdir}, "psl");
  &debug("Creating all event alignments file: ".$psl_file_path);
  open(my $psl_file, ">", $psl_file_path) ||
    &error("Can't open all event alignments file ".$psl_file_path." ".$!);
  return $psl_file;
}

sub AddAlignmentsToFile {
  # get function parameters
  my ($small_psl_path, $psl_file) = @_;
  open(my $small_psl_file, $small_psl_path) ||
    &error("Can't open intermediate event alignments file ".
      $small_psl_path." ".$!);
  while (<$small_psl_file>) {
    print $psl_file $_;
  }
  close($small_psl_file);
  return undef;
}

sub SetupInputPaths {
  # No parameters
  # Create an empty dictionary to hold the input directory names
  my %input_paths = ();
  for my $contig_type (sort keys %contig_types) {
    &status("Will examine ".&Desc($contig_type));
    %{$input_paths{$contig_type}} = ();
    for my $aligner (sort keys %aligner_types) {
      %{$input_paths{$contig_type}{aligners}{$aligner}} = ();
      if ("split" eq $aligner && $OPT{include_gap_events}) {
        %{$input_paths{$contig_type}{aligners}{gap}} = ();
      }
    }
  }
  return \%input_paths;
}

#Find the alignment files
sub Find_Files {
  # get function parameters
  my ($input_paths) = @_;
  # keep a count of the number of alignment files found
  my %num_align_files = ();
  my $total_num_align_files = 0;
  for my $contig_type (sort keys %{$input_paths}) {
    &debug("Finding $contig_type files");
    my $ctg_set = $contig_type;
    if ($contig_type eq 'adj' || $contig_type eq 'main') {
      $ctg_set = "merge.$contig_type";
    }
    my %aligner_paths = %{$input_paths->{$contig_type}{aligners}};
    for my $aligner (sort keys %aligner_types) {
      #if ((("gap"   ne $aligner && ! $OPT{no_split_check}) ||
      #     ("split" eq $aligner && ! $OPT{no_gap_check})) &&
      #    ($OPT{identify_candidates} || $OPT{cluster})) {
      if (("gap"   ne $aligner && ! $OPT{no_split_check}) ||
           ("split" eq $aligner && ! $OPT{no_gap_check})) {
        # Get the appropriate extension
        #my $ext = $aligner_types{$aligner}{align_ext};
        # Find the alignment files
        #&FindAlignmentFiles(\%aligner_paths, $contig_type, $aligner, $ext);
        &FindAlignmentFiles(\%aligner_paths, $contig_type, $aligner);
        $num_align_files{$contig_type}{$aligner} =
          @{$aligner_paths{$aligner}{aligns}};
        $total_num_align_files += $num_align_files{$contig_type}{$aligner};
      }
      # Get the candidate contig output paths
      #&GetCandidateContigPaths(\%aligner_paths, $aligner, $ctg_set);
    }
    #Find the source contigs
    $input_paths->{$contig_type}{contig_seq_file} = $CFG{sequences}{contigs};
  }
  &status(""); # a blank newline

  if ($OPT{debug}) {
    &DisplayInputFilesFound($input_paths, \%num_align_files);
  }
  if (0 == $total_num_align_files &&
      ($OPT{identify_candidates} || $OPT{cluster})) {
    &error("Cannot find alignment files.");
  }
  return undef;
}

sub DisplayInputFilesFound {
  # get function parameters
  my ($input_paths, $num_align_files) = @_;
  &status("Input Files:");
  for my $contig_type (keys %{$input_paths}) {
    &status("  Contig type: ".$contig_type);
    my $contig_seq_file = $input_paths->{$contig_type}{contig_seq_file};
    &status("    Contigs file: ".$contig_seq_file);
    my %aligner_paths = %{$input_paths->{$contig_type}{aligners}};
    for my $aligner (sort keys %aligner_paths) {
      &status("    Aligner: ".$aligner);
      if (exists $aligner_paths{$aligner}{aligns}) {
        if (0 == $num_align_files->{$contig_type}{$aligner}) {
          &status("      No $aligner alignment files found for ".
                  &Desc($contig_type));
        } elsif ($OPT{debug_high}) {
          for my $file (
            sort { my ($a_num) = $a =~ /seq\.(\d+)\.[a-z]+/;
                   my ($b_num) = $b =~ /seq\.(\d+)\.[a-z]+/;
                   int($a_num) <=> int($b_num)
                 } @{$aligner_paths{$aligner}{aligns}}) {
            &status("      Align file: ".$file);
          }
        } else {
          &status("      Num align files: ".
                  $num_align_files->{$contig_type}{$aligner});
        }
      }
      for my $candidate_path (@{$aligner_paths{$aligner}{candidates}}) {
        &status("      Candidate contigs file: ".$candidate_path);
      }
    }
  }
  &status(""); # a blank newline
}

sub FindAlignmentFiles {
  # get function parameters
  #my ($aligner_paths, $contig_type, $aligner, $ext) = @_;
  my ($aligner_paths, $contig_type, $aligner) = @_;
  &status("  Looking for $contig_type $aligner alignment files...");
  my @align_files = ();
  my $contig_group = $OPT{lib};
  unless ("contigs" eq $contig_type) {
    $contig_group .= substr($contig_type,0,1);
  }
  # first try using aligner for align group
  my $align_group = $contig_group."-".substr($aligner,0,4);
  #my $align_output_dir = join("/",
  #  $merge_dir, "cluster", $align_group, "output");
  #my $merge_dir = "";
  #my $align_output_dir =
  #  catdir($merge_dir, "cluster", $align_group, "output");
  #&debug("  Looking in $align_output_dir");
  #if (! -d $align_output_dir) {
  #  &debug("  ".$align_output_dir." is not a directory.");
  #  &status("  Using contig type to build alignment output directory path.");
  #  $align_group = $OPT{lib}."-".$contig_type;
  #  $align_output_dir =
  #    catdir($merge_dir, "cluster", $align_group, "output");
  #  &debug("  Looking in $align_output_dir");
  #}
  #if (-d $align_output_dir) {
  &debug("  Looking in ".$OPT{c2g_dir});
  if (-d $OPT{c2g_dir}) {
    #opendir(ALIGN_OUTPUT_DIR, $align_output_dir) ||
    #  &error("Can't open ".$contig_type." dir ".$align_output_dir." ".$!);
    #my @align_file_names =
    #  grep {/seq\.\d+\.$ext$/} grep {$_ !~ /^\./} readdir ALIGN_OUTPUT_DIR;
    #closedir(ALIGN_OUTPUT_DIR);
    my @raw_align_files = `ls $CFG{alignments}{c2g} 2>/dev/null`;
    my $num_align_files = @raw_align_files;
    if (0 == $num_align_files) {
      &error("Could not find any contig-to-genome alignment files: ".
        $CFG{alignments}{c2g});
    }
    for my $align_file (sort @raw_align_files) {
      #push @align_files, catfile($align_output_dir, $align_file_name);
      chomp($align_file);
      push @align_files, $align_file;
    }
    &debug("  Example c2g file: \"".$align_files[0]."\"");
  } else {
    #&error("Can't open ".$contig_type." dir ".
    #  $align_output_dir.": not a directory");
    &error("Contig-to-genome alignments directory does not exist: ".
      $OPT{c2g_dir});
  }
  # if no files were found yet
  if (0 == @align_files) {
    #&warn("could not find $contig_type $aligner alignments");
    &warn("could not find any contig-to-genome alignments");
  }
  if (@align_files) {
    my $num_align_files = @align_files;
    #&status("  Found $num_align_files $contig_type $aligner alignment files");
    &status("  Found ".$num_align_files." contig-to-genome alignment files");
    $aligner_paths->{$aligner}{aligns} = \@align_files;
  } else {
    #&status("  No $contig_type $aligner alignments found");
    &status("  No contig-to-genome alignments found");
  }
}

#sub GetCandidateContigPaths {
#  # get function parameters
#  my ($aligner_paths, $aligner, $ctg_set) = @_;
#  my $local_dir = catdir($OPT{results_dir}, "local_cid");
#  my $split_candidate_path = &GetOutputPath($local_dir, $ctg_set.".split");
#  #my $split_candidate_path = &GetOutputPath($OPT{results_dir}, $ctg_set.".split");
#  #$ctg_set.".".substr($aligner,0,4));
#  &debug("Adding $aligner split candidate path: $split_candidate_path");
#  $aligner_paths->{$aligner}{candidates} = ();
#  #if (! $OPT{cluster}) {
#  #  push @{$aligner_paths->{$aligner}{candidates}}, $split_candidate_path;
#  #}
#  if ("split" eq $aligner && $OPT{include_gap_events}) {
#    my $gap_candidate_path = &GetOutputPath($OPT{results_dir}, $ctg_set.".gap");
#    &debug("Adding gap candidate path: $gap_candidate_path");
#    $aligner_paths->{gap}{candidates} = ();
#    if (! $OPT{cluster}) {
#      push @{$aligner_paths->{gap}{candidates}}, $gap_candidate_path;
#    }
#  }
#  return undef;
#}

# Get the file names used for the gap filter
sub Get_Ctg_Seq_File {
  # get function parameters
  my ($align_out) = @_;
  # DEBUG # &print_full($logref, "align: ".$align_out."\n");
  my $ctg_seq_file = $align_out;
  $ctg_seq_file =~ s/output/input/;
  $ctg_seq_file =~ s/\.[^.]+$/.fa/;
  unless (-f $ctg_seq_file) {
    &error("Cannot find contig sequence file: ".$ctg_seq_file);
  }
  return $ctg_seq_file;
}

sub GetClusterOutputPath {
  # get function parameters
  my ($cluster_dir, $align_out_path, $align_ext, $aligner, $contig_type) = @_;
  &debug_high("Getting file name from: ".$align_out_path." Ext: ".$align_ext);
  my ($align_file_name,undef,undef) = fileparse($align_out_path, ($align_ext));
  &debug_high("  File name: ".$align_file_name);
  my ($align_file_num) = $align_file_name =~ /\.(\d+)\.$/;
  &debug_high("  File num: ".$align_file_num);
  &error("cannot get alignment output file number from: ".
    $align_out_path) unless defined $align_file_num;
  my $job_dir = catdir($cluster_dir, "job_".$align_file_num);
  # ensure that job directory exists
  make_path($job_dir) if (!-d $job_dir);
  #my $split_ext = "split.candidate.".substr($aligner,0,4);
  #my $gap_ext = "gap.candidate";
  my $split_ext = "split.candidates";
  my $gap_ext = "gap.candidates";
  unless ("contigs" eq $contig_type) {
    $split_ext = $contig_type.".".$split_ext;
    $gap_ext = $contig_type.".".$gap_ext;
  }
  my $split_candidate_path = catfile($job_dir, $align_file_name).$split_ext;
  my $gap_candidate_path = catfile($job_dir, $align_file_name).$gap_ext;
  return ($split_candidate_path, $gap_candidate_path);
}

# Get the arguments needed to use the gap realignment tool
sub GapRealignmentArgs {
  # get function parameters
  my ($id_candidates_args, $ctg_seq_file, $gap_candidate_path,
    $gap_tools_config) = @_;
  my $gap_realignment_args = join(" ",
    $id_candidates_args,
    "--gap-candidates",
    "--ctg-file",         $ctg_seq_file,
    "--gap-realigner",    $EXEC{gap_realigner},
    "--gap-min-size",     $OPT{min_gap_size},
    "--gap-min-identity", $OPT{min_gap_identity},
    "--gap-min-fraction", $OPT{min_gap_fraction},
    "--gap-max-len",      $OPT{gap_max_len},
    "--gap-config",       $gap_tools_config,
  );
  return $gap_realignment_args;
}

sub FindClusterResults {
  # get function parameters
  my ($input_paths, $num_gaps_found) = @_;
  my $results_dir = &GetLocalDir();
  if ($OPT{cluster}) {
    $results_dir = &GetClusterDir();
  }
  my $total_num_results_files = 0;
  for my $contig_type (sort keys %{$input_paths}) {
    &debug("Finding cluster results for ".&Desc($contig_type)."...");
    my $aligner_paths = $input_paths->{$contig_type}{aligners};
    for my $aligner (sort keys %aligner_types) {
      my $align_ext = $aligner_types{$aligner}{align_ext};
      &debug("  Aligner: ".$aligner." Ext: ".$align_ext);
      for my $align_out (
        sort { my ($a_num) = $a =~ /seq\.(\d+)\.$align_ext/;
               my ($b_num) = $b =~ /seq\.(\d+)\.$align_ext/;
               int($a_num) <=> int($b_num)
             } @{$aligner_paths->{$aligner}{aligns}}) {
        &debug("    AlignOut: ".$align_out);
        my ($split_candidate_path, $gap_candidate_path) =
          &GetClusterOutputPath($results_dir, $align_out,
          $align_ext, $aligner, $contig_type);
        &debug("    SPLIT: ".$split_candidate_path);
        &debug("    GAP:   ".$gap_candidate_path);
        if (not $OPT{no_split_check}) {
          push @{$aligner_paths->{$aligner}{candidates}},
            $split_candidate_path;
        }
        if ("split" eq $aligner && $OPT{include_gap_events}) {
          push @{$aligner_paths->{gap}{candidates}}, $gap_candidate_path;
        }
      }
      my $num_results_files = @{$aligner_paths->{$aligner}{candidates}};
      $total_num_results_files += $num_results_files;
      &status("Number of ".$aligner." candidate results files: ".
        $num_results_files."\n");
    }
  }
  &status("Total number of candidate results files: ".
    $total_num_results_files."\n");
  &status(""); # a blank newline
  return undef;
}

# Run the code to find candidate contigs on the alignment files
sub IdentifyCandidates {
  # get function parameters
  my ($input_paths, $num_gaps_found) = @_;
  &status("Searching for candidate contigs...");
  # Prepare for using the cluster
  my $cluster_jobs_file = &PrepareForCluster();
  # Create gap realignment config file, if checking gapped contigs
  my $gap_tools_config;
  unless ($OPT{no_gap_check}) {
    $gap_tools_config = &CreateGapToolsConfig();
  }
  # populate an @ignore_cmds array outside of the "for contig_type" loop
  my @ignore_cmds = ();
  #my $completed_cmds = "$outdir/ignore_cmds.txt";
  my $completed_cmds = catfile($OPT{outdir}, "ignore_cmds.txt");
  if (defined $OPT{ignore_cmds}) {
    &GetIgnoreCommandsList(\@ignore_cmds);
    $completed_cmds = $OPT{ignore_cmds};
  } elsif (-e $completed_cmds) {
    unlink($completed_cmds);
  }
  my $id_candidates_args = &CandidateIdentificationArgs();
  for my $contig_type (sort keys %{$input_paths}) {
    if ($OPT{identify_candidates}) {
      &status("Examining alignments of ".&Desc($contig_type));
    }
    my @commands = ();
    for my $aligner (sort keys %aligner_types) {
      &debug_high("Alignment file extension: ".
        $aligner_types{$aligner}{align_ext});
      my $new_commands = &CandidateIdentificationCommands
        ($input_paths->{$contig_type}{aligners}, $aligner,
         $aligner_types{$aligner}{align_ext}, $id_candidates_args,
         $cluster_jobs_file->{dir}, $contig_type, $gap_tools_config);
      push @commands, @{$new_commands};
    }
    if (defined $OPT{ignore_cmds}) {
      &IgnoreCommands(\@ignore_cmds, \@commands);
    }
    if ($OPT{cluster}) {
      &AddJobsToFile($cluster_jobs_file->{handle}, \@commands);
    } elsif ($OPT{identify_candidates}) {
      &RunCandidateIdentificationCommands
        ($completed_cmds, \@commands, $num_gaps_found, $contig_type);
    } else {
      &status("Candidate ".&Desc($contig_type)." files found");
      #&print_full($logref, join("\n\n",@commands));
      #&print_full($logref,"\n\n");
    }
  }
  if ($OPT{cluster}) {
    # close the cluster jobs file
    close $cluster_jobs_file->{handle};
  }
  #for my $contig_type (keys %files) {
  #&print_full($logref, "Num Gaps: $num_gaps_found->{$contig_type}\n");
  #}
  #if ($OPT{debug}) {
  #  &print_full($logref, "Files:\n");
  #  for my $contig_type (keys %input_paths) {
  #    &print_full($logref, "  Contig type: ".$contig_type."\n");
  #    if (exists $input_paths{$contig_type}) {
  #      my %aligners = %{$input_paths{$contig_type}};
  #      for my $aligner (grep {/candidate/} keys %aligners) {
  #        &print_full($logref, "    $aligner - $aligners{$aligner}\n");
  #      }
  #    } else {
  #      $msg = "No candidate files found for $contig_type contigs";
  #      &print_full($logref, "    $msg\n");
  #    }
  #  }
  #  &print_full($logref, "\n");
  #}
  return $cluster_jobs_file;
}

sub PrepareForCluster {
  # no parameters
  my %cluster_jobs_file;
  if ($OPT{cluster}) {
    # create the cluster directory
    $cluster_jobs_file{dir} = &GetClusterDir();
    make_path($cluster_jobs_file{dir});
    # construct the path to the jobs file
    $cluster_jobs_file{path} = catfile($cluster_jobs_file{dir}, "jobs");
    # open the jobs file for writing
    &debug("Creating cluster jobs file: ".$cluster_jobs_file{path});
    open($cluster_jobs_file{handle}, ">", $cluster_jobs_file{path}) ||
      &error("Can't open cluster jobs file ".$cluster_jobs_file{path}." ".$!);
  } else {
    $cluster_jobs_file{dir} = &GetLocalDir();
    make_path($cluster_jobs_file{dir});
  }
  return \%cluster_jobs_file;
}

sub GetClusterDir {
  return catdir($OPT{results_dir}, "cluster_cid");
}

sub GetLocalDir {
  return catdir($OPT{results_dir}, "local_cid");
}

sub CreateGapToolsConfig {
  # no function parameters
  my $gap_tools_config = &GetOutputPath($OPT{outdir}, "gap.cfg");
  # open the gap tools config file
  open(GAP_CFG, ">", $gap_tools_config) ||
    &error("Can't open gap tools config file ".$gap_tools_config." ".$!);
  # write the gap tools config file
  print GAP_CFG "BLATpath=".$CFG{commands}{blat}."\n";
  print GAP_CFG "BLAToptions=".$CFG{gap_realigner}{BLAToptions}."\n";
  print GAP_CFG "2bitGenPath=".$CFG{sequences}{genome_2bit}."\n";
  print GAP_CFG "relaxedVersion=1\n";
  print GAP_CFG "fullContigAlign=0\n";
  print GAP_CFG "dupAlignGapPenalty=10\n";
  #print GAP_CFG "NoGapSplitAlign=1\n";
  print GAP_CFG "BLATsource=exec\n";
  # close the gap tools config file
  close GAP_CFG;
  # return the file name
  return $gap_tools_config;
}

sub GetIgnoreCommandsList {
  # get function parameters
  my ($ignore_cmds) = @_;
  &status("Populating ignore commands array...");
  open(IGNORE, $OPT{ignore_cmds}) ||
    &error("Can't open ignore cmds file ".$OPT{ignore_cmds}." ".$!);
  # get a line from the file
  while (<IGNORE>) {
    chomp;
    my $ignore_cmd = $_;
    push @{$ignore_cmds}, $ignore_cmd;
  }
  close IGNORE;
  my $num_ignore_cmds = @{$ignore_cmds}; # DEBUG
  &debug("$num_ignore_cmds commands will be ignored");
  return undef;
}

sub CandidateIdentificationArgs {
  # no parameters
  my @arg_list = ();
  if ($OPT{no_split_check}) {
    push @arg_list, "--no-split-candidates";
  } else {
    push @arg_list, "--split-candidates";
  }
  push @arg_list, "--num-aligns", $OPT{num_aligns},
    "--min-identity", $OPT{min_identity};
  if ($OPT{add_gene_annotation}) {
    push @arg_list, "--genes", "--gene-coords",
      $CFG{coords}{gene_features};
  } else {
    push @arg_list, "--no-gene-info";
  }
  push @arg_list, "--single-align", $OPT{max_match_fraction},
    "--merge-overlap", $OPT{min_merge_overlap};
  if ($OPT{use_quick_chooser}) {
    push @arg_list, "--quick-chooser";
  } else {
    push @arg_list, "--smart-chooser";
  }
  push @arg_list, "--maintain-pared-groups";
  unless ($OPT{no_exon_preference}) {
    push @arg_list, "--prefer-exons";
  }
  push @arg_list, "--ctg-rep", $OPT{min_ctg_represented};
  unless ($OPT{mito_prefilter_off}) {
    push @arg_list, "--no-mito";
  }
  if (! $OPT{cluster}) {
    push @arg_list, "--log-file", $log_file_name;
  }
  if ($OPT{debug_high} || $OPT{cluster}) {
    push @arg_list, "--debug";
  }
  if ($OPT{disable_profiling_timer}) {
    push @arg_list, "--disable-profiling-timer";
  }
  my $id_candidates_args_standard = join(" ", @arg_list);
  my %id_candidates_args = (
    standard => $id_candidates_args_standard,
    append => $id_candidates_args_standard." --append",
  );
  return \%id_candidates_args;
}

sub CandidateIdentificationCommands {
  # get function parameters
  my ($aligner_paths, $aligner, $align_ext, $id_candidates_args,
      $cluster_dir, $contig_type, $gap_tools_config) = @_;
  my @commands = ();
  if (exists $aligner_paths->{$aligner}{aligns}) {
    #my $split_candidate_path = $aligner_paths->{$aligner}{candidates}[0];
    #if (! defined $split_candidate_path) {
    #  &error("Split candidate path not defined!");
    #}
    #my $gap_candidate_path = "";
    #if (exists $aligner_paths->{gap}{candidates}) {
    #  $gap_candidate_path = $aligner_paths->{gap}{candidates}[0];
    #  if (! defined $gap_candidate_path) {
    #    &error("Gap candidate path not defined!");
    #  }
    #}
    my $id_candidates_args_to_use = $id_candidates_args->{standard};
    if (1 < @{$aligner_paths->{$aligner}{aligns}} &&
        ! $OPT{cluster}) {
      #if (! defined $OPT{ignore_cmds}) {
      #  if (! $OPT{no_split_check} && -e $split_candidate_path) {
      #    push @commands, "rm ".$split_candidate_path;
      #  }
      #  if (! $OPT{no_gap_check} && -e $gap_candidate_path &&
      #      "split" eq $aligner) {
      #    push @commands, "rm ".$gap_candidate_path;
      #  }
      #}
      $id_candidates_args_to_use = $id_candidates_args->{append};
    }
    unless ("contigs" eq $contig_type) {
      $id_candidates_args_to_use .= " --contig-set ".$contig_type;
    }
    my $id_candidates_args_no_gap = $id_candidates_args_to_use;
    $id_candidates_args_to_use .= " --no-gap-candidates";
    for my $align_out (
      sort { my ($a_num) = $a =~ /seq\.(\d+)\.$align_ext/;
             my ($b_num) = $b =~ /seq\.(\d+)\.$align_ext/;
             int($a_num) <=> int($b_num)
           } @{$aligner_paths->{$aligner}{aligns}}) {
      #if ($OPT{cluster}) {
      #  ($split_candidate_path, $gap_candidate_path) = &GetClusterOutputPath
      #    ($cluster_dir, $align_out, $align_ext, $aligner, $contig_type);
      #}
      my ($split_candidate_path, $gap_candidate_path) = &GetClusterOutputPath
        ($cluster_dir, $align_out, $align_ext, $aligner, $contig_type);
      #&print_full($logref, "Adding $aligner file $align_out\n"); # DEBUG
      if ("split" eq $aligner && ! $OPT{no_gap_check}) {
        my $ctg_seq_file = &Get_Ctg_Seq_File($align_out);
        $id_candidates_args_to_use = &GapRealignmentArgs
          ($id_candidates_args_no_gap, $ctg_seq_file, $gap_candidate_path,
           $gap_tools_config);
      }
      #$id_candidates_args_to_use =~ s/CONTIG_SET/$contig_type/;
      my (undef, $candidate_dir) = splitpath($split_candidate_path);
      my $python = $CFG{commands}{python};
      if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
        $python = $CFG{commands}{python_cluster};
      }
      push @commands, join(" ",
        $python, $EXEC{identify_candidates},
        $align_out, $candidate_dir, $id_candidates_args_to_use,);
    }
  }
  return \@commands;
}

sub IgnoreCommands {
  # get function parameters
  my ($ignore_cmds, $commands) = @_;
  &status("Ignoring commands...");
  # make a copy of the commands list and clear the original at the same time
  my @all_commands = splice(@{$commands});
  # iterate through the copy
  for my $command (@all_commands) {
    &debug("CMD: $command");
    # assume the command should not be ignored
    my $ignore_this_cmd = 0;
    # iterate through the list of commands to ignore
    my $ignore_index = 0;
    for my $ignore_cmd (@{$ignore_cmds}) {
      my ($ignore_dir, $ignore_file) = split(/ /, $ignore_cmd);
      &debug("  IGNORE: $ignore_cmd");
      &debug("     DIR: $ignore_dir FILE: $ignore_file");
      if ($ignore_cmd eq $command or
          $command =~ /${ignore_dir}.*${ignore_file}/) {
        &debug("-- ignore it --");
        $ignore_this_cmd = 1;
        last;
      }
      $ignore_index++;
    }
    # if the command should be ignored
    if ($ignore_this_cmd) {
      # remove the ignore_cmd from the list
      splice(@{$ignore_cmds}, $ignore_index, 1);
    } else {
      # add the command back to the list of commands
      push @{$commands}, $command;
    }
  }
}

sub AddJobsToFile {
  # get function parameters
  my ($cluster_jobs_file, $commands) = @_;
  for my $command (@{$commands}) {
    &debug("Add job to list: ".$command);
    print $cluster_jobs_file &JobLine($command);
  }
  return undef;
}

sub RunCandidateIdentificationCommands {
  # get function parameters
  my ($completed_cmds, $commands, $num_gaps_found, $contig_type) = @_;
  my $total_id_candidates_start = new Benchmark;
  open(my $completed_cmds_file, ">>", $completed_cmds) ||
    &error("Can't open completed commands file ".$completed_cmds." ".$!);
  for my $command (@{$commands}) {
    # run the command
    &RunCandidateIdentificationCommand
      ($command, $num_gaps_found, $contig_type, $completed_cmds_file);
  }
  close $completed_cmds_file;
  &TimeSpent("identifying candidate contigs", $total_id_candidates_start, 1);
  return undef;
}

sub RunCandidateIdentificationCommand {
  # get function parameters
  my ($command, $num_gaps_found, $contig_type, $completed_cmds_file) = @_;
  #&debug("RUNNING CANDIDATE IDENTIFICATION COMMAND ".$command);
  my ($align_group, $align_file_name) = '';
  #&debug("INITIALIZING VARS...");
  if ($command =~ /$EXEC{identify_candidates}/) {
    # extract file name from command
    &debug("Extracting file name from command...");
    my @cmd_fields = split(/ /, $command);
    #my $num_cmd_fields = @cmd_fields;
    #&debug("# command fields: ".$num_cmd_fields);
    my $AFC_ARG_ALIGN_FILE = 2;
    &debug_high("Lib: ".$OPT{lib});
    &debug_high("File: ".$cmd_fields[$AFC_ARG_ALIGN_FILE]);
    ($align_group, $align_file_name) =
      $cmd_fields[$AFC_ARG_ALIGN_FILE] =~
      /\/($OPT{lib}\w?)-\w+\/.*\/([^\/]*)$/;
    if (!$align_group) {
      &warn("Could not get align group from \"".
        $cmd_fields[$AFC_ARG_ALIGN_FILE]."\"");
    }
    if ($align_file_name) {
      &status("Searching alignment file: ".$align_group."-".$align_file_name);
    } else {
      &error("Could not get file name from \"".
        $cmd_fields[$AFC_ARG_ALIGN_FILE]."\"");
    }
    &debug($command);
  } else {
    &status($command);
  }
  my $identify_candidates_start = new Benchmark;
  # temporarily close log file so that identify_candidates can use it
  close $logref;
  my $result = system($command);
  # reopen the log file
  open($logref,">>",$log_file_name) ||
    &error("Can't open log file $log_file_name $!");
  #my $output = `$command 2>&1`;
  #chomp $output;
  # need to shift by eight to get the actual return value
  my $shifted = $? >> 8;
  #my $result = $? >> 8;
  #&status($output);
  #&debug("RESULT: ".$result."\n"); # DEBUG
  #my $unshifted = $result << 8;
  #&debug("Unshifted: ".$unshifted."\n");
  #&debug("Success: $CI_SUCCESS ");
  #&debug("No Aligns: ".$CI_NO_ALIGNMENTS."\n");
  #&debug("Command: ".$command."\n");
  if ($CI_SUCCESS != $shifted) {
  #if ($CI_SUCCESS != $result) {}
    # check that command was "identify_candidate_contigs"
    if ($CI_NO_ALIGNMENTS == $shifted &&
    #if ($CI_NO_ALIGNMENTS == $result &&
        $command =~ /identify_candidate_contigs/) {
      &warn("no alignments were found in ".$align_file_name.".");
    } else {
      &error("while running ".$command.": ".$shifted);
    }
  }
  if ($command =~ /identify_candidate_contigs/ &&
      $command =~ /psl/ && $command =~ /gap/) {
    my @command_fields = split(/ /, $command);
    my $ac_input_path = $command_fields[2];
    my (undef, undef, $ac_file_name) = splitpath($ac_input_path);
    #my $new_tail = "counts.blat";
    my $new_tail = "counts";
    unless ("contigs" eq $contig_type) {
      $new_tail = $contig_type.".".$new_tail;
    }
    $ac_file_name =~ s/psl/$new_tail/;
    my $ac_out_dir = $command_fields[3];
    if ($OPT{cluster}) {
      $ac_out_dir = catdir($ac_out_dir, "cluster");
    }
    #my $count_file = $command_fields[3].".counts";
    my $count_file = catfile($ac_out_dir, $ac_file_name);
    &debug("COUNTS: ".$count_file);
    my $output = `tail -n2 $count_file`;
    my ($num_new_gaps) = $output =~ /Gapped:.*(\d+)/;
    if ($num_new_gaps !~ /./) {
      &error("could not get number of gaps from: \"".$output."\"");
    }
    #  $output =~ /Number of gapped alignments found: (\d+)/;
    &debug("New gaps found: \"$num_new_gaps\"\n");
    $num_gaps_found->{$contig_type} += $num_new_gaps;
  }
  unless ($command =~ /^rm/) {
    &debug("Finished searching!");
    &TimeSpent("processing contig alignment file",
               $identify_candidates_start);
    #my $cmd_dir = $OPT{lib}.substr($contig_type,0,1);
    # Add command identifier to completed commands list
    #&debug("COMMAND: ".$command);
    #&debug("ALIGN GROUP: ".$align_group);
    #&debug("FILE: ".$align_file_name);
    print $completed_cmds_file $align_group." ".$align_file_name."\n";
  }
}

sub SubmitJobsToCluster {
  # get function parameters
  my ($jobs_dir, $jobs_path, $job_type, $name_tail, $memory) = @_;
  &status("Submitting $job_type jobs to the cluster...");
  my $job_name = join("-", $OPT{lib}, $name_tail);
  #my %args_hash = (
  #            CLUSTER => $OPT{cluster},
  #            MQSUB => $CFG{commands}{mqsub},
  #            JOBS_DIR => $jobs_dir,
  #            JOB => $jobs_path,
  #            NAME => $job_name,
  #            MEM => '1G',
  #            WALL_TIME => '1:00:00',
  #           );
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $submit_cmd = join(" ",
    $python, $EXEC{submit_script}, $job_name, $jobs_path, $OPT{cluster},
    "--multiple", "--memory ".$memory, "--hostname ".$OPT{hostname}, "--queue ".$OPT{queue},
  );
  if ($OPT{use_wall_time}) {
    $submit_cmd .= " --wall-time 1:00:00";
  }
  if ($OPT{email}) {
    $submit_cmd .= " --email ".$OPT{email};
  }
  if ($OPT{debug_high}) {
    $submit_cmd .= " --debug";
  }
  &debug($submit_cmd);
  my $result = system($submit_cmd);
  if (0 != $result) {
    &error("Error while running ".$submit_cmd.": ".$result);
  }
  return undef;
}

#Parse the output of ac_align_check for generating summary file ->
# This gets called multiple times and serves to group the events by chr pairs
sub ParseCandidateContigs {
  # get function parameters
  my ($path, $aligner, $contig_type, $abs_member_id_ref,
      $ungrouped_members, $contigs_in_events) = @_;
  my $abs_member_id = ${$abs_member_id_ref};

  if ($OPT{cluster}) {
    if (0 == $abs_member_id % 100) {
      print ".";
    }
  } else {
    &status("Parsing candidate contigs in $path...");
  }
  open(CC_FILE, $path) ||
    &error("Can't open candidate contigs file ".$path." ".$!);
  # get a candidate contig line from the file
  while (<CC_FILE>) {
    chomp;
    if ("" eq $_) {
      next;
    }
    # split the candidate contig line into fields
    my @cc_fields = split;
    #&debug($_);
    if ($_ !~ /^CTG:/) {
      &error("Invalid candidate line: \"".$_."\"")
    }
    #next if ($cc_fields[$CC_TARGET] =~ /random/ ||
    #         $cc_fields[$CC_TARGET] =~ /hap/);
    # skip non-standard chromosomes
    #&debug("TARGET: \"".$cc_fields[$CC_TARGET]."\"");
    if ($cc_fields[$CC_TARGET] !~
        /:(chr)?(\d+|[XY]|MT?):.*,(chr)?(\d+|[XY]|MT?):/) {
      next;
    }

    my $candidate_contig = $_;
    # &print_full($logref, "$candidate_contig\n"); # DEBUG
    my ($ctg_id, $ctg_len) = $candidate_contig =~ /CTG:(\S+)\((\d+)bp\)/;

    # Mark the contig as found so we consider it when parsing the kalign
    # alignments later
    $contigs_in_events->{$ctg_id}{found} = 1;
    $contigs_in_events->{$ctg_id}{len} = $ctg_len;

    # Get the junction co-ordinates from the candidate contig
    my (undef, $junction_from, $junction_to) = $cc_fields[$CC_JUNCTION] =~
      /(JUNCTION|BREAKPOINTS):(chr[^-]+)-(chr.+)$/;
    my ($chr_from,$coord_from,$dir_from) =
      $junction_from =~ /chr(\d+|[XYM]|MT):(\d+)\((up|down)\)/;
    my ($chr_to,  $coord_to,  $dir_to) =
      $junction_to   =~ /chr(\d+|[XYM]|MT):(\d+)\((up|down)\)/;
    $chr_from = 'M' if $chr_from eq 'MT';
    $chr_to   = 'M' if $chr_to   eq 'MT';
    # &print_full($logref, "JUNCTION: 1) $chr_from 2) $chr_to\n"); # DEBUG
    if (not $chr_from or not $chr_to) {
      &status("Could not get junction chromosomes!");
      &status("JUNCTION FROM: ".$junction_from);
      &status("JUNCTION TO: ".$junction_to);
      &status("JUNCTION FIELD: ".$cc_fields[$CC_JUNCTION]);
      &status("LINE: ".$_);
      &error("parsing candidates file ".$path);
    }

    # KLUDGE! Should be handled upstream!!
    if ($chr_from eq $chr_to && 2 > abs($coord_from-$coord_to) &&
      (($coord_from <= $coord_to && "up" eq $dir_from && "down" eq $dir_to) ||
       ($coord_to <= $coord_from && "up" eq $dir_to && "down" eq $dir_from))) {
      &warn("\nSkipping \"contiguous\" alignment of contig ".$ctg_id.".");
      next;
    }

    # sort the junction co-ordinates
    my ($chr_from_s,$coord_from_s,$dir_from_s);
    my ($chr_to_s,$coord_to_s,$dir_to_s);
    if ((&ChrToInt($chr_from) < &ChrToInt($chr_to)) ||
        ($chr_from eq $chr_to && $coord_from  < $coord_to) ||
        ($chr_from eq $chr_to && $coord_from == $coord_to &&
         "up" eq $dir_from)) {
      $chr_from_s = $chr_from;
      $coord_from_s = $coord_from;
      $dir_from_s = $dir_from;
      $chr_to_s = $chr_to;
      $coord_to_s = $coord_to;
      $dir_to_s = $dir_to;
    } else {
      $chr_from_s = $chr_to;
      $coord_from_s = $coord_to;
      $dir_from_s = $dir_to;
      $chr_to_s = $chr_from;
      $coord_to_s = $coord_from;
      $dir_to_s = $dir_from;
    }
    my $junction_from_s = "$coord_from_s($dir_from_s)";
    my $junction_to_s   = "$coord_to_s($dir_to_s)";

    #&print_full($logref, "JUNCTION:\n"); # DEBUG
    #&print_full($logref, "    $junction_from\n"); # DEBUG
    #&print_full($logref, "  from ${dir_from_s}stream of ".
    #                       "chr${chr_from_s}:${coord_from_s}\n"); # DEBUG
    #&print_full($logref, "    $junction_to\n"); # DEBUG
    #&print_full($logref, "    to ${dir_to_s}stream of ".
    #                       "chr${chr_to_s}:${coord_to_s}\n"); # DEBUG

    my ($af1,$af2) =
      $cc_fields[$CC_ALIGN_METRICS] =~ /AF1:([^,]+),AF2:([^, ]+)/;
    my $min_af = $af1;
    if ($af2 < $af1) {
      $min_af = $af2;
    }

    # get the alignment topology
    my ($topology) = $cc_fields[$CC_TOPOLOGY] =~ /TOPOLOGY:(.*)/;
    #&debug ("Candidate topology: $topology");
    # Do not do this, because with the new format the contig coordinates are
    # coords1: from ctg2gen align
    # coords2: event region in ctg coords
    #if ($topology =~ /gap/) {
    #  # ensure that the contig co-ordinates are ordered
    #  my ($ctg_left1, $ctg_right1, $ctg_left2, $ctg_right2) =
    #    $cc_fields[$CC_CTG_COORDS] =~ /CONTIG:(\d+)-(\d+),(\d+)-(\d+)/;
    #  #&debug ("Contig Coords: $cc_fields[$CC_CTG_COORDS]");
    #  #&debug ("Contig 1: $ctg_left1-$ctg_right1");
    #  #&debug ("Contig 2: $ctg_left2-$ctg_right2");
    #  if ($ctg_left2 < $ctg_left1) {
    #    #&debug ("Flipping contig coordinates");
    #    my $new_ctg_coords =
    #      "CONTIG:$ctg_left2-$ctg_right2,$ctg_left1-$ctg_right1";
    #    $candidate_contig =~ s/$cc_fields[$CC_CTG_COORDS]/$new_ctg_coords/;
    #  }
    #}

    my $event_grouping_info = join(" ",
      $junction_from_s,
      $junction_to_s,
      $contig_type,
      $aligner,
      $candidate_contig);
    #&debug ("GROUP:\n  $event_grouping_info"); # DEBUG
    my ($jkey_from, $jkey_to);
    {
      use integer;
      my ($junction_from_coord) = $junction_from_s =~ /(\d+)\(/;
      my ($junction_to_coord)   = $junction_to_s   =~ /(\d+)\(/;
      $jkey_from = $junction_from_coord / 100;
      $jkey_to   = $junction_to_coord   / 100;
    }
    #&debug(join(" ", "CTG ID:", $ctg_id,
    #       "CHR FROM:", $chr_from_s, "CHR TO:", $chr_to_s,
    #       "JKEY FROM:", $jkey_from, "JKEY_TO:", $jkey_to));
    # Check for duplicates
    if (exists $ungrouped_members->{$chr_from_s} &&
        exists $ungrouped_members->{$chr_from_s}
                                   {$chr_to_s}   &&
        exists $ungrouped_members->{$chr_from_s}
                                   {$chr_to_s}
                                   {$jkey_from}  &&
        exists $ungrouped_members->{$chr_from_s}
                                   {$chr_to_s}
                                   {$jkey_from}
                                   {$jkey_to}
       ) {
      my %members_in_bin =
        %{$ungrouped_members->{$chr_from_s}{$chr_to_s}{$jkey_from}{$jkey_to}};
      for my $seen_member_id (keys %members_in_bin) {
        my $seen_event_info =
          $members_in_bin{$seen_member_id};
        if ($seen_event_info eq $event_grouping_info) {
          &warn("Duplicate event found: $event_grouping_info");
        }
      }
    }
    #$msg = "Member $abs_member_id: ".join(" ",
    #  $ctg_id, $ctg_len, $min_af, $chr_from_s, $chr_to_s);
    #&print_full($logref, "$msg\n");
    # Store the member grouping info by chr pairs for easier grouping later
    $ungrouped_members->{$chr_from_s}{$chr_to_s}
                        {$jkey_from}{$jkey_to}
                        {$abs_member_id} =
      $event_grouping_info;
    # increment the absolute member index
    $abs_member_id++;
  }
  close CC_FILE;
  ${$abs_member_id_ref} = $abs_member_id;
  return undef;
}

sub ChrToInt {
  # get function parameters
  my ($chr_str) = @_;
  if ("X" eq $chr_str) {
    return 0;
  }
  if ("Y" eq $chr_str) {
    return -1;
  }
  if ("M" eq $chr_str) {
    return -2;
  }
  return $chr_str;
}

#Group all the individual alignments
sub GroupEvents {
  # get function parameters
  my ($ungrouped_members, $contigs_in_events) = @_;
  &status("Grouping events...");
  my $group_start = new Benchmark;
  # Create an array for holding the groups
  my @unsorted_groups = ();
  # while there are ungrouped members
  &status("Creating groups");
  while (keys %{$ungrouped_members}) {
    # create a new group
    my $group =
      &CreateGroup($ungrouped_members);
    # add members to the group
    &AddMembersToGroup($group, $ungrouped_members);
    # check the "neighbour" status of the event
    &SetNeighbourStatus($group);
    # add the group to the array
    push(@unsorted_groups, $group);
    #&debug("-----------------------");
  }
  # sort the groups by max ctg len
  &status("Sorting groups by contig length");
  my @groups = sort(CompareGroups @unsorted_groups);
  # populate the contigs_in_events (ctg->grp mapping) hash
  &status("Mapping contigs to events");
  &MapContigsToEvents(\@groups, $contigs_in_events);
  my $num_groups = @groups;
  &status($num_groups." groups found");
  &TimeSpent("grouping events", $group_start);
  return \@groups;
}

sub CreateGroup {
  # get function parameters
  my ($ungrouped_members) = @_;
  &debug_high("Creating group...");
  my %group = (
    num_members => 0,
    max_ctg_len => 0,
  );
  # get the first member to be in the group
  my $member_id;
  ($group{member_keys}{chr_from},  $group{member_keys}{chr_to},
   $group{member_keys}{jkey_from}, $group{member_keys}{jkey_to},
   $member_id) =
    &GetFirstUngroupedMemberKeys($ungrouped_members);
  my $first_member_info = &PopUngroupedMember(
    $ungrouped_members,
    $group{member_keys}{chr_from}, $group{member_keys}{chr_to},
    $group{member_keys}{jkey_from}, $group{member_keys}{jkey_to},
    $member_id
  );
  &debug_high("FIRST MEMBER:\n  ".$first_member_info);
  # setup up the group based on the first member
  my @member_fields = split(/ /,$first_member_info);
  my ($junction_from, $junction_to) = @member_fields;
  &SetGroupingInfo(\%group, $junction_from, $junction_to);
  &AddMemberToGroup(\%group, $first_member_info, \@member_fields);
  return \%group;
}

sub AddMemberToGroup {
  # get function parameters
  my ($group, $member_info, $member_fields) = @_;
  &debug_high("Adding member to group...");
  $group->{num_members}++;
  # extract the basic info from the member info string
  my (undef, undef, $contig_type, $aligner) = @{$member_fields};
  # mark that the aligner is in the group
  $group->{aligners}{$aligner} = 1;
  my ($ctg_id) = $member_fields->[$EGI_CTG_ID] =~ /CTG:(\S+)\(/;
  if ($ctg_id !~ /./) {
    &error("could not get contig ID from: \"".
            $member_fields->[$EGI_CTG_ID]."\"");
  }
  # add the member info to the group
  &AddMemberInfo($group, $ctg_id, $aligner, $contig_type, $member_fields);
  # add the candidate contig's topology to the group
  # and return the candidate's topology
  my $topology = &AddTopologyToGroup($group, $member_fields->[$EGI_TOPOLOGY]);
  # process the target coordinates for the member
  my $flip_blocks = &ProcessTargetCoords
    ($group, $ctg_id, $topology,
     $member_fields->[$EGI_JUNCTION],
     $member_fields->[$EGI_TARGET]);
  # update the contig lengths and block info
  &UpdateCtgLengthsBlocksAndOverlap
    ($group, $flip_blocks, $topology,
     $member_fields->[$EGI_CTG_ID],
     $member_fields->[$EGI_BLOCKS],
     $member_fields->[$EGI_META]);
  return $ctg_id;
}

sub AddMemberInfo {
  # get function parameters
  my ($group, $ctg_id, $aligner, $contig_type, $member_fields) = @_;
  &debug_high("Adding member info for contig: \"".$ctg_id."\" ".
         "Aligner: \"".$aligner."\"");
  # remove the grouping info
  my $num_fields = @{$member_fields};
  my $member_info =
    join(" ", @{$member_fields}[$EGI_CTG_ID..($num_fields-1)]);
  #&debug("Member info: ".$member_info);
  my $full_info = join(" ", $member_info, $contig_type);
  if (! exists $group->{members}{$ctg_id} ||
      ! exists $group->{members}{$ctg_id}{info}{$aligner}) {
    @{$group->{members}{$ctg_id}{info}{$aligner}} = ();
  }
  push(@{$group->{members}{$ctg_id}{info}{$aligner}}, $full_info);
  return undef;
}

sub AddTopologyToGroup {
  # get function parameters
  my ($group, $new_topology_str) = @_;
  my ($new_topology) = $new_topology_str =~ /TOPOLOGY:(.*)/;
  $group->{topologies}{$new_topology} = 1;
  #&debug("NEW ALIGNMENT TOPOLOGY: ".$new_topology);
  return $new_topology;
}

sub ProcessTargetCoords {
  # get function parameters
  my ($group, $ctg_id, $topology, $junction_field, $target_field) = @_;
  &debug_high("Processing Target Coordinates...");
  # get the junction coordinates
  my (undef, $junction1, $junction2) =
    $junction_field =~ /(JUNCTION|BREAKPOINTS):(chr[^-]+)-(chr.+)$/;
  my ($junction_coord1) = $junction1 =~ /:(\d+)\(/;
  my ($junction_coord2) = $junction2 =~ /:(\d+)\(/;
  #&print_full($logref, "Junction1) $junction_coord1 ");
  #&print_full($logref, "Junction2) $junction_coord2\n");
  # get the sorted but unordered target coordinates
  my ($target_coords1, $target_coords2) =
    &GetSortedTargetCoords($target_field);
  # check whether this member would make the event local
  &CheckWhetherLocal($group, $ctg_id, $target_coords1, $target_coords2);
  my $flip_blocks;
  if (1 == $group->{num_members}) {
    # set the initial target coordinates for the event
    $flip_blocks = &SetTargetCoordinatesForGroup
        ($group, $topology,
         $target_coords1, $junction_coord1,
         $target_coords2, $junction_coord2);
  } else {
    $flip_blocks = &UpdateTargetCoordinatesForGroup
      ($group->{ordered_target_regions}, $topology,
       $target_coords1, $junction_coord1,
       $target_coords2, $junction_coord2);
  }
  return $flip_blocks;
}

sub CheckWhetherLocal {
  # get function parameters
  my ($group, $ctg_id, $target_coords1, $target_coords2) = @_;
  if ($group->{member_keys}{chr_from} eq
      $group->{member_keys}{chr_to} &&
      &AreCoordsLocal($target_coords1, $target_coords2)) {
    #&print_full($logref, "Group $group_id: $ctg_id1 = local\n");
    $group->{members}{$ctg_id}{islocal} = 1;
    $group->{islocal} = 1;
  }
  return undef;
}

sub UpdateCtgLengthsBlocksAndOverlap {
  # get function parameters
  my ($group, $flip_blocks, $topology,
      $ctg_field, $blocks_field, $meta_field) = @_;
  # extract the contig length from the contig field
  my ($contig_len) = $ctg_field =~ /\((\d+)bp\)/;
  # add the contig length for the group
  $group->{ctg_lens}{$contig_len} = 1;
  if ($contig_len > $group->{max_ctg_len}) {
    # use the alignment blocks from this contig
    #&print_full($logref, "Using block info for group $group->{id}\n");
    #&print_full($logref, "Ctg len: ".$contig_len."\n");
    $group->{max_ctg_len} = $contig_len;
    my ($blocks1, $blocks2) =
      $blocks_field =~ /BLOCKS:([^;]+);(\S+)/;
    if (! defined $blocks1 || ! defined $blocks2) {
      &error("Could not parse blocks from ".$blocks_field);
    }
    #&debug("Blocks: $blocks1; $blocks2\n");
    my %targets = %{$group->{ordered_target_regions}};
    if ($flip_blocks) {
      $targets{A}{blocks} = $blocks2;
      $targets{B}{blocks} = $blocks1;
    } else {
      $targets{A}{blocks} = $blocks1;
      $targets{B}{blocks} = $blocks2;
    }
    # update the contig overlap info
    $group->{ctg_overlap} = "N/A";
    if ($topology !~ /gap/) {
      ($group->{ctg_overlap}) = $meta_field =~ /CO:([^b]+)bp/;
    }
    if ($group->{ctg_overlap} !~ /./) {
      &error("Could not get contig overlap from ".$meta_field);
    }
  }
  return undef;
}

sub SetGroupingInfo {
  # get function parameters
  my ($group, $junction_from, $junction_to) = @_;
  #($group->{grouping_info}{set_from}, $group->{grouping_info}{dir_from}) =
  ($group->{grouping_info}{bl_from}, $group->{grouping_info}{br_from},
     $group->{grouping_info}{dir_from}) = &GetGroupingInfo($junction_from);
  #($group->{grouping_info}{set_to},   $group->{grouping_info}{dir_to}) =
  ($group->{grouping_info}{bl_to}, $group->{grouping_info}{br_to},
    $group->{grouping_info}{dir_to}) = &GetGroupingInfo($junction_to);
  return undef;
}

sub GetGroupingInfo {
  # get function parameters
  my ($junction_str) = @_;
  # Buffer to add to junction regions to determine whether to group pairs
  my $group_buffer = ceil($OPT{read_length} / 2);
  my ($coord, $dir) = $junction_str =~ /(\d+)\((up|down)\)/;
  if ($coord !~ /./ || $dir !~ /./) {
    &error("Could not get junction co-ordinates from: ".$junction_str);
  }
  my $buffer_left  = $coord-$group_buffer;
  my $buffer_right = $coord+$group_buffer;
  #my $buffer_range = join("-", $buffer_left, $buffer_right);
  #my $set = Set::IntSpan->new($buffer_range);
  #return ($set, $dir);
  return ($buffer_left, $buffer_right, $dir);
}

sub GetFirstUngroupedMemberKeys {
  # get function parameters
  my ($ungrouped_members) = @_;
  my $chr_from  =
    (sort CompareChrom keys %{$ungrouped_members})[0];
  my $chr_to    =
    (sort CompareChrom keys %{$ungrouped_members->{$chr_from}})[0];
  my $jkey_from =
    (sort keys %{$ungrouped_members->{$chr_from}{$chr_to}})[0];
  my $jkey_to   =
    (sort keys %{$ungrouped_members->{$chr_from}{$chr_to}{$jkey_from}})[0];
  my $member_id =
    (sort keys %{$ungrouped_members->{$chr_from}{$chr_to}
                                     {$jkey_from}{$jkey_to}})[0];
  return ($chr_from, $chr_to, $jkey_from, $jkey_to, $member_id);
}

sub CompareChrom {
  # parameters are automatically $a and $b
  return &ChrToInt($a) <=> &ChrToInt($b);
}

sub PopUngroupedMember {
  # get function parameters
  my ($ungrouped_members, $chr_from, $chr_to,
      $jkey_from, $jkey_to, $member_id) = @_;
  # remove the member from the hash
  my $member_info =
    delete $ungrouped_members->{$chr_from}{$chr_to}{$jkey_from}{$jkey_to}
                               {$member_id};
  # remove whole bins as they get empty
  if (! keys %{$ungrouped_members->{$chr_from}{$chr_to}
                                   {$jkey_from}{$jkey_to}}) {
    delete $ungrouped_members->{$chr_from}{$chr_to}{$jkey_from}{$jkey_to};
    if (! keys %{$ungrouped_members->{$chr_from}{$chr_to}{$jkey_from}}) {
      delete $ungrouped_members->{$chr_from}{$chr_to}{$jkey_from};
      if (! keys %{$ungrouped_members->{$chr_from}{$chr_to}}) {
        delete $ungrouped_members->{$chr_from}{$chr_to};
        if (! keys %{$ungrouped_members->{$chr_from}}) {
          delete $ungrouped_members->{$chr_from};
        }
      }
    }
  }
  return $member_info;
}

sub AddMembersToGroup {
  # get function parameters
  my ($group, $ungrouped_members) = @_;
  &debug_high("Adding members to group...");
  if (! exists $ungrouped_members->{$group->{member_keys}{chr_from}} ||
      ! exists $ungrouped_members->{$group->{member_keys}{chr_from}}
                                   {$group->{member_keys}{chr_to}}) {
    return undef;
  }
  #&debug("Adding additional members to group...");
  my %members_with_curr_chr_pair =
    %{$ungrouped_members->{$group->{member_keys}{chr_from}}
                          {$group->{member_keys}{chr_to}}};
  my @added_members = ();
  my $jkey_from = $group->{member_keys}{jkey_from};
  # add members with current junction from key
  push(@added_members, &AddMembersWithJunctionFromKey
    ($group, \%members_with_curr_chr_pair, $jkey_from));
  # add members with next junction from key
  push(@added_members, &AddMembersWithJunctionFromKey
    ($group, \%members_with_curr_chr_pair, $jkey_from+1));
  # remove added members from ungrouped_members
  &RemoveFromUngrouped
    (\@added_members, $ungrouped_members);
  return undef;
}

sub AddMembersWithJunctionFromKey {
  # get function parameters
  my ($group, $members_with_curr_chr_pair, $jkey_from) = @_;
  if (! exists $members_with_curr_chr_pair->{$jkey_from}) {
    return ();
  }
  #&debug("Looking for members with same from coord: ".$jkey_from."...");
  my %members_with_curr_jkey_from =
    %{$members_with_curr_chr_pair->{$jkey_from}};
  my @added_members = ();
  my $jkey_to = $group->{member_keys}{jkey_to};
  # add members with current junction to key
  push(@added_members, &AddMembersWithJunctionToKey
    ($group, \%members_with_curr_jkey_from, $jkey_from, $jkey_to));
  # add members with next junction to key
  push(@added_members, &AddMembersWithJunctionToKey
    ($group, \%members_with_curr_jkey_from, $jkey_from, $jkey_to+1));
  return @added_members;
}

sub AddMembersWithJunctionToKey {
  # get function parameters
  my ($group, $members_with_curr_jkey_from, $jkey_from, $jkey_to) = @_;
  if (! exists $members_with_curr_jkey_from->{$jkey_to}) {
    return ();
  }
  #&debug("Looking for members with same to coord: ".$jkey_to."...");
  my %members_with_curr_jkey_to =
    %{$members_with_curr_jkey_from->{$jkey_to}};
  my @added_members = ();
  for my $member_id (sort { $a <=> $b } keys %members_with_curr_jkey_to) {
    my $member_info = $members_with_curr_jkey_to{$member_id};
    my @member_fields = split(/ /,$member_info);
    my ($junction_from, $junction_to) = @member_fields;
    if (&MemberIsInGroup($group, $junction_from, $junction_to)) {
      # add the member to the group
      &AddMemberToGroup($group, $member_info, \@member_fields);
      # mark the member as being added
      push(@added_members,
        join(" ",
          $group->{member_keys}{chr_from}, $group->{member_keys}{chr_to},
          $jkey_from, $jkey_to, $member_id));
    }
  }
  return @added_members;
}

sub MemberIsInGroup {
  # get function parameters
  my ($group, $junction_from, $junction_to) = @_;
  #my ($set_from, $dir_from) = &GetGroupingInfo($junction_from);
  my ($bl_from, $br_from, $dir_from) = &GetGroupingInfo($junction_from);
  #my ($set_to,   $dir_to)   = &GetGroupingInfo($junction_to);
  my ($bl_to,   $br_to,   $dir_to)   = &GetGroupingInfo($junction_to);
  if ($group->{grouping_info}{dir_from} eq $dir_from &&
      $group->{grouping_info}{dir_to}   eq $dir_to   &&
      #$group->{grouping_info}{set_from}->intersect($set_from)->cardinality() &&
      &RangesOverlap($group->{grouping_info}{bl_from}, $group->{grouping_info}{br_from}, $bl_from, $br_from) &&
      #$group->{grouping_info}{set_to}->intersect($set_to)->cardinality()
      &RangesOverlap($group->{grouping_info}{bl_to}, $group->{grouping_info}{br_to}, $bl_to, $br_to)
     ) {
    # update the grouping coordinate sets
    ($group->{grouping_info}{bl_from}, $group->{grouping_info}{br_from}) = &MergeRanges(
      $group->{grouping_info}{bl_from}, $group->{grouping_info}{br_from}, $bl_from, $br_from);
    #$group->{grouping_info}{set_from} =
      #$group->{grouping_info}{set_from}->union($set_from);
    ($group->{grouping_info}{bl_to}, $group->{grouping_info}{br_to}) = &MergeRanges(
      $group->{grouping_info}{bl_to}, $group->{grouping_info}{br_to}, $bl_to, $br_to);
    #$group->{grouping_info}{set_to} =
      #$group->{grouping_info}{set_to}->union($set_to);
    return 1;
  }
  return 0;
}

sub RemoveFromUngrouped {
  # get function parameters
  my ($added_members, $ungrouped_members) = @_;
  for my $added_member (@{$added_members}) {
    my ($chr_from, $chr_to, $jkey_from, $jkey_to, $member_id) =
      split(/ /, $added_member);
    # remove the group from the ungrouped list
    &PopUngroupedMember
      ($ungrouped_members, $chr_from, $chr_to,
       $jkey_from, $jkey_to, $member_id);
  }
  return undef;
}

sub SetNeighbourStatus {
  # get function parameters
  my ($group) = @_;
  if ($group->{islocal}) {
    $group->{neighbourness} = 'neighbour';
  } elsif ($group->{member_keys}{chr_from} eq $group->{member_keys}{chr_to}) {
    $group->{neighbourness} = 'samechr';
  } else {
    $group->{neighbourness} = 'diffchr';
  }
  return undef;
}

sub CompareGroups {
  # parameters are automatically $a and $b
  return $b->{max_ctg_len} <=> $a->{max_ctg_len};
}

sub MapContigsToEvents {
  # get function parameters
  my ($groups, $contigs_in_events) = @_;
  my $group_id = 0;
  for my $group (@{$groups}) {
    $group->{id} = $group_id;
    for my $ctg_id (keys %{$group->{members}}) {
      # record that the contig is in the current event
      $contigs_in_events->{$ctg_id}{groups}{$group_id} = 1;
    }
    $group_id++;
  }
  return undef;
}

sub SetTargetCoordinatesForGroup {
  # get function parameters
  my ($group, $topology,
      $target_coords1, $junction_coord1,
      $target_coords2, $junction_coord2) = @_;
  &debug_high("Setting Target Coordinates for Group...");
  my ($regionA, $junctionA,
      $regionB, $junctionB,
      $flip_blocks) =
    &GetOrderedTargetCoords
    ($topology,
     $target_coords1, $junction_coord1,
     $target_coords2, $junction_coord2);
  %{$group->{ordered_target_regions}{A}} = (
    chrom    => $regionA->{chrom},
    left     => $regionA->{left},
    right    => $regionA->{right},
    strand   => $regionA->{strand},
    junction => $junctionA,
  );
  %{$group->{ordered_target_regions}{B}} = (
    chrom    => $regionB->{chrom},
    left     => $regionB->{left},
    right    => $regionB->{right},
    strand   => $regionB->{strand},
    junction => $junctionB,
  );
  &debug_high("Junction A: ".$junctionA." Junction B: ".$junctionB);
  return $flip_blocks;
}

sub UpdateTargetCoordinatesForGroup {
  # get function parameters
  my ($ordered_target_regions, $topology,
      $target_coords1, $junction_coord1,
      $target_coords2, $junction_coord2) = @_;
  &debug_high("Updating target coordinates for group");
  my ($regionA, $junctionA,
      $regionB, $junctionB,
      $flip_blocks) =
    &GetOrderedTargetCoords
    ($topology,
     $target_coords1, $junction_coord1,
     $target_coords2, $junction_coord2);
  # they should all have the same regionA->{chrom}
  my $check = $ordered_target_regions->{A}{chrom};
  if ($regionA->{chrom} ne $check) {
    &error("target chromosome ".$regionA->{chrom}.
      " is not the same as ".$check);
  }
  # they should all have the same junction side A
  if ("gap" eq $ordered_target_regions->{A}{junction}) {
    $ordered_target_regions->{A}{junction} = $junctionA;
  }
  $check = $ordered_target_regions->{A}{junction};
  if ("gap" ne $junctionA && $check ne $junctionA) {
    &error("junction side A: ".$junctionA." is not the same as ".$check);
  }
  if ($regionA->{left} < $ordered_target_regions->{A}{left}) {
    $ordered_target_regions->{A}{left} = $regionA->{left};
  }
  if ($regionA->{right} > $ordered_target_regions->{A}{right}) {
    $ordered_target_regions->{A}{right} = $regionA->{right};
  }
  # they should all have the same regionB->{chrom}
  $check = $ordered_target_regions->{B}{chrom};
  if ($regionB->{chrom} ne $check) {
    &error($regionB->{chrom}." is not the same as ".$check);
  }
  # they should all have the same junction side B
  if ("gap" eq $ordered_target_regions->{B}{junction}) {
    $ordered_target_regions->{B}{junction} = $junctionB;
  }
  $check = $ordered_target_regions->{B}{junction};
  if ("gap" ne $junctionB && $check ne $junctionB) {
    &error("junction side B: ".$junctionB." is not the same as ".$check);
  }
  if ($regionB->{left} < $ordered_target_regions->{B}{left}) {
    $ordered_target_regions->{B}{left} = $regionB->{left};
  }
  if ($regionB->{right} > $ordered_target_regions->{B}{right}) {
    $ordered_target_regions->{B}{right} = $regionB->{right};
  }
  return $flip_blocks;
}

sub GetOrderedTargetCoords {
  # get function parameters
  my ($topology,
      $target_coords1, $junction_coord1,
      $target_coords2, $junction_coord2) = @_;
  &debug_high("Getting ordered target coordinates...");
  my $flip_blocks = 0;
  my (%regionA, %regionB);
  #my ($chrA, $leftA, $rightA, $strandA,
  #    $chrB, $leftB, $rightB, $strandB);
  my ($junction_sideA, $junction_sideB) = ('gap', 'gap');
  &debug_high("Junction Side A: ".$junction_sideA." Junction Side B: ".$junction_sideB);
  if ($topology =~ /gap/) {
    &debug_high("Getting ordered t-coords for gap");
    $regionA{chrom} = $regionB{chrom} = $target_coords2->{chrom};
    $regionA{left}  = $regionA{right} = $target_coords2->{left};
    $regionB{left}  = $regionB{right} = $target_coords2->{right};
    # if it is a gap, do not set the junction
  # order target regions by the junction side, rather than the left-side
  } elsif ((&ChrToInt($target_coords1->{chrom}) <
            &ChrToInt($target_coords2->{chrom})) ||
           ($target_coords1->{chrom} eq $target_coords2->{chrom} &&
            $junction_coord1 < $junction_coord2)) {
    &debug_high("Maintaining target coord order");
    $regionA{chrom}  = $target_coords1->{chrom};
    $regionA{left}   = $target_coords1->{left};
    $regionA{right}  = $target_coords1->{right};
    $regionA{strand} = $target_coords1->{strand};
    if ("+" eq $regionA{strand}) {
      $junction_sideA = 'right';
    } else {
      $junction_sideA = 'left';
    }
    $regionB{chrom}  = $target_coords2->{chrom};
    $regionB{left}   = $target_coords2->{left};
    $regionB{right}  = $target_coords2->{right};
    $regionB{strand} = $target_coords2->{strand};
    if ("+" eq $regionB{strand}) {
      $junction_sideB = 'left';
    } else {
      $junction_sideB = 'right';
    }
  } else { # if (($target_coords1->{chrom} >
           #      $target_coords2->{chrom}) ||
           #     ($target_coords1->{chrom} eq $target_coords2->{chrom} &&
           #      $junction_coord1 >= $junction_coord2))
    #Here we flip the order
    &debug_high("Flipping target coord order");
    $flip_blocks = 1;
    $regionA{chrom}  = $target_coords2->{chrom};
    $regionA{left}   = $target_coords2->{left};
    $regionA{right}  = $target_coords2->{right};
    $regionA{strand} = $target_coords2->{strand};
    if ("+" eq $regionA{strand}) {
      $junction_sideA = 'left';
    } else {
      $junction_sideA = 'right';
    }
    $regionB{chrom}  = $target_coords1->{chrom};
    $regionB{left}   = $target_coords1->{left};
    $regionB{right}  = $target_coords1->{right};
    $regionB{strand} = $target_coords1->{strand};
    if ("+" eq $regionB{strand}) {
      $junction_sideB = 'right';
    } else {
      $junction_sideB = 'left';
    }
  }
  #&debug("1) &RegionCoordsStr($target_coords1) ".
  #       "2) &RegionCoordsStr($target_coords2)");
  #&debug("A) &RegionCoordsStr(\%regionA) ".
  #       "B) &RegionCoordsStr(\%regionB)");
  #&debug("Junction Side A: $junction_sideA ".
  #       "Junction Side B: $junction_sideB");
  #&debug("  Flip Blocks: ".$flip_blocks."\n");
  return (\%regionA, $junction_sideA,
          \%regionB, $junction_sideB,
          $flip_blocks);
}

sub GetSortedTargetCoords {
  # get function parameters
  my ($targets) = @_;
  my ($chr1, $range1, $chr2, $range2) =
    $targets =~ /TARGET:chr([^:]+):(\d+-\d+),chr([^:]+):(\d+-\d+)/;
  my ($left1,$right1) = &Sort_Range($range1);
  my $strand1 = &GetStrandFromRange($range1);
  my %target_coords1 = (
    chrom  => $chr1,
    left   => $left1,
    right  => $right1,
    strand => $strand1,
  );
  my ($left2,$right2) = &Sort_Range($range2);
  my $strand2 = &GetStrandFromRange($range2);
  my %target_coords2 = (
    chrom  => $chr2,
    left   => $left2,
    right  => $right2,
    strand => $strand2,
  );
  return (\%target_coords1, \%target_coords2);
}

sub GetStrandFromRange {
  # get function parameters
  my ($range) = @_;
  my ($start, $end) = $range =~ /(\d+)-(\d+)/;
  if ($start < $end) {
    return "+";
  } else {
    return "-";
  }
}

sub AreCoordsLocal {
  # get function parameters
  my ($coords1, $coords2) = @_;
  #Distance to determine whether we're dealing with local alignments
  my $local_buffer = 100000;
  # Create buffered sets and see if they overlap
  my $buffered_left1  = $coords1->{left}  - $local_buffer;
  my $buffered_right1 = $coords1->{right} + $local_buffer;
  my $buffered_left2  = $coords2->{left}  - $local_buffer;
  my $buffered_right2 = $coords2->{right} + $local_buffer;
  #my $set1 = Set::IntSpan->new("$buffered_left1-$buffered_right1");
  #my $set2 = Set::IntSpan->new("$buffered_left2-$buffered_right2");
  #if ($set1->intersect($set2)->cardinality()) {
  #  return 1;
  #}
  #return 0;
  return &RangesOverlap($buffered_left1, $buffered_right1,
    $buffered_left2, $buffered_right2);
}

sub FinishCandidateIdentification {
  # no parameters
  my $continue = $barnacle_cmd;
  my $barnacle_path = catfile($OPT{barnacle_src_dir}, "barnacle.pl");
  # make sure to use absolute paths in run_support script
  $continue =~ s/^.*barnacle\.pl/$barnacle_path/;
  $continue =~ s/-lib_dir [^ ]+/-lib_dir $OPT{lib_dir}/;
  $continue =~ s/-config [^ ]+/-config $OPT{config}/;
  $continue =~ s/-outdir [^ ]+/-outdir $OPT{outdir}/;
  $continue =~ s/identify_candidates/add_support/;
  if ($continue !~ / -assembly_ver /) {
    $continue .= " -assembly_ver ".$OPT{assembly_ver};
  }
  if ($continue !~ / -ver /) {
    $continue .= " -ver ".$OPT{ver};
  }
  if ($OPT{cluster}) {
    &status("All candidate identification jobs submitted to cluster. ".
      "Run the following command when jobs are complete:\n".$continue);
  } else {
    &status("All candidate identification jobs complete. ".
      "Run the following command to add support:\n".$continue);
  }
  #prepare job file for getting support
  my $support_job_path = catfile($OPT{outdir}, "run_support.sh");
  open(my $JOB, ">", $support_job_path) ||
    &error("cannot open ".$support_job_path.": ".$!);
  my $command = &JobLine($continue);
  $command =~ s/;/\n/g;
  print $JOB "#! /bin/sh\n#\$ -S /bin/sh\n\n".$command."\n";
  close $JOB;
  system("chmod u+x ".$support_job_path);
  return undef;
}

sub AddPairSupportToGroups {
  # get function parameters
  my ($groups, $previous_output_dir) = @_;
  if ($OPT{no_pair_to_genome}) {
    &status("Writing event coordinates...");
  } else {
    &status("Getting read pair to genome support...");
    &status("Using file: ".$CFG{alignments}{p2g});
    &debug("Should chromosome name have \"chr\"? ".
      $CFG{read_support}{p2g}{use_chr});
  }
  my $pair_start = new Benchmark;
  my %groups_with_no_pairs = ();
  # initialize the cluster info
  my %p2g_cluster_info = (
    subfile_counter => 1,
    group_counter   => 0,
  );
  if ($OPT{cluster_pair2gen}) {
    # prepare for creating pair to genome cluster jobs
    &PrepareForPairToGenomeCluster(\%p2g_cluster_info);
  }
  for my $group (@{$groups}) {
    # if ready to go on to the next file
    if ($OPT{pair2gen_split} <= $p2g_cluster_info{group_counter} &&
        $OPT{cluster_pair2gen}) {
      &GoToNextP2GSubfile(\%p2g_cluster_info);
    }
    # REMINDER: use all topologies
    my @topologies   = sort keys %{$group->{topologies}};
    my %targets = %{$group->{ordered_target_regions}};
    my %coordsA = (
      chrom    => $targets{A}{chrom},
      left     => $targets{A}{left},
      right    => $targets{A}{right},
      junction => $targets{A}{junction},
      blocks   => $targets{A}{blocks},
    );
    my %coordsB = (
      chrom    => $targets{B}{chrom},
      left     => $targets{B}{left},
      right    => $targets{B}{right},
      junction => $targets{B}{junction},
      blocks   => $targets{B}{blocks},
    );
    if ($OPT{debug_high}) {
      &status("   A) ".&RegionCoordsStr(\%coordsA).
               ", B) ".&RegionCoordsStr(\%coordsB)); # DEBUG
      &status("Blocks: $coordsA{blocks}; $coordsB{blocks}"); # DEBUG
    }
    my %pair_info = (
      group_id    => $group->{id},
      coords      => [\%coordsA, \%coordsB],
      # REMINDER: use all topologies
      topology    => $topologies[0],
      ctg_overlap => $group->{ctg_overlap},
      ctg_ids     => join(",", sort keys %{$group->{members}}),
      num_reads   => 0,
      p2g_error   => 0,
    );
    &debug_high("Pair info ctg_ids: $pair_info{ctg_ids}");
    # try not to use a gap topology
    my $num_topologies = @topologies;
    if ($pair_info{topology} =~ /gap/ && 1 < $num_topologies) {
      for my $topology (@topologies) {
        if ($topology !~ /gap/) {
          $pair_info{topology} = $topology;
          last;
        }
      }
    }
    if ($OPT{cluster_pair2gen}) {
      # write the pair info to the event sub-file
      &WritePairInfoToSubFile(\%pair_info, \%p2g_cluster_info);
      # increment the group counter
      $p2g_cluster_info{group_counter} += 1;
    }
    unless ($OPT{no_pair_to_genome} || $OPT{cluster_pair2gen}) {
      ($group->{pair_support},
       $group->{intronic_pair_support},
       $group->{filtered_pair_support},
       $group->{intronic_filtered_pair_support}) =
        &GetPairSupport(\%pair_info);
      $group->{p2g_error} = $pair_info{p2g_error};
      &debug("Pair support: $group->{pair_support} ".
             "Filtered: $group->{filtered_pair_support}.");
      &debug("Num Reads: ".$pair_info{num_reads}."\n");
      if (0 == $pair_info{num_reads}) {
        #&print_full($logref, "Found group with no sam reads: ".
        #                     $group_id."\n");
        $groups_with_no_pairs{$group->{id}} = 1;
      }
    }
  }
  if ($OPT{cluster_pair2gen}) {
    # close the current event sub-file
    close $p2g_cluster_info{subfile_handle};
  }
  if ($OPT{cluster_pair2gen}) {
    # add the last pair-to-genome job to the jobs file
    &AddP2GJobToFile(\%p2g_cluster_info);
    # submit the pair-to-genome jobs to the cluster
    &SubmitJobsToCluster(
      $p2g_cluster_info{dir},
      $p2g_cluster_info{jobs_path},
      "pair-to-genome",
      "p2g", $OPT{p2g_memory},
    );
  } else {
    my $num_groups_with_no_pairs = keys %groups_with_no_pairs;
    #&print_full($logref,
    #  "$num_groups_with_no_pairs groups with no pairs found\n");
    if (0 < $num_groups_with_no_pairs) {
      if ($OPT{debug_high}) {
        $warning = "Samtools view returned no reads for ".
                   "the following groups: ";
        $warning .= join(",", sort { $a <=> $b } keys %groups_with_no_pairs);
        &warn($warning);
      } else {
        &warn("Samtools view returned no reads for ".
              "$num_groups_with_no_pairs groups");
      }
    }
  }
  &TimeSpent("getting pair to genome support for events", $pair_start);
  return undef;
}

sub PrepareForPairToGenomeCluster {
  # get function parameters
  my ($p2g_cluster_info) = @_;
  # setup the pair2gen cluster directory
  $p2g_cluster_info->{dir} = &SetupPair2GenClusterDir();
  # open the first event sub-file
  &OpenEventSubFile($p2g_cluster_info);
  # open the pair-to-genome jobs file
  $p2g_cluster_info->{jobs_path} =
    catfile($p2g_cluster_info->{dir}, "pair2gen_jobs");
  open($p2g_cluster_info->{jobs_handle}, ">",
       $p2g_cluster_info->{jobs_path}) ||
    &error("Can't open pair-to-genome jobs file ".
        $p2g_cluster_info->{jobs_path}." ".$!);
  return undef;
}

sub SetupPair2GenClusterDir {
  # no function params
  my $p2g_cluster_dir = catdir($OPT{results_dir}, "cluster_p2g");
  # create the directory, if necessary
  make_path($p2g_cluster_dir) if (!-d $p2g_cluster_dir);
  return $p2g_cluster_dir;
}

sub GoToNextP2GSubfile {
  # get function parameters
  my ($p2g_cluster_info) = @_;
  # close the current event sub-file
  close $p2g_cluster_info->{subfile_handle};
  if ($OPT{cluster_pair2gen}) {
    # add a new pair-to-genome job to the jobs file
    &AddP2GJobToFile($p2g_cluster_info);
  }
  # open the next event sub-file
  $p2g_cluster_info->{subfile_counter} += 1;
  &OpenEventSubFile($p2g_cluster_info);
  # reset the group counter
  $p2g_cluster_info->{group_counter} = 0;
  return undef;
}

sub OpenEventSubFile {
  # get function parameters
  my ($p2g_cluster_info) = @_;
  $p2g_cluster_info->{job_dir} = catdir(
    $p2g_cluster_info->{dir},
    "job_".$p2g_cluster_info->{subfile_counter}
  );
  # create job_dir if necessary
  if (!-d $p2g_cluster_info->{job_dir}) {
    make_path($p2g_cluster_info->{job_dir});
  }
  my $subfile_name =
    $OPT{lib}.".".$p2g_cluster_info->{subfile_counter}.".p2g.data";
  $p2g_cluster_info->{subfile_path} =
    catfile($p2g_cluster_info->{job_dir}, $subfile_name);
  open($p2g_cluster_info->{subfile_handle}, ">",
       $p2g_cluster_info->{subfile_path}) ||
    &error("Can't open pair-to-genome events sub-file ".
      $p2g_cluster_info->{subfile_path}." ".$!);
  return undef;
}

sub AddP2GJobToFile {
  # get function parameters
  my ($p2g_cluster_info) = @_;
  # create pair-to-genome job command
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $p2g_command = join(" ",
    $python, $EXEC{p2g_script},
    $p2g_cluster_info->{subfile_path},
    $CFG{alignments}{p2g},
    # options
    "--read-length", $OPT{read_length},
    "--frag-len",    $OPT{frag_len},
    "--frag-fract",  $OPT{frag_diff},
    "--min-mapq",    $OPT{min_mapq},
    "--max-sam-retries", $OPT{max_sam_retries},
  );
  if (! $OPT{cluster_pair2gen}) {
    $p2g_command .= " --log-file ".$log_file_name;
  }
  if ($OPT{disable_profiling_timer}) {
    $p2g_command .= " --disable-profiling-timer";
  }
  # write command to pair-to-genome jobs file
  &debug("Add pair-to-genome job to list: ".$p2g_command);
  print { $p2g_cluster_info->{jobs_handle} } &JobLine($p2g_command);
  return undef;
}

sub WritePairInfoToSubFile {
  # get function parameters
  my ($pair_info, $p2g_cluster_info) = @_;
  my @pair_info_list = (
    $pair_info->{group_id},
    &PairInfoCoordsString($pair_info->{coords}[0]),
    &PairInfoCoordsString($pair_info->{coords}[1]),
    $pair_info->{topology},
    $pair_info->{ctg_overlap},
    $pair_info->{ctg_ids},
  );
  my $pair_info_string = join("\t", @pair_info_list);
  #&debug_high("PAIR INFO: ".$pair_info_string);
  print { $p2g_cluster_info->{subfile_handle} } $pair_info_string."\n";
  return undef;
}

sub PairInfoCoordsString {
  # get function parameters
  my ($coords_data) = @_;
  #&debug_high("PAIR_COORDS:\n  ".
  #  $coords_data->{chrom}."\n  ".
  #  $coords_data->{left}."\n  ".
  #  $coords_data->{right}."\n  ".
  #  $coords_data->{junction}."\n  ".
  #  $coords_data->{blocks}
  #);
  my $chrom = $coords_data->{chrom};
  if ($CFG{read_support}{p2g}{use_chr}) {
    $chrom = "chr".$chrom;
  }
  if ($CFG{read_support}{p2g}{use_MT}) {
    $chrom =~ s/M$/MT/;
  }
  my @coords_data_list = (
    $chrom,
    $coords_data->{left},
    $coords_data->{right},
    $coords_data->{junction},
    $coords_data->{blocks},
  );
  my $coords_string = join(" ", @coords_data_list);
  #&debug_high("PAIR_COORDS: ".$coords_string);
  return $coords_string;
}

#Parse the wtss snp lists from Ryan
sub Parse_RNA_Repeats {
  # get function parameters
  my ($groups) = @_;
  &status("Parsing RNA repeat file ".$CFG{coords}{structural_RNA}."...");
  my %rna_sets = ();
  my $parse_rna_start = new Benchmark;
  open(RNA_FILE, $CFG{coords}{structural_RNA}) ||
    &error("Can't open rna file ".$CFG{coords}{structural_RNA}." ".$!);
  while (<RNA_FILE>) {
    my ($chr,$left,$right,$rname) = split;
    if ($chr !~ /^chr/) {
      $chr = "chr$chr";
    }
    #my $rna_set = Set::IntSpan->new("$left-$right");
    if (exists $rna_sets{$chr}) {
      #$rna_sets{$chr} = $rna_sets{$chr}->union($rna_set);
      ($rna_sets{$chr}{left}, $rna_sets{$chr}{right}) = &MergeRanges(
        $rna_sets{$chr}{left}, $rna_sets{$chr}{right}, $left, $right);
    } else {
      #$rna_sets{$chr} = $rna_set;
      $rna_sets{$chr}{left} = $left;
      $rna_sets{$chr}{right} = $right;
    }
  }
  close RNA_FILE;

  # Check if any of the group members overlap the RNA regions
  #for my $group_id (sort {$a<=>$b} keys %{$groups}) {}
  for my $group (@{$groups}) {
    #&print_full($logref, "Group: ".$group_id."\n");
    for my $ctg_id (keys %{$group->{members}}) {
      #&print_full($logref, "  Ctg ID: ".$ctg_id."\n");
      for my $aligner (keys %{$group->{members}{$ctg_id}{info}}) {
        #&print_full($logref, "    Aligner: ".$aligner."\n");
        my @members = @{$group->{members}{$ctg_id}{info}{$aligner}};
        my $num_members = @members; # DEBUG
        #&debug("Num members: ".$num_members."\n"); # DEBUG
        for my $member_info (@members) {
          #Parse the relevant info
          my @member_fields = split(/ /,$member_info);
          my ($chr1,$start1,$end1) =
            $member_fields[$CC_TARGET] =~
            /TARGET:(chr[0-9XYMT]+):(\d+)-(\d+),/;
          my ($chr2,$start2,$end2) =
            $member_fields[$CC_TARGET] =~
            /TARGET:[^,]*,(chr[0-9XYMT]+):(\d+)-(\d+)/;
          $chr1 = 'chrM' if $chr1 =~ /MT/;
          $chr2 = 'chrM' if $chr2 =~ /MT/;

          #Check the edges against the RNA elements
          #my $set1 = Set::IntSpan->new("$start1-$start1");
          #my $set2 = Set::IntSpan->new("$end1-$end1");
          #my $set3 = Set::IntSpan->new("$start2-$start2");
          #my $set4 = Set::IntSpan->new("$end2-$end2");
          #if (!$set1->intersect($rna_sets{$chr1})->empty() ||
          #    !$set2->intersect($rna_sets{$chr1})->empty() ||
          #    !$set3->intersect($rna_sets{$chr2})->empty() ||
          #    !$set4->intersect($rna_sets{$chr2})->empty()) {}
          if (&CoordInRange($start1, $rna_sets{$chr1}{left}, $rna_sets{$chr1}{right}) ||
              &CoordInRange($end1,   $rna_sets{$chr1}{left}, $rna_sets{$chr1}{right}) ||
              &CoordInRange($start2, $rna_sets{$chr2}{left}, $rna_sets{$chr2}{right}) ||
              &CoordInRange($end2,   $rna_sets{$chr2}{left}, $rna_sets{$chr2}{right})) {
            &debug_high("  Coord ".$start1." within range ".$chr1.":".$rna_sets{$chr1}{left}."-".$rna_sets{$chr1}{right}." or...");
            &debug_high("  Coord ".$end1  ." within range ".$chr1.":".$rna_sets{$chr1}{left}."-".$rna_sets{$chr1}{right}." or...");
            &debug_high("  Coord ".$start2." within range ".$chr2.":".$rna_sets{$chr2}{left}."-".$rna_sets{$chr2}{right}." or...");
            &debug_high("  Coord ".$end2  ." within range ".$chr2.":".$rna_sets{$chr2}{left}."-".$rna_sets{$chr2}{right});
            $group->{rna} = 1;
          }
        }
      }
    }
  }
  &TimeSpent("parsing RNA repeat file", $parse_rna_start);
  return undef;
}

sub CalculateReadToCtgSupportExternal {
  # get function parameters
  my ($previous_output_dir) = @_;
  &status("Creating read-to-contig support jobs...");
  my $input_path = &GetOutputPath($previous_output_dir, "data");
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $r2c_cmd = join(" ",
    $python, $EXEC{r2c_submit}, $OPT{lib},
    $input_path, $CFG{alignments}{r2c},
    "--min-overlap", $OPT{r2c_min_overlap},
    "--ctgs-per-job", $OPT{r2c_ctgs_per_job},
  );
  if ($OPT{r2c_no_short_ctgs}) {
    $r2c_cmd .= " --no-short-ctgs";
    $r2c_cmd .= " --read-length ".$OPT{read_length};
  }
  if ($OPT{cluster}) {
    $r2c_cmd .= " --cluster-head ".$OPT{cluster};
    $r2c_cmd .= " --memory ".$OPT{r2c_memory};
    $r2c_cmd .= " --hostname ".$OPT{hostname};
    $r2c_cmd .= " --queue ".$OPT{queue};
    if ($OPT{email}) {
      $r2c_cmd .= " --email ".$OPT{email};
    }
  }
  if ($OPT{disable_profiling_timer}) {
    $r2c_cmd .= " --disable-profiling-timer";
  }
  if ($OPT{debug_high}) {
    $r2c_cmd .= " --debug";
  }
  &status($r2c_cmd);
  my $r2c_result = system($r2c_cmd) >> 8;
  if (0 != $r2c_result) {
    &error("while running ".$r2c_cmd.": ".$r2c_result);
  }
  &status("Finished submitting read-to-contig support jobs.\n");
}

sub OldGetReadToCtgSupport {
  # get function parameters
  my ($groups, $contigs_in_events, $previous_output_dir) = @_;
  # get the event region for each member of each group
  &GetEventRegions($groups);
  # Initialize read support to 0 and
  #  create ctg_id_map because ctg_ids are transformed in bam file
  my %ctg_id_map = ();
  for my $group (@{$groups}) {
    for my $ctg_id (keys %{$group->{members}}) {
      #&debug("Setting read count to 0 for group $group_id: $ctg_id");
      @{$group->{members}{$ctg_id}{num_reads}}        = (0, 0);
      @{$group->{members}{$ctg_id}{num_unique_reads}} = (0, 0);
      # transformed id has + and - replaced with _
      my $bam_ctg_id = $ctg_id;
      $bam_ctg_id =~ s/[+-],?/_/g;
      $bam_ctg_id =~ s/_$//;
      $ctg_id_map{$bam_ctg_id} = $ctg_id;
    }
  }
  &status("Checking read to contig support for events...");
  &ParseBamReadToCtg($groups, $contigs_in_events, \%ctg_id_map);
  &status(""); # a blank newline
  return undef;
}

sub ParseBamReadToCtg {
  # get function parameters
  my ($groups, $contigs_in_events, $ctg_id_map) = @_;
  &SortAndConvertBamToSam();
  #open(R2C_ALIGN, $CFG{read2contig_sam_align_path_rs}) ||
  #  &error("Can't open read-to-contig alignment file: ".
  #      $CFG{read2contig_sam_align_path_rs}." ".$!);
  #&status("Parsing $CFG{read2contig_sam_align_path_rs}...");
  my $convert_start = new Benchmark;
  &debug("Converting read-to-contig bam to sam with:\n  \"".
    $CFG{read2ctg_bam2sam_cmd}."\"");
  open(R2C_ALIGN, $CFG{read2ctg_bam2sam_cmd}."|") ||
    &error("Error while converting read-to-contig bam file to ".
      "sam format with:\n".$CFG{read2ctg_bam2sam_cmd});
  &TimeSpent("converting read to contig file", $convert_start);
  my $parse_r2c_start = new Benchmark;
  #&print_full($logref, "Looping through r2c file\n"); # DEBUG
  my @curr_read_aligns = ();
  my $line_num = 0;
  &status("Parsing read to contig file...");
  while (<R2C_ALIGN>) {
    chomp;
    my @errors_in_line = ();
    $line_num++;
    # track progress
    if (0 == ($line_num % 10000)) {
      my $mode = ">>";
      if (0 == ($line_num % 200000)) {
        $mode = ">";
      }
      my $tracking_path = &GetOutputPath($OPT{outdir}, "r2c.tracking");
      if (open(R2C_TRACKING, $mode, $tracking_path)) {
        my $date_str = `date`;
        print R2C_TRACKING $date_str;
        print R2C_TRACKING "Num lines processed: ".$line_num."\n";
        close R2C_TRACKING;
      }
    }
    # skip header and blank lines
    if ($_ =~ /^@/ || $_ !~ /[^ \t]/) {
      next;
    }
    #&print_full($logref, "  R2C_ALIGN: ".$_."\n"); # DEBUG
    my @align_fields = split(/\t/,$_);
    my $num_align_fields = @align_fields;
    if (4 > $num_align_fields) {
      push(@errors_in_line,
        "Too few fields when parsing sam line: ".$num_align_fields);
    }
    # skip reads that did not align anywhere (check flag)
    my $flag = $align_fields[$SAM_FLAG];
    # test if flag is numeric
    if ($flag !~ /^[0-9]+$/) {
      push(@errors_in_line, "Error parsing flag: \"".$flag."\"");
    }
    # make sure that the read is mapped
    my $read_is_mapped = 0;
    my $new_align;
    if (($flag & 4) == 0) {
      $read_is_mapped = 1;
      my $new_read_id = $align_fields[$SAM_READ_ID];
      # if this is the first alignment for the first read
      if (! @curr_read_aligns) {
        # initialize the alignment list with the read id
        @curr_read_aligns = ($new_read_id);
      }
      # if this is the first alignment for the next read
      if ($new_read_id ne $curr_read_aligns[0]) {
        # examine the alignments for the current read
        &ExamineReadToCtgAlignment
          ($groups, $contigs_in_events, \@curr_read_aligns, $ctg_id_map);
        # initialize the alignment list with the new read id
        @curr_read_aligns = ($new_read_id);
      }
      # add the current alignment to the list for the current read
      $new_align =
        &GetAlignInfo(\@align_fields, $line_num, \@errors_in_line);
    }
    my $num_errors = @errors_in_line;
    if (0 < $num_errors) {
      &warn("Errors in sam line #".$line_num."\n  ".$_);
      &status(join("\n", @errors_in_line));
      &status("Skipping read!");
      if (not $OPT{robust_r2c}) {
        &error("parsing sam file, see log for details");
      }
    } elsif ($read_is_mapped) {
      push @curr_read_aligns, $new_align;
    }
  }
  # examine the alignments for the last read
  &ExamineReadToCtgAlignment
    ($groups, $contigs_in_events, \@curr_read_aligns, $ctg_id_map);
  close R2C_ALIGN;
  &TimeSpent("checking read support", $parse_r2c_start);
  return undef;
}

sub GetAlignInfo {
  # get function parameters
  my ($sam_align_fields, $line_num, $errors_in_line) = @_;
  my $ctg_id = $sam_align_fields->[$SAM_CTG_ID];
  # undo any modifications that may have been done to the contig ID
  my $patt = "^".$OPT{lib}."[^k]*";
  $ctg_id =~ s/$patt//;
  $ctg_id =~ s/^(k[0-9]+)[^:0-9]+/$1:/;
  # test if left is numeric
  if ($sam_align_fields->[$SAM_LEFT] !~ /^[0-9]+$/) {
    #my $sam_line = join("\t", @{$sam_align_fields});
    push(@{$errors_in_line},
      "Error parsing left coordinate: \"".$sam_align_fields->[$SAM_LEFT]."\"");
  }
  # change the left position from 1-based (output by bowtie) to
  #  0-based (output by kaligner), because that is what
  #  ExamineReadToCtgAlignment expects
  my $left = $sam_align_fields->[$SAM_LEFT] - 1;
  my $align_info = $ctg_id." ".$left;
  return $align_info;
}

sub GetEventRegions {
  # get function parameters
  my ($groups) = @_;
  &debug("Looping through groups to get contigs and ".
         "event regions in each group...");
  #for my $group_id (sort {$a<=>$b} keys %{$groups}) {}
  for my $group (@{$groups}) {
    #&debug("Group: ".$group_id."\n"); # DEBUG
    #&print_full($logref, "  Looping through group members to get ".
    #                     "event info\n"); # DEBUG
    for my $ctg_id (keys %{$group->{members}}) {
      #&print_full($logref, "  Ctg ID: ".$ctg_id."\n");
      for my $aligner (keys %{$group->{members}{$ctg_id}{info}}) {
        #&print_full($logref, "    Aligner: ".$aligner."\n");
        my @members =
          @{$group->{members}{$ctg_id}{info}{$aligner}};
        my $num_members = @members; # DEBUG
        #&debug("Num members: ".$num_members."\n"); # DEBUG
        for my $member_info (@members) {
          #&print_full($logref, "    LINE: ".$member_info."\n");
          #Parse the relevant info
          my @member_fields = split(/ /,$member_info);
          #my ($member_ctg_id) = $member_fields[$CC_CTG_ID] =~ /CTG:(\S+)\(/;
          my ($ctg_coords1,$ctg_coords2) =
            $member_fields[$CC_CTG_COORDS] =~ /CONTIG:(\S+),(\S+)/;
          my ($topology) = $member_fields[$CC_TOPOLOGY] =~ /TOPOLOGY:(.*)/;
          my ($ctg_len) = $member_fields[$CC_CTG_ID] =~ /\((\d+)bp\)/;
          if ($topology =~ /gap/) {
            $group->{members}{$ctg_id}{event_coords} =
              &GetGapEventCoords($ctg_coords2, $ctg_len);
          } else {
            #&debug("Setting event coords for group $group_id: $ctg_id");
            $group->{members}{$ctg_id}{event_coords} =
              &GetSplitEventCoords($ctg_coords1, $ctg_coords2, $ctg_len);
              #&GetSplitEventCoords($ctg_coords1, $ctg_coords2, $topology);
          }
          if ($OPT{debug_high}) {
            my $debug_coords =
              $group->{members}{$ctg_id}{event_coords};
            &status("G$group->{id}-$ctg_id $aligner Event coords: ".
                    $debug_coords);
          }
        }
      }
    }
  }
  return undef;
}

sub ExamineReadToCtgAlignment {
  # get function parameters
  my ($groups, $contigs_in_events, $kalign_fields, $ctg_id_map) = @_;
  my $read_id = shift @{$kalign_fields};
  #&print_full($logref, "  READ: ".$read."\n"); # DEBUG
  my $num_alignments = @{$kalign_fields};

  # determine which contigs the read aligns to
  my $aligns_by_contig =
    &GetContigsReadAlignsTo($kalign_fields, $contigs_in_events, $ctg_id_map);

  # check each contig the read aligns to
  my %num_aligns_by_group = ();
  for my $ctg_id (keys %{$aligns_by_contig}) {
    #&print_full($logref, "Aligned Ctg: ".$ctg_id."\n"); # DEBUG
    # check each group the contig is in
    for my $group_id (keys %{$contigs_in_events->{$ctg_id}{groups}}) {
      #&print_full($logref, "  GROUP: ".$group_id."\n"); # DEBUG
      # Get all the contigs that are in the current group
      my %ctgs_in_group = %{$groups->[$group_id]{members}};
      # update the read support for events the current read supports
      &AddSupportToContigs
        ($aligns_by_contig, $num_alignments, $ctg_id,
         \%ctgs_in_group, \$num_aligns_by_group{$group_id});
    }
  }
}

sub GetContigsReadAlignsTo {
  # get function parameters
  my ($kalign_fields, $contigs_in_events, $ctg_id_map) = @_;
  my %aligns_by_contig = ();
  #&print_full($logref, "  Looping through KAlign fields\n"); # DEBUG
  for my $contig_str (@{$kalign_fields}) {
    #&print_full($logref, "    CTG: ".$contig_str."\n"); # DEBUG
    # Contig ids do not always start with k, especially during testing
    my ($ctg_id, $left) = $contig_str =~ /^(\S+)\s+(\d+)/;
    # adjust for left in kaligner being zero-based
    $left += 1;
    if (! defined $ctg_id || ! defined $left) {
      &error("parsing read info from: ".$contig_str." in:\n".
          join(" ", @{$kalign_fields}));
    }
    #&print_full($logref, "        \"$ctg_id\"\n"); # DEBUG
    # translate from bam contig id to real contig id
    if (exists $ctg_id_map->{$ctg_id}) {
      $ctg_id = $ctg_id_map->{$ctg_id};
    }
    #If the contig has an anomaly
    if (exists $contigs_in_events->{$ctg_id}) {
      #&print_full($logref, "    Adding $read to $ctg_id read aligns\n");
      my $right = $left + $OPT{read_length} - 1;
      push @{$aligns_by_contig{$ctg_id}}, "$left-$right";
    } #else { # DEBUG
      #&print_full($logref, "    No $ctg_id in contigs_in_events\n"); # DEBUG
      #} # DEBUG
  }
  return \%aligns_by_contig;
}

sub AddSupportToContigs {
  # get function parameters
  my ($aligns_by_contig, $num_alignments, $ctg_id,
      $ctgs_in_group, $num_aligned_ctgs_in_group) = @_;
  if ("N/A" eq $ctgs_in_group->{$ctg_id}{event_coords}) {
    &warn("Contig $ctg_id has no valid event region");
    return undef;
  }
  #&debug("  Looping through alignments to current contig"); # DEBUG
  #for my $key (sort keys %{$ctgs_in_group->{$ctg_id}}) {
  #  &debug("$key: $ctgs_in_group->{$ctg_id}{$key}");
  #}
  #&debug("Event coords: $ctgs_in_group->{$ctg_id}{event_coords}");
  # Get the event coords
  my ($event_left, $event_right) =
    $ctgs_in_group->{$ctg_id}{event_coords} =~ /^(\d+)-(\d+)/;
  if (! defined $event_left || ! defined $event_right) {
    &error("Error adding support to contig ".$ctg_id.": could not get event ".
      "position from ".$ctgs_in_group->{$ctg_id}{event_coords});
  }
  #&debug("Event left: $event_left Event right: $event_right");
  my ($event_left2, $event_right2);
  if ($ctgs_in_group->{$ctg_id}{event_coords} =~ /;/) {
    ($event_left2, $event_right2) =
      $ctgs_in_group->{$ctg_id}{event_coords} =~ /;(\d+)-(\d+)/;
    if (! defined $event_left2 || ! defined $event_right2) {
      &error("Error adding support to contig ".$ctg_id.": ".
        "could not get second event position from ".
        $ctgs_in_group->{$ctg_id}{event_coords});
    }
  }
  for my $read_align (@{$aligns_by_contig->{$ctg_id}}) {
    #&print_full($logref, "      READ: ".$read_align."\n"); # DEBUG
    my ($read_left, $read_right) = $read_align =~ /(\d+)-(\d+)/;
    if (! defined $read_left || ! defined $read_right) {
      &error("Error adding support to contig ".$ctg_id.
        ": could not get read position from ".$read_align);
    }
    #&debug("Read left: $read_left Read right: $read_right");
    # Check that the read coords cover the event region - if the event
    # region is shorter than the read length, or that the reads falls
    # within the event region - if the event region is longer than the
    # read length
    # REMINDER: this may need to change for events found by
    # the gap filter
    if (($read_left  < $event_left && $event_right < $read_right) ||
        ($event_left < $read_left  && $read_right  < $event_right)) {
      #&print_full($logref, "Read $read_id overlaps event $group_id: ".
      #                     $ctg_id."!\n"); # DEBUG
      $ctgs_in_group->{$ctg_id}{num_reads}[0]++;
      if (&AlignIsUnique
          ($num_alignments, $aligns_by_contig, $ctgs_in_group,
           $num_aligned_ctgs_in_group)) {
        #&print_full($logref, "Adding read to event $group_id: ".$ctg_id."\n");
        $ctgs_in_group->{$ctg_id}{num_unique_reads}[0]++;
      } else {
        #&print_full($logref, "Not unique enough: ".$num_alignments."\n");
      }
    }
    if ($ctgs_in_group->{$ctg_id}{event_coords} =~ /;/) {
      if ($read_left  < $event_left2 && $event_right2 < $read_right) {
        $ctgs_in_group->{$ctg_id}{num_reads}[1]++;
        if (&AlignIsUnique
            ($num_alignments, $aligns_by_contig, $ctgs_in_group,
             $num_aligned_ctgs_in_group)) {
          $ctgs_in_group->{$ctg_id}{num_unique_reads}[1]++;
        }
      }
    }
  }
  return undef;
}

sub AlignIsUnique {
  # get function parameters
  my ($num_alignments, $aligns_by_contig,
      $ctgs_in_group, $num_aligned_ctgs_in_group) = @_;
  # if there is only one alignment, the read is uniquely aligned
  if (1 == $num_alignments) {
    return 1;
  # if there are multiple alignments
  } else {
    # determine the number of contigs the read aligns to that are
    # in the current group
    if (not defined ${$num_aligned_ctgs_in_group}) {
      ${$num_aligned_ctgs_in_group} = 0;
      for my $aligned_ctg_id (keys %{$aligns_by_contig}) {
        if (exists $ctgs_in_group->{$aligned_ctg_id}) {
          ${$num_aligned_ctgs_in_group}++;
        }
      }
    }
    # check that the alignment is 'unique' enough
    # i.e. >50% of the contigs aligned to are in the current group
    my $group_fraction =
      ${$num_aligned_ctgs_in_group}/$num_alignments;
    $msg  = "# Aligned Contigs in group: ".
            ${$num_aligned_ctgs_in_group};
    #&print_full($logref, "$msg\n");
    #&print_full($logref, "# Alignments: ".$num_alignments."\n");
    #&print_full($logref, "Fraction: ".$group_fraction."\n");
    if ($group_fraction > 0.5) {
      return 1;
    }
  }
  return 0;
}

sub CompareMembers {
  # get function parameters
  my ($contigs_in_events, $ctg_id1, $ctg_id2) = @_;

  # Compare by contig lengths
  my ($ctg_len1) = $contigs_in_events->{$ctg_id1}{len};
  my ($ctg_len2) = $contigs_in_events->{$ctg_id2}{len};
  my $result = ($ctg_len2 <=> $ctg_len1);

  if (0 == $result) {
    # Compare by contig name
    $result = ($ctg_id1 cmp $ctg_id2);
  }
  return $result;
}

sub FlagExonBoundaryJunctions {
  # get function parameters
  my ($previous_output_dir) = @_;
  &status("Flagging junctions that match exon-boundaries...");
  my $flag_ebj_start = new Benchmark;
  my $input_path = &GetOutputPath($previous_output_dir, "data");
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $flag_ebj_cmd = join(" ",
    $python, $EXEC{exon_bounds_script}, $OPT{lib}, $input_path,
    $CFG{annotations}{genes}, "--buffer-size", $OPT{exon_bound_buffer},
  );
  if ($OPT{debug_high}) {
    $flag_ebj_cmd .= " --debug";
  }
  &status($flag_ebj_cmd);
  my $flag_ebj_result = system($flag_ebj_cmd) >> 8;
  if (0 != $flag_ebj_result) {
    &error("while running ".$flag_ebj_cmd.": ".$flag_ebj_result);
  }
  # construct two hashes to hold the exon coordinates
  #my $exons = &GetExonBoundaries();
  #for my $chrom (sort keys %{$exons}) {
  #  &print_full($logref, "$chrom\n");
  #}
  # iterate through the groups
  #&status("Checking groups for junctions matching exon-boundaries...");
  #for my $group_id (sort {$a<=>$b} keys %{$groups}) {}
  #for my $group (@{$groups}) {
  #  # initialize the number of exon-matching junctions to zero
  #  $group->{exon_junct} = 0;
  #  # use {ordered_target_regions} to check
  #  my $targets = $group->{ordered_target_regions};
  #  if (&DoesJunctionMatchExonBoundary ($targets->{A}, $exons)) {
  #    $group->{exon_junct}++;
  #  }
  #  if (&DoesJunctionMatchExonBoundary ($targets->{B}, $exons)) {
  #    $group->{exon_junct}++;
  #  }
  #}
  &TimeSpent("flagging exon-boundary junctions", $flag_ebj_start);
  return undef;
}

#sub DoesJunctionMatchExonBoundary {
#  # get function parameters
#  my ($region, $exons) = @_;
#  my $chrom = $region->{chrom};
#  if (! exists $exons->{$chrom}) {
#    &warn("no exons on chromosome $chrom");
#    return 0;
#  }
#  my $junction_side = $region->{junction};
#  #&debug("J-side: $junction_side");
#  if ("gap" ne $junction_side) {
#    my $junction_coord = $region->{$junction_side};
#    #&debug("J-coord: $junction_coord");
#    for my $offset (-3..3) {
#      my $test_coord = $junction_coord + $offset;
#      #&debug("Offset: $offset Test: $test_coord");
#      if (exists $exons->{$chrom}{$junction_side}{$test_coord}) {
#        #&debug("Exon boundary found");
#        return 1;
#      }
#    }
#  }
#  return 0;
#}

#sub GetExonBoundaries {
#  my %exons = ();
#  &status("Getting exon boundaries...");
#  # get exon boundaries from the primary exon coordinates file
#  &GetExonBoundariesFromFile($CFG{coords}{exons}, \%exons);
#  #if (exists $CFG{coords}{exons}{secondary}) {
#  #  # get exon boundaries from the secondary exon coordinates files
#  #  for my $coords_file (@{$CFG{coords}{exons}{secondary}}) {
#  #    &GetExonBoundariesFromFile($coords_file, \%exons);
#  #  }
#  #}
#  return \%exons;
#}

#sub GetExonBoundariesFromFile {
#  # get function parameters
#  my ($exon_coords_file, $exons) = @_;
#  &debug("Exon file: $exon_coords_file");
#  open(EXON_COORDS, $exon_coords_file) ||
#    &error("Can't open exon coordinates file ".$exon_coords_file." ".$!);
#  # get the coordinates of an exon from the file
#  while (<EXON_COORDS>) {
#    # skip blank lines
#    next if ($_ !~ /[^ \t\n\r]/);
#    chomp;
#    my ($chrom, $left, $right, $name) = split;
#    if ($chrom !~ /./) {
#      &error("Cannot get exon coordinates from: ".$_.
#        "\n  ".$exon_coords_file);
#    }
#    # remove the "chr" from the chromosome ids
#    $chrom =~ s/chr//;
#    #&debug("Chromosome: $chrom");
#    if ($chrom =~ /^([0-9]+|[XYM])$/) {
#      push @{$exons->{$chrom}{left}{$left}{$right}}, $name;
#      push @{$exons->{$chrom}{right}{$right}{$left}}, $name;
#    }
#  }
#  close EXON_COORDS;
#  return undef;
#}

sub RemoveUnparsableEvents {
  # no parameters
  &status("Removing unparsable events...");
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $raw_path = &GetOutputPath($OPT{outdir}, "data");
  my $check_events_cmd = join(" ",
    $python, $EXEC{check_format}, $OPT{lib}, $raw_path);
  if ($OPT{debug_high}) {
    $check_events_cmd .= " --debug";
  }
  &status($check_events_cmd);
  my $check_events_result = system($check_events_cmd) >> 8;
  if (0 != $check_events_result) {
    &error("while running ".$check_events_cmd.": ".$check_events_result);
  }
  &status("Finished removing unparsable events.\n");
}

sub CheckForRepeatsExternal {
  # get function parameters
  my ($previous_output_dir) = @_;
  &status("Checking for repeats with repeats script...");
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $input_path = &GetOutputPath($previous_output_dir, "data");
  my $flag_repeats_cmd = join(" ",
    $python, $EXEC{repeats_script}, $OPT{lib}, $input_path,
    join(",", @{$CFG{coords}{repeats}}),
    "--region-size", $OPT{repeat_search_size},
    "--repeat-filter-off"
  );
  if ($OPT{no_flag_RNA}) {
    $flag_repeats_cmd .= " --no-structural-RNA-check";
  }
  if ($OPT{disable_profiling_timer}) {
    $flag_repeats_cmd .= " --disable-profiling-timer";
  }
  if ($OPT{debug_high}) {
    $flag_repeats_cmd .= " --debug";
  }
  &status($flag_repeats_cmd);
  my $flag_repeats_result = system($flag_repeats_cmd) >> 8;
  if (0 != $flag_repeats_result) {
    &error("while running ".$flag_repeats_cmd.": ".$flag_repeats_result);
  }
}

sub MarkBreakpointGenes {
  # get function parameters
  my ($previous_output_dir) = @_;
  &status("Marking breakpoint genes...");
  my $python = $CFG{commands}{python};
  if ($OPT{cluster} && exists $CFG{commands}{python_cluster}) {
    $python = $CFG{commands}{python_cluster};
  }
  my $input_path = &GetOutputPath($previous_output_dir, "data");
  my $breakpoint_genes_cmd = join(" ",
    $python, $EXEC{breakpoint_genes_script}, $OPT{lib}, $input_path,
    $CFG{coords}{gene_features}
  );
  if ($OPT{disable_profiling_timer}) {
    $breakpoint_genes_cmd .= " --disable-profiling-timer";
  }
  if ($OPT{debug_high}) {
    $breakpoint_genes_cmd .= " --debug";
  }
  &status($breakpoint_genes_cmd);
  my $breakpoint_genes_result = system($breakpoint_genes_cmd) >> 8;
  if (0 != $breakpoint_genes_result) {
    &error("while running ".$breakpoint_genes_cmd.": ".$breakpoint_genes_result);
  }
}

# Summarize the results, this is where the filters are checked and the
# pass/fail status is assigned to each event
sub Summary {
  # get function parameters
  my ($groups, $contigs_in_events, $fasta_data) = @_;
  #my $num_groups = keys %{$groups};
  my $num_groups = @{$groups};
  &status("$num_groups potential events found.");
  &status("Applying filters...");
  my $out_all_events = &GetOutputPath($OPT{outdir}, "data");
  open(ALL_OUTPUT, ">$out_all_events") ||
    &error("Can't open final output file ".$out_all_events." ".$!);
  #my $out_passing_events = &GetOutputPath($outdir, "pass");
  #open(PASS_OUTPUT, ">$out_passing_events") ||
  #  &error("Can't open pass output file ".$out_passing_events." ".$!);
  #my $type_files = &SetupTypeFiles();
  my $num_passing_events = 0;
  #for my $group_id (sort {$a<=>$b} keys %{$groups}) {}
  for my $group (@{$groups}) {
    $num_passing_events +=
      &SummarizeGroup($group, $contigs_in_events, $fasta_data);
  }
  &status("$num_passing_events events passed all filters.");
  # close the files
  close ALL_OUTPUT;
  #close PASS_OUTPUT;
  #if ($OPT{split_out}) {
  #  for my $type (keys %{$type_files}) {
  #    close $type_files->{$type}{handle};
  #    if (not $type_files->{$type}{found}) {
  #      #&print_full($logref, "Removing $type_files->{$type}{file_name}\n");
  #      unlink($type_files->{$type}{file_name});
  #      #system("rm $type_files->{$type}{file_name}");
  #    }
  #  }
  #}
  return undef;
}

sub SummarizeGroup {
  # get function parameters
  my ($group, $contigs_in_events, $fasta_data) = @_;
  # dictionary for storing various "group_data"
  my $group_data = &InitializeGroupData();
  # summarize the group members of the current group
  &SummarizeGroupMembers
    ($group, $group_data, $contigs_in_events);
  #Check if we have RNA repeats
  my $rna_status = "N/A";
  #my $rna_status = "N";
  #if ($OPT{no_flag_RNA}) {
  #  $rna_status = "N/A";
  #} elsif (exists $group->{rna}) {
  #  $group_data->{fail_reasons}{'RNA_REPEAT'} = 1;
  #  $rna_status = "Y";
  #}
  # construct the genomic coordinates string
  my $targets = $group->{ordered_target_regions};
  # construct the full coordinates string
  my $coord_str = &GetTargetCoordsStr($targets, $group_data->{all_gap});
  # get the repeat overlap information
  #my $repeats_str = &GetRepeats($targets);
  # Apply the mitochondrial DNA filter
  &MitochondrialFilter($targets, \%{$group_data->{fail_reasons}});
  # Construct the event topologies string
  my $topologies = join(",", sort keys %{$group->{topologies}});
  # Construct the string containing the lengths of the contigs in the group
  my $ctg_len_str = &ContigLensString($group->{ctg_lens});
  # Apply the read-pair to genome filter
  my ($pair_support, $filtered_pair_support) = #("N/A", "N/A");
    &PairToGenomeFilter($group, \%{$group_data->{fail_reasons}});
  # Check if we have one or both aligners
  #my $num_aligners =
  #  &GetNumAligners($group->{aligners}, $topologies,
  #                  \%{$group_data->{fail_reasons}});
  # Count the number of junction coordinates matching exon boundaries
  my $exon_matching_junctions = "N/A";
  #unless ($OPT{no_flag_exon_boundary_junctions}) {
  #  $exon_matching_junctions = $group->{exon_junct}
  #}
  my $status_str = '';
  if (keys %{$group_data->{fail_reasons}}) {
    $status_str = 'FAIL('.join(",",keys %{$group_data->{fail_reasons}}).')';
  } else {
    $status_str = 'PASS';
  }
  # count the number of contigs appearing in the group
  #&print_full($logref, "Counting contigs for group $group_id\n"); # DEBUG
  my $num_contigs = 0;
  for my $k_value (keys %{$group_data->{k_values}}) {
    #&print_full($logref, "  K Value: ".$k_value."\n"); # DEBUG
    $num_contigs += keys %{$group_data->{k_values}{$k_value}};
    #for my $ctg (keys %{$group_data->{k_values}{$k_value}}) {
    #  &print_full($logref, "    Ctg: ".$ctg."\n"); # DEBUG
    #}
  }
  #&print_full($logref, "----\n"); # DEBUG
  my $num_group_members = @{$group_data->{members}};
  my $group_str = join(" ",
    "GRPNUM:"     .$group->{id},
    "TOPOLOGIES:" .$topologies,
    "COORDS:"     .$coord_str,
    "MEMBERS:"    .$num_group_members,
    "CONTIGS:"    .$num_contigs,
    #"ALIGNERS:"   .$num_aligners,
    "CLENS:"      .$ctg_len_str,
    "PAIR_TO_GENOME_0:".$pair_support,
    "PAIR_TO_GENOME_".$OPT{min_mapq}.":".$filtered_pair_support,
    #"EXON_BOUNDS:".$exon_matching_junctions,
    "RNA:"        .$rna_status,
    "STATUS:"     .$status_str
  );
  $group_str .= " LIB:".$OPT{lib};
  if (defined $OPT{ver}) {
    $group_str .=" VER:".$OPT{ver};
  }
  if ($group->{p2g_error}) {
    $group_str .= " NOTE:error_getting_pair_to_genome_support";
  }
  &PrintSummary($group_data, $group_str, $status_str, $fasta_data);
  #$group_data, $group_str, $status_str, $fasta_data, $group->{topologies});
  if ($status_str =~ /PASS/) {
    return 1;
  }
  return 0;
}

sub InitializeGroupData {
  # no parameters
  my %group_data = (
    #  group_data{members}     => @(),
    #  group_data{contigs}     => %(),
    # Represents the local status of the group ->
    # a single local contig changes the group status to local
    islocal => 'N',
    # assume that all the members of the group are gapped events
    all_gap => 1,
  );
  # record all the reasons the event may fail
  %{$group_data{fail_reasons}} = ();
  #Set this to fail by default until we find a group member that passes
  if ($OPT{max_num_groups}) {
    $group_data{fail_reasons}{'TOO_MANY_GROUPS'} = 1;
  }
  #Set this to fail by default until we find a group member that passes
  #unless ($OPT{no_read_to_contig} || $OPT{cluster}) {
  #  $group_data{fail_reasons}{'READS_TO_CONTIG'} = 1;
  #}
  return \%group_data;
}

sub GetTargetCoordsStr {
  # get function parameters
  my ($targets, $all_gap) = @_;
  if ($all_gap) {
    my %gap_coords = (
      chrom => $targets->{A}{chrom},
      left  => $targets->{A}{left},
      right => $targets->{B}{right},
    );
    return "chr".&RegionCoordsStr(\%gap_coords);
  } else {
    my $coordA = "chr".&RegionCoordsStr($targets->{A});
    my $coordB = "chr".&RegionCoordsStr($targets->{B});
    return "$coordA;$coordB";
  }
}

sub GetRepeats {
  # get function parameters
  my ($targets) = @_;
  my ($repeatsA, $repeatsB) =("N/A","N/A");
  if ($OPT{flag_repeats}) {
    my $num_repeatsA = keys %{$targets->{A}{repeats}};
    if (0 < $num_repeatsA) {
      $repeatsA = join(",", sort keys %{$targets->{A}{repeats}});
    } else {
      $repeatsA = "None";
    }
    my $num_repeatsB = keys %{$targets->{B}{repeats}};
    if (0 < $num_repeatsB) {
      $repeatsB = join(",", sort keys %{$targets->{B}{repeats}});
    } else {
      $repeatsB = "None";
    }
  }
  return join(";", $repeatsA, $repeatsB);
}

sub MitochondrialFilter {
  # get function parameters
  my ($targets, $fail_reasons) = @_;
  if ($OPT{noMito}) {
    if ("M" eq $targets->{A}{chrom} || "M" eq $targets->{B}{chrom}) {
      $fail_reasons->{'MITOCHONDRIAL'} = 1;
    }
  }
  return undef;
}

sub ContigLensString {
  # get function parameters
  my ($ctg_lens) = @_;
  return join(",", sort {
    my ($a_num) = $a =~ /(\d+)/;
    my ($b_num) = $b =~ /(\d+)/;
    $b_num <=> $a_num
  } map { $_."bp" } keys %{$ctg_lens});
}

#sub GetNumAligners {
#  # get function parameters
#  my ($aligners, $topologies, $fail_reasons) = @_;
#  my $num_aligners = 0;
#  if (exists $aligners->{exonerate}) {
#    $num_aligners++;
#  }
#  if (exists $aligners->{"split"} ||
#      exists $aligners->{"gap"}) {
#    $num_aligners++;
#  }
#  if (2 > $num_aligners && $OPT{require_both_aligners} && $topologies !~ /gap/) {
#    $fail_reasons->{'SINGLE_ALIGN'} = 1;
#  }
#  return $num_aligners;
#}

sub PairToGenomeFilter {
  # get function parameters
  my ($group_info, $fail_reasons) = @_;
  my $pair_support = "N/A";
  my $filtered_pair_support = "N/A";
  unless ($OPT{no_pair_to_genome} || $OPT{cluster_pair2gen}) {
    &debug_high("Setting PAIR SUPPORT...");
    $pair_support = $group_info->{pair_support};
    $filtered_pair_support = $group_info->{filtered_pair_support};
    &debug_high("  Unfiltered: ".$pair_support.", Filtered: ".
      $filtered_pair_support);
    #&debug("Pair support: $group_info->{pair_support} ".
    #       "Filtered: $group_info->{filtered_pair_support}.");
    #&debug("Pair support: $pair_support ".
    #       "Filtered: $filtered_pair_support.");
    #&debug("Min: ".$OPT{min_filt_pairs});
    #if ($OPT{min_pairs}      > $pair_support ||
    #    $OPT{min_filt_pairs} > $filtered_pair_support) {
    #  #&debug("Adding fail reason: PAIR_SUPPORT");
    #  $fail_reasons->{'PAIR_SUPPORT'} = 1;
    #}
    if (0 < $group_info->{intronic_pair_support}) {
      my $nointron_pair_support =
        $pair_support - $group_info->{intronic_pair_support};
      &debug_high("  Unfiltered intronic: ".
        $group_info->{intronic_pair_support}.", exonic: ".
        $nointron_pair_support);
      $pair_support = $nointron_pair_support."(".$pair_support.")";
    }
    if (0 < $group_info->{intronic_filtered_pair_support}) {
      my $nointron_filt_pair_support =
        $filtered_pair_support -
        $group_info->{intronic_filtered_pair_support};
      &debug_high("  Filtered intronic: ".
        $group_info->{intronic_filtered_pair_support}.", exonic: ".
        $nointron_filt_pair_support);
      $filtered_pair_support = $nointron_filt_pair_support.
        "(".$filtered_pair_support.")";
    }
  }
  &debug_high("    Unfiltered: ".$pair_support.", filtered: ".
    $filtered_pair_support);
  return ($pair_support, $filtered_pair_support);
}

sub SummarizeGroupMembers {
  # get function parameters
  my ($group, $group_data, $contigs_in_events) = @_;
  # iterate over the group's lines (each line represents
  # one candidate contig alignment pair)
  my $member_id = 'a';
  for my $ctg_id (sort { CompareMembers($contigs_in_events, $a, $b) }
                            keys %{$group->{members}}) {
    #&print_full($logref, "Ctg ID: ".$ctg_id."\n");
    if (! exists $group->{members}) {
      &status("Group{members} does not exist.");
    } elsif (! exists $group->{members}{$ctg_id}) {
      &status("Group{members}{$ctg_id} does not exist.");
    } elsif (! exists $group->{members}{$ctg_id}{info}) {
      &status("Group{members}{$ctg_id}{info} does not exist.");
    }
    my %curr_ctg = %{$group->{members}{$ctg_id}{info}};
    for my $aligner (sort { $a cmp $b } keys %curr_ctg) {
      for my $member_info (@{$curr_ctg{$aligner}}) {
        # Check if the member is part of a local alignment
        my $islocal = 'N';
        # Set the group local status if this member is local
        #if ($group->{members}{$ctg_id}{islocal}) {
        #  $islocal = 'Y';
        #  $group_data->{islocal} = 'Y';
        #  $group_data->{fail_reasons}{'LOCAL'} = 1 if $OPT{nolocal};
        #}
        my $member_summary =
          &SummarizeMemberInfo($member_info, $aligner, $islocal,
                               $group->{neighbourness}, $group->{members},
                               $contigs_in_events, $group_data);
        push @{$group_data->{members}},
          "$group->{id}$member_id) $member_summary";
        $member_id++;
      }
    }
  }
  return undef;
}

sub PrintSummary {
  # get function parameters
  my ($group_data, $group_str, $status_str, $fasta_data) = @_;
  #my ($group_data, $group_str, $status_str, $fasta_data, $topologies) = @_;

  my $group_info = "$group_str\n";
  $group_info .= join("\n",@{$group_data->{members}})."\n\n";
  print ALL_OUTPUT $group_info;
  if ($status_str =~ /PASS/) {
    #Store the contig for generating the fasta later
    for my $ctg_id (keys %{$group_data->{contigs}}) {
      my $contig_type = $group_data->{contigs}{$ctg_id};
      $fasta_data->{$contig_type}{ctgs}{$ctg_id} = 1;
    }
    #print PASS_OUTPUT $group_info;
    #if ($OPT{split_out}) {
    #  my $event_topology;
    #  my $num_topologies = keys %{$topologies};
    #  if (1 eq $num_topologies) {
    #    ($event_topology) = keys %{$topologies};
    #  } else {
    #    $event_topology = "multi";
    #  }
    #  #&print_full($logref, "Event topology: ".$event_topology."\n");
    #  if (exists $topology_files->{$event_topology}) {
    #    $topology_files->{$event_topology}{found} = 1;
    #    print {$topology_files->{$event_topology}{handle}} $group_info;
    #  } else {
    #    &warn("No file specified for event with topology: $event_topology");
    #  }
    #}
  }
  return undef;
}

sub SummarizeMemberInfo {
  # get function parameters
  my ($member_info, $aligner, $islocal, $neighbourness,
      $ctgs_in_group, $contigs_in_events, $group_data) = @_;
  #&print_full($logref, "$member_info\n"); # DEBUG
  #Parse the relevant info
  my @member_fields = split(/ /, $member_info);
  my ($ctg_str) = $member_fields[$CC_CTG_ID] =~ /CTG:(\S+)/;
  my ($ctg_id) = $member_fields[$CC_CTG_ID] =~ /CTG:(\S+)\(/;
  my ($topology) = $member_fields[$CC_TOPOLOGY] =~ /TOPOLOGY:(\S+)/;
  if ($topology !~ /gap/) {
    $group_data->{all_gap} = 0;
  }
  # Store the ids of the contigs in the current group
  my $contig_type = $member_fields[-1];
  $group_data->{contigs}{$ctg_id} = $contig_type;
  # get event co-ordinates
  my ($query1,$query2) =
    $member_fields[$CC_CTG_COORDS] =~ /CONTIG:(\S+),(\S+)/;
  my ($target1,$target2) =
    $member_fields[$CC_TARGET] =~ /TARGET:(\S+),(\S+)/;
  my ($pid1,$pid2,$af1,$af2) =
    $member_fields[$CC_ALIGN_METRICS] =~ /I1:(.*),I2:(.*),AF1:(.*),AF2:(.*)/;
  #Record the k value
  #&print_full($logref, "Getting K Value: $member_fields[$CC_CTG_ID]\n");
  my $k_value = 0;
  if ($contig_type =~ /k(\d+)/) {
    $k_value = $1;
  } elsif ($member_fields[$CC_CTG_ID] =~ /CTG:k(\d+)/) {
    ($k_value) = $1;
  }
  #&print_full($logref, "K Value: ".$k_value."\n");
  $group_data->{k_values}{$k_value}{$ctg_id}++;
  # Determine the contig length as a function of the form 2k+x
  my ($ctg_len) = $ctg_str =~ /\((\d+)bp\)/;
  my $two_k = 2 * $k_value;
  my $x_val = $ctg_len - $two_k;
  my $k_form;
  if (0 > $x_val) {
    $k_form = "2k$x_val";
  } else {
    $k_form = "2k+$x_val";
  }
  $ctg_str =~ s/bp\)/bp:$k_form\)/;
  unless ("contigs" eq $contig_type) {
    $ctg_str = $contig_type.":".$ctg_str;
  }
  #&print_full($logref, "K = $k_value, 2K = $two_k, X = $x_val\n"); # DEBUG
  # If we find a contig that does not appear in too many events,
  # remove the fail status
  my $num_groups = keys %{$contigs_in_events->{$ctg_id}{groups}};
  if ($OPT{max_num_groups} && $num_groups <= $OPT{max_num_groups}) {
    delete $group_data->{fail_reasons}{'TOO_MANY_GROUPS'};
  }
  # Check read to contig alignment support
  #my ($read_support, $unique_read_support) =
  #  ReadToContigFilter($ctgs_in_group->{$ctg_id},
  #                     \%{$group_data->{fail_reasons}});
  # check whether the alignments overlap or are near to any genes
  my $overlapping_genes = "N/A;N/A";
  my $nearby_genes = "N/A;N/A";
  if ($topology =~ /gap/) {
    $overlapping_genes = "N/A";
    $nearby_genes = "N/A";
  }
  if (defined $OPT{add_gene_annotation}) {
    #&debug($member_fields[$CC_GENES]);
    ($overlapping_genes) =
      $member_fields[$CC_GENES] =~ /GENES:(.*)/;
    # REMINDER: check for "intronic"
    # older versions of gap output do not have GENES field
    if ($topology =~ /gap/ && ! $overlapping_genes) {
      $overlapping_genes = "N/A";
    }
    ($nearby_genes) =
      $member_fields[$CC_NEARBY] =~ /NEARBY:(.*)/;
    # older versions do not have NEARBY field
    if (! $nearby_genes) {
      if ($topology =~ /gap/) {
        $nearby_genes = "N/A";
      } else {
        $nearby_genes = "N/A;N/A";
      }
    }
  }
  # construct the member summary string
  #&debug("TOPOLOGY: ".$topology);
  #&debug("CONTIG: ".$contig_type.":".$ctg_str);
  #&debug("TARGET: ".$event_coords);
  #&debug("ALIGNER:$aligner,$neighbourness");
  #&debug("READ_TO_CTG:".$read_support);
  #&debug("READ_TO_CTG_UNIQUE:".$unique_read_support);
  #&debug("NUM_GROUPS:".$num_groups);
  #&debug("OVERLAPPING_GENES:".$overlapping_genes);
  #&debug("META: ".$member_fields[$CC_META]);
  #&debug("JUNCTION: ".$member_fields[$CC_JUNCTION]);
  my $member_summary = join(" ",
    "TOPOLOGY:"          .$topology,
    "CONTIG:"            .$ctg_str,
    "ALIGN_A:"           .&GetMemberCoordStr($query1, $target1, $af1, $pid1),
    "ALIGN_B:"           .&GetMemberCoordStr($query2, $target2, $af2, $pid2),
    "ALIGNER:"           .$aligner,
    "READ_TO_CTG:N/A",       #.$read_support,
    "READ_TO_CTG_UNIQUE:N/A",#.$unique_read_support,
    "NUM_GROUPS:"        .$num_groups,
    "OVERLAPPING_GENES:" .$overlapping_genes,
    "NEARBY_GENES:"      .$nearby_genes,
    "REPEATS:N/A;N/A",
    $member_fields[$CC_META],
    $member_fields[$CC_JUNCTION],
    $member_fields[$CC_BLOCKS]
  );
  if ($topology =~ /gap/) {
    if ($member_fields[$CC_SEQ] =~ /EVENT_SEQ:/) {
      $member_summary .= " ".$member_fields[$CC_SEQ];
    }
    if ($member_fields[$CC_SEQ + 1] =~ /INS_SEQ:/) {
      $member_summary .= " ".$member_fields[$CC_SEQ + 1];
    }
  }
  return $member_summary;
}

sub GetMemberCoordStr {
  # get function parameters
  my ($query, $target, $af, $pid) = @_;
  my ($qstart,$qend) = $query =~ /(\d+)-(\d+)/;
  my $qspan = abs($qstart - $qend) + 1;
  my ($tstart,$tend) = $target =~ /:(\d+)-(\d+)/;
  my $tspan = abs($tstart - $tend) + 1;
  my $strand = "+";
  if ($tend < $tstart) {
    $strand = "-";
  }
  my $member_coord_str = "ctg:".$query."(".${qspan}."bp)=".$target.
    "(".${tspan}."bp,".$strand.");AF:".$af.",PID:".$pid;
  return $member_coord_str;
}

sub ReadToContigFilter {
  # get function parameters
  my ($ctg_info, $fail_reasons) = @_;
  my $read_support = "N/A";
  my $unique_read_support = "N/A";
  unless ($OPT{no_read_to_contig} || $OPT{cluster}) {
    my ($event_left, $event_right) =
      $ctg_info->{event_coords} =~ /^(\d+)-(\d+)/;
    if (! defined $event_left || ! defined $event_right) {
      &error("getting event coordinates for contig \"".$ctg_info->{id}.
        "\": \"".$ctg_info->{event_coords}."\"");
    }
    $read_support        = $ctg_info->{num_reads}[0];
    $unique_read_support = $ctg_info->{num_unique_reads}[0];
    if ($ctg_info->{event_coords} =~ /;/) {
      if ($ctg_info->{num_reads}[1] < $read_support) {
        $read_support = $ctg_info->{num_reads}[1];
      }
      if ($ctg_info->{num_unique_reads}[1] <
          $unique_read_support) {
        $unique_read_support = $ctg_info->{num_unique_reads}[1];
        ($event_left, $event_right) =
          $ctg_info->{event_coords} =~ /;(\d+)-(\d+)/;
        if (! defined $event_left || ! defined $event_right) {
          &error("getting second event coordinates for contig ".
            $ctg_info->{id}.": ".$ctg_info->{event_coords});
        }
      }
    }
    #if (($OPT{consider_all_reads} &&
    #     $read_support        >= $OPT{min_read_coverage}) ||
    #    ($unique_read_support >= $OPT{min_unique_read_coverage})) {
    #  delete $fail_reasons->{'READS_TO_CONTIG'};
    #}
    my $event_len = $event_right - $event_left + 1;
    my $norm_read_supp = "NaN";
    if (0 < $event_len) {
      $norm_read_supp = $unique_read_support / $event_len;
      $norm_read_supp = sprintf("%.2f", $norm_read_supp);
    }
    $unique_read_support .= "(${event_len}bp:$norm_read_supp)";
  }
  return ($read_support, $unique_read_support);
}

sub GetPairSupport {
  # get function parameters
  my ($pair_info) = @_;

  &debug_high("GRPNUM: $pair_info->{group_id} Contigs: $pair_info->{ctg_ids}");
  &debug_high("Topology: $pair_info->{topology}");
  # use 125% of the frag_len, since the fragment length is not exact
  my $region_len = 1.25 * $OPT{frag_len};
  if ("N/A" ne $pair_info->{ctg_overlap} && 0 < $pair_info->{ctg_overlap} ) {
    $region_len += $pair_info->{ctg_overlap};
  }
  &debug_high("Region Length: $region_len");

  # get the coordinates of the upstream region
  my %region1 = (
    chrom => $pair_info->{coords}[0]{chrom},
  );
  if ($CFG{read_support}{p2g}{use_chr}) {
    $region1{chrom} = "chr".$region1{chrom};
  }
  if ($CFG{read_support}{p2g}{use_MT}) {
    $region1{chrom} =~ s/M$/MT/;
  }
  #if ($pair_info->{topology} =~ /gap/) {
  #  &print_full($logref, "Not using blocks from gap event\n");
  #} else {
  #  &print_full($logref, "Blocks:  $pair_info->{coords}[0]{blocks}\n");
  #}
  ($region1{left}, $region1{right}) =
    &GetRegionToSearch
      ($pair_info->{coords}[0], $pair_info->{topology}, $region_len);
  &debug_high(join(" ",
    "Region1:",
    &RegionCoordsStr(\%region1),
    $pair_info->{coords}[0]{junction})
  );
  if ($region1{left} > $region1{right}) {
    &warn("Bad region1: \"".$region1{left}."-".$region1{right}.
          "\". For group ".$pair_info->{group_id}." with contigs: ".
          $pair_info->{ctg_ids});
    $pair_info->{p2g_error} = 1;
  }

  # get the coordinates of the downstream region
  my %region2 = (
    chrom => $pair_info->{coords}[1]{chrom},
  );
  if ($CFG{read_support}{p2g}{use_chr}) {
    $region2{chrom} = "chr".$region2{chrom};
  }
  if ($CFG{read_support}{p2g}{use_MT}) {
    $region2{chrom} =~ s/M$/MT/;
  }
  #&print_full($logref, "FC2:  $full_coords2\n");
  #if ($pair_info->{topology} =~ /gap/) {
  #  &print_full($logref, "Not using blocks from gap event\n");
  #} else {
  #  &print_full($logref, "Blocks:  $pair_info->{coords}[1]{blocks}\n");
  #}
  ($region2{left}, $region2{right}) =
    &GetRegionToSearch
      ($pair_info->{coords}[1], $pair_info->{topology}, $region_len);
  &debug_high(join(" ",
    "Region2:",
    &RegionCoordsStr(\%region2),
    $pair_info->{coords}[1]{junction})
  );
  if ($region2{left} > $region2{right}) {
    &warn("Bad region2: \"".$region2{left}-$region2{right}.
          "\". For group $pair_info->{group_id} with contigs: ".
          $pair_info->{ctg_ids});
    $pair_info->{p2g_error} = 1;
  }
  &debug_high("CO:      $pair_info->{ctg_overlap}");

  my $what_to_sam = $SAM_BOTH;
  my %regionA = (
    chrom => $region1{chrom},
    left  => $region1{left},
    right => $region1{right},
  );
  my %regionB = (
    chrom => $region2{chrom},
    left  => $region2{left},
    right => $region2{right},
  );
  if ($pair_info->{topology} =~ /gap/) {
    ##PYTHON: AdjustGapRegions()
    # region2 left and right are both the right side of the gap
    $regionA{right} = $region2{right};
    # region1 left and right are both the left side of the gap
    $regionB{left}  = $region1{left};
    # add the blocks information from the contig to genome alignment
    if ($pair_info->{topology} =~ /duplication/) {
      $regionB{blocks} = $pair_info->{coords}[0]{blocks};
    } else {
      $regionB{blocks} = &MergeBlocks
        ($pair_info->{coords}[0]{blocks},$pair_info->{coords}[1]{blocks});
    }
  } elsif ("left" eq $pair_info->{coords}[0]{junction}) {
    %regionA = (
      chrom => $region2{chrom},
      left  => $region2{left},
      right => $region2{right},
    );
    %regionB = (
      chrom => $region1{chrom},
      left  => $region1{left},
      right => $region1{right},
    );
    $what_to_sam = $SAM_DOWNSTREAM;
  } elsif ("left" eq $pair_info->{coords}[1]{junction}) {
    $what_to_sam = $SAM_UPSTREAM;
  }
  #&debug("What to SAM: $what_to_sam");

  # run the samtools view command to get the reads aligning to
  # the chosen alignment region
  my $sam_region = &RegionCoordsStr(\%regionA);
  &RunSamTools($sam_region, $pair_info->{ctg_ids});

  # Parse the results file from the samtools view command
  my ($pair_support, $intronic_pair_support,
      $filtered_pair_support, $intronic_filtered_pair_support);
  ($pair_support, $intronic_pair_support,
   $filtered_pair_support, $intronic_filtered_pair_support,
   $pair_info->{num_reads}) =
    &ParseSamFile($what_to_sam, \%regionB,
                  $pair_info->{ctg_ids}, $pair_info->{topology});
  &debug_high("# Reads: $pair_info->{num_reads}, # Support: $pair_support");
  return ($pair_support, $intronic_pair_support,
          $filtered_pair_support, $intronic_filtered_pair_support);
}

sub ParseSamFile {
  # get function parameters
  my ($what_to_sam, $region, $ctg_ids, $topology) = @_;
  my $num_reads = 0;
  my %potential_pairs = ();
  open(SAMFILE_PRIMARY, $CFG{read_support}{p2g}{sam_out});
  while (<SAMFILE_PRIMARY>) {
    # skip blank lines
    next if ($_ !~ /[^ \t\n\r]/);
    $num_reads++;
    chomp;
    my %pair = ();
    ($pair{id},$pair{flag},$pair{chrom},$pair{left},
     $pair{mapq},$pair{cigar},$pair{mchrom},$pair{mleft},
     $pair{frag_size}) = split (/\t/, $_);
    #my $mright = $pair{mleft} + $OPT{read_length};
    #&print_full($logref, "FULL: ".$_."\n");
    #&print_full($logref, "READ: $read_id ");
    #&print_full($logref, "$mapq ");
    #&print_full($logref, "$mchr ");
    #&print_full($logref, "$chr2 ");
    #&print_full($logref, "$mleft ");
    #&print_full($logref, "$region->{left} ");
    #&print_full($logref, "$mright ");
    #&print_full($logref, "$region->{right}\n");
    my $chr_check = $region->{chrom};
    if (not defined $chr_check   ||
        not defined $pair{chrom} ||
        not defined $pair{mchrom}) {
      &warn("Could not get mate-pair alignment info from:\n\t".$_);
    }
    if ($pair{left} !~ /^\d+$/) {
      &warn("Improperly formatted read left coordinate: ".
            "\"$pair{left}\" from \"$_\".");
    }
    if ($pair{mleft} !~ /^\d+$/) {
      &warn("Improperly formatted mate left coordinate: ".
            "\"$pair{mleft}\" from \"$_\".");
    }
    if ($chr_check eq $pair{chrom}) {
      $chr_check = "=";
    }
    # if the read has already been counted, do not count it again
    if (#(exists $potential_pairs{$pair{id}}) ||
      # if the mapping quality is too low, skip the read
        #($pair{mapq} < $OPT{min_mapq}) ||
      # if the mate's chromosome is not right, skip the read
        ($pair{mchrom} ne $chr_check)# ||
      # if the mate's position is not right, skip the read
        #($pair{mleft} < $region->{left} || $region->{right} < $pair{mleft})
       ) {
      # skip the read
      next;
    }
    # check the read flag to see whether the pair might be supporting
    if (&CheckReadFlag(\%pair, $ctg_ids, $topology)) {
      #&debug("Pair num: $pair{num}");
      #&debug("Potential read: $_");
      %{$potential_pairs{$pair{id}}{$pair{num}}} = %pair;
      #%{$potential_pairs{$pair{id}}{$pair{num}}} = (
      #  mleft     => $pair{mleft},
      #  frag_size => $pair{frag_size},
      #  mapq      => $pair{mapq},
      #  orient    => $pair{orient}
      #);
    }
  }
  close SAMFILE_PRIMARY;

  my %supporting_pairs = ();
  my $num_potential_pairs = keys %potential_pairs;
  if (0 < $num_potential_pairs) {
    if ($topology =~ /gap/) {
      &CheckGapRegion
        ($topology, $region, $ctg_ids, \%potential_pairs, \%supporting_pairs);
    } elsif ($SAM_BOTH == $what_to_sam) {
      # use samtools to get the reads in region B
      my $sam_region = &RegionCoordsStr($region);
      &CheckRegionWithSam
        ($sam_region, $ctg_ids, \%potential_pairs, \%supporting_pairs);
    } else {
      # looks for reads in the potential pairs hash with
      # the correct mleft value
      &CheckRegionWithMLeft
        ($region->{left}, $region->{right},
         \%potential_pairs, \%supporting_pairs);
    }
  }

  # filter the supporting reads
  my %filtered_supporting_pairs = ();
  for my $read_id (keys %supporting_pairs) {
    #if ($supporting_pairs{$read_id}{isintronic}) {
    #  &debug("Supporting pair is intronic");
    #}
    if ($OPT{min_mapq} < $supporting_pairs{$read_id}{mapq}) {
      $filtered_supporting_pairs{$read_id} =
        $supporting_pairs{$read_id}{isintronic};
    }
  }

  my $num_supporting_pairs = keys %supporting_pairs;
  my $num_intronic_pairs =
    grep { $supporting_pairs{$_}{isintronic} } keys %supporting_pairs;
  my $num_filtered_supporting_pairs = keys %filtered_supporting_pairs;
  my $num_filtered_intronic_pairs =
    grep { $filtered_supporting_pairs{$_} } keys %filtered_supporting_pairs;
  &debug_high("Support: ".$num_supporting_pairs.
    " (".$num_intronic_pairs." intronic), ".
    "Filtered: ".$num_filtered_supporting_pairs.
    " (".$num_filtered_intronic_pairs." intronic)");
  return ($num_supporting_pairs, $num_intronic_pairs,
          $num_filtered_supporting_pairs, $num_filtered_intronic_pairs,
          $num_reads);
}

sub CheckGapRegion {
  # get function parameters
  my ($topology, $region, $ctg_ids, $potential_pairs, $supporting_pairs) = @_;
  &debug_high("Checking samtools result for gap-supporting pairs");
  for my $read_id (keys %{$potential_pairs}) {
    for my $pair_num (keys %{$potential_pairs->{$read_id}}) {
      my %pair = %{$potential_pairs->{$read_id}{$pair_num}};
      #&debug("$region->{left}-$region->{right}");
      #&debug("$read_id: $pair{mleft}");
      # REMINDER: get "exonic frag length" somehow
      my ($frag_size, $isintronic) =
        &GetExonicFragmentLength(\%pair, $region->{blocks}, $ctg_ids);
      my $min_isize = (1 - $OPT{frag_diff}) * $OPT{frag_len};
      my $max_isize = (1 + $OPT{frag_diff}) * $OPT{frag_len};
      my $count_read = 0;
      if ("gap-tandem-duplication" eq $topology) {
        if (abs($pair{frag_size}) <= $min_isize) {
          $count_read = 1;
        }
      } elsif ("gap-nontandem-duplication" eq $topology) {
        if (abs($frag_size) <= $min_isize ||
            abs($frag_size) >  $max_isize ||
            "RF" eq $pair{orient}) {
          $count_read = 1;
        }
      } elsif ("gap-tandem-inverted_duplication"    eq $topology ||
               "gap-nontandem-inverted_duplication" eq $topology) {
        if (abs($frag_size) <= $min_isize ||
            "FF" eq $pair{orient} ||
            "RR" eq $pair{orient}) {
          $count_read = 1;
        }
      } elsif ("gap-internal_inversion" eq $topology) {
        # count read-pairs directed "in" to the inverted region
        if ("RR" eq $pair{orient}) {
          if ($region->{right} < ($pair{mleft} + $OPT{read_length})) {
            $count_read = 1;
          }
        }
        if ("FF" eq $pair{orient}) {
          if ($region->{left} > $pair{mleft}) {
            $count_read = 1;
          }
        }
      } else {
        &warn("unrecognized alignment topology: ".$topology);
      }
      # make sure to store the minimum mapq examined for each pair
      if (exists $supporting_pairs->{$read_id} &&
          $pair{mapq} >= $supporting_pairs->{$read_id}{mapq}) {
        $count_read = 0;
      }
      if ($count_read) {
        %{$supporting_pairs->{$read_id}} = (
          isintronic => $isintronic,
          mapq       => $pair{mapq},
        );
      }
    }
  }
}

sub GetExonicFragmentLength {
  # get function parameters
  my ($pair, $blocks_str, $ctg_ids) = @_;
  # assume that the first read in the pair is aligned between blocks
  my $isintronic = 1;
  my ($leftA, $leftB) = ($pair->{left}, $pair->{mleft});
  if ($leftB < $leftA) {
    ($leftA, $leftB) = ($pair->{mleft}, $pair->{left});
  }
  #&debug("Fragment: $leftA-$leftB");
  #&debug("BLOCKS: $blocks_str");
  my @blocks = split(/,/, $blocks_str);
  my $frag_size = 0;
  my $prev_right = $leftA;
  for my $block (@blocks) {
    #&debug("BLOCK: $block");
    my ($block_left, $block_right) = $block =~ /(\d+)-(\d+)/;
    #&debug("Left: $block_left, Right: $block_right");
    my $block_len;
    # skip blocks occurring before the fragment
    if ($block_right < $leftA) {
      next;
    }
    # skip blocks occurring after the fragment
    if ($leftB < $block_left) {
      $block_len = $leftB - $prev_right + 1;
      $frag_size += $block_len;
      #&debug("Adding $block_len, size = $frag_size");
      # the second read in the pair is aligned between blocks
      $isintronic = 1;
      last;
    }
    # make sure to start counting at the left side of the fragment
    if (0 == $frag_size) {
      # check whether the first read in the pair is aligned between blocks
      if ($block_left < $leftA) {
        $isintronic = 0;
      }
      $block_left = $leftA;
    }
    # if the block ends before the fragment,
    # add the whole block length
    if ($block_right < $leftB) {
      $block_len = $block_right - $block_left + 1;
      $frag_size += $block_len;
      #&debug("Adding $block_len, size = $frag_size");
    # if the fragment ends before the block,
    # add only to the end of the fragment
    } else { # $block_right >= $leftB
      $block_len = $leftB - $block_left + 1;
      $frag_size += $block_len;
      #&debug("Adding $block_len, size = $frag_size");
      last;
    }
    $prev_right = $block_right;
  }
  if ($prev_right <= $leftA && 0 == $frag_size) {
    $frag_size = $leftB - $leftA + 1;
  }
  if (1 > $frag_size) {
    &error("Could not get exonic fragment size for region ".$leftA.
      "-".$leftB." with blocks: ".$blocks_str."; contigs: ".$ctg_ids);
  }
  $frag_size += $OPT{read_length};
  return ($frag_size, $isintronic);
}

sub MergeBlocks {
  # get function parameters
  my ($blocks_str1, $blocks_str2) = @_;
  #&debug("MERGING: $blocks_str1; $blocks_str2");
  my @blocks1 = split(/,/, $blocks_str1);
  my @blocks2 = split(/,/, $blocks_str2);
  my $len1 = @blocks1;
  my $len2 = @blocks2;
  #&debug("Num Blocks: $len1; $len2");
  my @merged_blocks = ();
  my ($index1, $index2) = (0, 0);
  while ($index1 < $len1 && $index2 < $len2) {
    my ($left1) = $blocks1[$index1] =~ /(\d+)-\d+/;
    my ($left2) = $blocks2[$index2] =~ /(\d+)-\d+/;
    if ($left1 < $left2) {
      push @merged_blocks, $blocks1[$index1];
      $index1++;
    } else { # $left2 <= $left1
      push @merged_blocks, $blocks2[$index2];
      $index2++;
    }
  }
  while ($index1 < $len1) {
    push @merged_blocks, $blocks1[$index1];
    $index1++;
  }
  while ($index2 < $len2) {
    push @merged_blocks, $blocks2[$index2];
    $index2++;
  }
  my $merged_blocks_str = join(",", @merged_blocks);
  #&debug("Merged blocks: $merged_blocks_str");
  return $merged_blocks_str;
}

sub CheckRegionWithSam {
  # get function parameters
  my ($sam_region, $ctg_ids, $potential_pairs, $supporting_pairs) = @_;
  &debug_high("Using samtools to check region B");
  &RunSamTools($sam_region, $ctg_ids);
  # look for reads in region B with mates in the potential pairs hash
  open(SAMFILE_SECONDARY, $CFG{read_support}{p2g}{sam_out});
  while (<SAMFILE_SECONDARY>) {
    my ($read_id, $flag, $chrom, $left, $mapq) = split;
    my $mpair_num;
    # if the read is the first read in the pair
    if ($flag & 0x0040) {
      # its mate is the second
      $mpair_num = "second";
    # if the read is the second read in the pair
    } elsif ($flag & 0x0080) {
      # its mate is the first
      $mpair_num = "first";
    # otherwise
    } else {
      &error("read was neither the first nor second pair: ".$read_id);
    }
    if (exists $potential_pairs->{$read_id}{$mpair_num}) {
      my %mate = %{$potential_pairs->{$read_id}{$mpair_num}};
      #&debug("$read_id: $mate{mleft}");
      $supporting_pairs->{$read_id}{isintronic} = 0;
      # make sure to store the minimum mapq examined for each pair
      if ($mapq < $mate{mapq}) {
        $supporting_pairs->{$read_id}{mapq} = $mapq;
      } else {
        $supporting_pairs->{$read_id}{mapq} = $mate{mapq};
      }
    }
  }
  close SAMFILE_SECONDARY;
}

sub CheckRegionWithMLeft {
  # get function parameters
  my ($region_left, $region_right, $potential_pairs, $supporting_pairs) = @_;
  &debug_high("Checking samtools result for supporting pairs");
  for my $read_id (keys %{$potential_pairs}) {
    for my $pair_num (keys %{$potential_pairs->{$read_id}}) {
      my %pair = %{$potential_pairs->{$read_id}{$pair_num}};
      #&debug("$region_left-$region_right");
      #&debug("$read_id: $pair{mleft}");
      if ($region_left <= $pair{mleft} && $pair{mleft} <= $region_right) {
        # make sure to store the minimum mapq examined for each pair
        if (not exists $supporting_pairs->{$read_id}) {
          %{$supporting_pairs->{$read_id}} = (
            isintronic => 0,
            mapq       => $pair{mapq},
          );
        } elsif ($pair{mapq} < $supporting_pairs->{$read_id}{mapq}) {
          $supporting_pairs->{$read_id}{mapq} = $pair{mapq};
        }
      }
    }
  }
}

sub CheckReadFlag {
  # get function parameters
  my ($pair, $ctg_ids, $topology) = @_;
  # check that the paired-read flag is properly set
  if ($pair->{flag} & 0x0001) {
    # if the read or its mate is unmapped
    if ($pair->{flag} & 0x0004 || $pair->{flag} & 0x0008) {
      # skip the read
      #&print_full($logref, "Read or pair is unmapped: ".$read_id."\n");
      return 0;
    }
    # strand 0 is the forward strand, strand 1 is reverse
    my ($strand, $mstrand) = (0,0);
    if (($pair->{flag} & 0x0010) > 0) {
      $strand = 1;
    }
    if (($pair->{flag} & 0x0020) > 0) {
      $mstrand = 1;
    }
    #&print_full($logref,
    #   "Strand: $strand, MStrand: ".$mstrand."\n") if ($strand > 0);
    my ($read1, $read2);
    # if the read is the first read in the pair
    if ($pair->{flag} & 0x0040) {
      $pair->{num} = "first";
      $read1 = "$strand $pair->{left}";
      $read2 = "$mstrand $pair->{mleft}";
    # if the read is the second read in the pair
    } elsif ($pair->{flag} & 0x0080) {
      $pair->{num} = "second";
      $read1 = "$mstrand $pair->{mleft}";
      $read2 = "$strand $pair->{left}";
    # otherwise
    } else {
      &warn("read was neither the first nor second pair: $pair->{id}");
    }
    my $orientation;
    if (&CheckStrandAndOrientation
        ($read1, $read2, $ctg_ids, $topology, \$orientation)) {
      #&print_full($logref, "Good Pair: $pair->{id}\n");
      $pair->{orient} = $orientation;
      return 1;
    }
  } else {
    &warn("read was not paired in sequencing: $pair->{id}");
  }
  return 0;
}

sub RunSamTools {
  # get function parameters
  my ($region, $ctg_ids) = @_;
  my $sam_cmd = join(" ",
    $CFG{read_support}{p2g}{sam_cmd}." \"$region\"",
    "1> ".$CFG{read_support}{p2g}{sam_out},
    "2> ".$CFG{read_support}{p2g}{sam_out}.".err");
  &debug_high($sam_cmd);
  my $result = system($sam_cmd);
  my $exit_status = $?;
  my $err_num = $!; #$ERRNO;
  my $err_hash = $!{$err_num}; #$ERRNO{$err_num};
  my $compile_status = $@;
  my $extended_err = $^E;
  my $num_retries = 0;
  if (0 != $exit_status) {
    &status("EXIT: $exit_status ERRNO: $err_num ERRNO_HASH: $err_hash ".
            "COMPILE: $compile_status EXTENDED: $extended_err");
  }
  # if the sam_cmd cound not be started, try a again a few times
  while (-1 == $result && $num_retries < $OPT{max_sam_retries}) {
    $num_retries++;
    $warning = "Could not start samtools, will try again.\n".
      "Exit Status: ".$exit_status."\n".
      "Retry Attempt: ".$num_retries."\n";
    &warn($warning);
    # wait for ten minutes
    sleep(600);
    # try again
    $result = system($sam_cmd);
    $exit_status = $?;
    $err_num = $!; #$ERRNO;
    $err_hash = $!{$err_num}; #$ERRNO{$err_num};
    $compile_status = $@;
    $extended_err = $^E;
    if (0 != $exit_status) {
      &status("EXIT: $exit_status ERRNO: $err_num ERRNO_HASH: $err_hash ".
              "COMPILE: $compile_status EXTENDED: $extended_err");
    }
  }
  if (0 != $result) {
    my $shifted_result = $result >> 8;
    $warning = "Error while running $sam_cmd: $result = $shifted_result\n";
    $warning .= "Contigs: $ctg_ids";
    &warn($warning);
    system("cat $CFG{read_support}{p2g}{sam_out}.err");
    return 0;
  }
  #my $errors = `wc -l $CFG{pair2genome}{sam_out}.err`;
  my $errors = `cat $CFG{read_support}{p2g}{sam_out}.err`;
  #my ($num_errors) = $errors =~ /(\d+\) /;
  #if (0 < $num_errors) {}
  if ($errors) {
    &error("running sam command: ".$sam_cmd."\n".$errors);
  }
}

sub GetRegionToSearch {
  # get function parameters
  my ($full_coords, $topology, $region_len) = @_;
  #&debug("FC:  $full_coords");
  #my ($left,$right) = $full_coords =~ /chr[^:]+:(\d+)-(\d+)/;
  #&print_full($logref, "Left: $left, Right: ".$right."\n");
  my ($region_left,$region_right);
  if ($topology =~ /gap/) {
    # get the gapped alignment event region
    # just use the whole region
    ($region_left, $region_right) =
      ($full_coords->{left}, $full_coords->{right});
    #&warn("using whole region for gapped alignment.");
  } else {
    # get the split alignment event region
    if ((abs($full_coords->{left}-$full_coords->{right}) + 1) < $region_len) {
      # extend the full coords to the correct length
      #&print_full($logref, "Extending full coordinates\n");
      if ("left" eq $full_coords->{junction}) {
        $region_left  = $full_coords->{left};
        $region_right = $full_coords->{left} + $region_len;
      } else {
        $region_left  = $full_coords->{right} - $region_len;
        $region_right = $full_coords->{right};
      }
    } else {
      # walk along the blocks to get the region
      if ("left" eq $full_coords->{junction}) {
        #&print_full($logref, "Walking forwards from left\n");
        $region_left = $full_coords->{left};
        $region_right =
          &WalkForwards
            ($region_len, $full_coords->{left}, $full_coords->{blocks});
      } else {
        #&print_full($logref, "Walking backwards from right\n");
        $region_left =
          &WalkBackwards
            ($region_len, $full_coords->{right}, $full_coords->{blocks});
        $region_right = $full_coords->{right};
      }
    }
  }
  #&debug("RegLeft: $region_left, RegRight: ".$region_right."\n");
  return ($region_left, $region_right);
}

sub CheckStrandAndOrientation {
  # get function parameters
  my ($read1, $read2, $ctg_ids, $topology, $orientation) = @_;
  my ($strand1, $left1) = split(/ /, $read1);
  my ($strand2, $left2) = split(/ /, $read2);
  if ($left1 !~ /^\d+$/) {
    &warn("Improperly formatted read left coordinate: ".
          "\"$left1\" from \"$read1\".");
  }
  if ($left2 !~ /^\d+$/) {
    &warn("Improperly formatted read left coordinate: ".
          "\"$left2\" from \"$read2\".");
  }
  if ("interchr" eq $topology) {
    $$orientation = "NA";
    # count every pair
    return 1;
  }
  # determine the read-pair orientation
  if ($strand1 == $strand2) {
    # if the first read is on the forward strand
    if (0 == $strand1) {
      $$orientation = "FF";
    # if the first read is on the reverse strand
    } else { # 1 == $strand1
      $$orientation = "RR";
    }
  } else { # $strand1 != $strand2
    # assume the reads are oriented inwards
    $$orientation = "FR";
    # if the first read is on the forward strand
    if (0 == $strand1) {
      # if the second read is upstream
      if ($left1 - $left2 + $OPT{read_length} > 0) {
        # the reads are oriented outwards
        $$orientation = "RF";
      }
    # if the first read is on the reverse strand
    } else { # 1 == $strand1
      # if the second read is downstream
      if ($left2 - $left1 + $OPT{read_length} > 0) {
        # the reads are oriented outwards
        $$orientation = "RF";
      }
    }
  }
  if ("local-inversion"     eq $topology ||
      "intrachr-opp-strand" eq $topology) {
    # count read pairs with the same strand
    if ($strand1 == $strand2) {
      return 1;
    }
  } elsif ("junction-duplication"  eq $topology ||
           "intrachr-non-colinear" eq $topology ||
           "end-duplication"       eq $topology) {
    # count read pairs with outwards (RF) orientation
    if ("RF" eq $$orientation) {
      return 1;
    }
  } elsif ("read-through"              eq $topology ||
           "intrachr-same-strand"      eq $topology ||
           "gap-nontandem-duplication" eq $topology) {
    # count read pairs with opposite strands
    if ($strand1 != $strand2) {
      return 1;
    }
  } elsif ("gap-tandem-duplication" eq $topology) {
    # count read pairs with inwards (FR) orientation
    if ("FR" eq $$orientation) {
      return 1;
    }
  } elsif ("gap-tandem-inverted_duplication"    eq $topology ||
           "gap-nontandem-inverted_duplication" eq $topology) {
    # count pairs that do not have outwards (RF) orientation
    if ("RF" ne $$orientation) {
      return 1;
    }
  } elsif ("gap-internal_inversion" eq $topology) {
    # count left- or right-spooned reads
    if ("RR" eq $$orientation ||
        "FF" eq $$orientation) {
      return 1;
    }
  } elsif ("gap-genic_rearrangement" eq $topology) {
    # REMINDER: how to deal with gap genic-rearrangments?
    return 1; # JUST FOR NOW
  } else {
    &warn("in CheckStrandAndOrientation: contigs $ctg_ids ".
          "have unrecognized topology: ".$topology);
  }
  return 0;
}

sub WalkBackwards {
  # get function parameters
  my ($walk_len, $walk_start, $blocks_str) = @_;
  my $walk_end = -1;
  my @blocks = reverse split(/,/, $blocks_str);
  my $num_blocks = @blocks;
  my ($block_left, $block_right);
  if (1 == $num_blocks) {
    ($block_left, $block_right) = $blocks[0] =~ /(\d+)-(\d+)/;
    $walk_end = $block_right - $walk_len;
  } else {
    for my $block (@blocks) {
      ($block_left, $block_right) = $block =~ /(\d+)-(\d+)/;
      my $block_len = abs($block_left-$block_right) + 1;
      if ($block_len < $walk_len) {
        $walk_len -= $block_len;
      } else {
        $walk_end = $block_right - $walk_len;
        $walk_len = 0;
        last;
      }
    }
    if (0 < $walk_len) {
      ($walk_end) = $blocks[-1] =~ /(\d+)-\d+/;
    }
  }
  if (0 > $walk_end) {
    &error("Could not get walk end with blocks: ".$blocks_str);
  }
  return $walk_end;
}

sub WalkForwards {
  # get function parameters
  my ($walk_len, $walk_start, $blocks_str) = @_;
  my $walk_end = -1;
  my @blocks = split(/,/, $blocks_str);
  my $num_blocks = @blocks;
  my ($block_left, $block_right);
  if (1 == $num_blocks) {
    ($block_left, $block_right) = $blocks[0] =~ /(\d+)-(\d+)/;
    $walk_end = $block_left + $walk_len;
  } else {
    for my $block (@blocks) {
      ($block_left, $block_right) = $block =~ /(\d+)-(\d+)/;
      my $block_len = abs($block_left-$block_right) + 1;
      if ($block_len < $walk_len) {
        $walk_len -= $block_len;
      } else {
        $walk_end = $block_left + $walk_len;
        $walk_len = 0;
        last;
      }
    }
    if (0 < $walk_len) {
      ($walk_end) = $blocks[-1] =~ /\d+-(\d+)/;
    }
  }
  if (0 > $walk_end) {
    &error("Could not get walk end with blocks: ".$blocks_str);
  }
  return $walk_end;
}

sub GetGapEventCoords {
  # get function parameters
  my ($gap_coords, $contig_len) = @_;
  &debug_high("Getting gap event coords from: ".$gap_coords);
  my ($gap_left, $gap_right) = $gap_coords =~ /(\d+)-(\d+)/;
  my $gap_size = $gap_right - $gap_left + 1;
  # if the gap is small enough that a read can completely span it
  if ($gap_size + 2 * $OPT{r2c_min_overlap} < $OPT{read_length}) {
    &debug_high("using whole gap");
    my $event_left  = $gap_left  - $OPT{r2c_min_overlap};
    # ensure that the region does not extend past the start of the contig
    if (2 > $event_left) {
      $event_left = $gap_right - $OPT{r2c_min_overlap};
    }
    # really ensure that the region does not start too soon
    if (2 > $event_left) {
      $event_left = 2;
    }
    my $event_right = $gap_right + $OPT{r2c_min_overlap};
    # ensure that the region does not extend past the end of the contig
    if ($contig_len <= $event_right) {
      $event_right = $gap_left + $OPT{r2c_min_overlap};
    }
    # really ensure that the region does not end too late
    if ($contig_len <= $event_right) {
      $event_right = $contig_len - 1;
    }
    return join("-", $event_left, $event_right);
  # if the gap is too big for a read to completely span
  } else {
    &debug_high("using gap edges");
    my $event_leftA  = $gap_left  - $OPT{r2c_min_overlap};
    my $event_rightA = $gap_left  + $OPT{r2c_min_overlap};
    my $event_leftB  = $gap_right - $OPT{r2c_min_overlap};
    my $event_rightB = $gap_right + $OPT{r2c_min_overlap};
    my $event_coordsA = join("-", $event_leftA, $event_rightA);
    my $event_coordsB = join("-", $event_leftB, $event_rightB);
    # just in case a gap spans most of the contig (though this should
    # not happen)
    if (2 > $event_leftA && $contig_len <= $event_rightB) {
      return "N/A";
    # if the gap occurs at the very start of the contig
    } elsif (2 > $event_leftA) {
      # only use the right side of the gap
      return $event_coordsB;
    # if the gap occurs at the very end of the contig
    } elsif ($contig_len <= $event_rightB) {
      # only use the left side of the gap
      return $event_coordsA;
    } else {
      return join(";", $event_coordsA, $event_coordsB);
    }
  }
}

#Determine the estimated breakpoint coord(s) depending on the alignments
sub GetSplitEventCoords {
  # get function parameters
  my ($coords1, $coords2, $ctg_len) = @_;
  #my ($coords1, $coords2, $topology) = @_;
  my ($left1, $right1) = $coords1 =~ /(\d+)-(\d+)/;
  my ($left2, $right2) = $coords2 =~ /(\d+)-(\d+)/;
  my ($leftA, $rightA, $leftB, $rightB);
  # First we sort the range (do this because events found by the gapfilter
  # will not necessarily have these co-ordinates sorted already
  if ($left1 < $left2) {
    $leftA  = $left1;
    $rightA = $right1;
    $leftB  = $left2;
    $rightB = $right2;
  } else {
    $leftA  = $left2;
    $rightA = $right2;
    $leftB  = $left1;
    $rightB = $right1;
  }
  # determine whether there is an overlap or a space
  my ($event_left, $event_right);
  if ($rightA <= $leftB) {
    #here it's blunt or there's a space
    $event_left  = $rightA - $OPT{r2c_min_overlap};
    $event_right = $leftB  + $OPT{r2c_min_overlap};
  } else {
    #here they overlap
    $event_left  = $leftB  - $OPT{r2c_min_overlap};
    $event_right = $rightA + $OPT{r2c_min_overlap};
  }
  # ensure that the event does not extend past either end of the contig
  if (2 > $event_left) {
    $event_left = 2;
  }
  if ($ctg_len <= $event_right) {
    $event_right = $ctg_len - 1;
  }
  return join("-", $event_left, $event_right);
}

#Sorts the range for Set::IntSpan input
sub Sort_Range {
  # get function parameters
  my ($range) = @_;
  my ($coord1,$coord2) = $range =~ /(\d+)-(\d+)/;
  if ($coord1<$coord2) {
    return ($coord1,$coord2);
  } else {
    return ($coord2,$coord1);
  }
}

#Get the contigs to create a separate fasta file
sub Fasta {
  # get function parameters
  my ($fasta_data) = @_;
  &status("Generating contigs fasta file...");
  #Generate the fasta file
  my $barnacle_fasta_out = &GetOutputPath($OPT{outdir}, "contigs");
  open(FASTAOUT,">$barnacle_fasta_out") ||
    &error("Can't open contigs output file: ".$barnacle_fasta_out." ".$!);
  for my $contig_type (keys %{$fasta_data}) {
    #first create the fasta
    if (not open(FASTAIN,$fasta_data->{$contig_type}{seq_file})) {
      &warn("Can not open contig sequence input file: ".$!."\n  ".
            $fasta_data->{$contig_type}{seq_file});
      next;
    }
    my $ctg_id = '';
    my $cov = '';
    while (<FASTAIN>) {
      chomp;
      if ($_ =~ /^>/) {
        ($ctg_id,$cov) = $_ =~ /^>(\S+)(.*)$/;
        next;
      }
      if (exists $fasta_data->{$contig_type}{ctgs}{$ctg_id}) {
        print FASTAOUT ">$ctg_id$cov\n";
        print FASTAOUT "$_\n";
        $ctg_id = '';
        $cov = '';
      }
    }
    #my @lines = <FASTAIN>;
    #for (my $i=0;$i<@lines;$i++) {
    #  my $ctg_id = '';
    #  my $cov = '';
    #  if ($lines[$i] =~ /^>/) {
    #    ($ctg_id,$cov) = $lines[$i] =~ /^>(\S+)(.*)$/;
    #  }
    #  if (exists $fasta_data->{$contig_type}{ctgs}{$ctg_id}) {
    #    print FASTAOUT ">$ctg_id$cov\n";
    #    #print FASTAOUT ">$contig_type:$ctg_id$cov\n";
    #    print FASTAOUT "$lines[$i+1]\n";
    #    $ctg_id = '';
    #  }
    #}
    close FASTAIN;
  }
  close FASTAOUT;
  &status("Finished generating fasta sequences file.\n");
  return undef;
}

sub SetupOutputBaseDirectory {
  # no parameters
  my $out_base_dir = catdir($OPT{input_dir}, "barnacle");
  if (! -d $out_base_dir) {
    &debug_nolog("Creating output base directory: ".$out_base_dir."...");
    make_path($out_base_dir);
  }
  if (! -d $out_base_dir) {
    &error("output base directory does not exist and could not be created: ".
      $out_base_dir);
  }
  return $out_base_dir;
}

sub SetupOutputDirectory {
  # get function parameters
  my ($out_base_dir) = @_;
  my $output_directory;
  #if (defined $OPT{outdir}) {
  #  $output_directory = $OPT{outdir};
  #  chomp $output_directory;
  #  $output_directory =~ s/\/$//; # remove final slash(/) if necessary
  #} else {
  #  $output_directory =
  #    catdir($input_dir, "barnacle", $OPT{ver});
  #    #catdir($input_dir, "barnacle", $output_version);
  #}
  $output_directory = catdir($out_base_dir, $OPT{ver});
  # ensure that the output path is absolute (not relative)
  if (! file_name_is_absolute($output_directory)) {
    $output_directory = abs_path($output_directory);
    print("Absolute output path: \"$output_directory\"\n") if $OPT{debug};
    print("Current directory: ".getcwd()."\n") if $OPT{debug};
  }
  my $results_dir = $output_directory;
  if ($OPT{add_support}) {
    $output_directory = catdir($output_directory, "1_raw_candidates");
  }
  if (! -d $output_directory) {
    &debug_nolog("Creating output directory: $output_directory...");
    make_path($output_directory);
  }
  if (! -d $output_directory) {
    &error("output directory does not exist and could not be created: ".
      $output_directory);
  }
  return ($results_dir, $output_directory);
}

sub CheckDir {
  # get function parameters
  my ($path, $description) = @_;
  unless (-d $path) {
    &error("Cannot find ".$description." directory: ".$path);
  }
}

sub CheckPath {
  # get function parameters
  my ($path, $description) = @_;
  unless (-f $path) {
    &error("Cannot find ".$description.": ".$path);
  }
}

sub GetOutputVersion {
  # get function parameters
  my ($out_base_dir) = @_;
  # read in the version if a version file exists
  &debug_high_nolog("Getting output version...");
  my $version_path = catfile($OPT{barnacle_src_dir}, "version.py");
  &debug_high_nolog("  Version path: ".$version_path);
  if (-f $version_path) {
    &debug_high_nolog("    found version file!");
    open(VERSION, $version_path) ||
      &err("Can't open Barnacle version file ".$version_path." ".$!);
    # get a line from the file
    while (<VERSION>) {
      chomp; # remove any trailing newlines, etc.
      # skip blank lines and comments
      if ("" eq $_ or $_ =~ /^#/) {
        next;
      }
      if ($_ =~ /^VERSION = "/) {
        ($OPT{code_version}) = $_ =~ /^VERSION = "([0-9.]+)"/;
        print "Running Barnacle version ".$OPT{code_version}."\n";
        last;
      }
    }
    close VERSION;
  } else {
    &debug_nolog("Cannot find version file: ".$version_path);
  }
  if (! defined $OPT{code_version}) {
    &debug_nolog("Could not set code version!");
  }

  # get the proper output version
  if (defined $OPT{ver}) {
    &debug_high_nolog("  Output version parameter: ".$OPT{ver});
    if ($OPT{ver} !~ /^ver_/) {
      $OPT{ver} = "ver_".$OPT{ver};
    }
    if (defined $OPT{code_version} &&
      $OPT{ver} !~ /^ver_$OPT{code_version}\.[0-9]+$/) {
      &error("Invalid output version \"".$OPT{ver}."\": should start with \"".
        $OPT{code_version}."\" and end with an output specific number.");
    }
  } elsif (defined $OPT{code_version}) {
    # check output directory for existing versions and add 1
    if (! -d $out_base_dir) {
      &error("Can't find Barnacle output base directory ".$out_base_dir);
    }
    opendir(OUT_BASE_DIR, $out_base_dir) ||
      &error("Can't open Barnacle output base directory ".$out_base_dir." ".$!);
    my @barnacle_out_dirs =
      grep {/$OPT{code_version}\.[0-9]+$/} grep {$_ !~ /^\./} readdir OUT_BASE_DIR;
    closedir(OUT_BASE_DIR);
    my $outver = 0;
    # if no files were found
    if (0 == @barnacle_out_dirs) {
      &debug_nolog("  Found no existing Barnacle output directories");
    } else {
      my @outvers = ();
      for my $out_dir (@barnacle_out_dirs) {
        ($outver) = $out_dir =~ /\.([0-9]+)$/;
        &debug_high_nolog("  Output directory: ".$out_dir);
        unless (defined $outver) {
          &error("cannot extract output version from output directory");
        }
        &debug_high_nolog("   Output version: ".$outver);
        push @outvers, $outver;
      }
      $outver = (sort @outvers)[-1];
      unless (defined $outver) {
        &error("error selecting latest output version");
      }
      my $run_support_path = catfile($out_base_dir,
        "ver_".$OPT{code_version}.".".$outver, "run_support.sh");
      &debug_high_nolog("    Last output version: ".$outver.
        "\n  run_support.sh path: ".$run_support_path);
      if ($OPT{identify_candidates} && -f $run_support_path) {
        &debug_high_nolog("    run_support.sh exists: incrementing version.");
        $outver += 1;
      }
    }
    $OPT{ver} = "ver_".$OPT{code_version}.".".$outver;
    &debug_high_nolog("  No output version parameter. ".
      "Output version set to ".$OPT{ver});
  } else {
    &error("Output version not provided and cannot find code version to use ".
      "to construct output version.");
  }
}

sub ReadConfigFile {
  # no parameters
  my %config_options = (
    '-ConfigFile'      => $OPT{config},
    '-CComments'       => 'no',
    '-InterPolateVars' => 0
  );
  # create a config object from the file for the library-specific options
  my $config_obj = new Config::General(%config_options);
  # get the configuration options from the object
  %CFG = $config_obj->getall();

  # get the general BARNACLE configuration options
  &debug_nolog("Loading Barnacle configuration file...");
  my $barnacle_cfg_path = catfile($OPT{barnacle_src_dir}, "barnacle.cfg");
  open(BARNACLE_CFG, $barnacle_cfg_path) ||
    &error("Can't open BARNACLE configuration file ".
      $barnacle_cfg_path." ".$!);
  # get a line from the file
  my $section = "";
  while (<BARNACLE_CFG>) {
    chomp; # remove any trailing newlines, etc.
    # skip blank lines and comments
    if ("" eq $_ or $_ =~ /^#/) {
      next;
    }
    if ($_ =~ /^\[/) {
      ($section) = $_ =~ /^\[([^\]]+)\]/;
      &debug_nolog("Setting section to: ".$section);
      next;
    }
    # split the line into field name and value
    my ($field, $value) = split /=/;
    &debug_nolog("Setting ".$field." to: ".$value);
    if ("" eq $section) {
      $CFG{$field} = $value;
    } else {
      $CFG{$section}{$field} = $value;
    }
  }
  close BARNACLE_CFG;

  if ($OPT{debug}) {
    print "---CONFIG---\n";
    &PrintHash(\%CFG, "");
    print "------------\n";
  }

  # set up the input directory from the template
  &GetInputDir();

  # get contig-to-genome alignments template
  &SetupC2GAlignmentsTemplate();

  # check that the contig sequences path is valid
  &GetSequencesPath("contigs", "contig");

  # check that the Python path is valid
  &CheckPathFromConfig("commands", "python", "Python command");
  # if using the gap filter
  unless ($OPT{no_gap_check}) {
    # check that the BLAT path is valid
    &CheckPathFromConfig("commands", "blat", "BLAT command");
    # check that the 2bit genome path is valid
    #$CFG{sequences}{genome_2bit} =
    #  &ReplaceVars($CFG{sequences}{genome_2bit});
    &ReplaceCfgVars("sequences", "genome_2bit");
    &CheckPathFromConfig("sequences", "genome_2bit", "2bit genome file");
    # check that the blat options were given
    unless (exists $CFG{gap_realigner}{BLAToptions}) {
      &error("No value given for BLAToptions parameter in config file.");
    }
  }
  # check that the SAMtools path is valid
  &CheckPathFromConfig("commands", "samtools", "SAMtools command");
  # check that the mqsub path is valid
  #   mqsub script only exists on cluster...
  #&CheckPathFromConfig("commands", "mqsub", "mqsub script");

  # if flagging exon boundaries
  if ($OPT{flag_exon_boundary_junctions}){
    # check that the gene annotations file is valid
    &ReplaceCfgVars("annotations", "genes");
    &CheckPathFromConfig("annotations", "genes",
      "gene annotations file");
  }
  # if adding gene annotations
  if ($OPT{add_gene_annotation}) {
    # check that the gene feature coordinates file is valid
    &ReplaceCfgVars("coords", "gene_features");
    &CheckPathFromConfig("coords", "gene_features",
      "gene feature coordinates file");
  }
  # if filtering out repetitive structural RNA
  #unless ($OPT{no_flag_RNA}) {
  #  # check that the repetitive structural RNA coords file was found
  #  &ReplaceCfgVars("coords", "structural_RNA");
  #  &CheckPathFromConfig("coords", "structural_RNA",
  #    "structural RNA coordinates file");
  #}
  # if checking for repeat sequences
  if ($OPT{flag_repeats}) {
    # check that repeat file path(s) is/are given
    if (! exists $CFG{coords}{repeats}) {
      &error("No repeat coordinates file path(s) found ".
        "in config file: ".$OPT{config});
    }
    # ensure that the repeat file path entry is an array
    unless ("ARRAY" eq ref($CFG{coords}{repeats})) {
      my $single_repeat_path = $CFG{coords}{repeats};
      delete $CFG{coords}{repeats};
      @{$CFG{coords}{repeats}} = ($single_repeat_path);
    }
    my @replaced_repeat_paths = ();
    for my $repeat_path (@{$CFG{coords}{repeats}}) {
      # check that the repeat file is valid
      $repeat_path = &ReplaceVars($repeat_path);
      &CheckPath($repeat_path, "repeat coordinates");
      push @replaced_repeat_paths, $repeat_path;
    }
    @{$CFG{coords}{repeats}} = @replaced_repeat_paths;
    if (exists $CFG{repeats}{cols_to_use}) {
      # check the format of the cols_to_use
      if ($CFG{repeats}{cols_to_use} =~ /[^0-9\/]/) {
        &error("Improperly formatted cols_to_use value for repeat ".
          "coordinates file: \"".$CFG{repeats}{cols_to_use}."\"; ".
          "should only contain numbers and slashes.");
      }
    } else {
      # set the default cols_to_use if none was given
      $CFG{repeats}{cols_to_use} = "1";
    }
  }
  # if calculating pair-to-genome support
  unless ($OPT{no_pair_to_genome} || $OPT{identify_candidates}) {
    &StartPairToGenomeSetup();
  }
  # if calculating read-to-contig support
  unless ($OPT{no_read_to_contig} || $OPT{identify_candidates}) {
    &SetupReadToContigSupport();
  }
  # set the default cluster options, if present
  if (exists $CFG{cluster}) {
    &SetDefaultFromConfig("cluster", "hostname");
    &SetDefaultFromConfig("cluster", "queue");
  }
}

sub SetDefaultFromConfig {
  # get function parameters
  my ($section, $key) = @_;
  if (exists $CFG{$section}{$key}) {
    &SetDefaultValue($key, $CFG{$section}{$key});
  }
}

sub CheckPathFromConfig {
  # get function parameters
  my ($section, $key, $description) = @_;
  unless (exists $CFG{$section}{$key}) {
    &error("No ".$description." path found in config file: ".$OPT{config});
  }
  &CheckPath($CFG{$section}{$key}, $description);
}

sub StartPairToGenomeSetup {
  # no parameters
  # get pair2genome alignment info
  #my $pair2genome_align_dir=catdir($lib_dir, "Reads_to_genome");
  # any file that starts with the library id and ends with a .bam extension
  # in $pair2genome_align_dir will be used, may be multiples
  #my $pair2genome_align_file=$OPT{lib}."*.bam";
  #my $pair2genome_path =
  #  catfile($pair2genome_align_dir, $pair2genome_align_file);
  &GetAlignmentPath("p2g", "pair-to-genome");
}

sub FinishPairToGenomeSetup {
  # no parameters
  # put the temporary samtools output files in the out directory
  $CFG{read_support}{p2g}{sam_out} = catfile($OPT{outdir}, "sam_out_tmp");
  # construct the default samtools command
  $CFG{read_support}{p2g}{sam_cmd} =
    $CFG{commands}{samtools}." view ".$CFG{alignments}{p2g};
  # check whether to include "chr" in chromosome IDs
  # and whether mitochondrial DNA should be "M" or "MT"
  &SetupP2GChromosomeInfo();
}

sub SetupP2GChromosomeInfo {
  # no function params
  # get the first line of the pairs to genome alignment file
  my $sam_cmd =
    #"$CFG{pair2genome}{sam_cmd} | head -n1";
    $CFG{read_support}{p2g}{sam_cmd}." -H 2>".
    $CFG{read_support}{p2g}{sam_out}.".err";
  my $sam_result = `$sam_cmd`;
  my $exit_status = $?;
  my $err_num = $!; #$ERRNO;
  my $err_hash = $!{$err_num}; #$ERRNO{$err_num};
  my $compile_status = $@;
  my $extended_err = $^E;
  &debug_nolog("Checking whether chromosome should use \"chr\":\n".
         "  ".$sam_cmd);
       #"  ".$sam_cmd."\n  \"".$sam_result."\"");
  my $chr_line;
  my @sam_lines = split(/\n/, $sam_result);
  &debug_nolog(join("\n", @sam_lines[0..4]));
  for my $sam_line (@sam_lines) {
    if ($sam_line =~ /\@SQ/) {
      $chr_line = $sam_line;
      last;
    }
  }
  if (! defined $chr_line || 0 != $exit_status) {
    &status("ERRNO: $err_num ERRNO_HASH: $err_hash COMPILE: $compile_status ".
            "EXTENDED: $extended_err");
    system("cat ".$CFG{read_support}{p2g}{sam_out}.".err");
    &error("could not get chromosome info from sam header: \"".
            $sam_result."\"\nExit status: ".$exit_status."\nCommand used: ".
            $sam_cmd);
  }
  $CFG{read_support}{p2g}{use_chr} = &ShouldChromUseChr($chr_line, 1);
  $CFG{read_support}{p2g}{use_MT}  = $sam_result =~ /SN:\S*MT/;
  return undef;
}

sub SetupReadToContigSupport {
  # no parameters
  # check that the contig-sorted bam alignment file exists
  #if (exists $CFG{read2contig_bam_align_path_cs}) {
  #  &status("Using read2contig from config");
  #  $CFG{read2contig_bam_align_path_cs} =~ s/LIB/$OPT{lib}/g;
  #} else {
  #  &status("Using default read2contig");
  #  $CFG{read2contig_bam_align_path_cs} = catdir
  #    ($input_dir, "reads_to_contigs", $OPT{lib}."-contigs.bam");
  #}
  #unless (-f $CFG{read2contig_bam_align_path_cs}) {
  #  &error("Cannot find read-to-contig alignment file: ".
  #    $CFG{read2contig_bam_align_path_cs});
  #}
  #$CFG{read2contig_sort_by_read} = 1;
  &GetAlignmentPath("r2c", "read-to-contig");
}

sub GetInputDir {
  # no parameters
  unless (exists $CFG{path_templates}{input_dir}) {
    &error("No input directory path template found in config file: ".
      $OPT{config});
  }
  # av = assembly version, i.e. "${assembly_ver}" in template
  my ($before_av, $after_av) = $CFG{path_templates}{input_dir} =~
    /(.*)\${assembly_ver}(.*)/;
  &debug_high_nolog("Before: ".$before_av."\nAfter: ".$after_av);
  $before_av = &ReplaceVars($before_av);
  # get assembly version
  if ("current" eq $OPT{assembly_ver}) {
    &debug_nolog("Attempting to get assembly version from ".$before_av);
    my $current_dir = $before_av.$OPT{assembly_ver};
    &debug_nolog("  Current path: ".$current_dir);
    if (-f $current_dir) {
      my $current_target = `ls -ld $current_dir`;
      &debug_nolog("  Current target: ".$current_target);
      my $current_ver = "";
      ($current_ver) = $current_target =~
        /([^ ]+-[0-9]+\.[0-9]+(\.[0-9]+)?)\/?$/;
      &debug_nolog("  Current version: ".$current_ver);
      if (($OPT{assembler} && $current_ver =~ /^$OPT{assembler}/) ||
          (! $OPT{assembler} && $current_ver =~ /\w/)) {
        $OPT{assembly_ver} = $current_ver;
      } else {
        &debug_nolog("Could not get assembly version pointed to by \"current\".");
      }
    } else {
      &debug_nolog("  No \"current\" symlink");
      if (-d $before_av) {
        my %assembly_vers;
        my %assembly_seps;
        opendir(BEFORE_DIR, $before_av) ||
          &error("Can't open assemblies directory ".$before_av." ".$!);
        my @assembly_dirs =
          grep {/[0-9.]+$/} grep {$_ !~ /^\./} readdir BEFORE_DIR;
        closedir(BEFORE_DIR);
        my ($assembler, $sep, $av1, $av2, $av3, $av4);
        for my $assembly_dir (@assembly_dirs) {
          ($assembler, $sep, $av1, $av2, $av3, $av4) = $assembly_dir =~
            /^(.*)([-_]v?)([0-9]+)\.([0-9]+)\.([0-9]+)(\.[0-9]+)?$/;
          if (defined $av4) {
            $av4 =~ s/\.//g;
          } else {
            $av4 = -1;
          }
          $assembly_vers{$assembler}{$av1}{$av2}{$av3}{$av4} = 1;
          $assembly_seps{$assembler} = $sep;
        }
        # if no files were found yet
        if (0 == keys %assembly_vers) {
          &error("could not find any assembly directories");
        }
        #&PrintHash(\%assembly_vers, "");
        if (defined $OPT{assembler}) {
          unless (exists $assembly_vers{$OPT{assembler}}) {
            &error("cannot find any assembly directories for ".
              $OPT{assembler});
          }
        } elsif (1 == keys %assembly_vers) {
          ($OPT{assembler}) = keys %assembly_vers;
          &debug_high_nolog("  Setting assembler to ".$OPT{assembler});
        } else {
          &error("  assembly directories from multiple assemblers found in ".
            $before_av.". Please specify assembler to use with ".
            "\"-assembler\" option or specify specific assembly with ".
            "\"-assembly_ver\" option.");
        }
        my @avs = ();
        my %av1_opts = %{$assembly_vers{$OPT{assembler}}};
        push @avs, (sort keys %av1_opts)[-1];
        my %av2_opts = %{$av1_opts{$av1}};
        push @avs, (sort keys %av2_opts)[-1];
        my %av3_opts = %{$av2_opts{$av2}};
        push @avs, (sort keys %av3_opts)[-1];
        my %av4_opts = %{$av3_opts{$av3}};
        $av4 = (sort keys %av4_opts)[-1];
        if (-1 != $av4) {
          push @avs, $av4;
        }
        $OPT{assembly_ver} = $OPT{assembler}.$assembly_seps{$OPT{assembler}}.
          join(".", @avs);
        &debug_nolog("Assembly version set to: ".$OPT{assembly_ver});
      } else {
        &error("Can't find assemblies directory ".$before_av);
      }
    }
  }
  $OPT{input_dir} =
    &ReplaceVars($before_av.$OPT{assembly_ver}.$after_av);
  $OPT{input_dir} = &ReplaceWildcards($OPT{input_dir}, "input directory");
  &debug_nolog("Using input directory: ".$OPT{input_dir});
  &CheckDir($OPT{input_dir}, "input");
}

sub SetupC2GAlignmentsTemplate {
  &GetAlignmentPath("c2g", "contig-to-genome", 0);
  $OPT{c2g_dir} = dirname($CFG{alignments}{c2g});
  $OPT{c2g_template} = basename($CFG{alignments}{c2g});
  &CheckDir($OPT{c2g_dir}, "contig-to-genome alignments");
}

sub GetSequencesPath {
  # get function parameters
  my ($sequence_type, $description) = @_;
  unless (exists $CFG{sequences}{$sequence_type}) {
    &error("No ".$description." sequences path found in config file: ".
      $OPT{config});
  }
  #$CFG{sequences}{$sequence_type} =
  #  &ReplaceVars($CFG{sequences}{$sequence_type});
  &ReplaceCfgVars("sequences", $sequence_type);
  $CFG{sequences}{$sequence_type} = &ReplaceWildcards(
    $CFG{sequences}{$sequence_type}, $description." sequences");
  &debug_nolog("Using ".$description." sequences file: ".
    $CFG{sequences}{$sequence_type}."\n");
  &CheckPath($CFG{sequences}{$sequence_type}, $description." sequences");
}

sub GetAlignmentPath {
  # get function parameters
  my ($support_type, $description, $replace_wild_and_check) = @_;
  $replace_wild_and_check = 1 unless defined $replace_wild_and_check;
  unless (exists $CFG{alignments}{$support_type}) {
    &error("No ".$description." alignments path found in config file: ".
      $OPT{config});
  }
  #$CFG{alignments}{$support_type} =
  #  &ReplaceVars($CFG{alignments}{$support_type});
  &ReplaceCfgVars("alignments", $support_type);
  if ($replace_wild_and_check) {
    $CFG{alignments}{$support_type} = &ReplaceWildcards(
      $CFG{alignments}{$support_type}, $description." alignment");
  }
  &debug_nolog("Using ".$description." alignment file: ".
    $CFG{alignments}{$support_type}."\n");
  if ($replace_wild_and_check) {
    &CheckPath($CFG{alignments}{$support_type}, $description." alignment file");
  }
}

sub ReplaceVars {
  # get function parameters
  my ($in_path) = @_;
  &debug_nolog("Replacing variables in \"".$in_path."\"");
  my $out_path = $in_path;
  while ($out_path =~ /\${[^}]+}/) {
    my ($var_name) = $out_path =~ /\${([^}]+)}/;
    unless (exists $OPT{$var_name}) {
      &error("unrecognized variable in configuration path: ".$var_name."\n".
        $in_path);
    }
    $out_path =~ s/\${$var_name}/$OPT{$var_name}/g
  }
  &debug_nolog("After replacement: \"".$out_path."\"");
  return $out_path;
}

sub ReplaceCfgVars {
  # get function parameters
  my ($section, $key) = @_;
  $CFG{$section}{$key} = &ReplaceVars($CFG{$section}{$key});
}

sub ReplaceWildcards {
  # get function parameters
  my ($in_path, $description) = @_;
  &debug_nolog("Replacing wildcards in \"".$in_path."\"");
  my $out_path = $in_path;
  if ($in_path =~ /[\*\?]/) {
    &debug_nolog("Looking for ".$description." file: ".
      $in_path."\n");
    # "-t" sorts by modification time, newest first
    my @pair2genome_files =
      `ls -t $in_path 2>/dev/null`;
    my $num_pair2genome_files = @pair2genome_files;
    if (0 == $num_pair2genome_files) {
      &error("Could not find ".$description." file: ".$in_path);
    } elsif (1 < $num_pair2genome_files) {
      &warn("found multiple ".$description."files: ".$in_path);
    }
    #for my $p2g_path (@pair2genome_files) {
    #  &debug_nolog("Found p2g path: \"".$p2g_path."\"");
    #}
    #use the first (most recent) one only?
    $out_path = $pair2genome_files[0];
    chomp($out_path);
  }
  &debug_nolog("After replacement: \"".$out_path."\"");
  return $out_path;
}

sub PrintHash {
  # get function parameters
  my ($hash, $indent) = @_;
  for my $key (sort keys %{$hash}) {
    if ("HASH" eq ref($hash->{$key})) {
      print "$indent$key:\n";
      &PrintHash(\%{$hash->{$key}}, $indent."  ");
      #for my $key2 (sort keys %{$CFG{$key1}}) {
      #  if (ref $CFG{$key1}{$key2} eq 'HASH') {
      #    print "  $key2:\n";
      #    for my $key3 (sort keys %{$CFG{$key1}{$key2}}) {
      #      print "    $key3 => $CFG{$key1}{$key2}{$key3}\n";
      #    }
      #  } else {
      #    print "  $key2 => $CFG{$key1}{$key2}\n";
      #  }
      #}
    } elsif ("ARRAY" eq ref($hash->{$key})) {
      print $indent.$key." =>\n".$indent."  ".
        join(",\n".$indent."  ", @{$hash->{$key}})."\n";
    } else {
      print "$indent$key => $hash->{$key}\n";
    }
  }
}

sub PrintGroupCtgIDs {
  # get function parameters
  my ($groups) = @_;
  #for my $group_id (sort keys %{$groups}) {}
  for my $group (@{$groups}) {
    &status("Group $group->{id}");
    for my $ctg_id (sort keys %{$group->{members}}) {
      &status("  ".$ctg_id." ".
        join(" ", sort keys %{$group->{members}{$ctg_id}{info}})."."
      );
    }
  }
}

#Print to stdout and a filehandle at the same time
sub print_full {
  # get function parameters
  my ($handle,$string,$require_log) = @_;
  $require_log = 1 unless defined $require_log;
  if (defined $handle) {
    print $handle $string;
  } elsif ($require_log) {
    &error("File not opened, cannot print: ".$string);
  }
  print $string;
}

sub status {
  # get function parameters
  my ($msg, $require_log) = @_;
  &print_full($logref, $msg."\n", $require_log);
}

sub warn {
  # get function parameters
  my ($warning) = @_;
  &status("WARNING: ".$warning);
}

sub error {
  # get function parameters
  my ($err_msg) = @_;
  &status("ERROR: ".$err_msg, 0);
  die;
}

sub debug {
  # get function parameters
  my ($msg) = @_;
  &status($msg) if ($OPT{debug} || $OPT{debug_high});
}

sub debug_nolog {
  # get function parameters
  my ($msg) = @_;
  if ($OPT{debug} || $OPT{debug_high}) {
    print $msg."\n";
  }
}

sub debug_high {
  # get function parameters
  my ($msg) = @_;
  &status($msg) if $OPT{debug_high};
}

sub debug_high_nolog {
  # get function parameters
  my ($msg) = @_;
  if ($OPT{debug_high}) {
    print $msg."\n";
  }
}

sub TimeSpent {
  # get function parameters
  my ($process, $start, $total) = @_;
  my $end  = new Benchmark;
  my $time = timediff($end, $start);
  my $time_str = timestr($time);
  my $head = "Time";
  if ($total) { $head = "Total time"; }
  &status("$head spent $process: ".$time_str."\n");
}

sub GetOutputPath {
  # get function parameters
  my ($dir, $ext) = @_;
  return catfile($dir, $OPT{lib}.".barnacle.".$ext);
}

sub ShouldChromUseChr {
  # get function parameters
  my ($line, $col) = @_;
  chomp ($line);
  &debug("LINE: ".$line);
  # check whether the chromosome uses "chr" in its id
  my @fields;
  if ($line =~ /\t/) {
    @fields = split(/\t/, $line);
  } else {
    @fields = split(/ /, $line);
  }
  my $chrom = $fields[$col];
  &debug("CHROM:".$chrom);
  if ($chrom =~ /chr/) {
    return 1;
  }
  return 0;
}

sub ShouldChromUseChrFile {
  # get function parameters
  my ($file, $col) = @_;
  my $line = `head -n1 $file`;
  return &ShouldChromUseChr($line, $col);
}

sub RegionCoordsStr {
  # get function parameters
  my ($region) = @_;
  my $region_coords_str =
    "$region->{chrom}:$region->{left}-$region->{right}";
  return $region_coords_str;
}

sub Desc {
  # get function parameters
  my ($contig_type) = @_;
  if ($contig_type =~ /contigs/) {
    return $contig_type;
  }
  return $contig_type." contigs";
}

sub OpenLogFile {
  # no parameters
  my $log_file_path = &GetOutputPath($OPT{outdir}, "log");
  if ($OPT{add_support}) {
    $log_file_path .= ".continue";
  }
  print "Log file: ".$log_file_path."\n\n" if ($OPT{debug});
  open(my $logref,">",$log_file_path) ||
    &error("Can't open log file ".$log_file_path." ".$!);
  # turn off buffering on log file (make log file "hot")
  my $ofh = select $logref;
  $| = 1;
  select $ofh;
  return ($log_file_path, $logref);
}

sub CheckOptionsForConflicts {
  # no parameters
  my @conflicts = ();
  #&CheckConflict("only_blat_aligns", "only_exon_aligns", \@conflicts);
  #&CheckConflict("only_blat_aligns", "blat_and_exon_aligns", \@conflicts);
  #&CheckConflict("only_blat_aligns", "require_both_aligners", \@conflicts);
  #&CheckConflict("only_exon_aligns", "blat_and_exon_aligns", \@conflicts);
  #&CheckConflict("only_exon_aligns", "require_both_aligners", \@conflicts);
  &CheckConflict("use_smart_chooser", "use_quick_chooser", \@conflicts);
  &CheckConflict("mito_prefilter_on", "mito_prefilter_off", \@conflicts);
  &CheckConflict("flag_exon_boundary_junctions", "no_flag_exon_boundary_junctions", \@conflicts);
  &CheckConflict("add_gene_annotation", "no_gene_annotation", \@conflicts);
  &CheckConflict("prefer_exons", "no_exon_preference", \@conflicts);
  &CheckConflict("prefer_exons", "no_gene_annotation", \@conflicts);
  &CheckConflict("pair_to_genome", "no_pair_to_genome", \@conflicts);
  &CheckConflict("cluster_pair2gen", "no_pair_to_genome", \@conflicts);
  &CheckRequired("cluster_pair2gen", "cluster", \@conflicts);
  &CheckConflict("read_to_contig", "no_read_to_contig", \@conflicts);
  &CheckConflict("r2c_short_ctgs", "r2c_no_short_ctgs", \@conflicts);
  &CheckConflict("flag_RNA", "no_flag_RNA", \@conflicts);
  &CheckConflict("flag_repeats", "no_flag_repeats", \@conflicts);
  &CheckConflict("no_flag_repeats", "repeat_search_size", \@conflicts);
  &CheckConflict("breakpoint_genes", "no_breakpoint_genes", \@conflicts);
  &CheckConflict("check_split_and_gap", "no_split_check", \@conflicts);
  &CheckConflict("check_split_and_gap", "no_gap_check", \@conflicts);
  &CheckConflict("no_split_check", "no_gap_check", \@conflicts);
  &CheckConflict("no_gap_check", "include_gap_events", \@conflicts);
  if (@conflicts) {
    print "Conflicting options used:\n";
    for my $conflict (@conflicts) {
      print "  ".$conflict."\n";
    }
    die;
  }
  return undef;
}

sub CheckConflict {
  # get function parameters
  my ($opt1, $opt2, $conflicts) = @_;
  if ($OPT{$opt1} && $OPT{$opt2}) {
    push @{$conflicts}, "-".$opt1." and -".$opt2;
  }
}

sub CheckRequired {
  # get function parameters
  my ($opt1, $opt2, $conflicts) = @_;
  if ($OPT{$opt1} && ! $OPT{$opt2}) {
    push @{$conflicts}, "-".$opt1." without -".$opt2;
  }
}

sub SetDefaultValue {
  # get function parameters
  my ($key, $value) = @_;
  if (not defined $OPT{$key}) {
    $OPT{$key} = $value;
  }
}

sub SetDefaultFlag {
  # get function parameters
  my ($default_key, $other_key) = @_;
  $OPT{$default_key} = 1;
  if ($OPT{$other_key}) {
    $OPT{$default_key} = 0;
  }
}

sub SetOptionDefaults {
  # no parameters
  &SetDefaultValue("contig_types", "contigs");
  #&SetDefaultValue("only_blat_aligns", 1);
  &SetDefaultValue("require_both_aligners", 0);
  # alignment processing parameters
  &SetDefaultValue("cid_memory", "8G");
  &SetDefaultFlag("use_smart_chooser", "use_quick_chooser");
  &SetDefaultValue("max_match_fraction", 0.999);
  &SetDefaultValue("min_merge_overlap", 0.80);
  &SetDefaultValue("min_ctg_represented", 0.85);
  &SetDefaultValue("num_aligns", 500);
  &SetDefaultValue("min_identity", 40.0);
  &SetDefaultFlag("mito_prefilter_on", "mito_prefilter_off");
  &SetDefaultFlag("flag_exon_boundary_junctions",
    "no_flag_exon_boundary_junctions");
  &SetDefaultValue("exon_bound_buffer", 4);
  # check whether to annotate events with the genes involved
  &SetDefaultValue("add_gene_annotation", 1);
  &SetDefaultFlag("no_exon_preference", "prefer_exons");
  # Filter parameter: max number of events a contig can be found in to be
  # considered good evidence.
  #   adds fail reason TOO_MANY_GROUPS
  &SetDefaultValue("max_num_groups", 0);
  &SetDefaultFlag("no_pair_to_genome", "pair_to_genome");
  if ($OPT{cluster}) {
    &SetDefaultFlag("cluster_pair2gen", "no_pair_to_genome");
  } else {
    &SetDefaultValue("cluster_pair2gen", 0);
  }
  &SetDefaultValue("pair2gen_split", 5000);
  &SetDefaultValue("p2g_memory", "5G");
  &SetDefaultValue("frag_len", 200);
  &SetDefaultValue("frag_diff", 0.10);
  &SetDefaultValue("min_mapq", 10);
  &SetDefaultFlag("read_to_contig", "no_read_to_contig");
  &SetDefaultValue("r2c_memory", "5G");
  &SetDefaultValue("r2c_min_overlap", 5);
  &SetDefaultValue("r2c_ctgs_per_job", 1500);
  &SetDefaultFlag("r2c_no_short_ctgs", "r2c_short_ctgs");
  &SetDefaultFlag("flag_RNA", "no_flag_RNA");
  &SetDefaultFlag("flag_repeats", "no_flag_repeats");
  &SetDefaultValue("repeat_search_size", 100);
  &SetDefaultFlag("breakpoint_genes", "no_breakpoint_genes");
  &SetDefaultValue("read_length", 75);
  &SetDefaultFlag("include_gap_events", "no_gap_check");
  # gap filter parameters
  &SetDefaultValue("min_gap_size", 4);
  &SetDefaultValue("min_gap_identity", 0.95);
  &SetDefaultValue("min_gap_fraction", 0.3);
  &SetDefaultValue("gap_max_len", 50000);
  # use the default max number of samtools retries, if none is given
  &SetDefaultValue("max_sam_retries", 5);
  #&SetDefaultValue("hostname", "node*");
  #&SetDefaultValue("queue", "all.q");
  &SetDefaultValue("assembly_ver", "current");
  if ($OPT{debug_high}) {
    $OPT{debug} = 1;
  }
  return undef;
}

sub JobLine {
  # get function parameters
  my ($command) = @_;
  return join(";",
    "export BARNACLE_PATH=".$OPT{barnacle_src_dir},
    "export PATH=\$BARNACLE_PATH:\$PATH",
    "export PYTHONPATH=.:\$BARNACLE_PATH:\$PYTHONPATH",
    #"export PERL5LIB=".`echo -n \$PERL5LIB`,
    #"export PERL5LIB=".catdir($OPT{barnacle_src_dir},"perl_libs"),
    "time ".$command."\n",
  );
}

sub CoordInRange {
  # get function parameters
  my ($coord, $left, $right) = @_;
  if ($coord < $left || $right < $coord) {
    return 0;
  }
  return 1;
}

sub RangesOverlap {
  # get function parameters
  my ($left1, $right1, $left2, $right2) = @_;
  if ($right2 < $left1 || $right1 < $left2) {
    return 0;
  }
  return 1;
}

sub MergeRanges {
  # get function parameters
  my ($left1, $right1, $left2, $right2) = @_;
  my ($merged_left, $merged_right) = ($left1, $right1);
  if ($left2 < $left1) {
    $merged_left = $left2;
  }
  if ($right2 > $right1) {
    $merged_right = $right2;
  }
  return ($merged_left, $merged_right);
}
