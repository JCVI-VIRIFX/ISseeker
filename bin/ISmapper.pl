#!/usr/local/bin/perl

#############################################################################
#
#    ISmapper    - Finds portions of contigs flanking IS sequences
#                  and blasts them against annotated references 
#                  to infer IS insertion points
#
#    Written by Brian Bishop
# 
#
#    Copyright (C) 2015  JCVI
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    See example in doc directory
#
#
#############################################################################

use strict;      #require explicit scopes for all vars
use warnings;    #print debug warnings to console

use Cwd qw(abs_path);
use File::Basename;
use lib dirname(dirname abs_path $0)."/lib";
use lib dirname(dirname abs_path $0)."/perllib";

use Log::Log4perl;

#for command line argument parsing
use Getopt::Long;

#for CSV parsing
use Text::CSV;
use File::Basename;
use Bio::SeqIO;
use List::Util qw(min max);

use InsertionSeq::ContigBlastHit;
use InsertionSeq::AnnotationBlastHit;
use InsertionSeq::Annotation;
use InsertionSeq::Flank;

use InsertionSeq::Defaults qw(
    $BLASTALL_PATH
    $FORMATDB_PATH
    $EXTRACTSEQ_PATH
    $DEFAULT_REQ_IS_PERCENT_ID
    $DEFAULT_REQ_FLANK_PERCENT_ID
);


my $DEFAULT_CONFIG_FILE = (dirname abs_path $0)."/ISmapper.cfg";


my $DEFAULT_OUTPUT_PATH		= "output";
my $TMP_EXTRACTION_FILE 	= "seqTemp.fasta";
my $logfile_path;
my $output_path = $DEFAULT_OUTPUT_PATH;;
my $loglevel;
my $log;

my @annotation_file_extentions = (
	".coords.csv",
	".csv",
	"_coords.csv",
	".ptt",
	);
	
#global variables
use vars qw($is_path $genome_path $log_level $log_path $genome_name @plasmid_annot @annot_fastas @annotations %annotationHash $output_path %blast_hash);


my %LOG_LEVELS = 
(
	0 => "WARN",
	1 => "INFO",
	2 => "DEBUG",
	3 => "TRACE"
);




main();



sub init_logging
{
	my ($genome_name, $is_name, $loglevel_param) = @_;
	$logfile_path = "$output_path/ISmapper-$genome_name-$is_name-log.txt";
	
	$loglevel = "INFO";
	if (defined($loglevel_param))
	{
		$loglevel_param = 3 if ($loglevel_param > 3);
		$loglevel_param = 0 if ($loglevel_param < 0);
		$loglevel = $LOG_LEVELS{$loglevel_param};
	}

 	my $conf = 
		"log4perl.rootLogger=$loglevel,Logfile,Screen\n".
		
	
		"log4perl.appender.Logfile=Log::Log4perl::Appender::File\n".
		"log4perl.appender.Logfile.filename= $logfile_path\n".
		"log4perl.appender.Logfile.mode=write\n".
		"log4perl.appender.Logfile.layout=PatternLayout\n".
		"log4perl.appender.Logfile.layout.ConversionPattern=[%d] %m\n".
		"log4perl.appender.Logfile.Threshold=DEBUG\n".
	
		"log4perl.appender.Screen         = Log::Log4perl::Appender::Screen\n".
		"log4perl.appender.Screen.stderr  = 0\n".
		"log4perl.appender.Screen.layout=PatternLayout\n".
		"log4perl.appender.Screen.layout.ConversionPattern=%m\n".
		"log4perl.appender.Screen.Threshold=DEBUG\n";
		


  Log::Log4perl::init( \$conf );
        
  $log = Log::Log4perl->get_logger();
}



##
## Fail fast if we don't have the apps we need
##
sub checkForRequiredPrograms
{
	my $error = 0;

	if ( !-e $FORMATDB_PATH || !-x $FORMATDB_PATH )
	{
		print "Can't find formatdb at $FORMATDB_PATH\n";
		$error = 1;
	}
	if ( !-e $BLASTALL_PATH || !-x $BLASTALL_PATH )
	{
		print "Can't find blastall at $BLASTALL_PATH\n";
		$error = 1;
	}

	exit if ($error);
}


sub sizeFasta
{
    my $fileName = shift @_;

    my %result;
    
    my $reader=new Bio::SeqIO(-format=>'fasta',-file=>$fileName);
	while (my $seqRec=$reader->next_seq)
	{
		$result{$seqRec->id} = length($seqRec->seq());
    }

    return %result;
}

sub formatDB
{
	my $fileName = shift @_;

	my $formatDBOut = system("$FORMATDB_PATH -i $fileName -p F -o T");

	$log->logdie ("Error running $FORMATDB_PATH $fileName\n") unless ($formatDBOut == 0);

}

sub isReadableTextFile
{
	my $file = shift;
	return ( -e $file || -r $file || -T $file ) ;
}



##
## Read in annotations
## Verify matching Fasta exists
## Sort by location
## Remove identical dups
##
sub readAnnotations
{
	my $type = shift;
	my @files = @_;

	my @results;

	
	$log->info("Processing ".($#files+1)." annotation(s)......");
	
	for my $testFile (@files)
	{	
		my $fastaFile;
		my @annotationTargetNames;
		
		##
		## Test FASTA file
		##
		
		$log->info("Testing annotation fasta file name $testFile......");
		if (isReadableTextFile($testFile))
		{
			$fastaFile = $testFile;
			$log->info("OK\n");
			#last;
		}
		else
		{
			$log->info("Not found\n");
		}
	
		
		$log->logdie ("Could not find fasta file for annotatation fasta $testFile\n") unless ($fastaFile);
		
		
		formatDB($fastaFile);
		
		## Extract targets from file				
		my $reader=new Bio::SeqIO(-format=>'fasta',-file=>$fastaFile);
		while (my $seqRec=$reader->next_seq)
		{
			$log->info("Found annotation fasta defline ".$seqRec->id."\n");
			## Take 4the field if we see upbars
			my $targetName =  $seqRec->id;
			if ($seqRec->id =~ /^(.*)\|(.*)\|(.*)\|(.*)\|.*/)
			{
				$log->debug("($1) ($2) ($3) ($4)\n");
				$targetName = $4 ;
			}
			
		   $log->info("Using annotation target name: $targetName\n");
		   push @annotationTargetNames, $targetName;
		}
		

		##
		## Test annotation coords file
		##
		## Try foo.csv foo_coords.csv, etc
		##
		my $fastaDir = dirname($fastaFile);
		for my $targetName (@annotationTargetNames)
		{			
			my $annotationFile;
			for my $ext (@annotation_file_extentions)
			{		
				$testFile = $fastaDir."/".$targetName.$ext;	
				$log->info("Testing annotation file name $testFile......");
				if (isReadableTextFile($testFile))
				{
					$annotationFile = $testFile;
					$log->info("OK\n");
					last;
				}
				else
				{
					$log->info("Not found\n");
				}
			}
			
			$log->logdie ("Could not find coords annotation file for target $targetName from $fastaFile\n") unless ($annotationFile);
			##
			## Make the annotation and add it to the list
			##
			my $annotation = InsertionSeq::Annotation->new(					
						name => $targetName,
						type => $type,
						filename => $annotationFile,
						fasta_filename => $fastaFile
					);	
				
			## Print gene stats	
			my $annotSize = @{$annotation->{genelist}};
			$log->info("$annotation->{name} : $annotSize unique genes.\n");
			$log->logdie ("Error reading annotation file: no genes!\n") if ($annotSize <= 0);
			push @results, $annotation;
			
		}
		
		
		
		
	}
	
	@results;
}



#argument validation
sub main
{

	#local variables
	my $help;
	my $is_name;
	my $is_length;
	my $is_req_pct;
	my $flank_req_pct;
	my %genomeSizeHash;
	my %isSizeHash;
	my $config_file = $DEFAULT_CONFIG_FILE;
	$genome_path = "";
	$is_path     = "";

	#usage instructions for command line execution
	my $usage =
"USAGE: ./ISmapper.pl --genome <genome_fasta> --is <is_fstaile> --reference <reference_fasta> [--is_pct <float>] [--flank_pct <float>] [--out <output_dir>] [--config <config_file>] [--verbose] [--log <level>] [--help]\n"
	  . "\t--genome: path to query genome sequence in FASTA format\n"
	  . "\t--is: path to IS element sequence in FASTA format\n"
	  . "\t--reference: basename of reference FASTA file (with corresponding _coords files for each entry)\n"
	  . "\t--is_pct: percent ID cutoff for including an IS blast hit against a contig or the reference\n"
	  . "\t--flank_pct: percent ID cutoff for including a flank blast hit against the reference\n"
	  . "\t--out: path to output files\n"
	  . "\t--config: alternate config file with app locations\n"
	  . "\t--log: logging level (0-3)\n"
	  . "\t--help: print this usage statement\n";

	#read command line options for path to sequences and annotation behavior
	exit unless GetOptions(
		"help"     		=> \$help,
		"is=s"     		=> \$is_path,
		"genome=s" 		=> \$genome_path,
		"is_pct=f" 		=> \$is_req_pct,
		"flank_pct=f"  	=> \$flank_req_pct,
		"out=s" 		=> \$output_path,
		"reference=s@" 	=> \@annot_fastas,
		"log=i"  		=> \$log_level,
		"config=s"     	=> \$config_file
	);

	InsertionSeq::Defaults::Process( $config_file );

	$is_req_pct		= $DEFAULT_REQ_IS_PERCENT_ID unless(defined($is_req_pct));
	$flank_req_pct	= $DEFAULT_REQ_FLANK_PERCENT_ID unless (defined($flank_req_pct));

	#set path to genome info file

	checkForRequiredPrograms();

	if ( $is_path eq "" || $genome_path eq "" || $help )
	{
		print $usage;
		exit;
	}
	
	if ( !( -e $genome_path || -r $genome_path || -T $genome_path ) )
	{
		print $usage;
		die("ERROR: Could not open \"$genome_path\" - check if file exists, is readible, and is a FASTA file\n");
	}
	
	if ( !( -e $is_path || -r $is_path || -T $is_path ) )
	{
		print $usage;
		die("ERROR: Could not open \"$is_path\" - check if file exists, is readible, and is a FASTA file\n");
	}
	
	if ( !(-d $output_path) )
	{
   	 	mkdir $output_path || die "Cannot make output dir: $output_path\n";
	}
	
	my $param_string = "params:";
	$param_string .= " -is_pct $is_req_pct" if defined(($is_req_pct));
	$param_string .= " -flank_pct $flank_req_pct" if (defined($flank_req_pct));
	$param_string .= " -is $is_path";
	$param_string .= " -genome $genome_path";
	for my $outref (@annot_fastas)
	{
		$param_string .= " -reference $outref";
	}
	$param_string .= " -out $output_path";

	#get just the filename, and remove extension for printing purposes
	$genome_name = $genome_path;
	$genome_name =~ s/.*\///;
	$genome_name =~ s/\.fasta|\.fa|\.FA|\.FASTA//;

	$is_name = $is_path;
   	$is_name =~ s/.*\///;
   	$is_name =~ s/\.fasta|\.fa|\.FA|\.FASTA//;
    
    init_logging($genome_name, $is_name, $log_level);

	#set up log file

	$log->info("Writing log output to $logfile_path\n");
	
	$log->info("$param_string\n");
	$log->info("Writing log output to $logfile_path\n");
	$log->info("IS min percent    = $is_req_pct\n");
	$log->info("Flank min percent = $flank_req_pct\n");
	
	
	if (@annot_fastas)
	{
		push @annotations, readAnnotations( $InsertionSeq::Annotation::TYPE_CHROMOSOME, @annot_fastas )	;
	}
	else
	{
    	$log->logdie("Reference FASTA is required.\n");
	}

	##
	## build hash
	##
	for my $annotation (@annotations)
	{
    	$annotationHash{$annotation->{name}} = $annotation;
	}
	

	%isSizeHash = sizeFasta($is_path);

	$log->logdie("Multiple sequences in IS file $is_path\n\n") if ( scalar( keys %isSizeHash ) > 1 );

	$is_name   = ( keys %isSizeHash )[0];
	$is_length = $isSizeHash{$is_name};

	$log->debug("Running isfind on $genome_name for $is_path\n");

	#blast the IS element against the genome
	blast_is_genome( $is_name, $is_length, $is_path, $genome_path, $is_req_pct, $flank_req_pct );

}




#program subs
sub blast_is_genome()
{
	##
	## Params
	##
	my $is_name     	= shift(@_);
	my $is_length   	= shift(@_);
	my $is_path     	= shift(@_);
	my $genome_path 	= shift(@_);
	my $is_req_pct 		= shift(@_);
	my $flank_req_pct	= shift(@_);

	#local variable declare
	my $hit;

	my %genomeSizeHash = sizeFasta($genome_path);

	$log->info( "Processing IS $is_name SIZE $is_length\n" );

	my $annotationISSQL = "";
	for my $fasta (@annot_fastas)
    {
 		my @annotationISHits =  blast_is_against_annotation( $is_name, $is_length, $is_path, $fasta, $is_req_pct );

		for my $anBlastHit (@annotationISHits)
		{

			##
			## Match target name to which annotation we hit
			##
			my $annotation = $annotationHash{$anBlastHit->{name}};
			$log->logdie("Error - cannot match annotataion blast target name $anBlastHit->{name} in list of annotations!\n\n") unless ($annotation);
			##
			## Store entire IS seqs found
			##
			if ($anBlastHit->{type} eq $InsertionSeq::ContigBlastHit::TYPE_ENTIRE || $anBlastHit->{type} eq $InsertionSeq::ContigBlastHit::TYPE_EMBED)
			{
           		$annotationISSQL .= $anBlastHit->to_feat_mysql($is_name, $annotation->{name}, $annotation->{name});
			}
		}
    }


	$log->debug("\nBLASTing IS element against genome\n");

	formatDB($genome_path);
	
	#run blast of seq against genome
	my $command = "$BLASTALL_PATH -p blastn -d $genome_path -e .01 -m 8 -i $is_path";
	$log->info("$command\n");
	my $blastout = `$command`;
	$log->info("RAW BLAST OUTPUT:\n");
	$log->info("queryid\t\tsubjectid\t%id\tmatch\tmis\tgaps\tqstart\tqend\tsstart\tsend\te\tbit\n");
	$log->info("$blastout");


	#loop over each BLAST result and parse
	$log->info("\nIS ANALYSIS:\n");

	#pull each line of blast output into an array
	my @blastout = split( /\n/, $blastout );

	#loop by each blast result
	foreach $hit (@blastout)
	{

		#parse blast line by the tabs
		my @hit = split( /\t/, $hit );

		#extract data from the parse
		my $dir       = "";
		my $queryid   = $hit[0];
		my $subjectid = $hit[1];
		my $pct_id 	  = $hit[2];
		my $qstart    = $hit[6];
		my $qend      = $hit[7];
		my $sstart    = $hit[8];
		my $send      = $hit[9];

		#query length extracted from fasta file
		my $qlength = $is_length;

		##
		## Pull subject length from hash
		##
		my $slength = $genomeSizeHash{$subjectid};
		$log->logdie ("Error no size reference for $subjectid\n") unless defined($slength);

		my $hit = InsertionSeq::ContigBlastHit->new(
			q_name    => $is_name,
			is_name   => $is_name,
			pct_id    => $pct_id,
			req_pct_id  => $is_req_pct,
			queryid   => $queryid,
			qstart    => $qstart,
			qend      => $qend,
			qlength   => $qlength,
			subjectid => $subjectid,
			sstart    => $sstart,
			send      => $send,
			slength   => $slength,
			name 	  => $subjectid,
			filename => $genome_path
		);

		#$hit->initialize();

		push( @{ $blast_hash{$subjectid} }, $hit );
	

	}    #end BLAST loop foreach


	$log->info("*  = Rejected/Non-annotatable flank type\n");
	$log->info(" + = Rejected/Failed Pct ID threshold\n\n");
	
	my @flanks;

	
	## Calc max field widths
	my $maxSubWidth = 10;
	for my $subjectid ( sort keys %blast_hash )
    {
		for $hit ( sort { $a->{sstart} <=> $b->{sstart} }
                    @{ $blast_hash{$subjectid} }
          )
        {
        	$maxSubWidth = length($hit->{subjectid}) if (length($hit->{subjectid}) > $maxSubWidth);
		}
	}
	
	for my $subjectid ( sort keys %blast_hash )
	{
		## Sort by contig, then position
		for $hit ( sort { $a->{sstart} <=> $b->{sstart} } 
					@{ $blast_hash{$subjectid} }
		  )
		{		
			my $col1 = "";
			my $col2 = "";
			my $col3 = "";
			my $col4 = "";
			my $col5 = "";
			my $col6 = "";
			my $col7 = "";
			my $col8 = "";

			##
			##  Flag un-annotatable/poor quality hits
			##
			if ($hit->is_annotatable())
			{
				$col1 .= " ";
			}
			else
			{
				$col1 .= "*";
			}
			if ($hit->passes_pct_threshold())
			{
				$col1 .= " ";
			}
			else
			{
				$col1 .= "+";
			}

			my $tmplen = abs($hit->{sstart}-$hit->{send} )+1;
			$col2 = $hit->{subjectid};

			$col3 = $hit->get_pct_id() ; 
			$col4 = $hit->{type} ; 
			$col5 = "($hit->{sdir})";
			$col6 = "[$hit->{sstart} to $hit->{send} of $hit->{slength}]";
			$col7 = "(size=$tmplen)";
			

			if ($hit->is_annotatable() && $hit->passes_pct_threshold())
			{
				##
				## TODO = mate 2 flanks from the same contig automatically
				##
				my $beginFlank = $hit->get_flank($InsertionSeq::ContigBlastHit::BEGINNING_FLANK);
				if (defined($beginFlank))
				{
					$col8 .= "\tUpstream - $beginFlank->{contig_lower_coord}-$beginFlank->{contig_upper_coord}";
					push @flanks, $beginFlank;
				}
				
				my $endFlank = $hit->get_flank($InsertionSeq::ContigBlastHit::ENDING_FLANK);
				if (defined($endFlank))
				{
					$col8 .=  "\tDownstream - $endFlank->{contig_lower_coord}-$endFlank->{contig_upper_coord}";
					push @flanks, $endFlank;
				}	
				#if (defined($beginFlank) && defined($endFlank))
				#{
					#$logMessage .=" [Auto-pairing flanks from same IS in same contig.] ";
					#$beginFlank->mate_with($endFlank);
				#}	
		    }	
			my $logMessage = sprintf("%-3s %-".$maxSubWidth."s   %6s%%   %-11s   %6s   %-30s   %-10s   %-30s", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8);
			$log->info($logMessage."\n");
		}
	}
	$log->info("\n");
	
	

	$log->info("\nExtracting flanks and blasting references....\n");
	$log->info("RAW FLANK ANNOTATION BLAST OUTPUT:\n");
	if (@flanks)
	{
		$log->info("queryid\t\tsubjectid\t%id\tmatch\tmis\tgaps\tqstart\tqend\tsstart\tsend\te\tbit\n");
		for my $flank (@flanks)
		{
			my $command = "$EXTRACTSEQ_PATH $genome_path:$flank->{contig_name} -regions $flank->{contig_lower_coord}-$flank->{contig_upper_coord} -out $TMP_EXTRACTION_FILE -sid1 ".$flank->get_sid_name();
			$log->debug("$command\n");
			my $result = system("$command 2>/dev/null");
			$log->logdie ("Error $result while running: $command\n") if ($result);
			## Blast here
			#for my $annotation (@annotations)
			for my $annotation (@annot_fastas)
			{	
				blastAnnotation ($TMP_EXTRACTION_FILE, $flank, $annotation, $flank_req_pct  );
			}
		}
	}
	else
	{
    	$log->info("<NONE>\n");
	}
	
	$log->info("\nANNOTATION-BLASTED FLANKS:\n");
	$log->info("+ = Rejected/Failed Pct ID threshold\n\n");
	#$log->info("name contig_dir flank_loc annot_dir closest_base annotation\n");

	my $maxSidLen = 10;
	for my $flank (@flanks)
	{
		my $best = $flank->pick_best_blast();
       	$maxSidLen = length($flank->get_sid_name()) if (length($flank->get_sid_name()) > $maxSidLen);
	}

    $log->info(sprintf("%-3s %-".$maxSidLen."s   %6s%%   %-5s %-5s %-5s  %-20s   %-35s\n",
			"  ", "Contig Name", "PctID", "CDir", "FType", "FDir", "Closest Base", "Annotation"));
	for my $flank (@flanks)
	{
		my $best = $flank->pick_best_blast();
   		my $reject = " ";
	    $reject = "+" unless ($best && $best->passes_pct_threshold());

		my $logString = sprintf("%-3s %-".$maxSidLen."s   ",$reject, $flank->get_sid_name());
		if (defined($best))
		{
			$logString .= sprintf("%6s%%   %-5s %-5s %-5s  %-20s   %-35s\n",
				 $best->get_pct_id(), $flank->{contig_direction}, $flank->{location}, $best->{sdir},$best->closest_is_base(),$flank->{annotation_name});
		}
		else
		{
			$logString .= " -------------------- NO BLAST HITS --------------------\n";
		}
		$log->info($logString);
		#$log->info($flank->get_sid_name()." ".scalar(@{$flank->{blast_hits}})."\n");
		
	}

	my @bad_flanks;
	my @good_flanks;
	for my $flank (@flanks)
	{
		my $blast = $flank->pick_best_blast();
		if (!defined($blast) || ! $blast->passes_pct_threshold())
		{
			push @bad_flanks, $flank ;
		}
		else
		{
			push @good_flanks, $flank ;
		}
	}
	
	
	##
	## Sort by annotation and location and type
	## Mark lower blast scores as dups before higher
	## Mark PARTIAL matches as dups before WHOLEs
    ## Finally by contig name to keep things repeatable in case of ties
	##
	
	@good_flanks = sort {  
							$a->{annotation_name} cmp $b->{annotation_name} 
							|| $a->{closest_is_base} <=> $b->{closest_is_base} 
							|| $a->pick_best_blast()->get_pct_id() <=> $b->pick_best_blast()->get_pct_id()
							|| $a->pick_best_blast()->type() cmp $b->pick_best_blast()->type() 
							|| $b->pick_best_blast()->get_name() cmp $a->pick_best_blast()->get_name() 
				 		}  @good_flanks ;
	
	##
	## Calculate offsets from previous flank			 		
	##
	my $basesDiff = "";
	my $lastAnnot = "";
	my $lastLoc = -1;
   	my $last_flank;
	for my $flank (@good_flanks)
	{
		my $best = $flank->pick_best_blast();
		$basesDiff = 0;
		$lastLoc = -1 if ($flank->{annotation_name} ne $lastAnnot);
		$basesDiff = ($best->closest_is_base() - $lastLoc) if (($flank->{annotation_name} eq $lastAnnot) && $lastLoc >= 0);		
		$flank->{offset_from_last} =  $basesDiff;
		$lastLoc = $best->closest_is_base();
		$lastAnnot = $flank->{annotation_name};

		##
		## Flanks are sorted with Partial matches before Whole,
		## So mark the earlier one as the dup
	    $flank->{flank_dup} = 0;
		if (defined($last_flank) && $flank->is_dup_of($last_flank))
		{
			$last_flank->{flank_dup} = 1;
		}
		$last_flank = $flank;
	}

	##
	## Tally dups for "master flank"
	##
	for my $flank (@good_flanks)
	{
		if (!$flank->{flank_dup})
		{
			for my $test_dup_flank (@good_flanks)
			{
				if ($test_dup_flank->{flank_dup} && $flank->is_dup_of($test_dup_flank))
				{
                	push( @{$flank->{dup_list}}, $test_dup_flank);
				}
			}
		}
	}


	##
	## Tag flank pairs
	##
	my $lastFlank;
	for my $flank (@good_flanks)
	{
		if (!$flank->{flank_dup}
			&& defined($lastFlank) 
			&& !defined($flank->{mate}) 
			&& !defined($lastFlank->{mate}) 
			&& $flank->is_mate_of($lastFlank))
		{ 
			$flank->mate_with($lastFlank); 	
		}
	
		$lastFlank = $flank;	
	}
	


	my $current_annotation = "PLACEHOLDER_NO_ANNOTATION";
	## Sort flanks by annotation, then by location
	## (skip flanks with no best hit)
	
	
	
	for my $flank (@good_flanks)
	{
		##
		## Fresh Header, if this is a new annotation
		##
		if ($flank->{annotation_name} ne $current_annotation)
		{
			$current_annotation = $flank->{annotation_name};
			$log->info("\nANNOTATED FLANKS ($current_annotation):\n");
			
		}
		my $best = $flank->pick_best_blast();

		my $flank_desc_line;

		
		
	   	$flank_desc_line .= "*DUP* " if ($flank->{flank_dup});
		if ( defined($flank->{mate}) 
			&& (   ($flank->{direction} eq $InsertionSeq::Flank::DIRECTION_FWD && $flank->{is_location} eq $InsertionSeq::Flank::LOCATION_BEGIN)
				|| ($flank->{direction} eq $InsertionSeq::Flank::DIRECTION_REV && $flank->{is_location} eq $InsertionSeq::Flank::LOCATION_END)
				)
		)
		{
			$flank_desc_line .= "  _/--";
		}
		elsif (defined($flank->{mate}) 
			&& (($flank->{direction} eq $InsertionSeq::Flank::DIRECTION_REV && $flank->{is_location} eq $InsertionSeq::Flank::LOCATION_BEGIN)
				|| ($flank->{direction} eq $InsertionSeq::Flank::DIRECTION_FWD && $flank->{is_location} eq $InsertionSeq::Flank::LOCATION_END)
				)
				)
		{
			$flank_desc_line .= "   \\--";
		}
		else
		{
	   		$flank_desc_line .= "      " unless ($flank->{flank_dup});
		}
			
		$flank_desc_line .= "$flank->{annotation_name} $flank->{direction} $flank->{is_location} ".$best->type()." ";
		
		$flank_desc_line .= sprintf("%6s%%",$best->get_pct_id());
		$flank_desc_line .= sprintf("%10d",$best->closest_is_base());
		$flank_desc_line .= sprintf("%20s", " [".sprintf("%8s","+".$flank->{offset_from_last}) . "]    " );
		$flank_desc_line .= sprintf("%-75s", $best->{annotation}->locate_base_text($best->closest_is_base()) );
		$flank_desc_line .= sprintf("%-40s"," [".$flank->get_sid_name()."]" );

		$flank_desc_line .= " (".scalar(@{$flank->{dup_list}})." duplicates) " if ($flank->{dup_list});

		## Redundant with the left side ** DUP ** 
	    ##$flank_desc_line .= " ***** FLANK DUPLICATE ***** " if ($flank->{flank_dup});

		$log->info("$flank_desc_line\n");
		
	}
	
	$log->info("\nUNANNOTATED FLANKS:\n");
	for my $flank (sort {  $a->{contig_name} cmp $b->{contig_name} } @bad_flanks)
	{
		my $best = $flank->pick_best_blast();
		my $pct = "0";
		$pct = $best->get_pct_id() if ($best);
		$log->info(sprintf("%-40s %6s%% %-20s %-20s %-15s\n", $flank->{contig_name}, $pct, "($flank->{contig_lower_coord}-$flank->{contig_upper_coord})", "Orientation: $flank->{contig_direction}", "Begin/End: $flank->{location}"));
	}
	$log->info("<NONE>\n") if (-1 == $#bad_flanks);
	$log->info("\n");


	my $anno_file_ext = "";
	if (@annotations)
	{
		$anno_file_ext = "-".$annotations[0]->{name};
	}

	
	
	
	my $source_genome = basename($genome_path);
	$source_genome =~ s/\..*//;
	
	##
	## Delete previous data
	##
	my $cleanSQL;
	$cleanSQL .= "## \n## DELETE OLD ROWS\n##\n";
	$cleanSQL .= "START TRANSACTION;\n";

	my $isrunid_sql =  "select distinct is_run_id from is_run where is_run_element = '$is_name' and  is_run_genome = '$source_genome' and is_run_references = '".join(",",sort keys %annotationHash)."'";
	
	$cleanSQL .= "## Unlink is_flanks\n";
	$cleanSQL .= "UPDATE is_flank set mate_flank_id = NULL where is_run_id = ($isrunid_sql);\n";

	$cleanSQL .= "## DELETE is_flanks\n";
	$cleanSQL .= "DELETE FROM is_flank where is_run_id = ($isrunid_sql);\n";

    $cleanSQL .= "## DELETE is_query_feature\n";
   	$cleanSQL .= "DELETE FROM is_query_feature where is_run_id = ($isrunid_sql);\n";

   	$cleanSQL .= "## DELETE is_run \n";
   	$cleanSQL .= "DELETE FROM is_run where is_run_element = '$is_name' and  is_run_genome = '$source_genome' and is_run_references = '".join(",",sort keys %annotationHash)."';\n";
	$cleanSQL .= "COMMIT;\n\n\n";
	

	$cleanSQL .= "START TRANSACTION;\n";

	$cleanSQL .= "##\n## Inserting master project record\n##\n";
	$cleanSQL .= "INSERT INTO is_run(is_run_element,is_run_genome,is_run_references,is_req_pct_id,flank_req_pct_id)\n";
	$cleanSQL .= " VALUES ('$is_name','$source_genome','".join(",",sort keys %annotationHash)."',$is_req_pct,$flank_req_pct);\n";
	$cleanSQL .= 'SET @is_run_id = LAST_INSERT_ID();';
	$cleanSQL .= "\n";

	$cleanSQL .= "\n## IS sequences found in reference directly: \n";
	$cleanSQL .= $annotationISSQL;
	$cleanSQL .= "\n";
	
	$cleanSQL .= "## \n## Inserting paired flanks in is_flank with supporting is_query_feature rows\n##\n";
	#
	# Insert Ends before Begins
	##
	for my $flank (sort { $b->{is_location} cmp $a->{is_location} } @good_flanks)
	{

		if (!$flank->{flank_dup})
		{
			##
			## Only ENDS 
			##
			if ($flank->{is_location} eq $InsertionSeq::Flank::LOCATION_END)
			{
				my $flank_distance = "NULL";
				$flank_distance = abs ($flank->{closest_is_base} - $flank->{mate}->{closest_is_base} ) if (defined($flank->{mate}));
				$cleanSQL .= $flank->to_is_mysql($flank_distance,"NULL",$source_genome);
				if (defined($flank->{mate}))
				{
					##
					## BEGINS with matches are added here
					##
					$log->debug("Adding linked BEGIN IS...\n");
					## Add SQL for BEGIN mate using LAST_INSERT_ID for link
					$cleanSQL .= $flank->{mate}->to_is_mysql($flank_distance, $InsertionSeq::Flank::SQL_LAST_IS_FLANK_ID_VAR ,$source_genome);
				}
			}
			elsif (!defined($flank->{mate}))
			{
				##
				## Only BEGINS with no matches here
				##
				my ($distance, $id);
				$cleanSQL .= $flank->to_is_mysql($distance,$id,$source_genome);
			}
		}

	}
	
	$cleanSQL .= "## \n## Cross Link Mated Flanks\n##\n";
	$cleanSQL .=  "UPDATE is_flank a, is_flank b SET a.mate_flank_id = b.is_flank_id where b.mate_flank_id = a.is_flank_id and a.mate_flank_id IS NULL;\n\n";


  	##
   	## Write un-annotated flanks (once per annotation it was run against)
   	##
   	$cleanSQL .= "## \n## Save un-annotated flanks in is_query_feature table\n##\n";
	for my $annotation (@annotations)
	{
   		for my $flank (sort {  $a->{contig_name} cmp $b->{contig_name} } @bad_flanks)
   		{
   	  		#$cleanSQL .= $flank->to_feat_mysql($is_name, $source_genome, $annotation->{name} );
   	  		$cleanSQL .= $flank->to_feat_mysql($source_genome, $annotation->{name}, $annotation->{name} );
   		}
	}
   
   	## 
   	## Write contigs that were all IS and thus un-annotatable
   	##
   	$cleanSQL .= "## \n## Save all complete IS contigs (un-annotatable, no flanks available) in is_query_feature table\n##\n";
	for my $annotation (@annotations)
	{
   		for my $subjectid ( sort keys %blast_hash )
   		{
       		## Sort by contig, then position
       		for $hit ( sort { $a->{sstart} <=> $b->{sstart} }
                   		@{ $blast_hash{$subjectid} }
         		)
       		{
				##
				## Here we only care about contigs that are completely IS
				##
				if ($hit->{type} eq $InsertionSeq::ContigBlastHit::TYPE_ENTIRE)
				{
   					$cleanSQL .= $hit->to_feat_mysql($is_name, $source_genome, $annotation->{name});
				}
   			}
   		}
	}


	$cleanSQL .= "COMMIT;\n\n\n";
	$log->debug("SQL:\n");
	$log->debug($cleanSQL);

	my $sqlfile_path = "$output_path/ISmapper-$genome_name-$is_name$anno_file_ext.sql";
	$log->info("Writing SQL file $sqlfile_path\n");
	open SQL, "> $sqlfile_path" || $log->logdie ("Cannot open SQL file $sqlfile_path\n");
	print SQL $cleanSQL;
	close SQL;

	my $make_csv_file = 0;
	if ($make_csv_file)
	{
		my $csv_file_path = "$output_path/ISmapper-$genome_name-$is_name$anno_file_ext.csv";
    	$log->info("Writing CSV file $csv_file_path\n");
    	open CSV, "> $csv_file_path" || $log->logdie ("Cannot open CSV file $sqlfile_path\n");
		print CSV "$InsertionSeq::Flank::CSV_HEADER";

		for my $flank (@good_flanks)
		{
			my $flank_distance = "";
			$flank_distance = abs ($flank->{closest_is_base} - $flank->{mate}->{closest_is_base} ) if (defined($flank->{mate}));
			print CSV $flank->to_csv($flank_distance,"",$source_genome);
		}


    	for my $flank (@bad_flanks)
    	{
			print CSV $flank->to_csv("","",$source_genome);
   		}


		close CSV;

	}


	
}    #end sub


sub blastAnnotation
{
	my $flank_seq_file = shift;
	my $flank = shift;
	my $annotation_fasta = shift;
	my $flank_req_pct = shift;
	
	$log->debug("Blast $flank->{contig_direction} $flank->{contig_name} \n");
	
	my $command ="$BLASTALL_PATH  -p blastn -d $annotation_fasta  -e .0000000001 -m 8 -i $flank_seq_file";
	$log->debug("$command\n");
	
	my $blastout = `$command`;
	$log->info("$blastout");
	
	my @blastout = split( /\n/, $blastout );

	#loop by each blast result
	foreach my $hit (@blastout)
	{
		#parse blast line by the tabs
		my @hit = split( /\t/, $hit );

		#extract data from the parse
		my $dir       = "";
		my $queryid   = $hit[0];
		my $subjectid = $hit[1];
	    my $pct_id    = $hit[2];
		my $qstart    = $hit[6];
		my $qend      = $hit[7];
		my $sstart    = $hit[8];
		my $send      = $hit[9];
		

		 my $targetName =  $subjectid;
		##
		## Genbank style = grab accession
		##
         if ($targetName =~ /^(.*)\|(.*)\|(.*)\|(.*)\|.*/)
         {
             $targetName = $4 ;
         }

		my $annotation = $annotationHash{$targetName};
		$log->logdie ("Error matching blasted annotation $targetName to annotation list\n")  unless defined($annotation);

		if (!defined($annotation->{size}))
		{
			my %newSizes = sizeFasta($annotation_fasta);

			for my $key (keys %newSizes)
			{
				##
				## Genbank style = grab accession
				##
				my $name = $key;
         		if ($name =~ /^(.*)\|(.*)\|(.*)\|(.*)\|.*/)
         		{
             		$name = $4 ;
         		}
				my $annotationSizer = $annotationHash{$name};
				$log->logdie("Found unexpected annotation $name while looking up sizes!\n") unless ($annotationSizer);
            	$annotationSizer->{size} = $newSizes{$key};
			}
		}
		$log->logdie("Failed to find size for annotation : $annotation->{name}!\n") unless (defined($annotation->{size}));


		my $hit = InsertionSeq::AnnotationBlastHit->new(
			q_name    => $flank->get_sid_name(),
			is_name   => $flank->{is_name},
			pct_id    => $pct_id,
            req_pct_id  => $flank_req_pct,
			queryid   => $queryid,
			qstart    => $qstart,
			qend      => $qend,
			qlength   => $flank->{length},
			subjectid => $subjectid,
			sstart    => $sstart,
			send      => $send,
			slength   => $annotation->{size},
			name 	  => $subjectid,
			filename  => $annotation_fasta,
			annotation => $annotation
		);
		
		$flank->add_blast_hit($hit);
	}
	
}


sub blast_is_against_annotation()
{
    ##
    ## Params
    ##
    my $is_name     		= shift(@_);
    my $is_length   		= shift(@_);
    my $is_path    		 	= shift(@_);
    my $annot_fasta_path  	= shift(@_);
    my $is_req_pct			= shift(@_);

    my %genomeSizeHash = sizeFasta($annot_fasta_path);

    print "\nBLASTING IS $is_name SIZE $is_length against annotation FASTA $annot_fasta_path\n";

    #for the first BLAST of the IS against the annotation fasta
    formatDB($annot_fasta_path);

    #run blast of seq against annotation
    my $command = "$BLASTALL_PATH -p blastn -d $annot_fasta_path -e .0000000001 -m 8 -i $is_path";
    $log->info("$command\n");
    my $blastout = `$command`;
    $log->info("IS vs. ANNOTATION BLAST OUTPUT (pct id cutoff = $is_req_pct):\n");
    $log->info("queryid\t\tsubjectid\t%id\tmatch\tmis\tgaps\tqstart\tqend\tsstart\tsend\te\tbit\n");
    $log->info("$blastout");

    #pull each line of blast output into an array
    my @blastout = split( /\n/, $blastout );

	my @annotationList = ();
    #loop by each blast result
    foreach my $blast_line (@blastout)
    {

        #parse blast line by the tabs
        my @hit = split( /\t/, $blast_line );

        #extract data from the parse
        my $dir       = "";
        my $queryid   = $hit[0];
        my $subjectid = $hit[1];
        my $pct_id	  = $hit[2];
        my $qstart    = $hit[6];
        my $qend      = $hit[7];
        my $sstart    = $hit[8];
        my $send      = $hit[9];

        #query length extracted from fasta file
        my $qlength = $is_length;

        ##
        ## Pull subject length from hash
        ##
        my $slength = $genomeSizeHash{$subjectid};
        $log->logdie( "Error no size reference for $subjectid\n") unless defined($slength);

       	my $contig_blast_hit = InsertionSeq::ContigBlastHit->new(
           	q_name    	=> $is_name,
           	is_name   	=> $is_name,
   			pct_id    	=> $pct_id,
   			req_pct_id	=> $is_req_pct,
           	queryid   	=> $queryid,
           	qstart    	=> $qstart,
           	qend      	=> $qend,
           	qlength   	=> $qlength,
           	subjectid 	=> $subjectid,
           	sstart    	=> $sstart,
           	send      	=> $send,
           	slength   	=> $slength,
           	name      	=> $subjectid,
           	filename	=> $annot_fasta_path
       	);


       	push( @annotationList, $contig_blast_hit );


    }    #end BLAST loop foreach
	my $hitcount = scalar @annotationList ;

	my $valid_hit_count = 0;
	for my $hit_check (@annotationList)
	{
    	$valid_hit_count++ if ($hit_check->passes_pct_threshold());
	}

	my $hitInfo =       "IS annotation hits      : ";
	my $validHitInfo = "Valid IS annotation hits: ";

	$hitcount = "<NONE>" if ($hitcount == 0);
	$hitInfo .= $hitcount;

	$valid_hit_count = "<NONE>" if ($valid_hit_count == 0);
	$validHitInfo .= $valid_hit_count;

	$log->info("$hitInfo\n$validHitInfo\n\n");

	@annotationList;
}

