package InsertionSeq::Flank;

#############################################################################
#
#    ISseeker    - Finds portions of contigs flanking IS sequences
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

use strict;

use InsertionSeq::ContigBlastHit;
use Log::Log4perl;
 
#for CSV parsing
use Text::CSV;



##
## Is this adjacent to the beginnng or end of the IS seq?
##
our $LOCATION_BEGIN		= "B";
our $LOCATION_END   	= "E";

##
## What direction is the IS with respect to the contig or annotation
##
our $DIRECTION_FWD		= "F";
our $DIRECTION_REV		= "R";


##
##  A Whole match is within this many from each end:
##
our $ANNOT_BASE_WIGGLE_ROOM		= 10;

##
## Max distance we can consider to flanks mates
##
our $MAX_FLANK_PAIRING_DISTANCE = 15;

our $SQL_FLANK_FEAT_TYPE = 'F';

our $SQL_LAST_IS_FLANK_ID_VAR	= '@last_is_flank_id';
our $SQL_LAST_FEAT_ID_VAR		= '@last_is_feat_id';


my $log = Log::Log4perl->get_logger();


sub new 
{
    my $class = shift;
    my $self = {@_};

    bless($self, $class);

	die ("Flank missing flank location.\n") if ( !defined($self->{location}) );
	die ("Flank invalid flank location:($self->{location}).\n") if ( $self->{location} ne $LOCATION_BEGIN && $self->{location} ne $LOCATION_END);
	
	##
	## IS location is opposite of flank location
	## (an ending flank means IS is in front of it)
	##
	$self->{is_location} = $LOCATION_END if ($self->{location} eq $LOCATION_BEGIN);
	$self->{is_location} = $LOCATION_BEGIN if ($self->{location} eq $LOCATION_END);
	
	die ("Flank missing contig direction.\n") if ( !defined($self->{contig_direction}) );
	die ("Flank invalid contig direction:($self->{contig_direction}).\n") if ( $self->{contig_direction} ne $DIRECTION_FWD && $self->{contig_direction} ne $DIRECTION_REV);
	
	die ("Flank missing contig lower coord.\n") if ( !defined($self->{contig_lower_coord}) );
	die ("Flank missing contig upper coord.\n") if ( !defined($self->{contig_upper_coord}) );
	die ("Flank missing contig name.\n") if ( !defined($self->{contig_name}) );
	die ("Flank missing IS name.\n") if ( !defined($self->{is_name}) );
	
	die ("Flank missing contig fasta file name.\n") if ( !defined($self->{contig_file_name}) );
	
	$self->{length} = ($self->{contig_upper_coord} - $self->{contig_lower_coord}) +1;


	$self->{direction} = "";


	

	if ($self->{contig_name} =~ /\|/)
	{
		##
		## Remove pipes from contig name
		##
		if ($self->{contig_name} =~ /(gb|ref)\|(\w+\.\w+)/)
		{
			##
			## opt for full accession first
			##
#$log->info("Extracting short name $2 from $self->{contig_name} \n");
			$self->{contig_name} = $2;
		}
		elsif ($self->{contig_name} =~ /(gb|ref)\|(\w+)/)
		{
			##
			## opt for accession 
			##
#$log->info("Extracting short name $2 from $self->{contig_name} \n");
			$self->{contig_name} = $2;
		}
		elsif ($self->{contig_name} =~ /gi\|(\w+)/)
		{
			## 
			## Or take GI
			##
#$log->info("Extracting short name $1 from $self->{contig_name} \n");
			$self->{contig_name} = $1;
		}
		elsif ($self->{contig_name} =~ /^>(\w+|\.)\|/)
		{
			## Take what's first
#$log->info("Extracting short name $self->{contig_name} \n");
		}
		else
		{
			$log->logdie("Can't un-pipe this genome def line: $self->{contig_name}\n");
		}
	}

	return $self;
}



##
## Long name that includes source_type_contig coords
##
sub get_sid_name
{
	my $self = shift;	

	return $self->{contig_name}."_".$self->{location}.$self->{contig_direction}."_".$self->{contig_lower_coord}."_".$self->{contig_upper_coord};
}

##
## Determine the end of the blast hit closest to the IS
## Based on the original contig match orientation
## And the annotation match orientation
## And whether this is a beginning flank or an ending flank
##
sub get_location_closest_to_is
{
	my $self = shift;
	my $blast = shift;
	
	my $blast_hit_direction = $blast->{sdir};
	my $lower = $blast->{sstart};
	my $upper = $blast->{send};	
	
	
	die "Flank: unknown blast hit direction: $blast_hit_direction" if ( $blast_hit_direction ne $DIRECTION_FWD && $blast_hit_direction ne $DIRECTION_REV);
	
	if ($blast_hit_direction eq $DIRECTION_FWD)
	{
		if ($self->{contig_direction} eq $DIRECTION_FWD)
		{
			return $upper if ($self->{location} eq $LOCATION_BEGIN);
			return $lower if ($self->{location} eq $LOCATION_END);
		}
		elsif ($self->{contig_direction} eq $DIRECTION_REV)
		{
			## Reversed
			return $lower if ($self->{location} eq $LOCATION_BEGIN);
			return $upper if ($self->{location} eq $LOCATION_END);
		}
	}
	elsif ($blast_hit_direction eq $DIRECTION_REV)
	{
		## REV-REV == the same as FWD-FWD
		if ($self->{contig_direction} eq $DIRECTION_REV)
		{
			return $upper if ($self->{location} eq $LOCATION_BEGIN);
			return $lower if ($self->{location} eq $LOCATION_END);
		}
		elsif ($self->{contig_direction} eq $DIRECTION_FWD)
		{
			return $lower if ($self->{location} eq $LOCATION_BEGIN);
			return $upper if ($self->{location} eq $LOCATION_END);
		}
	}
}

sub add_blast_hit
{
	my $self = shift;
	my $blast_hit = shift;
	
	if ($blast_hit->{qstart} <= $ANNOT_BASE_WIGGLE_ROOM 
		&& $blast_hit->{qlength} - $blast_hit->{qend} <= $ANNOT_BASE_WIGGLE_ROOM
		&& $blast_hit->passes_thresholds())
	{
		$blast_hit->type($InsertionSeq::ContigBlastHit::TYPE_FLANK_WHOLE);
	}
	else
	{
		$blast_hit->type($InsertionSeq::ContigBlastHit::TYPE_FLANK_PARTIAL);
	}
	
	##
	## Record which base is closest to the IS
	##
	my $base = $self->get_location_closest_to_is($blast_hit);
	
	##
	## Save in in the blast_hit
	$blast_hit->closest_is_base($base);
	
	$log->debug("Closest_IS: ".$blast_hit->type()." $self->{contig_direction} $self->{location} $blast_hit->{sdir} ($blast_hit->{sstart} $blast_hit->{send}) = ".$blast_hit->closest_is_base()."\n");
	
	push @{$self->{blast_hits}}, $blast_hit;
	
	##
	## Reset the best, since we don't know anymore
	##
	undef($self->{best_blast_hit});
	$self->{direction} = "";
}


##
## Find the blast closest to the end of the flank adjacent to the IS
## There's probably a slicker way to do this, but it seems clearer to me
## to just break out the cases 
##
sub pick_best_blast
{
	my $self = shift;
	
	if (defined($self->{best_blast_hit}))
	{
		return $self->{best_blast_hit};
	}
		
	if (defined($self->{blast_hits}))
	{	
		my @sorted_list;
		
		if ($self->{contig_direction} eq $DIRECTION_FWD)
		{
			##
			## Original contig sequence direction is FWD
			##
			if ($self->{location} eq $LOCATION_BEGIN)
			{
				## Maximum query end
				@sorted_list = sort { $b->passes_thresholds() <=> $a->passes_thresholds() || $b->{qend} <=> $a->{qend} || $b->{qhitlen} <=> $a->{qhitlen} } @{$self->{blast_hits}}
			}
			elsif ($self->{location} eq $LOCATION_END)
			{
				# Minimum query start
				@sorted_list =  sort  { $b->passes_thresholds() <=> $a->passes_thresholds() || $a->{qstart} <=> $b->{qstart} || $b->{qhitlen} <=> $a->{qhitlen}  } @{$self->{blast_hits}}
			}
			else
			{
				die "Unexpected flank location: $self->{location}\n";
			}
		}	
		elsif ($self->{contig_direction} eq $DIRECTION_REV)
		{
			##
			## Original contig sequence direction is REV
			##
			if ($self->{location} eq $LOCATION_BEGIN)
			{
				# Minimum query start
				@sorted_list =  sort  { $b->passes_thresholds() <=> $a->passes_thresholds() || $a->{qstart} <=> $b->{qstart} || $b->{qhitlen} <=> $a->{qhitlen}  } @{$self->{blast_hits}}		
			}
			elsif ($self->{location} eq $LOCATION_END)
			{
				## Maximum query end
				@sorted_list =  sort  { $b->passes_thresholds() <=> $a->passes_thresholds() ||  $b->{qend} <=> $a->{qend} || $b->{qhitlen} <=> $a->{qhitlen}  } @{$self->{blast_hits}}
			}
			else
			{
				die "Unexpected flank location: $self->{location}\n";
			}	
		}
		else
		{
			die "Unexpected contig blast direction: $self->{contig_direction}\n";
		}
		

		##
		## Pick the best from the sort
		##
		my $best = $sorted_list[0];

		if ($best)
		{
			##
			## Copy up the relevant blast data into the flank
			##
			$self->{best_blast_hit} = $best;

			if ($self->{contig_direction} eq $DIRECTION_FWD and $best->{sdir} eq $DIRECTION_FWD)
			{
				$self->{direction} = $DIRECTION_FWD;
			}
			elsif ($self->{contig_direction} eq $DIRECTION_REV and $best->{sdir} eq $DIRECTION_REV)
			{
				$self->{direction} = $DIRECTION_FWD;
			}
			else
			{
				$self->{direction} = $DIRECTION_REV;
			}
			$self->{annotation_name} = $best->{annotation}->{name};
			$self->{closest_is_base} = $best->closest_is_base();
			$self->{flank_pct_id} = $best->{pct_id};
		}

	}
	

	return $self->{best_blast_hit};
}

sub closest_is_base
{
	my $self = shift;

	die "Can't find closest IS base to a flank with no blast hits"  if (!defined($self->{best_blast_hit}));

	$self->get_location_closest_to_is($self->{best_blast_hit});
}

sub is_dup_of
{
    my $self = shift;
    my $other = shift;


    ##
    ## Must be same IS, reference, query, location
	## Assume reference is the same....
    ##
    return ($self->{annotation_name} eq $other->{annotation_name}
		&& $self->{closest_is_base} == $other->{closest_is_base}
		&& $self->{is_name} eq $other->{is_name});


}

sub is_mate_of
{
	my $self = shift;
	my $other = shift;


	##
	## Must be same reference
	##
	return 0 if ($self->{annotation_name} ne $other->{annotation_name});

	
	if ($self->{direction} eq $DIRECTION_REV && $other->{direction} eq $DIRECTION_REV)
	{
		## Begin will be 15 bases toward 3' end from End
		#return 1 if ($self->{is_location} eq $LOCATION_BEGIN && $self->{offset_from_last} <= $MAX_FLANK_PAIRING_DISTANCE);
		#return 1 if ($other->{is_location} eq $LOCATION_BEGIN && $other->{offset_from_last} <= $MAX_FLANK_PAIRING_DISTANCE);

		return 1 if ($self->{is_location} eq $LOCATION_END 
			&& $other->{is_location} eq $LOCATION_BEGIN 
			&& $other->closest_is_base() > $self->closest_is_base() 
			&& ($other->closest_is_base() - $self->closest_is_base())  <= $MAX_FLANK_PAIRING_DISTANCE);

		return 1 if ($self->{is_location} eq $LOCATION_BEGIN 
			&& $other->{is_location} eq $LOCATION_END 
			&& $self->closest_is_base()  > $other->closest_is_base()
			&& ($self->closest_is_base() - $other->closest_is_base())  <= $MAX_FLANK_PAIRING_DISTANCE);
	}
	elsif ($self->{direction} eq $DIRECTION_FWD && $other->{direction} eq $DIRECTION_FWD)
	{
		## End will be 15 bases toward 3' end from Begin
		#return 1 if ($self->{is_location} eq $LOCATION_END && $self->{offset_from_last} <= $MAX_FLANK_PAIRING_DISTANCE);
		#return 1 if ($other->{is_location} eq $LOCATION_END && $other->{offset_from_last} <= $MAX_FLANK_PAIRING_DISTANCE);

		return 1 if ($self->{is_location} eq $LOCATION_END 
			&& $other->{is_location} eq $LOCATION_BEGIN 
			&& $self->closest_is_base()  > $other->closest_is_base()
			&& ($self->closest_is_base() - $other->closest_is_base())  <= $MAX_FLANK_PAIRING_DISTANCE);

		return 1 if ($self->{is_location} eq $LOCATION_BEGIN 
			&& $other->{is_location} eq $LOCATION_END 
			&& $other->closest_is_base() > $self->closest_is_base() 
			&& ($other->closest_is_base() - $self->closest_is_base())  <= $MAX_FLANK_PAIRING_DISTANCE);
	}


	return 0;
}

sub mate_with
{
	my $self = shift;
	my $other = shift;
	
	$other->{mate} = $self;
	$self->{mate} = $other;
}


sub get_valid_contig_count
{
	my $self = shift;

    my $contig_count = 1; ## 1 for me
	for my $flank ( @{$self->{dup_list}} )
	{
		my $blast = $flank->pick_best_blast();
    	++$contig_count if ($blast->passes_thresholds());
	}

   return $contig_count;
}

sub to_is_mysql
{
	my $self = shift;
	my $end_flank_distance = shift;
	my $end_flank_id = shift;
	my $source_genome = shift;
	
	$end_flank_distance  = "NULL" unless defined($end_flank_distance);
	$end_flank_id  = "NULL" unless defined($end_flank_id);
	
	my $sql;
	my $blast = $self->pick_best_blast();
	
	if (defined($blast))
	{
	    my  ($geneBefore,$gene,$geneAfter) = $blast->{annotation}->locate_base($blast->{closest_is_base}) ;
	    
	    ##
	    ## Quote them if they exist
	    ##
	    my $geneBeforeSql = "'$geneBefore->{accession}'" if defined($geneBefore);
	    my $geneInSql = "'$gene->{accession}'" if defined($gene);
	    my $geneAfterSql = "'$geneAfter->{accession}'" if defined($geneAfter);
	    
	    ##
	    ## NULL them if they don't
	    ##
	    $geneBeforeSql = "NULL" unless defined($geneBeforeSql);
	    $geneInSql = "NULL" unless defined($geneInSql);
	    $geneAfterSql = "NULL" unless defined($geneAfterSql);
	    

		##
		## Count up supporting contigs
		##
		my $contig_count = $self->get_valid_contig_count();
	    
		##
		## Write dependent record first:
		##
		$sql = $self->to_feat_mysql($source_genome,$self->{annotation_name},$self->{annotation_name} );

		my $flank_pct_id = $self->{flank_pct_id};
		$flank_pct_id = "NULL" unless $flank_pct_id;
	    
		$sql .= "INSERT INTO is_flank (is_run_id, is_qf_id, is_element, reference, q_genome, begin_end, orientation, nearest_base, after_gene, in_gene, before_gene, flank_pct_id, match_quality, flank_separation, mate_flank_id, contig_count ) VALUES (";
    
		$sql .= '@is_run_id,';
		$sql .= "$SQL_LAST_FEAT_ID_VAR,";
		$sql .= "'$self->{is_name}',";
		$sql .= "'$self->{annotation_name}',";
		$sql .= "'$source_genome',";
		$sql .= "'$self->{is_location}',";
		$sql .= "'$self->{direction}',";
		$sql .= "$self->{closest_is_base},";
		$sql .= "$geneBeforeSql,";
		$sql .= "$geneInSql,";
		$sql .= "$geneAfterSql,";
		$sql .= "$flank_pct_id,";
		$sql .= "'$blast->{type}',";
		$sql .= "$end_flank_distance,";
		$sql .= "$end_flank_id,";
		$sql .= "$contig_count";
		$sql .= ");\n" ;	

		$sql .= "SET $SQL_LAST_IS_FLANK_ID_VAR = LAST_INSERT_ID();\n";
	}	

	$sql;
}

sub to_feat_mysql
{
	my $self = shift;
	my $query_genome = shift;
	my $source_genome = shift;
	my $annotation = shift;


	my $annotated = 0;
	my $blast = $self->pick_best_blast();
	$annotated = 1 if (defined($blast) && $blast->passes_thresholds() );

	my $flank_pct_id = $self->{flank_pct_id} if ($self->{annotation_name} && $self->{annotation_name} eq $annotation);
	$flank_pct_id = "NULL" unless $flank_pct_id;

    my $sql = "INSERT INTO is_query_feature (is_run_id,is_element, reference, q_genome, s_genome, is_pct_id, flank_pct_id, contig_name, feat_type, flank_begin_end, orientation, begin_base, end_base, is_annotated ) VALUES (";
		$sql .= '@is_run_id,';
		$sql .= "'$self->{is_name}',";
		$sql .= "'$annotation',";
		$sql .= "'$query_genome',";
		$sql .= "'$source_genome',";
		$sql .= "$self->{is_pct_id},";
		$sql .= "$flank_pct_id,";
		$sql .= "'$self->{contig_name}',";
		$sql .= "'$SQL_FLANK_FEAT_TYPE',";
		$sql .= "'$self->{is_location}',";
		$sql .= "'$self->{contig_direction}',";
		$sql .= "'$self->{contig_lower_coord}',";
		$sql .= "'$self->{contig_upper_coord}',";
		$sql .= $annotated;
		$sql .= ");\n" ;	
		$sql .= "SET $SQL_LAST_FEAT_ID_VAR = LAST_INSERT_ID();\n";

	$sql;

}

##our $CSV_HEADER = "is_element, reference, q_genome, s_genome, contig_name, feat_type, flank_begin_end, orientation, begin_base, end_base, is_annotated, (flank) is_element, reference, q_genome, begin_end, orientation, nearest_base, after_gene, in_gene, before_gene, match_quality, flank_separation, mate_flank_id, contig_count\n";
our $CSV_HEADER = "is_element, genome, contig_name, feat_type, flank_begin_end, orientation, contig_flank_begin_base, contig_flank_end_base, is_annotated, reference, orientation, nearest_base, after_gene, in_gene, before_gene, match_quality, flank_separation, mate_flank_id, contig_count\n";

sub to_csv
{

	my $self = shift;
	my $end_flank_distance = shift;
	my $end_flank_id = shift;
	my $source_genome = shift;


	my $annotated = 0;
	$annotated = 1 if (defined($self->pick_best_blast()));


	my $blast = $self->pick_best_blast();
	
	my $csv = "";
    my $geneBeforeCsv = "";
    my $geneInCsv = "";
    my $geneAfterCsv = "";

	my $blastType = "";
		my $contig_count = 1; ## 1 for me
		$contig_count  += scalar(@{$self->{dup_list}}) if ($self->{dup_list});

   	my $annotation_name = $self->{annotation_name} if ($self->{annotation_name});

	if (defined($blast))
	{
    	$blastType = $blast->{type};
		if ($blast->{annotation})
		{
	    	my  ($geneBefore,$gene,$geneAfter) = $blast->{annotation}->locate_base($blast->{closest_is_base});
    		$geneBeforeCsv = $geneBefore->{accession} if defined($geneBefore);
    		$geneInCsv = $gene->{accession} if defined($gene);
    		$geneAfterCsv = $geneAfter->{accession} if defined($geneAfter);
		}
	}


		$csv .= "$self->{is_name},";
        $csv .= "$source_genome,";
        $csv .= "$self->{contig_name},";
        $csv .= "$SQL_FLANK_FEAT_TYPE,";
        $csv .= "$self->{is_location},";
        $csv .= "$self->{contig_direction},";
        $csv .= "$self->{contig_lower_coord},";
        $csv .= "$self->{contig_upper_coord},";
        $csv .= "$annotated,";

		## flank info 
		#$csv .= "$self->{is_name},";
		$csv .= "$annotation_name,";
		#$csv .= "$source_genome,";
		#$csv .= "$self->{is_location},";
		$csv .= "$self->{direction},";
		$csv .= "$self->{closest_is_base},";
		$csv .= "$geneBeforeCsv,",
		$csv .= "$geneInCsv,";
		$csv .= "$geneAfterCsv,";
		$csv .= "$blastType,";
		$csv .= "$end_flank_distance,";
		$csv .= "$end_flank_id,";
		$csv .= "$contig_count";
		$csv .= "\n";
	

	$csv;

}


1;
