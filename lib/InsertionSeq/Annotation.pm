package InsertionSeq::Annotation;

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

use strict;

#for CSV parsing
use Text::CSV;

use InsertionSeq::AnnotationGene;



our $TYPE_PLASMID		= "plasmid";
our $TYPE_CHROMOSOME   	= "chromosome";


my $log = Log::Log4perl->get_logger();


sub new 
{
    my $class = shift;
    my $self = {@_};

    bless($self, $class);

	die ("Missing annotation name.\n") if ( !defined($self->{name}) );
	die ("Missing annotation file name.\n") if ( !defined($self->{filename}) );	
	die ("Missing fasta file name.\n") if ( !defined($self->{fasta_filename}) );
	die ("Missing annotation type.\n") if ( !defined($self->{type}) );
	die ("Invalid annotation type.\n") if ( $self->{type} ne $TYPE_PLASMID && $self->{type} ne $TYPE_CHROMOSOME);

	$self->{filetype} = "csv" if ($self->{filename} =~ /.csv$/);
	$self->{filetype} = "ptt" if ($self->{filename} =~ /.ptt$/);
	die ("Unknown file type: $self->{filename}\n") if ( !defined($self->{filetype}) );



	$self->read_file();
	$log->debug( "GENES = \n");

	for my $gene (@{$self->{genelist}})
	{
    	$log->debug(  $gene->{start}." - ".$gene->{end}." [$gene->{strand}] ".$gene->{name}."\n" );
	}

	return $self;
}


sub read_csv_file
{
	my $self = shift;

	#51775,52908,10,"+","ACICU_00044","Tellurite resistance protein"

	my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();

 
    my %geneHash;
    my $geneParseCount = 0;
 
	open my $fh, "<:encoding(utf8)", $self->{filename} or die("Could not open file $self->{filename}.");
	while ( my $row = $csv->getline( $fh ) ) 
	{
		my $accession = $row->[4];
		$accession =~ s/^\s*//;
		$accession =~ s/\s*$//;
		if (length($accession) == 0)
		{
			$log->trace("Skipping gene with no accesssion $row->[5] start=$row->[0] end=$row->[1]\n");
		}
		else
		{
			$log->trace("Adding gene $row->[4] $row->[5] start=$row->[0] end=$row->[1]\n");
			my $gene = InsertionSeq::AnnotationGene->new(
				start		=> $row->[0],
				end 		=> $row->[1],
				name		=> $accession." ".$row->[5],
				accession   => $accession
			);
			$geneParseCount++;
	
			if (defined($geneHash{$gene->{accession}}))
			{
				my $oldGene = $geneHash{$gene->{accession}};
				die "Conflicting annotation information in file $self->{filename} for gene accession $gene->{accession} ($gene->{name})" unless ($oldGene->{start} == $gene->{start} && $oldGene->{end} == $gene->{end})
			}
			
			$geneHash{$gene->{accession}} = $gene;
		}

	}
	$csv->eof or $csv->error_diag();
	close $fh;

	$log->debug( "NAME:$self->{name} FILE:$self->{filename} GENES PARSED = $geneParseCount UNIQUE = ".keys( %geneHash )."\n" );
	

	@{$self->{genelist}} =  sort { $a->{start} <=> $b->{start} } values %geneHash;
}

sub read_ptt_file
{
	my $self = shift;

##
##Acinetobacter baumannii ACICU plasmid pACICU1, complete sequence - 1..28279
##28 proteins
##Location    Strand  Length  PID Gene    Synonym Code    COG Product
##375..1325   +   316 184159989   repAc1  ACICU_p0001 -   COG5527L    replicase
##1318..1893  +   191 184159990   -   ACICU_p0002 -   -   DNA replication protein

	my $csv = Text::CSV->new ( {	binary	=> 1,  
  								 	sep_char => '	'} )
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();

    my %geneHash;
    my $geneParseCount = 0;
 
	open my $fh, "<:encoding(utf8)", $self->{filename} or die("Could not open file $self->{filename}.");
	my $header = $csv->getline( $fh ) ;
	$header = $csv->getline( $fh ) ;
	$header = $csv->getline( $fh ) ;

	while ( my $row = $csv->getline( $fh ) ) 
	{
		
		my $location = $row->[0];
		my $accession = $row->[5];
		$accession =~ s/^\s*//;
		$accession =~ s/\s*$//;
		if (length($accession) == 0)
		{
			$log->trace("Skipping gene with no accesssion $row->[5] start=$row->[0] end=$row->[1]\n");
		}
		else
		{
			my ($start,$end);
		    if ($location =~ /(\d+)\.\.(\d+)/)
			 {
					$start 	= $1;
					$end     = $2;
			}
			else
			{
             	"Cannot parse start/end from .ptt entry: $location\n";
			}
			$log->trace("Adding gene $row->[4] $row->[5] start=$row->[0] end=$row->[1]\n");
			my $gene = InsertionSeq::AnnotationGene->new(
					start		=> $start,
					end 		=> $end,
					name		=> $accession." ".$row->[8],
					accession   => $accession
			);
			$geneParseCount++;
	
			if (defined($geneHash{$gene->{accession}}))
			{
				my $oldGene = $geneHash{$gene->{accession}};
				die "Conflicting annotation information in file $self->{filename} for gene accession $gene->{accession} ($gene->{name})" unless ($oldGene->{start} == $gene->{start} && $oldGene->{end} == $gene->{end})
			}
			
			$geneHash{$gene->{accession}} = $gene;
		}

	}
	$csv->eof or $csv->error_diag();
	close $fh;

	$log->debug( "NAME:$self->{name} FILE:$self->{filename} GENES PARSED = $geneParseCount UNIQUE = ".keys( %geneHash )."\n" );
	

	@{$self->{genelist}} =  sort { $a->{start} <=> $b->{start} } values %geneHash;

}

 

sub read_file
{
	my $self = shift;
	$self->read_csv_file() if ($self->{filetype} eq "csv");
	$self->read_ptt_file() if ($self->{filetype} eq "ptt");
}

sub locate_base_text
{
	my $self = shift;
	my $base = shift;
	
	my ($geneBefore,$geneIn, $geneAfter) = $self->locate_base($base);
	

	return "BEFORE $geneAfter->{accession} ($geneAfter->{start}-$geneAfter->{end})" if (defined($geneAfter)  && !defined($geneIn) && !defined($geneBefore));
	return "INSIDE $geneIn->{accession} ($geneIn->{start}-$geneIn->{end})" if (defined($geneIn));
	return "BETWEEN $geneBefore->{accession} ($geneBefore->{start}-$geneBefore->{end}) AND $geneAfter->{accession} ($geneAfter->{start}-$geneAfter->{end})" if (defined($geneBefore) && defined($geneAfter));
	return "AFTER $geneBefore->{accession} ($geneBefore->{start}-$geneBefore->{end})" if (defined($geneBefore) && !defined($geneIn) && !defined($geneAfter));
	
	return "NOT FOUND";
}

sub locate_base
{
	my $self = shift;
	my $base = shift;
	
	my $last_gene;
	my $gene;
	

	## Dummy placeholders:
	my $geneBefore;
	my $geneIn;
	my $geneAfter;
	
	for $gene (@{$self->{genelist}})
	{
		return ($geneBefore,$geneIn, $gene) if (!defined($last_gene)  && $gene->is_after_base($base));
		return ($geneBefore,$gene,$geneAfter) if ($gene->contains_base($base));
		return ($last_gene,$geneIn,$gene) if (defined($last_gene) && $last_gene->is_before_base($base) and $gene->is_after_base($base));
		
		$last_gene = $gene;		
	}	
	
	return($last_gene, $geneIn, $geneAfter) if ($last_gene->is_before_base($base));
	

}


1;
