package InsertionSeq::AnnotationGene;

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

my $PLUS_STRAND		= "+";
my $MINUS_STRAND	= "-";

##
## Figure out direction and make numbers go low -> high
## Add notation for +/-
##
sub new 
{
    my $class = shift;
    my $self = {@_};

    bless($self, $class);

	die ("AnnotationGene missing name.\n") if ( !defined($self->{name}) );
	die ("AnnotationGene $self->{name} missing accession name.\n") if ( !defined($self->{accession}) );
	die ("AnnotationGene $self->{name} missing start.\n") if ( !defined($self->{start}) );

	die ("AnnotationGene $self->{name} missing end.\n") if ( !defined($self->{end}) );
	#die ("AnnotationGene missing strand.\n") if ( !defined($self->{strand}) ;

	if ( $self->{start} > $self->{end} ) 
	{
		$self->{strand} = $MINUS_STRAND;

		## Swap
		my $tmp = $self->{end};
		$self->{end}   = $self->{start};
		$self->{start} = $tmp;
	}
	else 
	{
		$self->{strand} = $PLUS_STRAND;
	}
	
	die ("Bad gene annotation") if ( !defined($self->{start})  ||  !defined($self->{end}) ) ;

	return $self;
}

sub start {
	my $self = shift;
	if (@_) {
		$self->{start} = shift;
	}
	return $self->{start};
}

sub end {
	my $self = shift;
	if (@_) {
		$self->{end} = shift;
	}
	return $self->{end};
}

sub strand {
	my $self = shift;
	if (@_) {
		$self->{strand} = shift;
	}
	return $self->{strand};
}

sub contains_base
{
	my $self = shift;
	my $base = shift;
	
	return $base >= $self->{start}	&& $base <= $self->{end};
}

sub is_before_base
{
	my $self = shift;
	my $base = shift;
	
	return $base > $self->{end};
}

sub is_after_base
{
	my $self = shift;
	my $base = shift;
	
	return $base < $self->{start};
}
