package InsertionSeq::BlastHit;

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

use Log::Log4perl;

##
## Make getters and setters for our members
##
use base qw(Class::Accessor);
  InsertionSeq::BlastHit->follow_best_practice;
  InsertionSeq::BlastHit->mk_accessors(qw(
    q_name      
    is_name    
    pct_id     
    req_pct_id 
    queryid    
    qstart     
    qend       
    qlength    
    subjectid  
    sstart     
    send       
    slength    
    name       
    filename   
));



##
## What direction is the match with respect to the target seq
##
our $DIRECTION_FWD		= "F";
our $DIRECTION_REV		= "R";


my $log = Log::Log4perl->get_logger();

sub new 
{

    my $class = shift;
    my $self = {@_};

    bless($self, $class);

	die ("BlastHit missing s start.\n") if ( !defined($self->{sstart}) );
	die ("BlastHit missing s end.\n") if ( !defined($self->{send}) );
	die ("BlastHit missing s length.\n") if ( !defined($self->{slength}) );
	die ("BlastHit missing q start.\n") if ( !defined($self->{qstart}) );
	die ("BlastHit missing q end.\n") if ( !defined($self->{qend}) );
	die ("BlastHit missing q length.\n") if ( !defined($self->{qlength}) );
	die ("BlastHit missing target name.\n") if ( !defined($self->{name}) );
	die ("BlastHit missing query name.\n") if ( !defined($self->{q_name}) );
	die ("BlastHit missing file name.\n") if ( !defined($self->{filename}) );

	$self->initialize();

    return $self;
}

##
## Figure out direction and make Subject match numbers low -> high
##
sub initialize 
{
	my $self = shift;
	if ( $self->{sstart} > $self->{send} ) 
	{
		$self->{sdir} = $DIRECTION_REV;

		## Swap
		my $tmp = $self->{send};
		$self->{send}   = $self->{sstart};
		$self->{sstart} = $tmp;
	}
	else 
	{
		$self->{sdir} = $DIRECTION_FWD;
	}
	
	die ("Unexpected qstart > qend.") if ( $self->{qstart} > $self->{qend} ) ;

	$self->{shitlength} =  ($self->{send} - $self->{sstart} ) + 1;
	$self->{qhitlength} =  ($self->{qend} - $self->{qstart}) + 1;

}


sub type {
	my $self = shift;
	if (@_) {
		$self->{type} = shift;
	}
	return $self->{type};
}


sub passes_pct_threshold
{
	my $self = shift;
	return $self->get_pct_id() >= $self->get_req_pct_id();
}

sub passes_pct_threshold_sort
{
	my $self = shift;
	if ($self->get_pct_id() >= $self->get_req_pct_id())
	{
		return 1;
	}
	else
	{
    	return 0;
	}
}
1;
