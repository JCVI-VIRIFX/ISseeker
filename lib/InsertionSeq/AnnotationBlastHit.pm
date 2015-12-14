package InsertionSeq::AnnotationBlastHit;

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

use parent q(InsertionSeq::BlastHit);

use strict;

use Log::Log4perl;


my $log = Log::Log4perl->get_logger();

sub new 
{
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    die ("AnnotationBlastHit missing annotation.\n") if ( !defined($self->{annotation}) );
    die ("AnnotationBlastHit missing IS name.\n") if ( !defined($self->{is_name}) );
    
    return $self;
}


##
## Record which subject base is closest to the IS
##
sub closest_is_base {
	my $self = shift;
	if (@_) {
		$self->{closest_is_base} = shift;
	}
	return $self->{closest_is_base};
}




1;
