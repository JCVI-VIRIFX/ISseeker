package InsertionSeq;

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

use Log::Log4perl qw( :easy );



sub extract_gb_name
{
    my $base_sid = $self->{contig_name};

$log->info("Checking contig name: $self->{contig_name}\n");
    ##
    ## Remove pipes
    ##
    if ($self->{contig_name} =~ /gb\|(\w+|\.)+/)
    {
        ##
        ## opt for accession first
        ##
        $base_sid = $1;
$log->info("Extracting short name $base_sid\n");
    }
    elsif ($self->{contig_name} =~ /gi\|(\w+|\.)/)
    {
        ##
        ## Or take GI
        ##
        $base_sid = $1;
    }
    elsif ($self->{contig_name} =~ /^(\w+|\.)\|/)
    {
        ## Take what's first
        $base_sid = $1;
    }
    else
    {
        $log->logdie("Can't un-pipe this genome def line: $self->{contig_name}\n");
    }

}


1;
