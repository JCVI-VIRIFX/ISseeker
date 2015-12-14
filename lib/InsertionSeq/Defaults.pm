package InsertionSeq::Defaults;

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


use Config::General;
use strict;
use warnings;

BEGIN { 
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( 
	$BLASTALL_PATH
	$FORMATDB_PATH
	$EXTRACTSEQ_PATH
    $DEFAULT_REQ_IS_PERCENT_ID
	$DEFAULT_REQ_FLANK_PERCENT_ID
);
require Exporter;
}



our %config;

our $SIZEFASTA_PATH;
our $BLASTALL_PATH;
our $FORMATDB_PATH;
our $EXTRACTSEQ_PATH;
our $DEFAULT_REQ_IS_PERCENT_ID;
our $DEFAULT_REQ_FLANK_PERCENT_ID;


sub Process
{
	my $file = shift @_;

	%config = new Config::General(
        -ConfigFile => $file,
        )
        ->getall();

	

	##
	## Required
	##
	$BLASTALL_PATH           		= $config{"BLASTALL_PATH"} || die "BLASTALL_PATH not set in config file $file!\n";
	$FORMATDB_PATH           		= $config{"FORMATDB_PATH"} || die "FORMATDB_PATH not set in config file $file!\n";
	$EXTRACTSEQ_PATH         		= $config{"EXTRACTSEQ_PATH"} || die "EXTRACTSEQ_PATH not set in config file $file!\n";
	$DEFAULT_REQ_IS_PERCENT_ID		= $config{"DEFAULT_REQ_IS_PERCENT_ID"} || die "DEFAULT_REQ_IS_PERCENT_ID not set in config file $file!\n";
	$DEFAULT_REQ_FLANK_PERCENT_ID	= $config{"DEFAULT_REQ_FLANK_PERCENT_ID"} || die "DEFAULT_REQ_FLANK_PERCENT_ID not set in config file $file!\n";


	\%config;
}

1;
