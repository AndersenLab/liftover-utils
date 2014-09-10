#!/software/bin/perl -w 
#
#   remap_gff_between_releases.pl                 
# 
# by Gary Williams                         
#
# This remaps the clone positions and genomic sequences between two
# releases i.e. it remaps positions going forwards in time, 
# e.g. from release 140 to release 160
#
# Last updated by: $Author: gw3 $     


use strict;                                      
use Getopt::Long;
use Carp;
use warnings;


######################################
# variables and command-line options # 
######################################


my ($release1, $release2, $version, $gff, $output);

GetOptions (
	    "gff=s"      => \$gff,
	    "output=s"   => \$output,
	    "release1=i"  => \$release1,
	    "release2=i"  => \$release2,
	    );



if (! defined $release1 || ! defined $release2) {
  die "Specify the release numbers to use\n";
}

if (! defined $gff || ! defined $output) {
  die "Specify the input and output files\n";
}

##########################
# read in the mapping data
##########################

my %mapping_data = read_mapping_data($release1, $release2);


##########################
# MAIN BODY OF SCRIPT
##########################

my ($indel, $change);

open (OUT, "> $output") || die "Can't open $output";
open (GFF, "< $gff") || die "Can't open GFF file $gff\n";

while (my $line = <GFF>) {
  chomp $line;
  if ($line =~ /^\s*$/) {next;}
  if ($line =~ /^\#/) {print OUT $line,"\n"; next;}
  if ($line =~ /^>/) {print OUT $line,"\n"; while (<GFF>) {print OUT $_;} next;}
  my @f = split /\t/, $line;

  my ($chromosome, $start, $end, $sense) = ($f[0], $f[3], $f[4], $f[6]);

  # change III, chrIII, chr_III to the standard WormBase nomenclature CHROMOSOME_III
  $chromosome =~ s/^chr_//;
  $chromosome =~ s/^chr//;
  $chromosome =~ s/^CHROMOSOME_//;
  $chromosome = "CHROMOSOME_${chromosome}";

  # some checks for malformed GFF files
  if (!defined $chromosome || ! defined $start || ! defined $end || ! defined $sense) {die "Malformed line (Missing values? Spaces instead of TABs?): $line\n"}
  if ($chromosome !~ /(^CHROMOSOME_[IVX]+$)|(^CHROMOSOME_MtDNA$)/) {die "Malformed line (Invalid chromosome name? Chromosomes should look like: CHROMOSOME_IV): $line\n"}
  if ($start !~ /^\d+$/) {die "Malformed line (Non-numeric start?): $line\n"} 
  if ($end !~ /^\d+$/) {die "Malformed line (Non-numeric end?): $line\n"} 
  if ($sense !~ /^[\+|\-]$/) {die "Malformed line (Invalid sense?): $line\n"}

  #print "chrom, start, end=$chromosome, $start, $end\n";
  ($f[3], $f[4], $f[6], $indel, $change) = remap_gff($chromosome, $start, $end, $sense, $release1, $release2, %mapping_data);
  
  if ($indel) {
    print "$line\nThere is an indel in the sequence in CHROMOSOME $chromosome, $start, $end\n";
  } elsif ($change) {
    print "$line\nThere is a change in the sequence in CHROMOSOME $chromosome, $start, $end\n";
  }
  
  $line = join "\t", @f;
  print OUT $line,"\n";
}

close (GFF);
close (OUT);


exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################################
# 
# Name:      read_mapping_data
# Usage:     %mapping_data = &read_mapping_data($release1, $release2);
# Function:  reads the data used to remap across release
# Args:      $release1, $release2, the first and last wormbase release
#                  numbers to use e.g. 140, 150 to convert data made using wormbase
#                  release WS140 to the coordinates of release WS150
# Returns:   the mapping data, for use in remap_gff()
#

# the mismatch_start value is the start of the mismatch, it is the first position which doesn't match
# the mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
                                                                                                                                                            
#
# Input files look like:
# Chromosome: CHROMOSOME_I
# 4765780 4765794 14      4765780 4765794 14      0
                                                                                                                                                            
                                                                                                                                                            
sub read_mapping_data {

  my ($release1, $release2) = @_;
                                                                                                                                           
  # hash (keyed by release) of hashes (keyed by chromosome number) of list (one for each difference) of list (one for each field)
  # access like: $fields_hashref = @{ $mapping_data{$release}{$chrom}[$next_difference] }
  my %mapping_data;
                                                                                                                                                            
  foreach my $release (($release1+1) .. $release2) {
    my %chroms;
    my $infile = "CHROMOSOME_DIFFERENCES/sequence_differences.WS$release";
    open (IN, "< $infile") || die "Can't open $infile\n";
    my $chrom;
    while (my $line = <IN>) {
      chomp $line;
      if ($line =~ /Chromosome:\s+(\S+)/) {
        $chrom = $1;
      } else {
        my @fields = split /\t/, $line;
        #print "fields=@fields\n";
	push @{$mapping_data{$release}->{$chrom}}, [@fields];                                  
      }
    }
    close(IN);
  }
                                                                                                                                                            
  return %mapping_data;
}



##########################################################
# 
# Name:      remap_gff
# Usage:     ($new_start, $new_end, $new_sense, $indel, $change) = remap_gff($chromosome, $start, $end, $sense, $release1, $release2, %mapping_data);
# Function:  does the unmapping (remapping backwards to past versions) of a pair of location values for a GFF file
# Args:      $chromosome, the chromosome number, e.g. 'III'
#            $start, the start value of the chromosomal location coordinate
#            $end, the end value of the chromosomal location coordinate (always >= $start)
#            $sense, the sense "+" or "-" of the coordinate
#            $release2, $release1, the first and last wormbase release
#                  numbers to use e.g. 140, 150 to convert data made using wormbase
#                  release WS140 to the coordinates of release WS150
#            @mapping_data - data as returned by read_mapping_data
# Returns:   $new_start, $new_end, $new_sense - the updated chromosomal location coordinates
#            $indel - true if indels affect this location
#            $change - true if any sort of changes affect this location

sub remap_gff {
  my ($chromosome, $start, $end, $sense, $release1, $release2, %mapping_data) = @_;

  my $indel = 0;                # true if indels affect this location
  my $change = 0;               # true if non-indel base changes affect this location
  my $start_deleted = 0;
  my $end_deleted = 0;

  if ($chromosome eq 'CHROMOSOME_MtDNA') {return ($start, $end, $sense, 0, 0, 0, 0);} # the mitochondrion is invariant
  if ($chromosome !~ /^(CHROMOSOME_I|CHROMOSOME_II|CHROMOSOME_III|CHROMOSOME_IV|CHROMOSOME_V|CHROMOSOME_X)$/) {
    die "Unrecognised chromosome ID: '$chromosome'\nThe chromosome ID must be one of\nCHROMOSOME_I, CHROMOSOME_II, CHROMOSOME_III, CHROMOSOME_IV, CHROMOSOME_V, CHROMOSOME_X\n";
  }

  foreach my $release (($release1+1) .. $release2) {
    my $this_ver_start = $start;
    my $this_ver_end = $end;
    
    if (exists $mapping_data{$release}->{$chromosome}) {
      foreach  my $fields (@{$mapping_data{$release}->{$chromosome}}) {
	
        #print "$release $chromosome fields= @$fields \n";
	
	# The mismatch_start value is the start of the mismatch, it is the first position which doesn't match.
	# The mismatch_end value is the base past the end of the mismatch region, the first base which matches again
	# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
        my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$fields;
	
	# N.B. mismatch values are in the normal perl coordinate system starting at position 0.
	# Convert them to the GFF coordinates system starting at position 1.
        $mismatch_start1++;
        $mismatch_end1++;
	
        if ($flipped) {    # is the feature inside a flipped region?
	  
	  
	  if ($this_ver_start >= $mismatch_start1 && $this_ver_end < $mismatch_end1) {
	    # flip strand
	    $sense = ($sense eq '+') ? '-' : '+';
	    
	    # $start = $mismatch_end - ($start - $mismatch_start), which simplifies to:
	    $start = $mismatch_start1 + $mismatch_end1 - $start; # flip the start and end positions
	    
	    # $end = $mismatch_end - ($end - $mismatch_start), which simplified to:
	    $end = $mismatch_start1 + $mismatch_end1 - $end;
	    
	    if ($start > $end) {
	      ($start, $end) = ($end, $start);
	    }
	    
	    # does the edge of the flipped region overlap our location?
	  } elsif ($this_ver_start >= $mismatch_start1 && $this_ver_start < $mismatch_end1 ||
		   $this_ver_end >= $mismatch_start1 && $this_ver_end < $mismatch_end1) {
	    # don't change the location, but note that we have changes
	    $change = 1;
	  }
 	} else {
	  
	  # if there is a change inside our location, note it
	  if ($this_ver_start < $mismatch_end1 and
	      (($len1 > 0 and $this_ver_end >= $mismatch_start1) or
	       ($len1 == 0 and $this_ver_end > $mismatch_start1))) {
	    $change = 1;
	  }
	  
	  # note the length of our location so we can see any indels occurring
	  my $location_length = $this_ver_end - $this_ver_start;
	  
	  # if the start or end are beyond the start of the change region, apply any shift
	  if ($len1 > 0 and $this_ver_start >= $mismatch_end1 or
	      $len1 == 0 and $this_ver_start > $mismatch_end1) { # if past the end of the change region, shift it
	    $start += $len2 - $len1;
	  } elsif ($this_ver_start >= $mismatch_start1) {
	    # within a changed segment; if the position within the original segment
	    # maps to beyond the end of the replacement segment, then this position no
	    # longer exists, in wich case we barf
	    if ($this_ver_start - $mismatch_start1 + 1 > $len2 ) {
	      $start_deleted = 1;
	      $start = $mismatch_start2 + $len2;
	    }
	  }
	  
	  if ($len1 > 0 and $this_ver_end >= $mismatch_end1 or
	      $len1 == 0 and $this_ver_end > $mismatch_end1) {
	    # if past the end of the change region, shift it
	    $end += $len2 - $len1;
	  } elsif ($this_ver_end >= $mismatch_start1) {
	    # within a changed segment; if the position within the original segment
	    # maps to beyond the end of the replacement segment, then this position no
	    # longer exists, in wich case we barf
	    if ($this_ver_end - $mismatch_start1 + 1 > $len2) {
	      $end_deleted = 1;
	      $end = $mismatch_start2 + $len2;
	    }
	  }
	  
	  # see if we have any indels in our location
	  if ($location_length != $end - $start) {
	    $indel = 1;
	  }
	  
 	}
      }
    } else {
      #print "no change: doesn't exist: $release $chromosome\n";
    }
  }
  
  # print " Remapped to $start $end $sense $indel $change\n";
  
  return ($start, $end, $sense, $indel, $change, $start_deleted, $end_deleted);
}

                           


##########################################

# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - remap_gff_between_releases.pl

=head1 USAGE

=over 4

=item remap_gff_between_releases.pl [options]

=back

This script reads in a GFF file and the numbers of two releases of
wormbase and maps the chromosomal locations of the first release to
the second release.


script_template.pl MANDATORY arguments:

=over 4

=item -release1 The first (earlier) database to convert from e.g. 140

=back

=item -release2 The second (later) database to convert to e.g. 155

=back

=item -gff the name of the GFF file to read in and convert

=back

=item -outfile the name of the converted GFF file to write out

=back



=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
