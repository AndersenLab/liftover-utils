#!/usr/bin/perl -w 
use strict;

# Lift over coordinates from GFF, WIG/BED or WIG/variable step from one 
# WS build to a newer WS version.
# This is a wrapper for Gary Williams' liftover script, which works.

# IMPORTANT! edit this constant to specify where you have installed this script
# (in the remap_gff directory).  So, for example, if you specify the location
# /home/username, then it is assumed the script and remapping data reside in
# /home/username/remap_gff
use constant LOCATION => "/home/smckay/tracks";

# the ftp URL to grab updated genome assembly data, if available
use constant FTP      => "ftp://ftp.sanger.ac.uk/pub2/wormbase/gary/remapping-data/";


# Do a sanity check here: Do we have the correct directory?
unless (-d LOCATION . "/remap_gff") {
  my $path = LOCATION . "/remap_gff";
  die <<END;
I count not find the path '$path'
Did you remember to edit the LOCATION constant?

# IMPORTANT! edit this constant to specify where you have installed this script
# (in the remap_gff directory).  So, for example, if you specify the location
# $ENV{HOME}, then it is assumed the script and remapping data reside in
# $ENV{HOME}/remap_gff
use constant LOCATION => "/home/smckay/tracks";
END

}


my ($ref,$idx,$data,$last_ref,$header);

# Call with three args: filename WSXXX WSXXX
# Arg1 must be > Arg2
if (@ARGV < 3) {
  die "I was expecting at least three arguments:\n  Usage: WS_liftover.pl WSXXX WSXXX filename(s)\n";
}

(my $from = shift) =~ s/WS//i;
(my $to   = shift) =~ s/WS//i;

my $pwd = `pwd`;
chomp $pwd;

# make sure we have the necessary map files 
check_for_data($from,$to);

while (my $file = shift) {

  # is file compressed?
  my $gz;
  if ($file =~ /\.gz$/) {
    $gz++;
    print STDERR "unpacking $file...";
    system "gunzip -f $file";
    print STDERR " Done!\n\n";
    $file =~ s/\.gz//;
  }


  # disposition of output file
  my $outfile = $file;
  if ($outfile =~ /WS\d+/) {
    $outfile =~ s/WS\d+/WS$to/;
  }
  elsif ($outfile =~ /(\.[^\.]+)$/) {
    my $ext = $1;
    $outfile =~ s/$ext$/_WS$to$ext/;
  }
  else {
    $outfile .= "WS$to";
  }

  if (-e $outfile) {
    next;
  }
  
  print STDERR "Converting $file from WS$from to WS$to.  The outputfile is named $outfile\n";
  
  open IN, $file or die $!;
  open RAW, ">$pwd/tmpin$$";
  
  my $type = get_type($file);
  
  
  # first pass, convert input to GFF
  while (<IN>) {
    if ( /^#|wiggle|track/ ) {
      print RAW "#COMMENT$_";
      next;
    }
    if (/chrom=(\S+)/) {
      my $old_ref = $ref = $1;
      $ref =~ s/[^IVXM]//g;
      $ref = "CHROMOSOME_$ref";
      s/chrom=$old_ref/chrom=$ref/;
      print RAW "#COMMENT$_";
      next;
    }
    
    my $line;
    chomp;
    
    
    if ($type eq 'gff' || $type eq 'bed') {
      $ref = my $old_ref = (split)[0];
      if ($ref !~ /^CHROMOSOME_/) {
	$ref =~ s/[^IVXM]//g;
	$ref = "CHROMOSOME_$ref";
	s/$old_ref/$ref/;
      }    
    }
    if ($type eq 'gff') {
      $line = "$_\n";
    }
    if ($type eq 'bed') {
      chomp;
      my ($ref,$start,$end,$score) = split;
      $line = join("\t",$ref,qw/bla bla/,$start,$end,$score,qw/. ./,'bla bla') . "\n";
    }
    if ($type eq 'wig') {
      chomp;
      my ($start,$score) = split;
      $line = join("\t",$ref,qw/bla bla/,$start,$start+50,$score,qw/. ./,'bla bla') . "\n";
    }
    
    print RAW $line;
  }
  

  convert_gff("$pwd/tmpin$$",$outfile,$type);

  cleanup();


  if ($gz) {
    print STDERR "Recompressing $file...";
    system "gzip -f $file";
    print STDERR " Done!\n\n";
    print STDERR "compressing $outfile...";
    system "gzip -f $outfile";
    print STDERR " Done!\n\n";
  }

  print STDERR "Done with $file!\n\n";

}

sub convert_gff {
  my ($infile,$outfile,$type) = @_;

  my $path = LOCATION;
  my $remap_path = "$path/remap_gff";

  system "$remap_path/remap_gff_between_releases.pl --location=$path --gff=$infile --out=$pwd/tmpout$$ --release1=$from --release2=$to";

  open TMP, "$pwd/tmpout$$" or die $!;
  open OUT, ">$outfile"     or die $!;  

  while (<TMP>) {
    if (/#COMMENT/) {
      s/^#COMMENT//;
      print OUT $_;
      next;
    }
    s/CHROMOSOME_//;
    
    if ($type eq 'gff') {
      print OUT $_;
    }
    elsif ($type eq 'bed') {
      my ($ref,$start,$end,$score) = (split)[0,3,4,5];
      print OUT join("\t",$ref,$start,$end,$score), "\n";
    }
    elsif ($type eq 'wig') {
      my ($start,$score) = (split)[3,5];
      print OUT join("\t",$start,$score), "\n";
    }
  }
}

sub get_type {
  my $f = shift;
  open TMP, $f or die $!;
  my $type;
  while (<TMP>) {
    next if /wiggle|chrom=|^\#/;
    my @cols = split;
    $type = 'wig' if @cols == 2 && $cols[1] =~ /^[-+.Ee0-9]+$/;
    $type = 'bed' if @cols == 4 && $cols[1] =~ /^[-+0-9]+$/ && $cols[2] =~ /^[-+0-9]+$/;
    $type = 'gff' if @cols >= 7 && $cols[3] =~ /^[-+0-9]+$/ && $cols[4] =~ /^[-+0-9]+$/;
    last if $type;
  }
  close TMP;
  return $type;
}

sub cleanup {
  system "rm -f $pwd/tmp*$$";
}

sub check_for_data {
  my ($from,$to) = @_;
  my $path = LOCATION . "/remap_gff/CHROMOSOME_DIFFERENCES";
  chdir $path or die $!;
  my $not_found;
  # we need all WS versions between the start and the end
  for ($from .. $to) {
    unless (-e "sequence_differences.WS$_") {
      print "Uh oh! The data for WS$_ are missing, I will try to download them... ";
      system "wget ftp://ftp.sanger.ac.uk/pub2/wormbase/gary/remapping-data/sequence_differences.WS$_ 2>/dev/null >/dev/null";
      if (-e "sequence_differences.WS$_") {
	print "Success!\n";
      }
      else {
	die "Failure!\nSorry, I could not find the remapping data for WS$_, I have to give up now\n";
      }
    } 
    
  }
  chdir $pwd or die $!;
}
