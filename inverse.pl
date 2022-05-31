#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my $char_file = "data/character.json";
my $frame_time = 1/120;
my $cutoff_freq = 2;
my $ground_offset = 0;
my $outdir = "output";
my $pose_file;

GetOptions(
    "char_file=s" => \$char_file,
    "frame_time=i" => \$frame_time,
    "cutoff_freq=i" => \$cutoff_freq,
    "ground_offset=i" => \$ground_offset,
    "outdir=s"=> \$outdir,
);

if (@ARGV != 1) {
    die <<"USAGE";
usage: $0 [options] pose_file

options:
    --char_file=string
    --frame_time=double
    --cutoff_freq=double
    --ground_offset=double
    --outdir=string
USAGE
} else {
    $pose_file = shift @ARGV;
}

mkdir $outdir unless -e $outdir;

system "python find-contact-node.py '$pose_file' '$outdir'";
system "./inverse.out --char_file='$char_file' --frame_time=$frame_time --cutoff_freq=$cutoff_freq --ground_offset=$ground_offset --outdir='$outdir' '$pose_file' '$outdir/contact_nodes.txt'";
