#!/usr/bin/env perl
use Getopt::Long;
use strict;
use warnings;

my $char_file = "data/character.json";
my $frame_time = 1/120;
my $cutoff_freq = 2;
my $ground_offset = 0;
my $reg = 10;
my $outdir = "output";
my $pose_file;

GetOptions(
    "j|char_file=s" => \$char_file,
    "f|frame_time=i" => \$frame_time,
    "c|cutoff_freq=i" => \$cutoff_freq,
    "g|ground_offset=i" => \$ground_offset,
    "r|regularization=i" => \$reg,
    "o|outdir=s"=> \$outdir,
);

if (@ARGV != 1) {
    die <<"USAGE";
usage: $0 [options] pose_file

options:
    -j, --char_file=string
    -f, --frame_time=double
    -c, --cutoff_freq=double
    -g, --ground_offset=double
    -r, --regularization=double
    -o, --outdir=string
USAGE
} else {
    $pose_file = shift @ARGV;
}

mkdir $outdir unless -e $outdir;

my $cmd;
$cmd = "python find-contact-node.py '$pose_file' '$outdir'";
print "$cmd\n";
system $cmd;
$cmd = "./inverse.out --char_file='$char_file' --frame_time=$frame_time --cutoff_freq=$cutoff_freq --ground_offset=$ground_offset --regularization=$reg --outdir='$outdir' '$pose_file' '$outdir/contact_nodes.txt'";
print "$cmd\n";
system $cmd;
