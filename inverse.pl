#!/usr/bin/env perl
use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;

my $char_file = "data/character.json";
my $frame_time = 1/120;
my $cutoff_freq = 2;
my $ground_offset = 0;
my $reg = 0.0001;
my $outdir = "output";
my $pose_file;
my @opts;

sub add_opts {
    push @opts, "-$_[0]", $_[1];
}

GetOptions(
    "j|char_file=s" => \&add_opts,
    "f|frame_time=i" => \&add_opts,
    "c|cutoff_freq=i" => \&add_opts,
    "u|mu=f" => \&add_opts,
    "s|step_length=i" => \&add_opts,
    "r|regularization=i" => \&add_opts,
    "o|outdir=s"=> sub { push @opts, "-$_[0]", $_[1]; $outdir = $_[1] },
    "A|start_frame=i" => \&add_opts,
    "E|end_frame=i" => \&add_opts,
    "F|filter_type=s" => \&add_opts,
    "S|use_sim_state" => sub { push @opts, "-$_[0]" },
);

if (@ARGV != 1) {
    die <<"USAGE";
usage: $0 [options] pose_file

options:
    -j, --char_file=string
    -f, --frame_time=double
    -c, --cutoff_freq=double
    -u, --mu=double
    -s, --step_length=int
    -r, --regularization=double
    -o, --outdir=string
    -A, --start_frame=int
    -E, --end_frame=int
    -F, --filter_type=none|position|velocity
    -S, --use_sim_state
USAGE
} else {
    $pose_file = shift @ARGV;
}

mkdir $outdir unless -e $outdir;

my @cmd;
@cmd = ('python', 'find-contact-node.py', $pose_file, $outdir);
print "@cmd\n";
system @cmd;
@cmd = ('./inverse.out', @opts, $pose_file, "$outdir/contact_nodes.txt");
print "@cmd\n";
system @cmd;
