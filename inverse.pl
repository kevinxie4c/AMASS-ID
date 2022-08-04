#!/usr/bin/env perl
use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;

my $outdir = "output";
my $pose_file;
my @contact_opts;
my @inverse_opts;

sub add_contact_opts {
    push @contact_opts, "-$_[0]", $_[1];
}

sub add_inverse_opts {
    push @inverse_opts, "-$_[0]", $_[1];
}

GetOptions(
    "j|char_file=s" => \&add_inverse_opts,
    "f|frame_time=f" => \&add_inverse_opts,
    "c|cutoff_freq=f" => \&add_inverse_opts,
    "u|mu=f" => \&add_inverse_opts,
    "s|step_length=i" => \&add_inverse_opts,
    "r|regularization=f" => \&add_inverse_opts,
    "o|outdir=s"=> sub { push @inverse_opts, "-$_[0]", $_[1]; push @contact_opts, "-$_[0]", $_[1]; $outdir = $_[1] },
    "A|start_frame=i" => \&add_inverse_opts,
    "E|end_frame=i" => \&add_inverse_opts,
    "F|filter_type=s" => \&add_inverse_opts,
    "S|use_sim_state" => sub { push @inverse_opts, "-$_[0]" },
    "t|threshold=f" => \&add_contact_opts,
    "g|ground=f"    => \&add_contact_opts,
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
    -t, --threshold=double
    -g, --ground=double
USAGE
} else {
    $pose_file = shift @ARGV;
}

mkdir $outdir unless -e $outdir;

my @cmd;
@cmd = ('python', 'find-contact-node.py', @contact_opts, $pose_file);
print "@cmd\n";
system @cmd;
@cmd = ('./inverse.out', @inverse_opts, $pose_file, "$outdir/contact_nodes.txt");
print "@cmd\n";
system @cmd;
