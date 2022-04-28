#!/usr/bin/env perl
use Mocap::BVH;
use Mocap::BVH::Joint;
use strict;

die "usage: $0 input.bvh output.json\n" unless @ARGV == 2;


my ($in_fn, $out_fn) = @ARGV;
my $bvh = Mocap::BVH->load($in_fn);
open my $fout, '>', $out_fn;

sub print_joint {
    my $joint = shift;
    print $fout "{\n";
    print $fout "\"name\": \"", $joint->name, "\",\n";
    print $fout "\"pos\": [", join(', ', $joint->offset), "],\n";
    if ($joint->parent) {
	print $fout "\"type\": \"ball\"";
    } else {
	print $fout "\"type\": \"free\"";
    }

    my @ends;
    if ($joint->end_site) {
	push @ends, [$joint->end_site];
    }
    for my $jt($joint->children) {
	push @ends, [$jt->offset];
    }

    if (@ends) {
	print $fout ",\n";
	print $fout "\"shape\": [\n";
	my $first = 1;
	for my $it (@ends) {
	    print $fout ",\n" unless $first;
	    $first = 0;
	    print $fout "{\n";
	    print $fout "\"type\": \"box\",\n";
	    my @size = (0.05) x 3;
	    my $idx = 0;
	    for my $i (1 .. 2) {
		$idx = $i if abs($it->[$i]) > abs($it->[$idx]);
	    }
	    $size[$idx] = abs($it->[$idx]);
	    print $fout "\"size\": [", join(', ', @size), "],\n";
	    
	    my @pos;
	    for (@$it) {
		push @pos, $_ / 2;
	    }
	    print $fout "\"pos\": [", join(', ', @pos), "]\n";

	    print $fout "}";
	}
	print $fout "\n]";
    }

    if ($joint->children) {
	print $fout ",\n";
	print $fout "\"children\": [\n";
	my $first = 1;
	for my $jt ($joint->children) {
	    print $fout ",\n" unless $first;
	    $first = 0;
	    print_joint($jt);
	}
	print $fout "]\n";
    }

    print $fout "\n}";
}

print_joint($bvh->root);
close $fout;
