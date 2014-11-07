#!/usr/local/bin/perl

use warnings;
use strict;
my $tab = "\t";
my $annotation = ">tr|F7UQ11|F7UQ11_SYNYG Spore coat polysaccharide biosynthesis protein OS=Synechocystis sp. (strain PCC 6803 / GT-S) GN=spsC PE=3 SV=1";
$annotation =~ /^>\w{2}\|(.*?)\|.*?\s(.*?)\sOS=.*?GN=(.*?)\s.*/;
        print $1, $tab, $2, $tab, $3, $tab;
