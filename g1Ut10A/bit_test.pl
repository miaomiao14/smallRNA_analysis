#!/usr/bin/env perl
$a = 0b00110;
$a |= 1<<3;
printf ("%b\n", $a & 1<<2);
printf ("%b\n", $a | 1<<10);
printf ("%b\n", ~($a));
