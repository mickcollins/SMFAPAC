#!/usr/bin/perl
#****************************************************************************
# install_hints.pl
#
# Function: Help user with termcap/terminfo tweeks for perlmenu.pm file.
#
# Version:  4.0
#
# Date:     February 1997
#
# Author:   Steven L. Kunz
#           Networked Applications
#           Iowa State University Computation Center
#           Ames, IA  50011
#****************************************************************************
use Config;
use lib join('/',$ENV{"COLUMBUS"},"CPAN") ;
use Curses;
&initscr;

$method1 = $method2 = $method3 = $sim_getcap = 0;

#
# Check for "getcap"
#
eval { $ku = &getcap('ku'); };
if ($@) {
  if ($] >= 5.001) {
    $method1 = 1;
    $sim_getcap = 1;
  } else {
    $method1 = 0;
    $sim_getcap = 0;
  }
} else {
  $method1 = 1;
}

#
# Check for "tigetstr"
#
eval { $ku = &tigetstr('kcuu1'); };
if ($@) {
  $method2 = 0;
} else {
  $method2 = 1;
}

#
# Check for "tput"
#
eval { $ku = `tput kcuu1`; };
if ($@||$?) {
  $method3 = 0;
} else {
  $method3 = 1;
}

&endwin;

#
# Summarize what we found
#
$total = $method1+$method2+$method3;
if ($total) {
  print "You should try one of the following methods:\n\n";
  if ($method1) {
    print "- Method 1 (getcap)\n";
    print "  (If the demo did not work you probably have a buggy getcap)\n";
    if ($sim_getcap) {
      print "  (You need to simulate the getcap function)\n";
    }
  }
  if ($method2) { print "- Method 2 (tigetstr)\n"; }
  if ($method3) { print "- Method 3 (tput)\n"; }
} else { print "Nothing works!  You probably need to install \"tput\".\n"; }
print "\nRefer to the INSTALLATION document for what to do next.\n";
