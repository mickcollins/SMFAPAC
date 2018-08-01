#!/usr/bin/perl
#use strict;


#use warnings;

#**************************************************************************
# SMFA control script
#**************************************************************************
use Config;
use Shell qw(date pwd);
use File::Copy 'cp';
use Cwd;


BEGIN {
$place=`which SMFA`;
chomp($place);
@arr=split("/bin/SMFA",$place);
$CODEDIR=join("",$arr[0],"/bin");
$EXEDIR=join("",$arr[0],"/exe");
system("echo $CODEDIR > CODENAME");
system("echo $EXEDIR > EXENAME");
}

# Perl5+Curses ONLY!
# Comment these lines for use with Perl4/curseperl
use lib "$CODEDIR/CPAN";
BEGIN { $Curses::OldCurses = 1; }
use Curses;                     # PerlMenu needs "Curses"
use perlmenu;                   # Main menu package (Perl5 only)

$cont="RUN";

$| = 1;                                # Flush after every write to stdout


  $JDIR = pwd();             # job directory
  $JDIR =~ s/\n//;
  $WORK = "WORK";            # calculations are carried out in WORK
  $MODIR = "MOCOEFS";         # MO coefficient files
  $LIST = "LISTINGS";        # listing files
  $REST = "RESTART";         # MCSCF restart files
  $GEOMDIR = "GEOMS";
  $OUTSMFA = "OUT_SMFA";

# global variables for SMFA
 
  $fraglvl = 0;
  $cutoff = 0;
  $nchgs = 0;
  $deriv = 0;

#
# Required global variables are $window, $row, and $col.
# These variables are used by the menuutil.pl routines.
#
   $window=initscr();
   &menu_curses_application($window);

# Default prefs active at start
  $numbered_flag = 1;       # Numbered menus
  $num_pref = "numbered";
  $gopher_pref = "default"; # Non-gopherlike arrows/scrolling
  $gopher_flag = 0;
  $mult_pref = "single";    # Single column menus
  $mult_flag = 0;
  $arrow_pref = "arrow";    # Arrow selection indicator menus
  $arrow_flag = 0;

  $menu_default_top = 0;    # Storage for mainline top item number.
  $menu_default_row = 0;    # Storage for mainline arrow location.
  $menu_default_col = 0;    # Storage for mainline arrow location.
  $row = $col = 0;      # Storage for row/col for menuutil.pl
  $title_cnt = 0;       # To trigger different subtitles/bottom titles

   while (1) {  # always true

     &menu_init($numbered_flag, "SMFA CONTROL",0,
                ,"\n-main menu options");
     &menu_item ("Set up input","inp");
     &menu_item ("SMFA examines the input and prints comments or queries to OUT_SMFA","viewinp");
#     &menu_item ("\(Optional\) Include embedded charges for explicit solvents","solventcharges");
     &menu_item ("Fragment the molecule","fragonly");
     &menu_item ("Time limit, memory, queue, etc","makepbsfile");
     &menu_item ("Run all the electronic structure calculations","runall");
     &menu_item ("Restart the electronic structure calculations","restart");
#     &menu_item ("Analyse the results","anal_MAC");
     &menu_item ("Utilities","utilities");
     &menu_item ("Exit SMFA","Exit");

     $sel = &menu_display("",$menu_default_row,$menu_default_top,
                             $menu_default_col);
     if ($sel eq "Exit") { last; }
     if ($sel eq "%EMPTY%" ) { die "not enough screen lines \n";}
    if ($sel ne "%UP%" )
      { # call the corresponding subprogram (assuming existence)
         $ret=&$sel();
         if ($ret) {$window=&initscr();
                    &menu_curses_application($window);
                    }}

   }

  &endwin();
  system("stty sane");
  exit(0);
#
######################### input control ##################################
  sub inp {
#
#
    while (1) {  # always true
     &menu_init($numbered_flag, "SMFA Control ",0,
                ,"\n-submenu: Input control");
     &menu_item ("Specify the Level of fragmentation","get_frag");
     &menu_item ("Specify the cartesian coordinates file","get_coords");
     &menu_item ("Specify the electronic structure calculation","get_abdata");
     &menu_item ("Specify hydrogen bonding and any unusual bonding (optional)","get_bonds");
     &menu_item ("Specify charges for metals & other atoms (optional)","get_charges");
     &menu_item ("Exit input control","Exit");

     $sel = &menu_display("",$menu_default_row,$menu_default_top,
                             $menu_default_col);

     if ($sel eq "%EMPTY%" ) { die "not enough screen lines \n";}
     if ($sel eq "Exit") {
      last;
     }

       if ($sel ne "%UP%" )
        { # call the corresponding subprogram (assuming existence)
         $ret=&$sel();
           if ($ret)
            {
            $window=&initscr();
            &menu_curses_application($window);
            } # end of: if ($ret)
        } # end of: if ($sel ne "%UP%" )
    } # while (1)
  } # end of: sub inp

