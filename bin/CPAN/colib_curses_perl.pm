#fp: separate the part of colib that needs Curses from the remaining part
#    then Columbus can still be run even if CPAN is not installed.
	package  colib_curses_perl;
	require Exporter;
        use Curses;
        require "menuutil.pl";
	@ISA       =qw(Exporter);
	@EXPORT    =qw(displaytext mvstr);
	


#---------------------------------------------------------------------------
#
  sub displaytext
  { my ($txt,@txt,@xline) = (@_);
    my ($currline,$txt1,$txt2,$i);
    # width of window: $COLS
    # text is given as a single line with \n
    # and need to be broken into separate lines at word boundaries
    # if necessary
    # strip off leading blanks

    @txt=split(/\n/,$txt);
    $currline=2;
    foreach $txt1 ( @txt )
     {  $txt1=~s/  */ /g; if (length($txt1) < ($COLS-2) ) { $txt1=~s/^ *//;
                                         &mvstr($currline,1,$txt1);
                                         $currline++;}
        else { @xline = split(/\s+/,$txt1);
               $txt2='';
               for ( $i=0;$i <= $#xline;$i++)
                 { if (length($txt2)+length($xline[$i]) > ($COLS-3) )
                      { $txt2=~s/^ *//;
                        mvstr($currline,1,$txt2);
                        $currline++;
                        $txt2=$xline[$i];}
                    else { $txt2 = join(' ',$txt2,$xline[$i]);}
                 }
                if ( length($txt2) > 1 ) {$txt2=~s/^ *//;
                                          mvstr($currline,1,$txt2);
                                          $currline++;}
              }
     }
     $row=$currline++;
     &pause("press return to continue");
    return;
  }

 
  sub mvstr
  { my ( $r,$c,$str) = (@_);
    &move ( $r, $c);
    &addstr( $str);
    return;
  }

#===========================================================
#