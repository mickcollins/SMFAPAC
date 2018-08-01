	package  colib_perl;
	require Exporter;
#        use Curses;
#        require "menuutil.pl";
	@ISA       =qw(Exporter);
	@EXPORT    =qw(appendtonl keyinfile mcscfconv ciconv callprog2 frootci intc changekeyword  add_mcscf_flag analyse_installconfig extractmolcas get_geom_info get_geom_extremes write_cigrdin);
	@EXPORT_OK =qw($found  $nconv  $niter $nintcoor $intcoor);


   sub analyse_installconfig {
    local ($filename,%installdata,$a1,$a2);
    ($filename ) =@_;
#
#  puts the contents of install.config KEY value into $installdata{KEY}=value
#
    open INSTCONFIG, "<$filename" || die "cannot find $filename exiting" ;
    $/="\n";
     while ( <INSTCONFIG> )
      { chop; 
        s/ *$//;  # remove trailing blanks 
        ($a1,$a2) = split(/\s+/,$_,2);
        if ($a2) { $installdata{"$a1"}=$a2;}
      }
    close INSTCONFIG;
    return \%installdata;
    }

    sub extractmolcas {
    local ($filename, $outmodule,$mcmodule,$modulename,@input,$i,$j,$inputline,@mymodule);
    ($filename,$modulename) =@_;
     $/="\n";
     open INP,"<$filename" ||  die "cannot find $filename exiting" ;
     @input=<INP>;
     # remove \n
     for ($i=0; $i<$#input; $i++) { chop $input[$i];
     # eliminate comment lines (*)
     # eliminate lines refering to shell execution (!)
                                   if (grep (/^\*/,$input[$i]) ) {$input[$i]="";}}
     $inputline=join(':',@input);
     #remove empty lines
     $inputline=~s/::/:/g;
     #split into separate sections for each module
     @input = split (/End of input/i,$inputline);
     $notfound=1;
     $outmodule="";
     for ($j=0;$j<$#input; $j++)
      { # identify module name
        @mymodule=split(/:/,$input[$j]);
        for ($i=0; $i<$#mymodule; $i++) { $_=$mymodule[$i];if (/&(.*) &END/) { $mcmodule=$1; $mcmodule=~s/ //g; last;}}
# also consider the case that several runs of the same type are to be executed
# e.g. CASPT2 with several type of fock matrices
        if ($modulename eq $mcmodule ) 
         { $outmodule=$outmodule.":".$input[$j]."End of input".":";  # $mcmodule=~s/:/\n/g; 
                                      $notfound=0;} 
      }
     if ($notfound) { $outmodule=-1;}
     return $outmodule;
    }

#
#==========================================================
#
   sub keyinfile {
#
#==========================================================
#
# returns the value of the namelist variable $key in file $filename
# namelist variables are separated by , or newline
# the routine stops parsing, when the key is found for the first
# time and all values corresponding to this key are returned.
#==========================================================
#
     local ($i,$j,$k,@keysperline,$found,@keyval,$filename,$key);
     ($filename,$key)=@_;
     local ( $finished);
#
      open (ANYFILE, "$filename");
      $/="\n";
      @keyval=<ANYFILE>;
      close ANYFILE;
      $found="";
      foreach $i ( @keyval ) {
      chop $i ;
      $i =~ s/ *$//g;
      @keysperline = split ',',$i ;
      $k=0;
      $finished=0;
        for ( $j= 0 ; $j <= $#keysperline; $j++)
         {  $_=$keysperline[$j];
           if ( /$key/i ) {$found=$keysperline[$j]; $k=1;next;}
           if ( ( ! /[\$&A-Z\/]/i ) && $k )
                { $found = join ':',$found,$keysperline[$j] ; }
           if (( /[\$&A-Z\/]/i ) && $k ) {$finished=1;}
         } # end of: for ( $j= 0 ; $j <= $#keysperline; $j++)
       if ( $finished) {last;}
       }
       $found =~ s/^.*= *//g;
       return $found;
      }
#
#======================================================================
#
  sub mcscfconv
{
  ($filein)=@_;
  my ($dum);
#
    $nconv=0;
      open(FILE,"$filein") or die" cannot open file: $filein\n";
      $/="\n";
         while (<FILE>)
         {
            if (/final mcscf/)
             {
              $_=<FILE> ; chop ; s/^ *//;
              ($niter)=split(/\s+/,$_,3);
              $nconv = /not conv/;
             }
         }
    close FILE;
#
    return $nconv,$niter;
#
} # end of: sub mcscfconv
#
#======================================================================
#
#
  sub ciconv
{
  ($filename)=@_;
  my ($dum,$niter,$nconv);
#
    $nconv=0;
      open(FILE,"$filename") or die" cannot open file: $filename\n";
      $/="\n";
         while (<FILE>)
         {
            if (/iterations/)
             {
              chop ; s/^ *//;
              (@dum[0..4],$niter)=split(/\s+/,$_,7);
              $nconv = /not reached/;
             }
         }
    close FILE;
#
    return $nconv,$niter;
#
} # end of: sub ciconv
#
#
#
  sub callprog2 {
  #======================================================================
  # call a fortran program , check for successful execution
  #  used for scripts outside of runc
  #======================================================================
   use Shell qw(date );
   local ($cbus,$errno,$x,$prog);
   $cbus = $ENV{"COLUMBUS"};
   $prog = $_[0] ;
   $errno=0;
   print "calling $prog ...\n";
   $errno = system ("$cbus/$prog 2>> runcerror");
#
   if (!-s "bummer"){$errno = $errno + 128 ;
    print "file bummer missing\n"; system "pwd;";}
    open (BUMMER, "bummer");
     while ( <BUMMER> )
      {
       if ( ! /normal|warning/) {$errno = $errno + 256 ;
       print "content of bummer:\n $_";}
      }
    close (BUMMER);
   if ($errno)
    {
     $!=11;
     print "Error occurred in $prog errno=$errno \n";
     die "Error occurred in $prog errno=$errno \n";
    }
   unlink ("bummer");
   }
#
#==========================================================
#
# extracts the root followed in CI from ciudgls*
#==========================================================
#
   sub frootci {

     local ($energy,$found);
     ($filename)=@_;

      open (ANYFILE, $filename) or die " Error in open file $filename\n";
      $/="\n";
       while ( <ANYFILE> )
        {if ( /vector at position/ )
          { $found=$_; last;}
        }
        chop $found;
       $energy=$found;
       $energy=~s/^.*energy=//g;
       $energy=~s/ //g;
       $found =~ s/energy.*$//g;
       $found =~ s/[ a-zA-Z]//g;
#      print "$energy\n";
       return ($found,$energy);
      }
#
#==========================================================
#
  sub intc {
      my (@coortype,$i,@field,$string);
      @coortype=qw(stre bend out tors lin1 lin2 tor2);
      @hessdiag=();
      $/="\n";
      $intcoor="";
      open(INTC,"intcfl") or die(" File: intcfl missing; performe gradient input first!");
      $title=<INTC>;
      $i=1;
      $nintcoor=0;
      $icount=0;
      while (<INTC>) {
      $line=$_;
      $line=~tr/A-Z/a-z/;
      @field=split('');
      $string=join'',@field[20..23];
      $kstring=$field[0];
      $string=~tr/A-Z/a-z/;
      $string=~s/\W.*//;
      $kstring=~s/\W.*//;
      $string=~tr/A-Z/a-z/;
      $kstring=~tr/A-Z/a-z/;
        for ($i=0; $i<=7;$i++)
         {
          if ($i==7)
           {
# read the force constat diagonals
            $line=~s/^ *//;
            (@tmp)=split(/\s+/,$line,8);
            @hessdiag= (@hessdiag,@tmp);
           }
          if ($string eq $coortype[$i] && $kstring eq "k"){$intcoor=$intcoor.$line;$nintcoor++;last}
          if ($string eq $coortype[$i]){$intcoor=$intcoor.$line;last}
         } # end of for ($i=0; $i<=6;$i++)
      } # end of: while (<INTC>)
      close(INTC);
     return $nintcoor,$intcoor,\@hessdiag ;
  } # end of: sub intc
#
#===========================================================
#

     sub changekeyword {
       local (@tmpv,$file_old,$file_new,$keyword,$value,$found);
       ($file_old,$file_new,$keyword,$value)=@_;
#
       $/="\n";
       open OLD, "<$file_old";
       while(<OLD>) {if (! /^\s+$/ ) {push @tmpv,$_;}
                     if ( /[&\$]end\s+$/ ) { last;}}
       @tmpvv=();
       while (<OLD>) { push @tmpvv,$_;}
       close OLD;
#      print "changekyword old file $file_old \n";
#      system ("cat $file_old");
        
#      print "in tmpv=",join(':',@tmpv),"NNNN\n";
#
       open NEW, ">$file_new";
       $found=0;
       while ( $#tmpv > 0 )
        { $_ = shift @tmpv;
         if (/$keyword/i) { print NEW "  $keyword=$value\n";$found=1; }
#                           print "Substituting $keyword=$value\n"; }
         else         { print NEW $_;}
        }
#      if (! $found == 1){ die "keyword $keyword not found in file: $file_old\n";}
       if (! $found == 1) { # print "Adding $keyword=$value\n";
                            print NEW " $keyword=$value\n"; }
       print NEW shift @tmpv;
       print NEW @tmpvv;
       close NEW;
#      print "changekyword new file $file_new \n";
#      system ("cat $file_new");
       return 0;
      }

#===========================================================
sub add_mcscf_flag {
   local($flag);
   ($flag)=@_;
#fp: central routine for modifying mcscfin - this could be called from makmc.pl or runc
#   maybe this routine should be made more consistent by checking if the flag is already there
           print "adding flag $flag to mcscfin\n";
           $/="\n"; my @inp;  if (! -f "mcscfin") {die ("could not find mcscfin\n");}
                        open MCSCFIN,"<mcscfin";
           @inp=<MCSCFIN>;
           close MCSCFIN;
           for ($i=1;$i<$#inp; $i++)
           {
              if (grep(/NPATH/i,$inp[$i]))
              {
                 chop $inp[$i]; $inp[$i]=~s/, *$//; $inp[$i]=~s/ *$//; $inp[$i].=",$flag\n";last;
              }
           }
           open MCSCFIN,">mcscfin";
           print MCSCFIN @inp;
           close MCSCFIN;
}
#===========================================================
#
     sub appendtonl {
       local (@tmpv,$file_old,$file_new,$string);
       ($file_old,$file_new,$string)=@_;
#
       $/="\n";
       open OLD, "<$file_old";
       while(<OLD>) {if (! /^\s+$/ ) {push @tmpv,$_;}
                     if ( /[&\$]end\s+$/ ) { last;}}
       open NEW, ">$file_new";
       print NEW @tmpv;
       print NEW $string;
       close NEW;
       return 0;
      }
#

  sub get_geom_info
  {
    my ($filgeo)=@_;
    $num_at = 0;
    $charge = 0.;
    $mass = 0.;
    %sumformula = ();
    open(GEOM,"<$filgeo");
    while(<GEOM>)
    {
      $num_at++;
      @line = split(/\s+/,$_);
      $symb = $line[1];
      #$sumformula{"$symb"}+=1;
      if ( exists $sumformula{"$symb"} ) {$sumformula{"$symb"}+=1;}
      else {$sumformula{"$symb"}=1;}
      $charge += $line[2];
      #print "@line";
      #print "line2: $charge $line[2]\n";
      $mass += $line[6];
    }
    close(GEOM);
    return ($num_at, $charge, $mass, %sumformula);
  }
  
  sub get_geom_extremes
  {
    # get the values of the outer most atoms
    $xmin =  1000.;
    $xmax = -1000.;
    $ymin =  1000.;
    $ymax = -1000.;
    $zmin =  1000.;
    $zmax = -1000.;
    
    my ($filgeo)=@_;
    open(GEOM,"<$filgeo");
    while(<GEOM>)
    {
        @line = split(/\s+/,$_);
        $x = $line[3];
        $y = $line[4];
        $z = $line[5];
        
        if ($x < $xmin) {$xmin = $x;}
        if ($x > $xmax) {$xmax = $x;}
        if ($y < $ymin) {$ymin = $y;}
        if ($y > $ymax) {$ymax = $y;}
        if ($z < $zmin) {$zmin = $z;}
        if ($z > $zmax) {$zmax = $z;}
    }
    return ($xmin, $xmax, $ymin, $ymax, $zmin, $zmax);
  }
  
  sub write_cigrdin
  {
  #fp: write cigrdin (called from makmc.pl and cidrtms.pl
    open CIGRDIN, ">cigrdin";
    print CIGRDIN " \&input \n nmiter= 100, print=0, fresdd=1, \n";
    print CIGRDIN " fresaa=1, fresvv=1, \n"; 
    print CIGRDIN " mdir=1, \n"; # if flag 11 is set in mcscfin 
    print CIGRDIN " cdir=1, \n"; # if the mcscf calc. was run with linear
                                # convergence only.
    print CIGRDIN " rtol=1e-6, dtol=1e-6,\n";
    print CIGRDIN " wndtol=1e-7,wnatol=1e-7,wnvtol=1e-7 \n \&end\n";
    close CIGRDIN;
}
