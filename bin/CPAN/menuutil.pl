#****************************************************************************
# menuutil.pl -- Curses Utility Functions
#
# Version: 1.1
#
# Author:  Steven L. Kunz
#          Networked Applications
#          Iowa State University Computation Center
#          Ames, IA  50011
#
# Official PerlMenu WWW home page:
#          http://www.cc.iastate.edu/perlmenu/
#
# Official Package Distributions:
#          ftp://ftp.iastate.edu/pub/perl
#
# Bugs:    skunz@iastate.edu
# Cheers:  skunz@iastate.edu
#
# Date:    Version 1.0 -- January, 1995 -- Original version 
#          Version 1.1 -- February, 1996 -- Added "clear_screen"
#
# Notes:   Perl4 - Requires curseperl
#          Perl5 - Requires William Setzer's "Curses" extension
#
# Routines:	&top_title      -- Clear screen, center title
#		&print_nl	-- Print text on screen
#		&new_line	-- Skip to new line on screen
#		&query		-- Ask question, validate response key
#		&pause		-- Pause (with prompt)
#		&popup_ask	-- Ask question in pop-up overlay box
#		&clear_screen   -- Clear the screen
#
# PerlMenu - Perl library module for curses-based menus & data-entry templates
# Copyright (C) 1992-97  Iowa State University Computation Center
#                        Ames, Iowa  (USA)
#
#    This Perl library module is free software; you can redistribute it
#    and/or  modify it under the terms of the GNU Library General Public 
#    License (as published by the Free Software Foundation) or the
#    Artistic License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of 
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Library General Public License for more details.
#
#    You should have received a copy of the GNU Library General Public
#    License along with this library; if not, write to the Free
#    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#****************************************************************************

#***********
#  CLEAR_SCREEN -- Clear screen (and row, col values)
#
#  Arguments:	Boolean flag (0=don't do refresh, 1=do refresh)
#
#  Notes:	
#***********
sub clear_screen {
  local($do_refresh) = @_;
  &clear();
  $row = $col = 0;
  if ($do_refresh) { &refresh(); }
}

#***********
#  TOP_TITLE -- Center title at top of screen
#
#  Arguments:	Title string
#
#  Notes:	Title string is normally placed on screen centered in
#		"standout" rendition.  If you do NOT want standout
#		rendition start the string with a "-" (it will not appear
#		in the title).
#***********
sub top_title {
  local($title) = @_;
  local($title_col) = 0;
  local($title_standout) = 1;

# Check for title format character.
  if (substr($title,0,1) eq '-') {
    $title = substr($title,1);
    $title_standout = 0;
  }

# Truncate or center title as necessary
  if (length($title) > $COLS) {
    $title = substr($_[1],0,$COLS);
  }
  else {
    $title_col = int($COLS/2) - int(length($title)/2);
  }

# Clear screen and put title at the top.
  &clear();  $row = $col = 0;
  &move(0,$title_col);
  if ($title_standout) { &standout(); }
  &addstr($title); $col += length($title);
  if ($title_standout) { &standend(); }

# Setup row,col for those that follow.
  $row = 2; $col = 0;
  &move($row,$col);
  &refresh();
}

#**********
# PRINT_NL -- Print text followed by optional new-line(s) on screen.
#
# Arguments:	- Text string
#		- New-line count. Optional
#**********
sub print_nl {
  local($text,$skip) = @_;

# Check for end-of-line wrap
  while ($col+length($text) > $COLS) {
    &addstr(substr($text,0,$COLS));
    $row++; $col = 0;
    &move($row,$col);
    $text = substr($text,$COLS);
  }

# Put last (or only) chunk on the screen
  &addstr($text);
  $col += length($text);

# Add any new-lines specified
  if ($skip) { &new_line($skip); }
}

#**********
# NEW_LINE -- Skip to new line (with end-of-page clear)
#
# Arguments:	New-line count
#**********
sub new_line {
  local($skip) = @_;

  $row += $skip;
  $col = 0;
  if ($row > $LINES - 1) {
    print "\007";
    &refresh();
    &clear(); $row = 0; $col = 0;
    &refresh();
  }
  &move($row,$col);
  &refresh();
}

#**********
#  QUERY -- Perform query of user
#
#  Arguments:	- Prompt text
#		- List of single-char values allowed ("ynq")
#
#  Returns:	One of the single characters from the "allowed" list.
#**********
sub query {
  local($prompt,$allowed) = @_;
  local($ch,$sel_prompt);
  local($nsl,$i) = 0;

# Construct "allowed selections" part of the prompt.
  $nsl = length($allowed);
  if ($nsl < 1) { return; }
  $sel_prompt = "(";
  while ($i < $nsl) {
    $sel_prompt .= substr($allowed,$i,1)."|";
    $i++;
  }
  chop($sel_prompt);
  $sel_prompt .= ")";

# Prompt until we get a good answer.
  while (1) {
    &getyx($window,$row,$col);
    $ch = &menu_getstr($row,$col,"$prompt $sel_prompt ",0); &new_line(1);
    if ($ch ne "") {
      $ch = substr($ch,0,1);
      $ch =~ tr/A-Z/a-z/;
      if (index($allowed,$ch) >= 0) { return($ch); }
    }
    print "\007";
    &new_line(1);
    &print_nl("Invalid answer - must be one of $sel_prompt.",1);
  }
}

#**********
#  PAUSE -- Pause until key is pressed
#
#  Arguments:	Optional message text string.
#		Defaults to "[Press any key to continue]"
#
#  Returns:	Nothing.
#
#  Note:	This routine always skips to a new line first.
#**********
sub pause {
  local($msg) = @_;

# Skip a row and put our message.
  $row += 1; $col=0;
  &move($row,$col);
  if ($msg eq "") { $msg = "[Press any key to continue]"; }
  &addstr($msg);

# Leave cursor on a fresh line.
  $row++;
  &move($row,$col);

# Refresh screen and wait for a keypress.
  &refresh();
  $ch = &getch();
}

#**********
#  POPUP_ASK -- Pop-up box to ask a question
#
#  Function:	Provides an outlined box, centered on the screen, to allow
#		edited data-entry by the user.
#
#  Arguments:	Note: All parms optional
#		- Prompt string
#		- Maximum data input length
#		- Default value string
#		- Boolean flag (0=show, 1=hidden)
#		- Boolean flag (0=alpha-numeric, 1=numeric)
#
#  Returns:	Value supplied by the user (may be null).
#**********

sub popup_ask {
  local($prompt,$maxlen,$default,$noshow,$numeric) = @_;
  local($box_row,$box_col,$box_win,$data_win,$val);

# Compute some key values
  if (!$maxlen) {
    $maxlen = $COLS - length($prompt) - 2;
  }
  $box_row = int($LINES/2) - 2;
  $box_col = int($COLS/2) - int((length($prompt)+$maxlen)/2) - 1;
  if ($box_col < 0) { $box_col = 0; } 

# Create the main pop-up box (with the prompt in it)
  $box_win = &newwin(3,length($prompt)+$maxlen+2,$box_row,$box_col);
  &wmove($box_win,1,1);
  &waddstr($box_win,$prompt);
  &box($box_win,ord('|'),ord('-'));
  &wrefresh($box_win);

# Create a second window for data entry within the main pop-up box
  $data_win = &newwin(1,$maxlen,$box_row+1,$box_col+length($prompt)+1);
  $val = &main'menu_getstr(0,0,"",0,$default,$maxlen,
			   $noshow,$numeric,$data_win);

# Clean up the data-entry window
  &wclear($data_win);
  &wrefresh($data_win);
  &delwin($data_win);

# Clean up the pop-up window
  &wclear($box_win);
  &wrefresh($box_win);
  &delwin($box_win);

# Clean up the main window and return value
  &touchwin($window);
  &refresh();

  $val;
}

1;
