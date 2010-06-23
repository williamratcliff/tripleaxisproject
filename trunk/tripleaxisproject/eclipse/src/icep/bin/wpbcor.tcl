#!/usr/bin/wish

# Tk-gui for polarized beam correction
# using pipe to pbscript executable
# and plotting with BLT(tk extension) or gnuplot.
# Can run without plotting. Then just get textfile results.
# Can save your session.

# this version has been tested with Linux and Windows (using cygwin)

# Unices could start this file with  #!wishprogpath
# For Windows must run: wishprogpath wishscriptpath args
# where wishscriptpath is this file
# and args can locate the pbscript executable and/or the gnuplot executable

# If the pbscript executable is not named on the command line,
# try to find pbscript.exe using recursive FindFile on env(PATH) dirs

######## get PATH
if { [string match windows $tcl_platform(platform)] } {
    set pbprogtail pbscript.exe    
    # pipe enabled version of gnuplot for windows
    set gnuprogtail pgnuplot.exe
    # switch to / in pathes for tcl inistead of the escape backslash
    set pathes [regsub -all {\\} $env(PATH) /]
    set pathlist [split $pathes \;]
    # never search the WINDOWS system folders
    while { [set i [lsearch $pathlist C:/WINDOWS*]] >= 0 } {
	set pathlist [lreplace $pathlist $i $i]
    }
} else {
    set pbprogtail pbscript
    set gnuprogtail gnuplot
    # I havent yet checked that this works for BSD Mac
    set pathlist [split $env(PATH) :]
}

set pbprogfound 0
set gnuprogfound 0


############## startup utility procs
# NB array set aname list only appends elements in list to current array
# need a proc to clear all elements or to set the array = values from list

proc arrayClear {aname} {
    upvar #0 $aname v
    foreach nam [array names v] {
	unset v($nam)
    }
}
proc arraySet {aname lst} {
    arrayClear $aname
    upvar #0 $aname v
    array set v $lst
}

# recursive file finder. Remember to cd to orig working dir when done.
proc FindFile { startDir namePat } {
    set pwd [pwd]
    if [catch "cd $startDir" err] { 
	    Dialog_Warn "Can't cd to $startDir.\n"
	    return ""
    }
    #puts "FindFile: looking in $startDir for $namePat"
    set match [glob -nocomplain -- $namePat]
    if { [llength $match] > 0 } {
	return [file join $startDir [lindex $match 0]]
    }
    foreach fil [glob -nocomplain *] {
	if [file isdirectory $fil] {
	    set match [FindFile [file join $startDir $fil] $namePat]
	    if { [llength $match] > 0 } {
		return [file join $startDir $fil [lindex $match 0]]
	    }
	}
    }
    cd $pwd
    return ""
}

############ dialog procedures

set warnStack {}
set msgfont ""
array set dialog {}
proc Dialog_Warn {string} {
    global warn warnStack msgfont

    if { [llength [after info]] > 0 } {
	after idle [list Dialog_Warn $string]
	return
    }

    set f .warn
    set c $f.f
    set b $c.buttons
    set m $c.msg
    if [Dialog_Create $f "Warning"] {
	frame $c
	pack $c -side top -expand true -fill both
	frame $b
	frame $m
	pack $m -side top -expand true -fill both
	pack $b -side top

	button $b.ok -text OK -command {set warn(ok) 1}
	pack $b.ok -side top
	message $m.msg -text $string -fg red -padx 10
	pack $m.msg -side top -expand true -fill both


	set msgfont [$m.msg cget -font]
    } else {
	if { [string match .warn [grab current]] } {
	    lappend warnStack $string
	    return
	}
    }

    set lines [split $string \n]

    set pixwid 0
    foreach line $lines {
	set ll [font measure $msgfont $line]
	if { $ll > $pixwid } { set pixwid $ll }
    }
    set pixwid [expr 10 + $pixwid]

    # sometimes dialog comes up with text and or buttons off widget??
    # try not setting the width
    $m.msg configure -text $string
    #-width $pixwid
    set warn(ok) 0
    Dialog_Wait $f warn ""
    Dialog_Dismiss $f
}


proc Dialog_Create {top title args} {
    global dialog
    if [winfo exists $top] {
	switch -- [wm state $top] {
	    normal { raise $top }
	    withdrawn -
	    iconified {
		wm deiconify $top
		catch {wm geometry $top $dialog(geo,$top)}
	    }
	}
	return 0
    } else {
	eval {toplevel $top} $args
	wm title $top $title
	return 1
    }
}
proc Dialog_Wait {top varName {focus {}} } {
    global pausing $varName
    if { ![winfo exists $top] } { return }
    set varname $varName\(ok\)
    eval set var $$varname
    bind $top <Destroy> [list set $varname $var]
    raise $top
    update idletasks
    if {[string length $focus] == 0} { set focus $top }
    set old [focus -displayof $top]
    focus $focus
    # getting hung up on tkwait visibility
    tkwait visibility $top
    grab $top

    # continually raise the dialog window (some wm's may allow hiding it)

    set ntry 0
    set pausing 0
    set var 0
    set i 0
    while { 1 } {
	set id [after 20 "setPausing $varname 1"]
	tkwait variable $varname
	# after 20 if varname not set by user then after event happens
	# which sets pausing 1 and varname 1
	eval set var $$varname
	if { ! $pausing } {
	    after cancel $id
	    break
	}
	incr i
	if { $i > 100 } { raise $top; set i 0 }
	set pausing 0
	#incr ntry
    }
    catch {grab release $top}
    focus $old
}
proc Dialog_Dismiss {top} {
    global dialog warnStack
    if { [string match .warn $top] && [llength $warnStack] > 0 } {
	set string [lindex $warnStack 0]
	set warnStack [lreplace $warnStack 0 0]
	Dialog_Warn $string
	return
    }
    catch {
	set dialog(geo,$top) [winfo geometry $top]
	wm withdraw $top
    }
}
proc setPausing {varname value} {
    upvar $varname vname
    global pausing
    set pausing 1
    set vname $value
}

########################################################################
# program startup with command line argument processing

set DEBUG 0
# could set pbprog = pbscript location manually here. cmndline would override.
# cmndline options
# debug
# pbs*=fullpath_to_pbscript
# gnu*=fullpath_to_gnuplot
# else fullpath_to_pbscript

set wantGNU 0
foreach arg $argv {
    if { [string match deb* $arg] || [string match DEB* $arg] || \
	     [string match dbg* $arg] || [string match DBG* $arg] } {
	set DEBUG 1
	# this will open a window to monitor io to pbscript
	# also for windows, open the console to get puts output
	if { [string match windows $tcl_platform(platform)] } {
	    console show
	}
    } elseif { [string match gnu* $arg] } {
	set wantGNU 1
	if { $DEBUG } { puts "checking arg=$arg for gnu" }
	if { [regexp {gnu[^=]*=(.*)$} $arg match gnuprog] && \
	  [string length $gnuprog] > 0 } {
	    if { $DEBUG } { puts "regexp got gnuprog=$gnuprog" }
	    if { ! [file exists "$gnuprog"] } {
		Dialog_Warn "File \"$gnuprog\" doesnt exist"
		exit
	    }
	    if { ! [file executable "$gnuprog"] } {
		Dialog_Warn \
		    "File named on command line=\"$gnuprog\" not executable!!"
		exit
	    }
	    set gnuprogfound 1
	}
    } else {
	if { ! [regexp {pbs[^=]*=(.*)$} $arg match pbprog] } {
	    set pbprog $arg
	}
	if { ! [file exists "$pbprog"] } {
	    Dialog_Warn "File \"$pbprog\" doesnt exist"
	    exit
	}
	if { ! [file executable "$pbprog"] } {
	  Dialog_Warn "File named on command line=\"$pbprog\" not executable!!"
	    exit
	}
	set pbprogfound 1
    }
}

# record the startup directory, which should be users DATA directory
set currentDir [pwd]

# last chance to locate pbprogtail with FindFile on pathlist
if { ! $pbprogfound } {
    foreach trydir $pathlist {
	if { $DEBUG } { puts "trying $trydir for $pbprogtail" }
	set pbprog [FindFile $trydir $pbprogtail]
	if { [string length $pbprog] > 0 } { break }
    }
    # puts "pbprog=$pbprog"
    # NB .bat or .lnk files show up as not executable so make sure
    # to use a real executable
    if { ! [file executable "$pbprog"] } {
	Dialog_Warn "required program=$pbprog not executable!!"
	#Dialog_Warn "program=$pbprog doesnt exist!!"
	exit
    }
    # make sure to go back to the users startup directory
    cd $currentDir
}

######## check for plotting capability
# N.B. Windows BLT is still only supported for TclTk8.4

set GNU 0
set BLT 0
set PLT 0

# one gnuplot pipe per group plotting window
set gpipe ""
array set gpipes {}
proc findAndCheckGnu {} {
    # try to find gnuplot with FindFile if it hasnt already been found
    # then try to open it
    # return 1 on success
    global gpipe gnuprog gnuprogfound pathlist gnuprogtail gpipes DEBUG
    if { ! $gnuprogfound } {
	set currentDir [pwd]
	foreach trydir $pathlist {
	    if { $DEBUG } { puts "trying $trydir for $gnuprogtail" }
	    set gnuprog [FindFile $trydir $gnuprogtail]
	    if { [string length $gnuprog] > 0 } { break }
	}
	if { [string length $gnuprog] < 1 } { return 0 }
	#puts "pbprog=$pbprog"
	# NB .bat or .lnk files show up as not executable so make sure
	# to use the executable
	if { ! [file executable "$gnuprog"] } {
	    return 0
	} else {
	    set gnuprogfound 1
	}
	cd $currentDir
    }
    # make sure the program can be opened on a pipe
    
    #puts "about to open gnuplot=$gnuprog"
    if [catch {open "|\"$gnuprog\"" r+} gpipe] {
	return 0
    }
    set gpipes(1) $gpipe
    return 1
}

if { $wantGNU } {
    if { [findAndCheckGnu] } {
	set GNU 1
	set PLT 1
    } else { 
        if [catch {package require BLT}] {
	    Dialog_Warn "Couldn't find BLT or gnuplot. Plotting disabled."
	} else {
	    set BLT 1
	    set PLT 1
	}
    }
} else {
    if [catch {package require BLT}] {
        if { [findAndCheckGnu] } {
	    set GNU 1
	    set PLT 1
        } else {
           Dialog_Warn "Couldn't find BLT or gnuplot. Plotting disabled."
        }
    } else {
	set BLT 1
	set PLT 1
    }
}

if { $BLT } { namespace import blt::* }

if { $GNU } {
  # gnuplot terminal type is platform dependent
  if { [string match windows $tcl_platform(platform)] } {
    set screen windows
  } else {
    set screen wxt
  }
}

############# debug interface includes gui access to pbscript and debugserver

if { $DEBUG } {
    # setup gui for debugging io to pbscript pipe.
    # use the intxt text widget to record sent commands
    # and the outtxt widget to record replies.
    # uses globals readfinished readlines from pbcor

    set db [toplevel .debug]
    wm title $db "pbscript command debugger"
    wm iconname $db "debug"
    wm deiconify $db


    proc Scrolled_Text { f args } {
	frame $f
	eval {text $f.text \
		  -xscrollcommand [list $f.xscroll set] \
		  -yscrollcommand [list $f.yscroll set]} $args
	scrollbar $f.xscroll -orient horizontal \
	    -command [list $f.text xview]
	scrollbar $f.yscroll -orient vertical \
	    -command [list $f.text yview]
	grid $f.text $f.yscroll -sticky news
	grid $f.xscroll -sticky news
	grid rowconfigure $f 0 -weight 1
	grid columnconfigure $f 0 -weight 1
	return $f.text
    }
    
    proc sendCur {} {
	global intxt outtxt
	global readfinished readlines
	
	if { [catch {set line \
			 [$intxt get "insert linestart" "insert lineend"]}] } {
	    return
	}
	set line [string trim $line]
	if { [string length $line] < 1 } { return }
	pipeWrite $line
	$outtxt insert end "pipeWrite: $line\n"
	$outtxt insert end "pipeRead:\n"
	vwait readfinished
	# for debugging pipeRead dumps to outtxt
    }

    proc sendClr {} {
	global outtxt
	$outtxt delete 1.0 end
    }
    
    
    set sendmb [frame $db.mb]
    set sendin [frame $db.in]
    set sendout [frame $db.out]
    
    grid $sendmb -sticky news -pady 5
    grid $sendin -sticky news
    grid $sendout -sticky news
    grid rowconfigure . 1 -weight 1
    grid rowconfigure . 2 -weight 1
    grid columnconfigure . 0 -weight 1


    set sendCur [button $sendmb.ex -text "PIPE current line" \
		     -command sendCur]
    
    set sendClr [button $sendmb.clr -text "Clear Output" \
		     -command sendClr]


    grid $sendCur $sendClr -sticky news


    set inlbl [label $sendin.lbl -text "PIPE write text editor" -justify center]
    set intxt [Scrolled_Text $sendin.ftxt]
    $intxt configure -height 6 -width 60

    set outlbl [label $sendout.lbl -text "PIPE read" -justify center]
    set outtxt [Scrolled_Text $sendout.ftxt]
    $outtxt configure -height 20 -width 60
    
    grid $inlbl -sticky news
    grid $sendin.ftxt -sticky news
    
    grid $outlbl -sticky news
    grid $sendout.ftxt -sticky news
    
    update idletasks

    # will save all input to pbscript pipe in a file
    # for later debugging by running pbscript through ddd
    set debugFPT [open debugCMDS w]

    # also try connection to the debugServer (send replacement via sockets)

    # this is debugClient code using socket
    # to be used by a tcl app that wishes to use the
    # debugServer to communicate with it while running.
    # source this code into a tcl script for debugging
    # start with: connectDebugServer localhost 
    # (or whatever host the server is on)
    
    # standard socket port for debugServer is 1500
    # connect to the debug Server
    # using the connectDebugServer proc
    # Once the connection is established debugStart is called to
    # send the apps name (toplevel window name) e.g. winfo name .
    # which should be echoed back by the server.
    # Then starts listening for commands. 
    # Exec any received cmnds in global context
    # and return the result followed by an END-RESULT line
    # to signal end of result
    # use catch for errors and send back the tcl stack trace as the result.
    
    # This is more flexible than sendTk as it works for a plane tcl script.
    # It is more secure, since the tcl script must agree to receive commands
    # from the server by connecting to it.
    # In fact tk send is broken since X uses xhost list

    if { [string match windows $tcl_platform(platform)] } { console show }

    set debugsock ""
    set debugServerPort 1500
    set debugconnect 0

    proc connectDebugServer {host} {
	global debugsock debugServerPort debugconnect
	# connect to the debugServer, and return the socket channel
	if [catch {socket $host $debugServerPort} debugsock] {
	    puts "connectDebugServer $host $debugServerPort failed: $debugsock"
	    return
	} else {
	    puts "debugServer connected."
	    set debugconnect 1
	}
	fconfigure $debugsock -buffering line
	debugStart
    }
    
    proc debugStart {} {
	global debugsock
	# send the appname and get one line back from the server
	# to confirm the connection
	puts $debugsock "[winfo name .]"
	flush $debugsock
	if { [eof $debugsock] || [gets $debugsock cmd] < 0 } {
	    catch {close $debugsock}
	    puts "no initial reply from debugServer"
	} else {
	    if { [string match *OK* $cmd] } {
		puts "debugServer replied OK"
	    } else {
		puts "debugServer replied without OK"
	    }
	    # init the client cmd reader
	    flush $debugsock
	    fileevent $debugsock readable debugExec
	}
    }
    
    proc debugExec {} {
	global debugsock errorInfo
	if { [eof $debugsock] || [gets $debugsock cmd] < 0 } {
	    catch {close $debugsock}
	    puts "debugServer may have died. Socket closed."
	} else {
	    # eval the cmd in global context
	    # and return the result
	    puts "debugClient socket=$debugsock received cmd: $cmd"
	    # let the tcl app catch any errors
	    if [catch {uplevel \#0 eval $cmd} rslt] {
		#puts "cmd result: $rslt"
		#if server sends debugClose, debugsock will have been closed
		if { [string match debugClose $cmd] } { return }
		puts $debugsock $rslt
		puts $debugsock "*** Tcl stack TRACE ***"
		puts $debugsock $errorInfo
	    } else {
		if { [string match debugClose $cmd] } { return }
		puts $debugsock $rslt
	    }
	    puts $debugsock "END-RESULT"
	    # do we need flush for line buffering?
	    flush $debugsock
	}
    }

    proc debugClose {} {
	# user app can close the socket or the server will
	# request debugClose before it shutsdown
	global debugsock debugconnect
	if { ! $debugconnect } { return }
	catch {fileevent $debugsock readable {}}
	catch {close $debugsock}
    }

    # try the debugServer socket connection
    connectDebugServer localhost
}


############## open the pbscript pipe and setup fileevents

# set to 1 if you want a scrollbar on the group-selector
set gselSB 1

# open a rw pipe to the pbprog program. turn off prompting via -p
# The double quoting for open, works for Windows XP with
# the prog from command line in quotes allowing spaces.

if [catch {open "|\"$pbprog\" -p" r+} pipe] {
    Dialog_Warn "Couldn't open pbscript program: $pipe"
    exit
}


# global pipe io variables
set readfinished 0
set readlines {}
set errlines ""
set warnlines ""
set nread 0
set nerr 0
set nwarn 0

set greadfinished 0
set greadlines {}
set nread 0

set sellectlist ""

# make the pipe non-blocking, and when pipe is ready for reading call pipeRead
# NB the pipe is line buffered by default

fileevent $pipe readable [list pipeRead]
fconfigure $pipe -blocking 0

# startup plot window
set iwin 1

if { $GNU } {
    fileevent $gpipe readable [list gpipeRead]
    fconfigure $gpipe -blocking 0
}


############### arrow bitmaps

set downarrowdata "#define downarrow9x7_width 7
#define downarrow9x7_height 9
static unsigned char downarrow9x7_bits[] = \{
   0x00, 0x08, 0x08, 0x08, 0x49, 0x3e, 0x1c, 0x08, 0x00
};"


set uparrowdata "#define uparrow9x7_width 7
#define uparrow9x7_height 9
static unsigned char uparrow9x7_bits[] = \{
   0x00, 0x08, 0x1c, 0x3e, 0x49, 0x08, 0x08, 0x08, 0x00
};"

set rightarrowdata "#define rightarrow9x6_width 6
#define rightarrow9x6_height 9
static unsigned char rightarrow9x6_bits[] = \{
    0x02, 0x04, 0x0c, 0x18, 0x3f, 0x18, 0x0c, 0x04, 0x02,
};"

set leftarrowdata "#define leftarrow9x6_width 6
#define leftarrow9x6_height 9
static unsigned char leftarrow9x6_bits[] = \{
    0x10, 0x08, 0x0c, 0x06, 0x3f, 0x06, 0x0c, 0x08, 0x10,
};"


image create bitmap rightarrow -data $rightarrowdata
image create bitmap leftarrow  -data $leftarrowdata

image create bitmap uparrow -data $uparrowdata
image create bitmap downarrow  -data $downarrowdata

####################### pipe fileevent callbacks

set gnuplotFPT ""

# gnuplot pipe write for group plotting
proc gpipeWrite {cmd} {
    global gpipe
    global greadfinished
    global gnread greadlines
    global saveGnuCmnds gnuplotFPT
    if { $saveGnuCmnds && [string length $gnuplotFPT] > 0 } {
       puts $gnuplotFPT $cmd
    }
    set greadlines {}
    set gnread 0
    set greadfinished 0
    puts $gpipe $cmd
    flush $gpipe
}

# gnuplot pipe write for dataset plotting only created if dsmode enabled
proc gdspipeWrite {cmd} {
    global gdspipe
    global gdsreadfinished
    global gdsnread gdsreadlines
    global saveGnuCmnds gnuplotFPT
    if { $saveGnuCmnds && [string length $gnuplotFPT] > 0 } {
       puts $gnuplotFPT $cmd
    }
    set gdsreadlines {}
    set gdsnread 0
    set gdsreadfinished 0
    puts $gdspipe $cmd
    flush $gdspipe
}


# pbscript pipe write
proc pipeWrite {cmd} {
    global pipe
    global readfinished
    global nread nerr nwarn
    global readlines errlines warnlines
    global DEBUG debugFPT

    # clear the readlines list and the readfinished flag
    set readlines {}
    set errlines ""
    set warnlines ""
    set nread 0
    set nerr 0
    set nwarn 0
    set readfinished 0
    puts $pipe $cmd
    flush $pipe
    if { $DEBUG } {
	global outtxt debugFPT
	$outtxt insert end "$cmd\n"
	puts $debugFPT $cmd
	flush $debugFPT
    }
}


proc waitRead {maxms} {
    # could replace the vwait after pipeWrite with this proc; waitRead maxms
    # if problems with pipe death.
    global readfinished
    if {$maxms > 0} {
	after $maxms set readfinished -1
    }
    vwait readfinished
    if { $readfinished < 0 } {
	Dialog_Warn "read time limit exceded. pbscript pipe may have died"
    } elseif {$maxms > 0} {
	after cancel set readfinished -1
    }
}

# gnuplot pipeWrite callback for reading result for group plotting
proc gpipeRead {} {
    global gpipe DEBUG
    global greadfinished
    global greadlines gnread

    if { [eof $gpipe] } {
	puts "eof"
	catch {close $gpipe}
	saveSession 1
	Dialog_Warn "gnuplot pipe died. Session saved."
	exit
	return
    }
    # gets discards nl
    set nb [gets $gpipe line]
    #puts "got $nb bytes back from pipe"
    if { $nb < 0 } {
	# this should have been caught as eof above
	#puts $line
    } else {
	set line [string trim $line]
	# we are not waiting for any return from gnuplot cmnds
	# if there is an error write it to the debug window
	if { $DEBUG } {
	    global outtxt
	    $outtxt insert end "FROM gnuplot: $line\n"
	}
	if { [string length $line] > 0 } {
	    lappend greadlines $line
	    incr gnread
	}
    }
    # hope to run gnuplot without ever reading back anything
    set greadfinished 1
}

# gnuplot pipeWrite callback for reading result for dataset plotting
proc gdspipeRead {} {
    global gdspipe DEBUG
    global gdsreadfinished
    global gdsreadlines gdsnread

    if { [eof $gdspipe] } {
	puts "eof"
	catch {close $gdspipe}
	saveSession 1
	Dialog_Warn "gnuplot pipe died. Session saved."
	exit
	return
    }
    # gets discards nl
    set nb [gets $gdspipe line]
    #puts "got $nb bytes back from pipe"
    if { $nb < 0 } {
	# this should have been caught as eof above
	#puts $line
    } else {
	set line [string trim $line]
	# we are not waiting for any return from gnuplot cmnds
	# if there is an error write it to the debug window
	if { $DEBUG } {
	    global outtxt
	    $outtxt insert end "FROM gnuplot: $line\n"
	}
	if { [string length $line] > 0 } {
	    lappend gdsreadlines $line
	    incr gdsnread
	}
    }
    # hope to run gnuplot without ever reading back anything
    set gdsreadfinished 1
}

# pbscript callback for reading
proc pipeRead {} {
    # called via fileevent when pipe ready for reading a line
    global pipe
    global readfinished
    global readlines errlines warnlines
    global nread nerr nwarn
    global DEBUG

    if { [eof $pipe] } {
	puts "eof"
	catch {close $pipe}
	saveSession 1
	Dialog_Warn "pbscript pipe died. Session saved. Please restart."
	exit
	return
    }
    # gets discards nl
    set nb [gets $pipe line]
    #puts "got $nb bytes back from pipe"
    if { $nb < 0 } {
	# this should have been caught as eof above
	#puts $line
    } else {
	set line [string trim $line]
	if { $DEBUG } {
	    global outtxt
	    $outtxt insert end "$line\n"
	}
	# end of output=END is specific to pbscript
	# collect readlines until END
	if [string match $line "END"] {
	    set readfinished 1
	} elseif { [string length $line] > 0 } {
	    if { [string match ERROR* $line] } {
		append errlines "$line\n"
		incr nerr
	    } elseif { [string match WARN* $line] } {
		append warnlines "$line\n"
		incr nwarn
	    } else {
		lappend readlines $line
		incr nread
	    }
	}
    }
}

###################################################################
# start building the main GUI window

wm title . "PolarizedBeam Corection for 3-axis    $currentDir"
wm iconname . "PB cor"
wm geometry . +10+10

# main frame with a left(files/options) and right(group/read) subframe
set a [frame .a -height 500 -width 500]
grid $a
set left [frame $a.left -bd 2 -relief raised]
set right [frame $a.right -bd 2 -relief raised]
grid $left $right -padx 5
set top [frame $left.top -bd 2 -relief flat -bg tan]
set mid [frame $left.mid -bd 2 -relief flat]
set bot [frame $left.bot -bd 2 -relief flat -bg tan]
grid $top -sticky new
grid $mid -sticky ew
grid $bot -sticky ews
set heightlines 1
array set bframes [list 1 $top 2 $mid 3 $bot]

# make 3 frames for left: top mid bot.
# these are used for combinations of
# fileselector, prep-solve and auto-prep (f/g/a)
# depending on the users height request for the application
# through heightlines variable.
# heightlines = 1 10L only only one of top/mid/bot is viewed and selectable
# heigthlines = 2 26L use  two of three frames
# heightlines = 3 36L use  all 3 left frames for all 3 (f/g/a)

# storage for group info
set ngroup 0
array set group {}
# the libPB flipper-state representation is 1=offoff 2=onon  3=onoff 4=offon
# in pbgui we use the ICP-ICE    1=offoff 2=onoff 3=offon 4=onon
# fss are the labels in the pbgui constraint matrix
# fs is used as a vector variable for BLT so it will create it as array
# so use fss instead of fs to store the standard flipper-state order
array set fss {1 offoff 2 ONoff 3 offON 4 ONON}
array set cm {}

# monochormator+filter options
set monoOptions [list "PGbt7 no filter" "PGbt7 2cm filter"]
# beam monitor options
set monitorOptions [list "none" "before He3 pol" "after He3 pol"]
# init option values
set monofilter [lindex $monoOptions 0]
set monitor [lindex $monitorOptions 2]

set acszero 0
set cellFiles {}
# standard constraint labels
set nsflbls {none nsf=0 ONON=offoff offoff=ONON}
set sflbls {none sf=0 offON=ONoff ONoff=offON}
set nsfcons [lindex $nsflbls 0]
set sfcons [lindex $sflbls 0]
set npts 0
set nmsr 0
set calc 0
set isfop 0
set textf(1) ""
set textf(2) ""

proc buildGroupPrep {} {
    # this goes in mid

    global cellFiles
    global monoOptions monitorOptions
    global BLT PLT

    global group gnum
    global fss
    
    global nsf sf
    global conenw conlbw conmmw

    # these are the setup variables/arrays
    global cellFile monofilter monitor
    global nsflbls sflbls
    global nsfcons sfcons
    #global nsfv sfv
    global acszero
    global cm cen

    # store the setup for each group in group(i)
    global bframes

    set f $bframes(2)
    clearGrid $f

    set exe $f.exec
    set optac $f.ac
    set con $f.cm

    set solvact $exe.solvact
    set solvsel $exe.solvsel
    set solvuns $exe.solvuns
    set textsln $exe.textsln
    set plotsln $exe.plotsln

    set cellmen $exe.cells
    set monocmen $exe.monoc
    set monitmen $exe.monit
    set nsfmen $exe.nsf
    set sfmen $exe.sf

    set cellmenlbl $exe.cellslbl
    set monocmenlbl $exe.monoclbl
    set monitmenlbl $exe.monitlbl
    set nsfmenlbl $exe.nsflbl
    set sfmenlbl $exe.sflbl

    set acszerow $optac.acszero

    set conen $con.en
    set conlb $con.lb
    set conmm $con.mm
    set i 0
    set conenw($i) $conen$i
    for {set i 1} {$i <= 4} {incr i} {
	set conenw($i) $conen$i
	set conlbw($i) $conlb$i
	for {set j 1} {$j <= 4} {incr j} {
	    set conmmw($i$j) $conmm$i$j
	}
    }

    if { ! [winfo exists $exe] } {

	frame $exe
	frame $optac
	frame $con

	button $solvact -text "SOLVEcur" -width 9 \
	    -command "solveGroup 0" \
	    -background green -activebackground yellow
	button $solvsel -text "SOLVEsel" -width 9 \
	    -command "solveGroup 1" \
	    -background green -activebackground yellow
	button $solvuns -text "SOLVEunsel" -width 9 \
	    -command "solveGroup 2" \
	    -background green -activebackground yellow
	button $textsln -text "TEXTsoln" -width 9 \
	    -command "showText" \
	    -background green -activebackground yellow
	button $plotsln -text "PLOTsoln" -width 9 \
	    -command "showPlot" \
	    -background green -activebackground yellow
	if { $PLT < 1 } { $plotsln configure -state disabled }

	label $cellmenlbl -text cellFile
	label $monocmenlbl -text filter
	label $monitmenlbl -text monitor
	eval tk_optionMenu $cellmen cellFile $cellFiles
	eval tk_optionMenu $monocmen monofilter $monoOptions
	eval tk_optionMenu $monitmen monitor $monitorOptions

	$cellmen configure -width 11
	$monocmen configure -width 11
	$monitmen configure -width 11

	label $nsfmenlbl -text nsfCons
	label $sfmenlbl -text sfCons
	eval tk_optionMenu $nsfmen nsfcons $nsflbls
	eval tk_optionMenu $sfmen sfcons $sflbls

	$nsfmen configure -width 11
	$sfmen configure -width 11

	trace variable nsfcons w setNSF
	trace variable sfcons w setSF

	checkbutton $acszerow -variable acszero \
	    -text "CS=0 if no data"

	set i 0
	label $conenw($i) -text enable
	for {set i 1} {$i <= 4} {incr i} {
	    checkbutton $conenw($i) -variable cen($i) \
			     -text $fss($i)=
	    label $conlbw($i) -text $fss($i)
	    for {set j 1} {$j <= 4} {incr j} {
		entry $conmmw($i$j) -textvariable cm($i$j) \
			 -width 5 -relief sunken
		set cm($i$j) 0
		if { $i == $j } { $conmmw($i$j) configure -state disabled }
	    }
	}	
    }

    grid $exe -sticky w
    grid $optac -sticky w
    grid $con -sticky w

    grid $acszerow -sticky w

    grid $cellmenlbl $cellmen -sticky w
    grid $solvact -sticky e -padx 10 -row 0 -column 2
    grid $monocmenlbl $monocmen -sticky w
    grid $solvsel -sticky e -padx 10 -row 1 -column 2
    grid $monitmenlbl $monitmen -sticky w
    grid $solvuns -sticky e -padx 10 -row 2 -column 2
    grid $nsfmenlbl $nsfmen -sticky w
    grid $textsln -sticky e -padx 10 -row 3 -column 2
    grid $sfmenlbl $sfmen -sticky w
    grid $plotsln -sticky e -padx 10 -row 4 -column 2

    grid $conenw(0) -row 0 -column 0
    for {set i 1} {$i <= 4} {incr i} {
	grid $conenw($i) -row $i -column 0 
	grid $conlbw($i) -row 0 -column $i
	for {set j 1} {$j <= 4} {incr j} {
	    grid $conmmw($i$j) -row $i -column $j
	}
    }
}

# procs for prep-solve

proc setCSrow {row args} {
    global conenw conmmw
    if {[llength $args] < 4} { return }
    for {set i 1; set j 0} {$j < 4} {incr i; incr j} {
	$conmmw($row$i) delete 0 end
	$conmmw($row$i) insert 0 [lindex $args $j]
    }
    $conenw($row) select
}
proc setNSF {varName indx op} {
    #upvar $varName var
    #global nsf nsfv
    global nsfcons
    global conenw
    global nsflbls

    if { [set i [lsearch -exact $nsflbls $nsfcons]] < 0 } { return }

    if { $i == 0 } {
	$conenw(1) deselect
	$conenw(4) deselect
    } elseif { $i == 1 } {
	setCSrow 1 0 0 0 0
	setCSrow 4 0 0 0 0
    } elseif { $i == 2 } {
	setCSrow 4 1 0 0 0
	$conenw(1) deselect
    } elseif { $i == 3 } {
	setCSrow 1 1 0 0 0
	$conenw(4) deselect
    }
}

proc setSF {varName indx op} {
    global conenw
    global sflbls
    global sfcons
    if { [set i [lsearch -exact $sflbls $sfcons]] < 0 } { return }

    if { $i == 0 } {
	$conenw(2) deselect
	$conenw(3) deselect
    } elseif { $i == 1 } {
	setCSrow 2 0 0 0 0
	setCSrow 3 0 0 0 0
    } elseif { $i == 2 } {
	setCSrow 3 0 1 0 0
	$conenw(2) deselect
    } elseif { $i == 3 } {
	setCSrow 2 0 0 1 0
	$conenw(3) deselect
    }
}

array set stndTolArr {}
set notallowed {y e c ce scl fs ts Mon Time S* pt}
set stndTolNames {}
proc makeTolArray {} {
    #grab all the standard column tolerances into array
    global stndTolArr stndTolNames notallowed 
    global readfinished readlines
    pipeWrite "c S"
    vwait readfinished
    for {set i 1} {$i < [llength $readlines]} {incr i} {
	set ln [lindex $readlines $i]
	if { [llength $ln] < 3 } { continue }
	set stndName [lindex $ln 1]
	set ok 1
	foreach na $notallowed {
	    if { [string match $na $stndName] } { set ok 0 ; break }
	}
	if { $ok < 1 } { continue }
	set stndTol [lindex $ln 2]
	set stndTolArr($stndName) $stndTol
	lappend stndTolNames $stndName
    }
}

array set autoStndTol {}
proc setStndTolerances {} {
    global stndTolNames stndTolArr
    global readfinished
    # set the stndTol in libDF from the users stndTolArr
    foreach snam $stndTolNames {
	pipeWrite "c St $snam $stndTolArr($snam)"
	vwait readfinished
    }
}

proc setStndTol {name} {
    # binding to choosing a tol col type from menu
    global stndTolArr tolerance
    if { [info exists stndTolArr($name)] } {
	set tolerance $stndTolArr($name)
    }
}
proc getStndTol {} {
    # binding to Leave tolerance entry
    global stndTolArr tolerance toltype
    set stndTolArr($toltype) $tolerance
}

proc nextEmptyGroup {ig} {
    global group
    if { $ig < 0 } { set ig 0 }
    while {1} {
	incr ig
	if { ! [info exists group($ig)] } { return $ig }
	set fls [getFilesInGroup $ig]
	if { [llength $fls] < 1 } { return $ig }
    }
}
proc getEntriesForGroup {ig} {
    # return list of entries for group ig
    global lbxs lbxi
    set igs [$lbxs(1) get 0 end]
    set lbxi {}
    for {set i 0} {$i < [llength $igs]} {incr i} {
	set igi [lindex $igs $i]
	if { $igi == $ig } { lappend lbxi $i }
    }
    return $lbxi
}

# lbxi listbox indices when selecting or creating file lists
set lbxi {}
proc getFilesInGroup {ig} {
    global lbxs lbxi
    global group

    if { $ig > 0 && ! [info exists group($ig)] } { return {} }
    # if ig < 1 return all files
    # else just files with ig group number
    set fls [$lbxs(0) get 0 end]
    if { $ig < 1 } { return $fls }

    set lbxi [getEntriesForGroup $ig]
    set ret {}
    for {set i 0} {$i < [llength $lbxi]} {incr i} {
	set igi [lindex $lbxi $i]
	lappend ret [lindex $fls $igi]
    }
    return $ret
}

proc appendLbxData {data} {
    global lbxs nlbxs lbxlabels lbxfmts
    for {set i 0} {$i < $nlbxs} {incr i} {
	set idata [lindex $data $i]
	for {set j 0} {$j < [llength $idata]} {incr j} {
	    set jdata [lindex $idata $j]
	    if { [string length $lbxfmts($i)] > 0 } {
		for {set k 0} {$k < [llength $jdata]} {incr k} {
		    set kdata [format $lbxfmts($i) [lindex $jdata $k]]
		    set jdata [lreplace $jdata $k $k $kdata]
		}
	    }
	    $lbxs($i) insert end $jdata
	}
	setLbxWidth $lbxs($i) [lindex $lbxlabels $i]
    }
}

proc getLbxFiles {} {
    # get short filenames from lbxs(0)
    global lbxs lbxfiles
    set lbxfiles [$lbxs(0) get 0 end]
}
set lbxfiles {}
proc getLbxData {fil} {
    # data should be grp ind except for grp num
    global lbxs nlbxs lbxfiles
    #set fils [$lbxs(0) get 0 end]
    set data {}
    if { [set i [lsearch -exact $lbxfiles $fil]] < 0 } { return $data }
    lappend data [lindex $lbxfiles $i]
    for {set j 1} {$j < $nlbxs} {incr j} {
	lappend data [$lbxs($j) get $i]
    }
    return $data
}

# flipper-state translations gui -> libPB and libPB -> gui
array set fst {1 1 2 3 3 4 4 2}
array set fsg {1 1 2 4 3 2 4 3}
proc writeGroupOpts {} {
    # write the current group options to pbscript
    # to setup for solving
    global gnum
    global cm cen acszero monitor monofilter 
    global cellFile
    global monoOptions monitorOptions
    global fss fst fsg
    global errlines nerr readfinished

    pipeWrite "p h $cellFile"
    vwait readfinished

    # all cellFile errors are fatal
    if { $nerr > 0 } {
	Dialog_Warn "Failed cellFile for group $gnum:\n$errlines"
	return 0
    }


    # check independent S (non-zero coef in equ) are not used in cnstrnt eq
    for {set i 1} {$i <= 4} {incr i} {
	if { ! $cen($i) } { continue }
	for {set j 1} {$j <= 4} {incr j} {
	    if { $i == $j } { continue }
	    if { [scan $cm($i$j) {%g} cmval] < 1 } {
		Dialog_Warn "In group $gnum\nConstraint equ for $fss($i) has non-numeric entry for $fss($j).\nPlease fix this before solving group."
	    }
	    if { [expr abs($cm($i$j))] != 0 } {
		if { $cen($j) } {
		    Dialog_Warn "In group $gnum\nConstraint equ for $fss($i) uses ind var $fss($j)\nwhich is then constrained.\nPlease fix this before solving group."
		    return 0
		}
	    }
	}
    }

    #   gui1-4 libPB1-4 fs index translation
    # offoff 1 -> 1
    # ONoff  2 -> 3
    # offON  3 -> 4
    # ONON   4 -> 2

    #   libPB1-4  gui1-4
    # offoff 1 -> 1
    # ONON   2 -> 4
    # ONoff  3 -> 2
    # offON  4 -> 3

    # p f c e i offoff ONON ONoff offON
    # p f c e i i1     i4   i2    i3
    #set fst(1) 1
    #set fst(2) 3
    #set fst(3) 4
    #set fst(4) 2

    for {set i 1} {$i <= 4} {incr i} {
	set fi $fst($i)
	set cmvs ""
	for {set j 1} {$j <= 4} {incr j} {
	    set fj $fsg($j)
	    append cmvs " $cm($i$fj)"
	}
	if { $cen($i) } {
	    pipeWrite "p f c e $fi $cmvs"
	    vwait readfinished
	} else {
	    pipeWrite "p f c e -$fi"
	    vwait readfinished
	}
    }

    set imon [lsearch -exact $monitorOptions $monitor]
    set iflt [lsearch -exact $monoOptions $monofilter]
    if { $imon < 0 } { set imon 0 }
    if { $iflt < 0 } { set iflt 0 }

    pipeWrite "p f m c $imon"
    vwait readfinished
    pipeWrite "p f m t $iflt"
    vwait readfinished

    pipeWrite "p f c a $acszero"
    vwait readfinished
    return 1
}

array set fopgroup {}

proc solveGroup {is} {
    # is==0 solve the current group
    # is==1 solve selected
    # is==2 solve unselected
    global readfinished readlines
    global nread nerr errlines nwarn warnlines
    global gnum calcgroup
    global cellFile
    global npts nmsr calc isfop fopgroup
    global textf
    global BLT PLT

    set grps {}
    if { $is == 0 } {
	set grps [list $gnum]
    } elseif { $is == 1 } {
	# all grps selected
	set grps [groupsSelected 1]
    } elseif { $is == 2 } {
	# all grps un-selected
	set grps [groupsSelected 0]
    }
    if { [llength $grps] < 1 } {
	Dialog_Warn "No Groups to Solve"
    }
    # save current group
    saveGroup
    set emptygroups ""
    foreach grp $grps {

	# make grp the current pbgroup
	pipeWrite "g$grp"
	vwait readfinished

	setGnum $grp

	set textf(1) ""
	set textf(2) ""
	set calc 0
	set calcgroup($grp) 0
	set fopgroup($grp) $isfop
	if { $isfop } { continue }
	hiliteCalcGroup

	# apply all constraints and options to this group
	# set the cellfile and then solve
	
	# p f m imonitor imono
	# imonitor 0 none 1 before-pol 2 after-pol
	# imono    0 BT7PG-no-filter 1 w/2cmPGfilter  2 BT4PGnofilter
	# p f c a autoconstrain
	# p f c e index 4coefs
	# p f c e -index disable

	# send group options to pbscript
	if { ![writeGroupOpts] } { continue }

	# send the stndTolArr for this group
	setStndTolerances
	# get the files for this group listed in lbxs(0)
	set fls [getFilesInGroup $grp]
	if { [llength $fls] < 1 } {
	    append emptygroups " $grp"
	    continue
	}

	fileListToSetList $fls
	# move this setlist to pbgroup $gnum
	# would like to keep correspondence between pbgroups and gnums

	# move setlist 1 to this group. This clears any previous list
	pipeWrite "g sl 1"
	vwait readfinished
	
	# prep the group
	# prep should probably only fail on memory problems
	# as it does little error checking
	# most of the error checking is done when solution requested (p cmd)
	pipeWrite "g p"
	vwait readfinished
	if { $nerr > 0 } {
	    Dialog_Warn "group $grp prep ERRORS:\n$errlines"
	    continue
	}
	if { $nwarn > 0 } {
	    Dialog_Warn "group $grp prep WARNINGS:\n$warnlines"
	}

	# after prep, check group basics to see that nmsr and pts>0
	# NB g l prints all groups.
	pipeWrite "g$grp"
	vwait readfinished
	set nscan \
	    [scan [lindex $readlines 1] {%[*]%d %d %d %d %d} \
		 star igrp calc nset npts nmsr]
	if { $nscan < 6 || $npts < 1 || $nmsr < 1 } {
	    Dialog_Warn "group $grp solution failed to get npts>0 nmsr>0"
	    continue
	}

	# set the group output control flags for extended output level 1
	# this is for the .m .p datasets and files
	# the other defaults are:
	# f=1 flags to header
	# h=0 no data headers output
	# c=# comment char
	# x=1 x-col sel output
	# m=0 no extended output, which we are changing here
	# e.g. this makes srcfile last column of .m output
	pipeWrite "g oc m 1"
	vwait readfinished

	# solve
	pipeWrite "p"
	vwait readfinished
	if { $nerr > 0 } {
	    Dialog_Warn "group $grp calc ERRORS:\n$errlines"
	    continue
	}
	if { $nwarn > 0 } {
	    Dialog_Warn "group $grp prep WARNINGS:\n$warnlines"
	}
	# readlines only line should have
	# group # correctOK ic datasetsOK id
	if { $nread < 1 } {
	    Dialog_Warn "group $grp calc didnt report success"
	    continue
	}
	set nscan \
	    [scan [lindex $readlines 0] {%s %d %s %d %s %d} \
		 glbl igrp clbl cok slbl sok]
	if { $nscan < 4 || ! $cok } {
	    Dialog_Warn "group $grp calc failed"
	    continue
	}
	if { $nscan < 6 || ! $sok } {
	    Dialog_Warn "group $grp calc OK but didnt produce result datasets"
	    continue
	}
	
	pipeWrite "g$grp"
	vwait readfinished
	set nscan \
	    [scan [lindex $readlines 1] {%[*]%d %d %d %d %d} \
		 star igrp calc nset npts nmsr]
	if { $nscan < 6 || $npts < 1 || $nmsr < 1 } {
	    append errs "group $grp ERROR: could not get nmsr/npts after calc\n"
	    continue
	}

	# write the group results
	pipeWrite "g w"
	vwait readfinished

	if { $nerr > 0 } {
	    Dialog_Warn "group $grp write ERRORS:\n$errlines"
	    continue
	}
	if { $nwarn > 0 } {
	    Dialog_Warn "group $grp write WARNINGS:\n$warnlines"
	}


	# should have made group datasets .p and .m
	# check the group root-file-name to find out what files will be written
	pipeWrite "g f"
	vwait readfinished

	set nscan \
	    [scan [lindex $readlines 0] {%s %d %s %s} \
		 glbl igrp rlbl root]
	if { $nscan < 4 || [string match NA $root] } {
	    append errs "group $grp ERROR: couldn't get root filename\n"
	    continue
	}

	# check that root.p and root.m exist
	set textf(1) $root
	append textf(1) .p
	set textf(2) $root
	append textf(2) .m

	foreach ix [array names textf] {
	    set fil $textf($ix)
	    if { ![istextfile $fil] } {
		Dialog_Warn "group $grp ERROR: failed to write $fil"
		set textf($ix) ""
	    }
	}
	# save the textfile names
	set calc 1
	set calcgroup($grp) 1
	set isfop 0
	set fopgroup($grp) 0
	# save the src and .m and .p file names
	groupSrcFiles
	# save this group
	saveGroup
	# hilite calc group number
	hiliteCalcGroup
	# nolonger need trace on gnum for groupSrcFiles

	if { $PLT } {
	    groupPlotSetup
	    savePlot
	    showPlot
	}
    }
    if { [string length $emptygroups] > 0 } {
	Dialog_Warn "Empty groups: $emptygroups"
    }
}

proc subtractBG {} {
    global bg
    if { ![winfo exists $bg] } { buildFBG }
    if { ![winfo ismapped $bg] } { wm deiconify $bg }
    raise $bg
}

# fast background subtraction is a separate toplevel window
set bg ""
proc buildFBG {} {
    global bg
    global FBGcps FBGcpsErr

    set bg [toplevel .bg]
    wm title $bg "fast background subtraction"
    wm iconname $bg "PB bg"

    # entry FBGcps  FBGcpsErr
    set bge $bg.e
    set bgcpsL $bge.cpsL
    set bgcps $bge.cps
    set bgerrL $bge.errL
    set bgerr $bge.err

    set bgb $bg.b
    set bgcurg $bgb.curg
    set bgselg $bgb.selg
    set bgself $bgb.self

    frame $bge
    frame $bgb
    grid $bge
    grid $bgb

    label $bgcpsL -text "fast bg (cps)"
    entry $bgcps -textvariable FBGcps -relief sunken -bd 2
    label $bgerrL -text "+-"
    entry $bgerr -textvariable FBGcpsErr -relief sunken -bd 2

    button $bgself -text "sub bg sel files" -command "subtractFBG 0" \
	    -background green -activebackground yellow
    button $bgcurg -text "sub bg cur group" -command "subtractFBG 1" \
	    -background green -activebackground yellow
    button $bgselg -text "sub bg sel groups" -command "subtractFBG 2" \
	    -background green -activebackground yellow

    grid $bgcpsL $bgcps $bgerrL $bgerr
    grid $bgself $bgcurg $bgselg -sticky ew
}

set FBGcps 0
set FBGcpsErr 0
proc subtractFBG {in} {
    global lbxs lbxlabels lbxi fbglbx gnum
    global dsFile srcfiles
    global readlines readfinished
    global fbgFile FBGcps FBGcpsErr

    # subtract fast background from
    # selected files in=0  lbxs(0) or lbxs(1)
    # current group  in=1
    # selected groupsin=2  lbxs(0) or lbxs(1)
    # FBG is here just proportional to Time
    # although eventually may add proportional to angle ttm tts tta ?
    # the FBG toplevel lets user enter FBGcps +- cps and
    # then command buttons to subtract selfiles or groups

    # use libDF clone each source file, returns the dsnumber
    # d nm to set its new name= src.fbg.x
    # put src and bg file into setlist via fileListToSetlist
    # make sure we have the fop for this setlist (1) with columns -1 -2 -3 2
    # Qx Qy Qz En via: s1 o c 2 -1 -2 -3
    # which should distinguish the datapoints for subtraction
    # purposes
    # prepare the data for the background set using
    # c df commands: c dfn 1 to init to one ind var
    # c dfx Time to make function value proportional to Time
    # c dfc1 FBGcps   c dfe1 FBGcpsErr
    # c df DetectorCol,DetectorErrCol
    # set the subtractor to the FBGsetnumber for this setlist
    # s ss FBGsetnumber
    # do the subtraction
    # s s
    # output should be src.s
    # can also find the result set from s1 l
    # which reports in second line op lnr opc lns
    # lnr is setlist link to result
    # that setlist should have the result file

    if { $in == 1 } {
	set fls [getFilesInGroup $gnum]
    } elseif { $in == 0 } {
	set fls [filesSelected 1]
    } elseif { $in == 2 } {
	set fls [filesInGroupsSelected 1]
    }
    if { [llength $fls] < 1 } {
	Dialog_Warn "No files to subtract FBG from.\n"
	return
    }

    for {set i 0} {$i < [llength $fls]} {incr i} {
	# clone the source file
	set fil [lindex $fls $i]
	if { [string match *.s $fil] || [string match *.fbg* $fil] } {
	    if { ![Dialog_Query "File $fil looks like it has already been used in subtraction. Continue anyway?"] } { continue }
	}
	set sset $dsFile($fil)
	pipeWrite "d$sset c"
	vwait readfinished
	# get the cloned setnumber
	if { [scan [lindex $readlines 0] {%s %d} lbl bset] < 2 } {
	    Dialog_Warn "failed to get cloned src file setnumber.\n"
	    return
	}
	# name the clone src.fbg
	pipeWrite "d$bset nm $fil\.fbg.x"
	vwait readfinished
	# set the clones y1 and e1 using the FBGcps and FBGcpsErr
	pipeWrite "c dfn 1"
	vwait readfinished
	pipeWrite "c dfx Time"
	vwait readfinished
	pipeWrite "c dfc1 $FBGcps"
	vwait readfinished
	pipeWrite "c dfe1 $FBGcpsErr"
	vwait readfinished
	# find the y1 and e1 column numbers
	pipeWrite "c al y1"
	vwait readfinished
	if { [scan [lindex $readlines 1] {%s %d %d %s} typ inx ycol lbl] < 4 } {
	    Dialog_Warn "failed to get cloned src file y1 col"
	    return
	}
	pipeWrite "c al e1"
	vwait readfinished
	if { [scan [lindex $readlines 1] {%s %d %d %s} typ inx ecol lbl] < 4 } {
	    Dialog_Warn "failed to get cloned src file e1 col"
	    return
	}
	# calc and set the FSB and its error
	pipeWrite "c df $ycol,$ecol"
	vwait readfinished
	# ready to put these two sets into setlist 1 for the subtraction
	pipeWrite "s1 lc"
	vwait readfinished
	pipeWrite "s l+ $sset $bset"
	vwait readfinished
	# define the subtractor
	pipeWrite "s ss $bset"
	vwait readfinished
	# make sure the fop has Qx Qy Qe E columns to compare
	pipeWrite "s o c 2 -1 -2 -3"
	vwait readfinished
	# do the subtraction
	pipeWrite "s s"
	vwait readfinished
	# list the setlist to find out which setlist it wrote the result to
	pipeWrite "s1 l"
	vwait readfinished
	# the second line of output should have: op lnr opc lns, lnr is link
	# to result
	if { [scan [lindex $readlines 1] {%s %s %s %s %d %d %d %d} \
		  lop llnr lopc llns op lnr opc lns] < 8 } {
	    Dialog_Warn "failed to get subtraction result setlist.\n"
	    return
	}
	# list the result setlist to find the result dataset number
	pipeWrite "s$lnr l"
	vwait readfinished
	# the 4th line should have: ilst match iset filename
	if { [scan [lindex $readlines 3] {%d %d %d %s} \
		  ilst mtch rset rfil] < 4 } {
	    Dialog_Warn "failed to get subtraction result dataset.\n"
	    return
	}
	# certainly want to write the result to disk and then
	# replace the src in lbxs(0) with the bgsub file
	# should we also write the background file, I guess
	pipeWrite "d$bset w"
	vwait readfinished
	pipeWrite "d$rset w"
	vwait readfinished
	# can get the lbx0 index for the processed file from lbxi
	# which was set when the fls was constructed.
	set indx [lindex $lbxi $i]
	$lbxs(0) delete $indx
	set fshort [file tail $rfil]
	$lbxs(0) insert $indx $fshort
	set dsFile($fshort) $rset
	set srcfiles($rset) $fshort
	set fbgFile($fshort) "$FBGcps $FBGcpsErr"
	$lbxs($fbglbx) delete $indx
	$lbxs($fbglbx) insert $indx $fbgFile($fshort)
    }
    setLbxWidth $lbxs(0) [lindex $lbxlabels 0]
    setLbxWidth $lbxs($fbglbx) [lindex $lbxlabels $fbglbx]
    hiliteActiveGroup
}

########### text file viewer procs ###################

proc istextfile {fil} {
    return [expr [file exists $fil] && [file isfile $fil] && \
	    [file readable $fil]]
}

set tv ""
set txtw ""
set mp 1

proc dsFromSn {n} {
    # n setnmuber 1 - ngrpsets
    # returns pbscript setnumber
    global dsns ngrpsets
    if { $n < 1 || $n > $ngrpsets } { return 0 }
    return [lindex $dsns [expr $n - 1]]
}
proc snFromDs {ds} {
    # ds pbscript datasetnumber >= 1
    # returns setnumber in group list 1 - ngrpsets
    global dsns ngrpsets
    for {set i 0; set j 1} {$i < $ngrpsets} {incr i; incr j} {
	set dsni [lindex $dsns $i]
	if { $dsni == $n } { return $j }
    }
    return 0
}
proc fileFromSn {sn} {
    global srcfiles ngrpsets
    if {$sn < 1 } {set sn 1}
    if {$sn > $ngrpsets} {set sn $ngrpsets}
    set dsni [dsFromSn $sn]
    return $srcfiles($dsni)
}
proc snFromFile {fnam} {
    global solnfiles
    if { [set i [lsearch $solnfiles $fnam]] >= 0 } { return [expr $i + 1] }
    #for {set i 0; set j 1} {$i < $ngrpsets} {incr i; incr j}
    #	set dsni [lindex $dsns $i]
    #	if { [string match $fnam $srcfiles($dsni)] } { return $j }
    return 0
}

set viewfile ""
# put file text in toplevel text-viewer tv
proc showText {args} {
    global tv txtw tvtitle viewfile
    global textf mp snum
    global srcfiles solnfiles dsns ngrpsets
    global sellectlist lbxs currentDir
    # deiconify the scrolledText window
    # and load the .p file
    # srcfiles dsns ngrpsets should already be available from grpSrcFiles
    if { ! [winfo exists $tv] } {
	buildShowText -width 60 -height 20 -wrap none
    }
    wm deiconify $tv
    raise $tv
    $txtw delete 1.0 end
    if { $mp < 3 } {
	set viewfile $textf($mp)
    } elseif { $mp == 3 } {
	# print the group sets, will incl p m
	# groupSrcFiles will set array srcfiles(setnumber) = filename
	if { $ngrpsets < 1 } {
	    Dialog_Warn "NO src files listed for current group"
	    return
	}
	if { $snum < 1 } { set snum 1 }	
	if { $snum > $ngrpsets } { set snum $ngrpsets }
	set viewfile [file tail [lindex $solnfiles [expr $snum - 1]]]
    } else {
	# sel can be in sellectlist or lbxs(0)
	set sel [selectionOwned 1 1 1]
	if { [llength $sel] < 1 } { return }
	set viewfile [lindex $sel 0]

	#set fil [file join $currentDir $fil]
    }
    if { ! [istextfile $viewfile] } {
	#Dialog_Warn "$viewfile is not a readable textfile"
	return
    }
    set in [open $viewfile]
    $txtw insert end [read $in]
    close $in
    wm title $tv "$tvtitle   $viewfile"
}



array set srcfiles {}
set ngrpsets 0
set dsns {}
array set isrcplot {}
set pset 0
set mset 0
set solnfiles {}

proc groupSrcFiles {args} {
    # only called from solveGroup
    # this is not session robust as it was
    # designed for call after solveGroup, but
    # so have session recovery resolve any groups with calc=1.

    # what is needed is the list of solnfiles
    # solnfiles should be identical to the groupL grpF

    # so retask srcfiles as simply the inverse of dsFile
    # not to be session saved,
    # instead save solnfiles foreach group
    # along with ngrpsets which counts solnfiles
    # also get mset pset from textf
    # after textf is session recovered,
    # regenerate mset pset at the end.
    # isrcplot can be saved per group.

    # call this whenever group is solved
    # then save on gnum change via setGroup/saveGroup
    # then srcfiles, pset, mset are available
    # isrcplot($sn) flags srcdata viewed
    # init to all on

    global gnum
    global srcfiles solnfiles ngrpsets dsns isrcplot
    global pset mset
    global nread readlines
    # applies to current group

    set solnfiles {}
    set pset 0
    set mset 0

    pipeWrite "g$gnum"
    vwait readfinished

    pipeWrite "g ls"
    vwait readfinished
    # DATASETS for group #
    # gset iset npts ncol srcfile
    # NB this wont recognize the m and p sets (not integer sn)
    # which are handled separately in showText

    for {set i 2} {$i < $nread} {incr i} {
	set rl [lindex $readlines $i]
	set nscan [scan $rl {%d %d %d %d %s} igs ids npt nc srcf]
	if { $nscan < 5 } {
	    if { [scan $rl {%s %d} igs ids] == 2 } {
		if { [string match p $igs] } {
		    set pset $ids
		} elseif { [string match m $igs] } {
		    set mset $ids
		}
	    }
	    continue
	}
	lappend solnfiles [file tail $srcf]
    }
    #set dsns [array names srcfiles]
    set ngrpsets [llength $solnfiles]
    for {set i 1} {$i <= $ngrpsets} {incr i} { set isrcplot($i) 1 }
    setSn 1
    updateSrcFileMenus
}

proc setSn {sn} {
    global snum mp isrcon isrcplot srcfile srcfiles
    set snum $sn
    if { $mp == 3 } { showText }
    if { [info exists isrcplot($snum)] } {
	set isrcon $isrcplot($snum)
    }
    set ids [dsFromSn $snum]
    if { [info exists srcfiles($ids)] } {
	set srcfile [file tail $srcfiles($ids)]
    }
}
proc incrSn {} {
    global snum ngrpsets
    if { $snum >= $ngrpsets } { return }
    incr snum
    setSn $snum
}
proc decrSn {} {
    global snum
    if { $snum <= 1 } { return }
    incr snum -1
    setSn $snum
}

array set viewcollabels {}
set viewncol 0
set saveset 0
proc saveColumnsFromFileView {} {
    global viewfile viewcollabels viewncol nread readlines
    global saveset readfinished
    # check that viewfile belongs to a dataset
    pipeWrite "d ne $viewfile"
    vwait readfinished
    # returns: dataset # filename= filename  #=0 if not found
    set rl [lindex $readlines 0]
    if { [scan $rl {%s %d} dum iset] < 1 || $iset < 1 } {
	# try to read the file
	set ret [readFile $viewfile]
	if { [llength $ret] < 1 } {
	    Dialog_Warn "Current file in viewer is not in a dataset and read failed."
	    return
	}
	set iset [lindex $ret 0]
    }

    # open a toplevel that lets user select columns from this dataset
    # make it the current dataset and get the list of columnlabels
    # for the user to select from
    pipeWrite "d$iset"
    vwait readfinished
    set rl [lindex $readlines 1]
    #puts "saveColumnsFromFileView: d$iset returned nread=$nread\n[lindex $readlines 0]\n$rl"
    set nscan [scan $rl {%[*]%d %d %d %d %d %d %s} \
		   star iset npts ncol ny nas df src]
    if { $nscan < 8 } {
	Dialog_Warn "failed to get dataset info for set $iset.\n"
	return
    }
    pipeWrite "c l"
    vwait readfinished
    if { $nread <= $ncol } {
	Dialog_Warn "failed to return all column labels for dataset $iset.\n"
	return
    }
    # icol columnlabel followed by data lines
    arrayClear viewcollabels
    for {set i 1} {$i <= $ncol} {incr i} {
	set rl [lindex $readlines $i]
	set nscan [scan $rl {%d %s} icol labl]
	if { $nscan < 2 } {
	    Dialog_Warn "failed to read column label $i"
	    return
	}
	set viewcollabels($icol) $labl
    }
    set viewncol $ncol

    # make sure this toplevel is built, mapped and raised
    # and has menu of current columns list
    set saveset $iset
    openSaveColumns
}

proc updateSaveColumns {} {
    global lbsrc lbout viewfile viewsave
    global viewncol viewcollabels
    set lbh $viewncol
    if { $lbh > 16 } { set lbh 16 }
    $lbsrc configure -height $lbh
    $lbout configure -height $lbh
    $lbsrc delete 0 end
    $lbout delete 0 end
    for {set i 1} {$i <= $viewncol} {incr i} {
	$lbsrc insert end $viewcollabels($i)
    }
    set viewfile [file tail $viewfile]
    set viewsave $viewfile
}
proc toOut {} {
    global lbsrc lbout
    set sel [$lbsrc curselection]
    if { [llength $sel] < 1 } { return }
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	set lbl [$lbsrc get $ix]
	$lbout insert end $lbl
    }
    set sel [lsort -decreasing -integer $sel]
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	$lbsrc delete $ix
    }
}
proc toSrc {} {
    global lbsrc lbout
    set sel [$lbout curselection]
    if { [llength $sel] < 1 } { return }
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	set lbl [$lbout get $ix]
	$lbsrc insert end $lbl
    }
    set sel [lsort -decreasing -integer $sel]
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	$lbout delete $ix
    }
}

# utility for saving user selected columns to a new file for plotting progs
proc saveColumns {} {
    global viewsave inclhead saveset currentDir
    global lbout viewcollabels viewncol
    global readlines readfinished
    # make a list of the output columnindices
    set outlbls [$lbout get 0 end]
    if { [set nout [llength $outlbls]] < 1 } {
	Dialog_Warn "No columns selected for save file"
	return
    }
    foreach {ix lbl} [array get viewcollabels] { set inverse($lbl) $ix }
    set outix {}
    for {set i 0} {$i < $nout} {incr i} {
	lappend outix $inverse([lindex $outlbls $i])
    }
    set outi [join $outix ,]
    # now ready to copy the source dataset using just the spec columns
    pipeWrite "d$saveset c $outi"
    vwait readfinished
    # this will return the output setnumber
    if { [scan [lindex $readlines 0] {%s %d} dum outset] < 2 } {
	Dialog_Warn "Set column extraction failed."
	return
    }
    # delete the header if inclhead == 0
    if { $inclhead == 0 } {
	pipeWrite "d$outset hd"
	vwait readfinished
    }
    # set the filename for outset
    pipeWrite "d$outset nm [file join $currentDir $viewsave]"
    vwait readfinished
    # now write the outset to disk in raw mode (no DFheader stuff)
    pipeWrite "d$outset wr"
    vwait readfinished
    pipeWrite "d$outset"
    vwait readfinished
    # finally delete the outset
    pipeWrite "d-$outset"
    vwait readfinished
}



set inclhead 1
set lbout ""
set lbsrc ""
proc buildSaveColumns {} {
    global viewfile viewsave inclhead
    global lbout lbsrc
    set vc [toplevel .vcol]
    wm title $vc "Save viewfile columns"
    wm iconname $vc SaveCol
    # mb has viewfileName  saveAsName  saveButton  inclHeaderCheckBox
    set mb [frame $vc.mb]
    set vfNameLbl [label $mb.vfl -text viewfile]
    set vfName [label $mb.vf -textvariable viewfile -width 24 \
		    -relief sunken -bd 1 -bg white]
    set saNameLbl [label $mb.sal -text saveAs]
    set saName [entry $mb.sa -textvariable viewsave -width 24]
    set inch [checkbutton $mb.inch -variable inclhead -text inclHead]
    set save [button $mb.save -text "SAVE" \
		  -command saveColumns \
		  -background green -activebackground yellow]


    set f [frame $vc.f]
    set lbsrc [listbox $f.src -relief raised -bd 2 \
		  -yscrollcommand "$f.srcscroll set" -selectmode extended]
    set srcscroll [scrollbar $f.srcscroll -command "$f.src yview"]
    set move [frame $f.move]
    set lbout [listbox $f.out -relief raised -bd 2 \
		  -yscrollcommand "$f.outscroll set" -selectmode extended]
    set outscroll [scrollbar $f.outscroll -command "$f.out yview"]

    set toOut $move.out
    set toSrc $move.src
    button $toOut -command toOut -image rightarrow \
	-background green -activebackground yellow
    button $toSrc -command toSrc -image leftarrow \
	-background green -activebackground yellow


    grid $mb
    grid $f

    grid $vfNameLbl $vfName $saNameLbl $saName $inch $save
    grid $srcscroll $lbsrc -sticky news
    grid $move -row 0 -column 2
    grid $lbout -row 0 -column 3 -sticky news
    grid $outscroll -row 0 -column 4 -sticky news
    grid $toOut
    grid $toSrc
}
proc openSaveColumns {} {
    if { ! [winfo exists .vcol] } { buildSaveColumns }
    updateSaveColumns
    wm deiconify .vcol
    raise .vcol
}

####################  plot a dataset #########################

array set setplot {}

proc plotDS {} {
    # binding to Button1 in sellectlist or lbxs(0)
    # if dsmode, plot a dataset. NOT PB-plot although can handle multiple y
    # so can plot FOP result(.p file) from PB
    # if selected fil is already in a dataset, use that
    # else readFile

    # this proc finds the dataset, and determines dataset column-labels
    # and column-data indices for plotting, then calls showDSplot

    # store the column labels/stats in splabelsarr and splabelslst
    # try to figure out what X Y Err are for the plot
    # e.g. should recognize Suu Sdd .. pb-datapts .p file data
    # or use mostDistinct to find an X if x assigned (and E stnd)
    # if no y assigned chose a random vary column

    global readlines nread readfinished
    global spfile spsn
    global setplot spcfg
    global dsmode dseditm dseditY fsSymb gfsSymb

    if { ! $dsmode } { return }

    # save any previous DS file plot configuration
    if { [info exists spfile] && [string length $spfile] > 0 } {
	set setplot($spfile) [array get spcfg]
    }

    set sel [selectionOwned 1 1 1]
    if { [llength $sel] < 1 } {
	Dialog_Warn "Selection is not in one of the file selector listboxes.\n"
	return {}
    }
    set fil [lindex $sel 0]
    set spsn 0
    set spfile [file tail $fil]
    # look for the dataset by input src filename in case
    # user is using sellectlist for the filename
    pipeWrite "d ne $spfile"
    vwait readfinished
    # returns: dataset # filename= filename  #=0 if not found
    set rl [lindex $readlines 0]
    if { [scan $rl {%s %d} dum spsn] < 1 || $spsn < 1 } {
	# try to read the file
	set ret [readFile $spfile]
	if { [llength $ret] < 1 } {
	    Dialog_Warn "file $spfile is not in a dataset and read failed.\n"
	    return
	}
	set spsn [lindex $ret 0]
    }
    pipeWrite "d$spsn"
    vwait readfinished
    set rl [lindex $readlines 1]
    if { [scan $rl {%[*]%d %d} dum iset npt] < 3 || $npt < 1 } {
	Dialog_Warn "file $spfile has no data.\n"
	return
    }


    # store the column stats
    pipeWrite "c s"
    vwait readfinished
    arrayClear splabelsarr
    set splbls {}
    for {set i 1} {$i <= $nread} {incr i} {
	set rl [lindex $readlines $i]
	if { [scan $rl {%d %s %g %g %d %d} icol lbl mn mx nch var] < 6 } {
	    continue
	}
	set splabelsarr($lbl) [list $var $icol]
	lappend splbls [list $lbl $var $icol]
    }
    set splabelslst [array names splabelsarr]

    # if file previously read, should be able to recover columns from setplot
    set colsOK 0
    if { [info exists setplot($spfile)] } {
	array set spcfg $setplot($spfile)
	# try to determine plot columns from spcfg if file prev read
	set colsOK [getSPcfgCols]
	if { [info exists spcfg(y1-lbl)] && \
		 [string match $spcfg(y1-lbl) Suu] } {
	    configSxy
	    set setplot($spfile) [array get spcfg]
	}
	if { $colsOK } {
	    showDSplot
	    return
	}
    }

    # try assigned cols
    pipeWrite "c al"
    vwait readfinished
    if [info exists stndcol] { unset stndcol }
    for {set i 1} {$i <= $nread} {incr i} {
	set rl [lindex $readlines $i]
	set stdlbl [lindex $rl 0]
	set indx [lindex $rl 1]
	set icol [lindex $rl 2]
	set lbl [lindex $rl 3]
	if { $indx > 0 } { append stdlbl $indx }
	set stndcol($stdlbl) [list $lbl $icol]
    }
    set ya [array names stndcol y*]
    set ea [array names stndcol e*]
    set xa [array names stndcol x*]
    # need at least 1 y and 1 x or En or other col-var=1 
    # if En is assigned and varies use that for x
    set xok 0
    if { [info exists stndcol(En)] && \
	     [lindex $splabelsarr([lindex $stndcol(En) 0]) 0] > 0 } {
	set spcfg(xlbl) [lindex $stndcol(En) 0]
	set spcfg(xcol) [lindex $stndcol(En) 1]
	set xok 1
    } else {
	foreach x $xa {
	    set lbl [lindex $stndcol($x) 0]
	    if { [lindex $splabelsarr($lbl) 0] > 0 } {
		set spcfg(xlbl) $lbl
		set spcfg(xcol) [lindex $stndcol($x) 1]
		set xok 1
		break
	    }
	}
    }

    set ny [llength $ya]
    set isSxy 0
    foreach y $ya {
	set spcfg(${y}-lbl) [lindex $stndcol($y) 0]
	set spcfg(${y}-col) [lindex $stndcol($y) 1]
	if { [lsearch {Suu Sdd Sdu Sud} $spcfg(${y}-lbl)] >= 0 } { set isSxy 1}
    }
    set ne [llength $ea]
    foreach e $ea {
	set spcfg(${e}-lbl) [lindex $stndcol($e) 0]
	set spcfg(${e}-col) [lindex $stndcol($e) 1]
    }
    if { $ne != $ny || $ny == 0 } { set ne 0 }
    # if at this point ny == 0 
    # first check for Suu Sdd Sdu Sud labels and Exx for errors
    # if there are no Sxx labels find a vary==2,1,0 column

    if { $ny == 0 } {
	foreach sxx {Suu Sdd Sdu Sud} exx {Euu Edd Edu Eud} {
	    if { [info exists splabelsarr($sxx)] } {
		incr ny
		set spcfg(y${ny}-lbl) $sxx
		set spcfg(y${ny}-col) [lindex $splabelsarr($sxx) 1]
	    }
	    if { [info exists splabelsarr($exx)] } {
		incr ne
		set spcfg(e${ny}-lbl) $exx
		set spcfg(e${ny}-col) [lindex $splabelsarr($exx) 1]
	    }
	}
	if { $ny > 0 } { set isSxy 1 }
    }
    set splbls [lsort -decreasing -integer -index 1 $splbls]
    set nsplbs [llength $splbls]
    if { $ny == 0 } {
	set ne 0
	set splb [lindex $splbls 0]
	set spcfg(y1-lbl) [lindex $splb 0]
	set spcfg(y1-col) [lindex $splb 2]
	set ny 1
    }
    if { ! $xok } {
	for {set i [expr $nsplbls - 1]} {$i > 0} {incr i -1} {
	    set splbl [lindex $splbls $i]
	    if { [lindex $lbl 1] > 0 } {
		set spcfg(xlbl) [lindex $splb 0]
		set spcfg(xcol) [lindex $splb 2]
		set xok 1
		break
	    }
	}
    }
    if { ! $xok || $ny < 1 } {
	Dialog_Warn "couldn't find x-y columns in $spfile.\n"
	return
    }
    set spcfg(ny) $ny
    set spcfg(ylbl) ""
    $dseditm delete 0 end
    $dseditm add radio -variable dseditY -label all -value all
    for {set i 1} {$i <= $ny} {incr i} {
	append spcfg(ylbl) "$spcfg(y${i}-lbl) "
	# put the ylabels as radio into dseditm
	$dseditm add radio -variable dseditY \
	    -label $spcfg(y${i}-lbl) -value $spcfg(y${i}-lbl)
	# make sure y#-gsym pixels color exist
	if { ![info exists spcfg(y${i}-color)] } { set spcfg(y${i}-color) blue }
	if { ![info exists spcfg(y${i}-pixels)] } { set spcfg(y${i}-pixels) 6 }
	if { ![info exists spcfg(y${i}-fill)] } {set spcfg(y${i}-fill) defcolor}
	if { ![info exists spcfg(y${i}-sym)] } { set spcfg(y${i}-sym) plus }
	if { ![info exists spcfg(y${i}-gsym)] } { set spcfg(y${i}-gsym) 1 }
    }
    set spcfg(sym) circle-fill
    if { $isSxy } { configSxy }
    set spcfg(color) blue
    set spcfg(pixels) 6
    set spcfg(fill) defcolor
    set spcfg(logscale) 0
    set dseditY all
    set spcfg(errs) 1
    set spcfg(syms) 1
    if { $ne != $ny } { set spcfg(errs) 0 }

    # save spcfg
    set setplot($spfile) [array get spcfg]
    # ready to showDSplot
    showDSplot
}

proc configSxy {} {
    global spcfg fsSymb gfsSymb
    for {set i 1} {$i <= 4} {incr i} {
	set spcfg(y${i}-sym) $fsSymb($i)
	set spcfg(y${i}-gsym) $gfsSymb($i)
    }
}

proc getSPcfgCols {} {
    global spcfg
    if { ! [info exists spcfg(ny)] } { return 0 }
    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	if { ! [info exists spcfg(y${i}-col)] } { return 0 }
	if { ! [info exists spcfg(y${i}-lbl)] } { return 0 }
    }
    return 1
}

proc showDSplot {} {
    # called to actually show the DS plot after columns determined
    # or for reconfig plot

    global spfile spsn spcfg
    global BLT GNU
    global gvds gwds gvtitle
    global gdspipe gnuprog

    if { $BLT } {
	if { ! [info exists gvds] || ! [winfo exists $gvds] } {
	    # create a new top level and put a graph widget in it
	    set gvds [toplevel .gds]
	    raise $gvds
	    set gwds [graph $gvds.gds]
	    pack $gwds -side top -fill both -expand true
	    wm iconname $gvds
	}
	wm title $gvds "$gvtitle  file: $spfile"
	bltDSGraphSetup
	# apply the per group graph options. Wont need any updating for BLT
	# Do NOT call applyGraphOptions 
	#applyGraphOptions 0
	return
    }

    # rest of this is for GNU

    if { ! [info exists gdspipe] || \
	     [string length $gdspipe] < 1 } {
	if [catch {open "|\"$gnuprog\"" r+} gdspipe] {
	    Dialog_Warn "failed to open new gnuplot window for DSmode.\n"
	    return 0
	}
    }


    if { $saveGnuCmnds } {
	set gnuplotFPT [open gnuplotCMNDS w]
    }

    if { $toscreen } {
	gdspipeWrite "set terminal $screen title \"$gvtitle file: $spfile\""
	gdspipeWrite "set output"
    } else {
	set gnuopts ""
	if { [string match postscript $gnuftyp] } {
	    set gnuopts "$gnuPSopt(landscape) $gnuPSopt(color)"
	}
	gdspipeWrite "set terminal $gnuftyp $gnuopts"
	gdspipeWrite "set ouput $spfile.$gnuftyp"
    }
    gdspipeWrite "set xlabel \"$spcfg(xlbl)\""
    gdspipeWrite "set ylabel \"$spcfg(ylbl)\" font \"Helvetica,16\" tc rgb \"blue\""
    gdspipeWrite "set ytics mirror font \"Helvetica,14\" tc rgb \"blue\""
    gdspipeWrite "set xlabel font \"Helvetica,16\""
    gdspipeWrite "set xtics mirror font \"Helvetica,14\""
    gdspipeWrite "set key outside top"

    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	set gnudsstyl($i) "set style line $i pt $spcfg(y${i}-gsym) lc rgb \"$spcfg(y${i}-color)\""
	set gnudscol($i) "$spcfg(xcol):$spcfg(y${i}-col):$spcfg(e${i}-col)"
    }
    
    set pltcmds {}

    set dswith points
    if { $spcfg(errs) } { set dswith yerrorbars }

    title
    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	if { [info exists spcfg(y${i}-show)] && \
		 [oneC $spcfg(y${i}-show)] } { continue }
	set pss [expr $spcfg(y${i}-pixels)/2]
	gdspipeWrite "$gnudsstyl($i) ps $pss"
	lappend pltcmds "\'$spfile\' using $gnudscol($i) with $dswith ls $i axes x1y1 title \"$spcfg(y${i}-lbl)\""
    }

    # apply current group graph options, no replot
    #applyGraphOptions 0
    gdspipeWrite "plot [join $pltcmds ,]"
    if { $saveGnuCmnds } {
	close $gnuplotFPT
	set gnuplotFPT ""
    }
}



proc buildPlotColumns {} {
    # build toplevel .pc plotcolumns gui to select plot columns
    # to plot any dataset
    # store by short-filename
    # when this is called a file(spfile) has been selected and is
    # already in a dataset(spsn)
    # the master set labels list is in splabelsarr(lbl) which are
    #  [list vary icol]
    # the previous pcxye gets saved as pcxyelast

    # the plotcolumns(labels) symbols symbol-sizes colors (multiple y)
    # axis-scaling axis-labels plot-title
    # array set spconfig $setplot(short-fil)
    # set setplot(short-fil) [array get spconfig]
    # spconfig array elements for:
    #   ny  xlbl
    #   y1-lbl y1-elbl y1-symbolname y1-symbolsize y1-colora y1-colorb
    #   y1-showsym   y1-showerr ...
    #   xaxlbl yaxlbl xaxautomin xaxautomax xaxmin xaxmax
    #                 yax same-but-incl log
    # load all into spconfig

    # this toplevel will set spconfig(xlbl) spconfig(y1-lbl) spconfig(y1-elbl)
    # etc for multiple y

    # for x-selection use only columns that vary (vary>0) with non-rand pref(1)
    # for y-selection put columns that vary rand (vary==2) at top
    # to do this put collabel vary as list foreach collabel
    # then first sort -dictionary on -index 0 then -integer


    global spconfig spfile spsn
    # the move label procs will use the listbox widgets
    global lbspx lbspy lbspe lbspsrc

    set pc [toplevel .pc]
    wm title $pc "plot columns"
    wm iconname $pc PlotCol
    # mb has: Apply(c1)  X-radio(c2)  Y-radio(c4) E-radio(c5)
    # but I think this will build OK all on one grid (single frame)
    set apply [button $pc.apply -text Apply \
		   -command applySetColumns \
		   -background green -activebackground yellow]

    # xye radio select buttons
    set xr [radiobutton $pc.x -variable pcxye -value 1 -text X -command pcradio]
    set yr [radiobutton $pc.y -variable pcxye -value 2 -text Y -command pcradio]
    set er [radiobutton $pc.e -variable pcxye -value 3 -text Err \
	       -command pcradio]

    # 4 lbxs
    set lbspsrc [listbox $pc.lbs -relief raised -bd 2 \
		  -yscrollcommand "$pc.sscroll set" -selectmode extended]
    set sscroll [scrollbar $pc.sscroll -orient vertical \
		     -command "$pc.lbs yview"]

    set lbspx [listbox $pc.lbx -relief raised -bd 2 \
		    -selectmode single]
    set lbspy [listbox $pc.lby -relief raised -bd 2 \
		  -yscrollcommand "$pc.yscroll set" -selectmode extended]
    set yscroll [scrollbar $pc.yscroll -orient vertical \
		     -command "$pc.lby yview"]
    set lbspe [listbox $pc.lbe -relief raised -bd 2 \
		  -yscrollcommand "$pc.escroll set" -selectmode extended]
    set escroll [scrollbar $pc.escroll -orient vertical \
		     -command "$pc.lbe yview"]

    
    # 4 lbx labels
    set xlbl [label $pc.xlbl -text X-column]
    set ylxl [label $pc.ylbl -text Y-columns]
    set elbl [label $pc.elbl -text Err-columns]
    set slbl [label $pb.slbl -text Src-columns]

    # move col-labels from src to selected lbox and vice versa
    set toOut [button $pc.out -command pctoOut -image rightarrow \
		   -background green -activebackground yellow]
    set toSrc [button $pc.src -command pctoSrc -image leftarrow \
		   -background green -activebackground yellow]


    grid $apply -row 1 -column 1
    grid $xr -row 1 -column 2
    grid $yr -row 1 -column 4
    grid $er -row 1 -column 5

    grid $sscroll -row 4 -rowspan 17 -column 0 -sticky ns
    grid $lbspsrc -row 4 -rowspan 17 -column 1 -columnspan 2 -sticky news
    grid $slbl    -row 3 -column 1 -columnspan 2 -sticky ew

    grid $toOut   -row 8 -column 3 -padx 10
    grid $toSrc   -row 11 -column 3 -padx 10

    grid $lbspx   -row 4 -column 4 -columnspan 2 -sticky news
    grid $xlbl    -row 3 -column 4 -columnspan 2 -sticky ew

    grid $lbspy   -row 7 -rowspan 6 -column 4 -columnspan 2 -sticky news
    grid $ylbl    -row 6 -column 4 -columnspan 2 -sticky ew
    grid $yscroll -row 7 -rowspan 6 -column 6 -sticky ns

    grid $lbspe   -row 15 -rowspan 6 -column 4 -columnspan 2 -sticky news
    grid $elbl    -row 14 -column 4 -columnspan 2 -sticky ew
    grid $escroll -row 15 -rowspan 6 -column 6 -sticky ns
}

proc pctoOut {} {
    global lbspsrc lbspx lbspy lbspe
    global pcxye
    set sel [$lbspsrc curselection]
    if { [llength $sel] < 1 } { return }
    switch -exact -- $pcxye {
	1 {set lbout $lbspx}
	2 {set lbout $lbspy}
	3 {set lbout $lbspe}
    }
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	set lbl [$lbsrc get $ix]
	$lbout insert end $lbl
    }
    set sel [lsort -decreasing -integer $sel]
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	$lbsrc delete $ix
    }
}
proc pctoSrc {} {
    global lbspsrc lbspx lbspy lpspe
    global pcxye
    switch -exact -- $pcxye {
	1 {set lbout $lbspx}
	2 {set lbout $lbspy}
	3 {set lbout $lbspe}
    }
    set sel [$lbout curselection]
    if { [llength $sel] < 1 } { return }
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	set lbl [$lbout get $ix]
	$lbspsrc insert end $lbl
    }
    set sel [lsort -decreasing -integer $sel]
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set ix [lindex $sel $i]
	$lbout delete $ix
    }
}

set pcxye 0
set pcxyelast 0
array set splabelsarr {}
set splabelslst {}
proc pcradio {} {
    # choose and sort from master list according to pcxye
    # the master list splabelsarr is array
    # splabelsarr(collabel) = [list vary icol]
    # so we can find the vary or icol for any label quickly
    # and the master list of all labels is just [array names splabelsin]
    # which is in splabelslst
    # the invariant columns are in splabelsinv

    # the current list in lbspsrc is the truncated list splabels
    # because some labels may have been selected (moved to x y OR e lbx)
    # in order to select/sort the current splabels
    # make a list of lists on elem foreach label appending the vary value.
    # then do -index sorts and truncations
    # the result goes into the lbspsrc listbox

    # for x-selection use only columns that vary (vary>0) with non-rand pref(1)
    # for y-selection put columns that vary rand (vary==2) at top
    # to do this put collabel vary as list foreach collabel
    # then first sort -dictionary on -index 0 then -integer

    global splabelsarr splabelslst
    global lbspsrc pcxye pcxyelast
    # this gets called on radiobutton pcxye change
    # so if pcxyelast > 1 & pcxye > 1 just return
    if { $pcxyelast > 1 && $pcxye > 1 } { return }

    # take the master and remove every label in x y e lbxs
    # if pcxye == 1 also remove every label that has vary = 0
    # do sort on icol then on vary (incr for x) (decr for ye)

    $lbspsrc delete 0 end
    set splabels $splabelslst
    foreach lbl [$lbspx get 0 end] {
	if { [set i [lsearch $splabels $lbl]] < 0 } { continue }
	set splabels [lreplace $splabels $i $i]
    }
    foreach lbl [$lbspy get 0 end] {
	if { [set i [lsearch $splabels $lbl]] < 0 } { continue }
	set splabels [lreplace $splabels $i $i]
    }
    foreach lbl [$lbspe get 0 end] {
	if { [set i [lsearch $splabels $lbl]] < 0 } { continue }
	set splabels [lreplace $splabels $i $i]
    }
    if { $pcxye == 1 } {
	foreach lbl $splabels {
	    if { [lindex $splabelsarr($lbl) 0] == 0 } {
		if { [set i [lsearch $splabels $lbl]] < 0 } { continue }
		set splabels [lreplace $splabels $i $i]
	    }
	}
    }
    if { [llength $splabels] < 1 } { return }
    set nlst {}
    foreach lbl $splabels {
	lappend nlst [list $lbl [lindex $splabelsarr($lbl) 0] \
			  [lindex $splabelsarr($lbl) 1]]
    }
    set nlst [lsort -integer -increasing -index 2 $nlst]
    if { $pcxye == 1 } {
	set sd -increasing
    } else {
	set sd -decreasing
    }
    set nlst [lsort -integer $sd -index 1 $nlst]
    foreach n $nlst { $lbspsrc insert end $n }
}

proc applySetColumns {} {
}

set tvtitle "text file viewer"
set sfvm ""
proc buildShowText {args} {
    # top frame for group and file type selection
    # bottom is scrolledText
    global tv txtw
    global gnum
    global mp snum sfvm
    global tvtitle

    set tv [toplevel .txtvw]
    wm title $tv $tvtitle
    wm iconname $tv "PB txt"
    wm iconify $tv

    set mb [frame $tv.mb]

    set inc $mb.inc
    set dec $mb.dec
    set cur $mb.cur
    set curlbl $mb.curlbl

    set dptf $mb.dptf
    set msrf $mb.msrf
    set srcf $mb.srcf
    set self $mb.self

    #set sninc $mb.sninc
    #set sndec $mb.sndec
    #set sncur $mb.sncur
    #set sncurlbl $mb.sncurlbl
    set sfv $mb.sfv

    set fl $mb.file
    menubutton $fl -text File -menu $fl.menu -bg #00aa00 -relief raised \
	    -activebackground #ffaa00
    set fm [menu $fl.menu -tearoff 0]
    $fm add command -label "save columns from current file to new file" \
	-command saveColumnsFromFileView

    button $inc -command incrGroup -image uparrow \
	-background green -activebackground yellow
    button $dec -command decrGroup -image downarrow \
	-background green -activebackground yellow
    label  $cur -textvariable gnum -width 3 -relief sunken -bd 1 -bg white
    label  $curlbl -text group

    radiobutton $dptf -variable mp -text datapts -value 1
    radiobutton $msrf -variable mp -text msrmnts -value 2
    radiobutton $srcf -variable mp -text srcfiles -value 3
    radiobutton $self -variable mp -text select -value 4
    trace variable mp w showText

    # replace these with just a menu of the src files with radio
    #button $sninc -command incrSn -image uparrow \
    #	-background green -activebackground yellow
    #button $sndec -command decrSn -image downarrow \
    #	-background green -activebackground yellow
    #label  $sncur -textvariable snum -width 3 -relief sunken -bd 1 -bg white
    #label  $sncurlbl -text dataset

    menubutton $sfv -text srcFile -menu $sfv.menu -bg #00aa00 -relief raised \
	    -activebackground #ffaa00
    set sfvm [menu $sfv.menu -tearoff 0]



    grid $fl $curlbl $cur $inc $dec \
	$dptf $msrf $srcf $sfv $self
    #$sninc $sndec $sncur $sncurlbl $self


    set tf $tv.tf
    set txtw [eval scrolledText $tf $args]
    pack $mb -side top -fill x
    pack $tf -side top -fill both -expand true
}

proc Scroll_Set {sb geo offset size} {
    if { $offset != 0.0 || $size != 1.0 } {
	eval $geo
	$sb set $offset $size
    } else {
	set manager [lindex $geo 0]
	$manager forget $sb
    }
}

proc scrolledText {f args} {
    frame $f
    text $f.text \
	-xscrollcommand [list Scroll_Set $f.xscroll \
			     [list grid $f.xscroll -row 1 -column 0 \
				  -sticky we]] \
	-yscrollcommand [list Scroll_Set $f.yscroll \
			     [list grid $f.yscroll -row 0 -column 1 \
				  -sticky ns]]
    eval {$f.text configure} $args
    scrollbar $f.xscroll -orient horizontal \
	-command [list $f.text xview]
    scrollbar $f.yscroll -orient vertical \
	-command [list $f.text yview]
    grid $f.text $f.yscroll -sticky news
    grid $f.xscroll -sticky news
    grid rowconfigure $f 0 -weight 1
    grid columnconfigure $f 0 -weight 1
    return $f.text
}

################# group plot procs ###########################

proc vectorExists {vecName} {
    return [expr [llength [vector names ::$vecName]] \
	    || [llength [vector names $vecName]]]
}

set maxgrpsets 0
set pVecNames {suu sdd sdu sud euu edd edu eud qx qy qz en}
set pDatNames {Suu Sdd Sdu Sud Euu Edd Edu Eud Qx Qy Qz En}
set pStdNames {y1  y2  y3  y4  e1  e2  e3  e4  x1 x2 x3 En}

set mVecNames {cc ccerr qxm qym qzm enm fs}
set mDatNames {Cc Ccerr Qxm Qym Qzm Enm Fs}
set mStdNames {y1  e1   x1  x2  x3  En  fs}
foreach pdat $pDatNames { set $pdat {} }
foreach mdat $mDatNames { set $mdat {} }
set Xpts {}
set Xmsr {}
set iX 1
set Src {}
set iSrc {}
set Rf {}
set iWts {}
set stylst {}
set rstylst {}
array set pen {}
array set pencolor {}
set icolor ff0000
set fcolor 00ff00
# symbols for uu dd du ud
array set fsSymb {1 circle 2 square 3 arrow 4 triangle}
# gnu equivalent filled pointstyle numbers from wxt & windows term (pt index)
array set gfsSymb {1 7     2 5      3 11    4 9}
array set bltSymbols {
    1 plus 2 cross
    3 square 4 circle
    5 triangle 6 arrow
    7 diamond
    8 square-fill 9 circle-fill
    10 triangle-fill 11 arrow-fill
    12 diamond-fill
}
set bltSymbolsLst {
    plus cross square circle triangle arrow diamond
    square-fill circle-fill triangle-fill arrow-fill diamond-fill
}
foreach {i val} [array get bltSymbols] { set bltSymbolIndex($val) $i }

# gnu symbols (pt) corresponding to BLT
array set gnuSymbols {
    1 1  2 2
    3 4  4 6
    5 8  6 10
    7 12
    8 5  9 7
    10 9 11 11
    12 13
}

proc mostDistinct {args} {
    # args is list of lists of data
    # sort increasing to easily calc number of distinct values
    # return index of list with most distinct values
    # on tie, earliest list wins
    set nl [llength $args]
    if { $nl < 1 } { return 0 }
    if { $nl == 1 } { return 1 }
    set ibest 1
    set mostdist 0
    for {set i 0} {$i < $nl} {incr i} {
	set lst [lsort -increasing -real [lindex $args $i]]
	set ll [llength $lst]
	if { $ll < 1 } { continue }
	set ndist 1 ;
	set a [lindex $lst 0]
	for {set j 1} {$j < $ll} {incr j} {
	    set b [lindex $lst $j]
	    if { $b > $a } { incr ndist }
	    set a $b
	}
	if { $ndist > $mostdist } {
	    set mostdist $ndist
	    set ibest [expr $i + 1]
	}
    }
    return $ibest
}

proc checkPens {} {
    global g maxgrpsets ngrpsets pen
    # create any needed pens for BLT
    set is [expr $maxgrpsets + 1]
    set penlst [$g pen names]
    #for {set i $is} {$i <= $ngrpsets} {incr i}
    for {set i 1} {$i <= $ngrpsets} {incr i} {
	for {set j 1} {$j <= 4} {incr j} {
	    set pen($j$i) p$j$i
	    if { [lsearch $penlst p$j$i] < 0 } { $g pen create p$j$i }
	    if { [lsearch $penlst r$j$i] < 0 } { $g pen create r$j$i }
	}
    }
    if { $ngrpsets > $maxgrpsets } { set maxgrpsets $ngrpsets }    
}

array set nsetfs {}
array set pcolumns {}
array set isrcplot {}
array set xaxlbl {1 E 2 Qx 3 Qy 4 Qz}
array set gnufil {}
array set gnustyl {}
array set gnucol {}

proc groupPlotSetup {} {
    # call this after a group is calc and PLT
    # this setup will work for gnuplot as well
    # this will init plot data and vectors
    # showPlot to display the results for cur grp is a separate call
    global g gnum
    global mset pset
    global nmsr npts
    global ngrpsets maxgrpsets
    # following are lists of data list names
    global pDatNames mDatNames
    global pStdNames mStdNames
    global Xpts Xmsr
    global Src Rf iSrc iWts stylst rstylst iX
    global isrcplot
    global readlines readfinished
    global pen pencolor
    global icolor fcolor
    global fsSymb gfsSymb
    global titlefont tickfont
    global BLT GNU
    global pcolumns nsetfs
    global textf
    global gnufil gnustyl gnucol

    # the corresponding vectors have same name l.c.
    # mset and pset set in groupSrcFiles when write gnum

    # call after solveGroup
    # get all the vector data lists
    # decide what to vary in the plot
    # already have ndatasets and srcfiles listings, pset mset
    # from groupSrcFiles

    # would like to set it up as for groupOptions
    # where there is a list of commands that gets run
    # when group changes
    # need vector names:

    # get vector data using
    # c al stndCol(x# y# En -> ctyp indx dscol label
    # c dl label -> lbl data ...
    # from .d need y1-4 = Suu Sdd Sdu Sud e1-4 x1-3 = Qxyz En
    # from .m need y1 = Cc e1 = Ccerr x1-3 = Qxyz En fs srcfile
    # symbols uu=+ du=arrow ud=triangle dd=circle
    # colors S-blue sets red->green ff0000->00ff00
    # use element -hide
    # ? can I use .m as a single element
    # could turn off a pen by setting symbol blank
    # so could have a pen for each set but then
    # how do change symbols for spin state
    # would need 4 pens per set
    # p1-4 weigths 1-4 p5-8 weights 5-8 etc
    # so weight for each msr = 4(setn-1)+fs
    # pens p{fs}{sn}
    # the hack to turn off the view for a set is to
    # set the color for its pens to the plotbackground color

    if { $ngrpsets < 1 } {
	Dialog_Warn "No group sets"
	return
    }
    if { $mset < 1 || $pset < 1 } {
	Dialog_Warn "NO msr/pts datasets for plotting current group"
	return
    }

    # nothing in this proc should depend on BLT graph existing yet
    # if { $BLT } { checkPens }

    arrayClear pcolumns

    # get the pset(.p)  data lists
    pipeWrite "d$pset"
    vwait readfinished
    foreach pdat $pDatNames pstd $pStdNames {
	global $pdat
	pipeWrite "c al $pstd"
	vwait readfinished
	set rl [lindex $readlines 1]
	# ctyp indx dscol label
	if { [scan $rl {%s %d %d} ct inx dsc] < 3 } {
	    Dialog_Warn "failed to get p data"
	    return
	}
	pipeWrite "c dl $dsc"
	vwait readfinished
	# get rid of the leading label
	set $pdat [lreplace [lindex $readlines 0] 0 0]
	if { [eval llength $$pdat] < $npts } {
	    Dialog_Warn "not enough pts read for p data"
	    return
	}
	# save the column numbers for GNU
	set pcolumns($pdat) $dsc
    }
    # get the mset(.m)  data lists
    pipeWrite "d$mset"
    vwait readfinished
    foreach mdat $mDatNames mstd $mStdNames {
	global $mdat
	pipeWrite "c al $mstd"
	vwait readfinished
	set rl [lindex $readlines 1]
	# ctyp indx dscol label
	if { [scan $rl {%s %d %d} ct inx dsc] < 3 } {
	    Dialog_Warn "failed to get m data"
	    return
	}
	pipeWrite "c dl $dsc"
	vwait readfinished
	# get rid of the leading label
	set $mdat [lreplace [lindex $readlines 0] 0 0]
	if { [eval llength $$mdat] < $nmsr } {
	    Dialog_Warn "not enough msr read for m data"
	    return
	}
    }
    # also get the src column from mset
    pipeWrite "c dl src"
    vwait readfinished
    set Src [lreplace [lindex $readlines 0] 0 0]
    if { [llength $Src] < $nmsr } {
	Dialog_Warn "not enough msrs read for m data"
	return
    }
    # also get the R flipping ratio column from mset
    pipeWrite "c dl R"
    vwait readfinished
    set Rf [lreplace [lindex $readlines 0] 0 0]
    if { [llength $Rf] < $nmsr } {
	Dialog_Warn "not enough msrs read for m data"
	return
    }

    # get the snum indices for each srcfile in .m, and set the weights
    set iSrc {}
    set iWts {}
    arrayClear nsetfs
    for {set i 1} {$i <= $ngrpsets } {incr i} {
	for {set j 1} {$j <= 4} {incr j} {
	    set nsetfs($i,$j) 0
	}
    }
    for {set i 0} {$i < $nmsr} {incr i} {
	set fil [file tail [lindex $Src $i]]
	set is [snFromFile $fil]
	if { $is < 1 || $is > $ngrpsets } {
	    Dialog_Warn "failed to index src files"
	    return
	}
	lappend iSrc $is
	# get msr pt fs
	set ifs [lindex $Fs $i]
	lappend iWts [expr 4*($is - 1) + $ifs]
	incr nsetfs($is,$ifs)
    }
    # software configure the pens and make the styles list
    # NB initially all sets are viewed (pens have color)
    # so should init the iSrc list

    set dcolor 1
    set stylst {}
    set rstylst {}
    if { $ngrpsets > 1 } { set dcolor [expr 1./($ngrpsets-1)] }
    for {set i 1} {$i <= $ngrpsets} {incr i} {
	set pcolor [hexColor $icolor $fcolor [expr ($i-1)*$dcolor]]
	for {set j 1} {$j <= 4} {incr j} {
	    set pen($j$i) p$j$i
	    set pencolor($j$i) $pcolor
	    set pw [expr 4*($i-1)+$j]
	    lappend stylst [list p$j$i $pw $pw]
	    lappend rstylst [list r$j$i $pw $pw]
	}
	# record the view state in isrcplot array already done
	set isrcplot($i) 1
    }

    # determine the bestX from datapt En Qx Qy Qz
    set iX [mostDistinct $En $Qx $Qy $Qz]
    if { $iX == 1 } {
	set pcolumns(X) $pcolumns(En)
	set Xpts $En
	set Xmsr $Enm
    } elseif { $iX == 2 } {
	set pcolumns(X) $pcolumns(Qx)
	set Xpts $Qx
	set Xmsr $Qxm
    } elseif { $iX == 3 } {
	set pcolumns(X) $pcolumns(Qy)
	set Xpts $Qy
	set Xmsr $Qym
    } elseif { $iX == 4 } {
	set pcolumns(X) $pcolumns(Qz)
	set Xpts $Qz
	set Xmsr $Qzm
    }

    # setting the graph vectors from the data lists for a group
    # is how we display each group
    # so we must store the data lists
    # store pencolor array to view and unview each msr set
    # store the iviewsrc array
    # store stylst iWts lists

    # for GNU
    # we still read Sab Eab from the .p file, using pcolumns

    # for msr and R write separate files
    # for msr: file for each set and fs (this is easier than plot '-' stdin)
    #           group#_set#_fs#M.dat
    # for R: file for each set symbol pt 2  group#_set#R.dat

    if { $GNU } {
	# write the files for msr and R, create the styles
	# store the filenames and style strings by style index
	# append the pointsize to the style string when showPlot
	# setnum*5 + (fs-1) for msr and setnum*5 + 4 for R

	# styles for Suu ... 
	for {set i 1;set j 0;set k 4} {$i <= 4} {incr i;incr j; incr k} {
	    set gnufil($i) $textf(1)
	    set gnustyl($i) "set style line $i pt $gfsSymb($i) lc rgb \"blue\""
	    set gnucol($i) "$pcolumns(X):$pcolumns([lindex $pDatNames $j]):$pcolumns([lindex $pDatNames $k])"
	}
	for {set i 1} {$i <= $ngrpsets} {incr i} {
	    set si [expr 5*$i]
	    set ri [expr $si + 4]
	    set gnufil($ri) "group${gnum}_set${i}R.dat"
	    set gnustyl($ri) "set style line $ri pt 2 lc rgb \"\#$pencolor(1$i)\""
	    set gnucol($ri) "1:2"
	    set f [open $gnufil($ri) w]
	    for {set j 0} {$j < $nmsr} {incr j} {
		set is [lindex $iSrc $j]
		if { $is != $i } { continue }
		set x [lindex $Xmsr $j]
		set r [lindex $Rf $j]
		puts $f "$x $r"
	    }
	    close $f
	    for {set k 1} {$k <= 4} {incr k} {
		if { $nsetfs($i,$k) == 0 } { continue }
		set sk [expr $si + $k - 1]
		set gnufil($sk) "group${gnum}_set${i}_fs${k}M.dat"
		set gnustyl($sk) "set style line $sk pt $gfsSymb($k) lc rgb \"\#$pencolor(1$i)\""
		set gnucol($sk) "1:2:3"
		set f [open $gnufil($sk) w]
		for {set j 0} {$j < $nmsr} {incr j} {
		    set is [lindex $iSrc $j]
		    if { $is != $i } { continue }
		    set ifs [lindex $Fs $j]
		    if { $ifs != $k } { continue }
		    set x [lindex $Xmsr $j]
		    set c [lindex $Cc $j]
		    set ce [lindex $Ccerr $j]
		    puts $f "$x $c $ce"
		}
		close $f
	    }
	}	
    }
}

proc setPlotGeo {} {
    # put gv under the main window
    global gv tv
    if { ![info exists gv] || ![winfo exists $gv] || ![winfo ismapped $gv] } {
	 return
    }
    set gg [wm geometry $gv]
    set mg [wm geometry .]
    #set tg ""
    #if { [winfo exists $tv] && [winfo ismapped $tv] } {
    #set tg [wm geometry $tv]
    #}
    if {[scan $mg {%dx%d%d%d} xl yl x y] < 4} {return}
    # try to place gv below main
    if {[scan $gg {%dx%d%d%d} xlg ylg xg yg] < 4} {return}
    set yg [expr $y + $yl + 5]
    set xg $x
    set gg [format %dx%d%+d%+d $xlg $ylg $xg $yg]
    wm geometry $gv $gg
}

set plotbg white
set toscreen 1
array set groupiwin {}
array set viewS {1 suu 2 sdd 3 sdu 4 sud}
proc showPlot {} {
    # showPlot plots the current group via all local variables
    # into the current iwin
    # PLT flag must be set and group must be calc==1

    # call this when gnum changes and PLT, or plot request
    # groupSrcfile
    # groupPlotSetup gets called before this
    # to set up all the data lists

    global PLT BLT GNU gnuplotFPT screen toscreen gnuftyp gnuPSopt gnuprog
    global groupiwin
    global saveGnuCmnds
    global g gv gnum calc snum srcfiles gvtitle gpipe gpipes gws gvs iwin
    global ngrpsets
    global pDatNames pVecNames
    global mDatNames mVecNames
    global Xpts Xmsr iX Rf xaxlbl titlefont tickfont
    global stylst rstylst iWts
    global plotbg pencolor fsSymb
    global isrcplot eview viewS textf pcolumns nsetfs
    global eview viewS plotsymbolsize plotconfig
    global gnufil gnustyl gnucol
    global yctslbl y2ctslbl y2Rlbl
    global tcl_platform

    if { ! $PLT || ! $calc } { return }


    # only build one plot control window = gv
    # Apr 2010 changed to put plot-selector below group-selector in main window
    if { ![winfo exists $gv.title] } { buildPlotSelector }
    #if { ![winfo ismapped $gv] } { wm deiconify $gv }
    #wm title $gv "$gvtitle   group $gnum  plotwindow $iwin"
    #setPlotGeo
    #raise $gv

    set groupiwin($iwin) $gnum

    # may need to open a new gnuplot channel, or blt toplevel for graph
    if { $BLT } {
	if { ! [info exists gws($iwin)] || ! [winfo exists $gws($iwin)] } {
	    # create a new top level and put a graph widget in it
	    set gvs($iwin) [toplevel .g$iwin]
	    raise $gvs($iwin)
	    set gws($iwin) [graph $gvs($iwin).g$iwin]
	    pack $gws($iwin) -side top -fill both -expand true
	}
	# make this graph current g
	set g $gws($iwin)
	wm title $gvs($iwin) "$gvtitle  group $gnum  plotwindow $iwin"
	bltGraphSetup
	# apply the per group graph options. Wont need any updating for BLT
	applyGraphOptions 0
	return
    }

    # rest of this is for GNU

    if { ! [info exists gpipes($iwin)] || \
	     [string length $gpipes($iwin)] < 1 } {
	if [catch {open "|\"$gnuprog\"" r+} gpipe] {
	    Dialog_Warn "failed to open new gnuplot window $iwin"
	    return 0
	}
	set gpipes($iwin) $gpipe
    } else {
	set gpipe $gpipes($iwin)
    }


    if { $saveGnuCmnds } {
	set gnuplotFPT [open gnuplotCMNDS w]
    }
    set gpipe $gpipes($iwin)
    # configure the xaxis
    # set Xlabel <xlabel> 
    # set style line <index> pt <#> ps <#> lc rgb #rgbcolor
    # for M files index = 4*(setnum-1) + fs
    #
    if { $toscreen } {
	gpipeWrite "set terminal $screen title \"$gvtitle group $gnum  plotwindow $iwin\""
	gpipeWrite "set output"
    } else {
	set gnuopts ""
	if { [string match postscript $gnuftyp] } {
	    set gnuopts "$gnuPSopt(landscape) $gnuPSopt(color)"
	}
	gpipeWrite "set terminal $gnuftyp $gnuopts"
	gpipeWrite "set ouput grp$gnum.$gnuftyp"
    }
    gpipeWrite "set xlabel \"$xaxlbl($iX)\""
    
    gpipeWrite "set ylabel \"$yctslbl\" font \"Helvetica,16\" tc rgb \"blue\""
    gpipeWrite "set y2label \"$y2ctslbl\" font \"Helvetica,16\" tc rgb \"red\""
    gpipeWrite "set ytics nomirror font \"Helvetica,14\" tc rgb \"blue\""
    gpipeWrite "set y2tics autofreq font \"Helvetica,14\" tc rgb \"red\""
    gpipeWrite "set xlabel font \"Helvetica,16\""
    gpipeWrite "set xtics font \"Helvetica,14\""
    gpipeWrite "set key outside top"
    # when R selected change to flip-ratio
    # set style line <index> for the msr and R files
    # must be done when ngrpsets is known in showPlot
    
    
    # construct the plot cmd
    # suu sdd sdu sud styles 1-4
    # set1 msrfs  5- 8  R  9
    # set2 msrfs 10-13  R 14
    # set3 msrfs 15-18  R 19
    # msrfsIndex=setnum*5 + (fs-1)
    # Rindex    =setnum*5 + 4
    # send the style definitions with current plotsymbolsize(elem) appended
    
    set pltcmds {}
    
    if { $eview(R) } {
	set msraxes x1y1
	gpipeWrite "set y2label \"$y2Rlbl\""
    } else {
	set msraxes x1y2
	gpipeWrite "set y2label \"$y2ctslbl\""
    }

    set pss [expr $plotsymbolsize(pbcor)/2]
    set msrwith points
    set pbwith points
    if { $plotconfig(errmsr) } { set msrwith yerrorbars }
    if { $plotconfig(errpb) } { set pbwith yerrorbars }
    for {set i 1} {$i <= 4} {incr i} {
	if { ! $eview($viewS($i)) } { continue }
	gpipeWrite "$gnustyl($i) ps $pss"
	lappend pltcmds "\'$gnufil($i)\' using $gnucol($i) with $pbwith ls $i axes x1y1 title \"$viewS($i)\""
    }
    for {set i 1} {$i <= $ngrpsets} {incr i} {
	set si [expr 5*$i]
	if { $i == 1 } {
	    set rtitle "title \"flip-ratio\""
	    set stitle "title \"msr\""
	} else {
	    set rtitle notitle
	    set stitle notitle
	}
	if { $eview(R) } {
	    set ri [expr $si + 4]
	    set pss [expr $plotsymbolsize(R)/2]
	    gpipeWrite "$gnustyl($ri) ps $pss"
	    lappend pltcmds "\'$gnufil($ri)\' using $gnucol($ri) with points ls $ri axes x1y2 $rtitle"
	}
	if { $eview(msr) } {
	    set pss [expr $plotsymbolsize(msr)/2]
	    if { ! $isrcplot($i) } { continue }
	    for {set j 1} {$j <= 4} {incr j} {
		if { $nsetfs($i,$j) < 1 } { continue }
		if { $j < 3 && ! $eview(nsfm) } { continue }
		if { $j > 2 && ! $eview(sfm) } { continue }
		set sj [expr $si + $j - 1]
		gpipeWrite "$gnustyl($sj) ps $pss"
		lappend pltcmds "\'$gnufil($sj)\' using $gnucol($sj) with $msrwith ls $sj axes $msraxes $stitle"
	    }
	}
    }
    # apply current group graph options, no replot
    applyGraphOptions 0
    gpipeWrite "plot [join $pltcmds ,]"
    if { $saveGnuCmnds } {
	close $gnuplotFPT
	set gnuplotFPT ""
    }
    return
}


proc applyGraphOptions {update} {
    global GNU dsmode
    # update just calls showPlot if GNU
    applyAxScl 0
    plotLogScale 0
    plotSymbolSize 0
    plotErrbars 0
    plotConfig
    if { ! $dsmode } {
	viewMsrSet 0
	viewMsrFS 0
    } else {
	plotSymbol 0
	plotSymbolColor 0
    }
    viewElem 0
    if { $GNU && $update } { if { $dsmode } { showDSplot } else { showPlot } }
}

proc oneC {on} {
    return [expr 1 - $on]
}

proc viewElem {update} {
    global g gwds eview BLT GNU spcfg dsmode dseditY
    if { $BLT } {
	if { $dsmode } {
	    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
		if { ![string match $dseditY all] && \
			 ![string match $dseditY $spcfg(y${i}-lbl)] } { continue }
		if { [info exists spcfg(y${i}-show)] } {
		    set hidestate [oneC $spcfg(y${i}-show)]
		    $gwds element configure y$i -hide $hidestate
		}
	    }
	} else {
	    foreach el {suu sdd sdu sud msr R} {
		set hidestate [oneC $eview($el)]
		$g element configure $el -hide $hidestate
		if { $el == "R" } {
		    if { $hidestate == 1 } {
			$g y2axis use y2
			$g axis configure y2 -hide 0
			$g axis configure rf -hide 1
		    } else {
			$g y2axis use rf
			$g axis configure y2 -hide 1
			$g axis configure rf -hide 0
		    }
		}
	    }
	}
    }
    if { $GNU && $update } {
	if { $dsmode } { showDSplot } else { showPlot }
    }
}

proc viewMsrFS {update} {
    global eview ngrpsets
    global g
    global plotbg
    global pencolor fsSymb
    global BLT GNU
    if { $BLT } {
	for {set i 1} {$i <= $ngrpsets} {incr i} {
	    if { $eview(sfm) } {
		$g pen configure p3$i \
		    -errorbarcolor #$pencolor(3$i) -symbol $fsSymb(3)
		$g pen configure p4$i \
		    -errorbarcolor #$pencolor(4$i) -symbol $fsSymb(4)
	    } else {
		$g pen configure p3$i -errorbarcolor $plotbg -symbol ""
		$g pen configure p4$i -errorbarcolor $plotbg -symbol ""
	    }
	    if { $eview(nsfm) } {
		$g pen configure p1$i \
		    -errorbarcolor #$pencolor(1$i) -symbol $fsSymb(1)
		$g pen configure p2$i \
		    -errorbarcolor #$pencolor(2$i) -symbol $fsSymb(2)
	    } else {
		$g pen configure p1$i -errorbarcolor $plotbg -symbol ""
		$g pen configure p2$i -errorbarcolor $plotbg -symbol ""
	    }
	}
    }
    if { $GNU && $update } { showPlot }
}

proc viewMsrSet {update} {
    # command for selecting which src sets are viewed in plot
    # isrcplot isrcon button
    # April 2010 isrcplot is controlled directly from menu sfm
    # so isrcon variable no longer used
    # this proc is now called when msrSet radio is selected.

    global isrcon ngrpsets
    global g
    global plotbg
    global pencolor fsSymb
    global isrcplot plotconfig
    global BLT GNU
    if { $BLT } {
	# this is a hack to turn off data associated with a pen
	# set its color to plotbg
	# or better set it symbol to ""
	# the pens for snum are 4*(snum-1)+j  j=1,4
	#set snum $i
	#if { $snum < 1 || $snum > $ngrpsets } { return }
	for {set j 1} {$j <= $ngrpsets} {incr j} {
	    for {set i 1} {$i <= 4} {incr i} {
		if { $isrcplot($j) } {
		    if { $plotconfig(errmsr) } {
			$g pen configure p$i$j -errorbarcolor #$pencolor($i$j)
		    } else {
			$g pen configure p$i$j -errorbarcolor $plotbg
		    }
		    $g pen configure p$i$j -symbol $fsSymb($i)
		} else {
		    $g pen configure p$i$j -errorbarcolor $plotbg
		    $g pen configure p$i$j -symbol ""
		}
	    }
	}
    }

    if { $GNU && $update } { showPlot }
}

# hexColorLinear implements linear color interpolation
proc hexColorLinear {icolor fcolor frac} {
    set r ""
    set nh [string length $icolor]
    for {set i 0} {$i < $nh} {incr i} {
	set di [scan [string index $icolor $i] %x]
	set df [scan [string index $fcolor $i] %x]
	set dx [format %1x [expr int((1-$frac)*$di + $frac*$df)]]
	append r $dx
    }
    return $r
}
# implement color vector rotation with hexColor
# hexout = cos(frac*Pi/2)icolor + sin(frac*Pi/2)fcolor
# colorVector(theta) = R(thet) (cos(thet)ich + sin(thet)(unit(ichxfch))xich)
# for thet=(0,thetaf) cos(thetaf)=ich.fch
# R(thet) = Ri + (Rf-Ri)(theta/thetaf)

proc unitvector {vecarray} {
    upvar $vecarray a
    set mag 0.
    foreach {i val} [array get a] { set mag [expr $mag + $val*$val] }
    if { $mag <= 0. } { return 0. }
    set mag [expr sqrt($mag)]
    foreach {i val} [array get a] { set a($i) [expr $a($i)/$mag] }
    return $mag
}
proc dotpro {avec bvec} {
    upvar $avec a $bvec b
    set dot 0.
    foreach {i av} [array get a] {j bv} [array get b] {
	set dot [expr $dot + $av*$bv]
    }
    return $dot
}
proc scalevector {scl avec ascl} {
    upvar $avec a $ascl s
    foreach {i val} [array get a] { set s($i) [expr $scl*$a($i)] }
}
proc addvectors {avec bvec cvec} {
    upvar $avec a $bvec b $cvec c
    foreach {i av} [array get a] {j bv} [array get b] {
	set c($i) [expr $av + $bv]
    }
}

proc crosspro {avec bvec cvec} {
    # 3-vecs
    upvar $avec a $bvec b $cvec c
    set c(1) [expr $a(2)*$b(3) - $a(3)*$b(2)]
    set c(2) [expr $a(3)*$b(1) - $a(1)*$b(3)]
    set c(3) [expr $a(1)*$b(2) - $a(2)*$b(1)]
}

proc hexColor {icolor fcolor frac} {
    array set hexfmts {
	1 {%01x%01x%01x}
	2 {%02x%02x%02x}
	3 {%03x%03x%03x}
	4 {%04x%04x%04x}
    }
    # break the colors into rgb components
    set nh [string length $icolor]
    set nc [expr $nh/3]
    scan $icolor $hexfmts($nc) irgb(1) irgb(2) irgb(3)
    scan $fcolor $hexfmts($nc) frgb(1) frgb(2) frgb(3)
    # get the magnitude and unit vector for icolor and fcolor
    set ri [unitvector irgb]
    set rf [unitvector frgb]
    if { $ri <= 0. || $rf <= 0. } { return [format {%06x} 0] }
    array set p {1 0. 2 0. 3 0.}
    crosspro irgb frgb p
    set pmag [unitvector p]
    if { $pmag <= 0. } { return [format {%06x} 0] }
    array set ip {1 0. 2 0. 3 0.}
    crosspro p irgb ip
    set tf [expr acos([dotpro irgb frgb])]
    set t [expr $frac*$tf]
    set ct [expr cos($t)]
    set st [expr sin($t)]
    set rt [expr $ri + ($rf - $ri)*$frac]
    array set si {1 0. 2 0. 3 0.}
    scalevector [expr $rt*$ct] irgb si
    array set sf {1 0. 2 0. 3 0.}
    scalevector [expr $rt*$st] ip sf
    array set rcolor {1 0. 2 0. 3 0.}
    addvectors si sf rcolor
    return [format $hexfmts($nc) \
		[expr round($rcolor(1))] \
		[expr round($rcolor(2))] \
		[expr round($rcolor(3))]]
}

set g ""
array set eview {suu 1 sdd 1 sdu 1 sud 1 msr 1 R 0 sfm 1 nsfm 1}
set gvtitle "PB plotting"
set titlefont [list Helvetica 20 bold]
set tickfont [list Helvetica 16]
set symbolsizes {0 2 4 6 8 10 12 14 16 18 20}
set imagefmts {gif bmp eps eps2 jpg miff pdf png ppm ps ps2 rgb sgi tiff xpm}
set plotimagetype "gif"

set yctslbl "pb corrected counts"
set y2ctslbl "mon cor msr counts"
set y2Rlbl "flip-ratio"

array set bltPSopt {
    -landscape 0
    -center 1
    -maxpect 0
    -decorations 0
    -colormode color
}
proc bltPS {} {
    global bltPSopt g gwds gnum iwin spfile dsmode
    set gw $g
    if { $dsmode } { set gw $gwds }
    foreach opt $bltPSopt {
	$gw postscript configure $opt $bltPSopt($opt)
    }
    if { $dsmode } {
	set outfil $spfile
    } else {
	set outfil pbcor_group$gnum
    }
    $gw postscript output ${outfil}.ps
}
proc bltGIF {} {
    global g gwds gnum iwin spfile dsmode
    set image [image create photo]
    set gw $g
    if { $dsmode } { set gw $gwds }
    $gw snap $image
    if { $iwin < 1 } {
	set outfil $spfile
    } else {
	set outfil pbcor_group$gnum
    }
    $image write ${outfil}.gif -format GIF
    image delete $image
}

set gnuftyp gif
array set gnuPSopt {
    landscape portrait
    color color
}

proc killPlotWindow {} {
    global iwin GNU BLT gpipe gpipes g gws gvs
    if { $iwin == 1 } {
	Dialog_Warn "Don't destroy plot window 1"
	return
    }
    if { $GNU } {
	if { [info exists gpipes($iwin)] && \
		 [string length $gpipes($iwin)] > 0 } {
	    close $gpipes($iwin)
	    unset gpipes($iwin)
	    set gpipe $gpipes(1)
	}
    } elseif { $BLT } {
	# user can just close window to destroy
	if { [info exists gvs($iwin)] && \
		 [winfo exists $gvs($iwin)] } {
	    destroy $gvs($iwin)
	    if { [winfo exists $gvs($iwin)] } { tkwait window $gvs($iwin) }
	    set g $gws(1)
	}
    }
    set iwin 1
    updatePlotWindowsMenu 0
}

proc selectPlotWindow {args} {
    # trace on writing iwin
    global iwin BLT GNU gpipes gpipe gws g gvs gvtitle gnum iwin
    if { $iwin < 1 } { return }
    if { $GNU } {
	if {[info exists gpipes($iwin)] && \
		[string length $gpipes($iwin)] > 0} {
	    set gpipe $gpipes($iwin)
	}
	if {[info exists gvs($iwin)] && [winfo exists $gvs($iwin)] } {
	    wm title $gvs($iwin) "$gvtitle  group $gnum  plotwindow $iwin"
	}
    } elseif { $BLT } {
	if { [info exists gws($iwin)] && [winfo exists $gws($iwin)] } {
	    set g $gws($iwin)
	}
    }
}

proc updatePlotWindowsMenu {new} {
    global iwin miwin
    global gpipes gws gvs
    global GNU BLT

    if { $GNU } {
	set pws [array names gpipes]
    } elseif { $BLT } {
	set pws {}
	foreach pw [array names gws] {
	    if { [winfo exists $gws($pw)] } {
		lappend pws $pw
	    } else {
		unset gws($pw)
		unset gvs($pw)
	    }
	}
    }
    if { [string match new $new] } {
	set iwin 1
	while { [lsearch $pws $iwin] >= 0 } { incr iwin }
	lappend pws $iwin
    } elseif { $new > 0 } {
	set iwin $new
	if { [lsearch $pws $new] < 0 } { lappend pws $new }
    }
    set pws [lsort -integer -increasing $pws]
    $miwin delete 0 end
    #puts "updatePlotWindowsMenu pws=$pws"
    foreach pw $pws {
	$miwin add radio -label $pw -variable iwin -value $pw
    }
}

proc updatePlotSelector {} {
    # call this when switching between PB and DS modes
    global pbWidgets dsWidgets
    global psEdit psPBmenu psDSmenu
    global dsmode
    if { $dsmode } {
	# DS mode
	foreach w $pbWidgets { grid remove $w }
	foreach w $dsWidgets { grid $w }
	$psEdit configure -menu $psDSmenu
    } else {
	foreach w $dsWidgets { grid remove $w }
	foreach w $pbWidgets { grid $w }
	$psEdit configure -menu $psPBmenu
    }
}



# gws array of BLT graph widgets
array set gws {}
array set gvs {}
set psEdit ""
set psPBmenu ""
set psDSmenu ""
set pbWidgets {}
set dsWidgets {}
set dseditY ""
set dseditm ""
set dsmode 0
set miwin ""
set onePlotPerGrp 1
set rgbColors {
    red
    green
    blue
    orange
    yellow
    violet
    purple
    grey
    black
    brown
}
proc buildPlotSelector {} {
    # call this to setup plot-selector only if PLT

    # I want to be able to reconfigure the plot selector
    # for plotDS showDSplot
    # for which there are no msr
    # so build a separate edit menu that can be attached to the same
    # edit menu button
    # need to save the widgets so we can use grid remove and re grid later

    # this build preps for PB as initial state of the plot selector

    global g gv gvs gvtitle gws iwin miwin onePlotPerGrp
    global gnum snum iwin
    global eview
    global srcon srcfile isrcon sfm
    global fsSymb gfsSymb
    global pVecNames mVecNames xpts xmsr
    global titlefont tickfont iX
    global plotsymbolsize symbolsizes
    global plotconfig plotimagetype imagefmts
    global BLT GNU PLT gnuftyp gnuPSopt bltPSopt toscreen
    global saveGnuCmnds
    global yctslbl y2ctslbl y2Rlbl 
    global psEdit psPBmenu psDSmenu pbWidgets dsWidgets
    global dseditY dseditm dsmode
    global spfile spcfg
    global bltSymbolsLst rgbColors

    # gv is the frame in group-selector which will hold the plot-selector
    #set gv [toplevel .g]
    #wm title $gv $gvtitle
    #wm iconname $gv "PB plt"
    #wm iconify $gv
    #set gvs(1) $gv

    #-font *-helvetica-bold-r-normal--*-240-*
    set title [label $gv.title -text "plot selector" \
	       -font {helvetica 24 bold} -foreground #f00]

    set mb $gv

    # throw away group number adjustment in favor of Edit menu
    #set inc $mb.inc
    #set dec $mb.dec
    #set cur $mb.cur
    #set curlbl $mb.curlbl

    set edit $mb.edit
    menubutton $edit -text Edit -menu $edit.pbMenu -bg #00aa00 -relief raised \
	    -activebackground #ffaa00

    # make 2 versions of the edit menu. One for PB and one for DS
    set psEdit $edit
    set psPBmenu $edit.pbMenu
    set psDSmenu $edit.dsMenu
    foreach m [list $psPBmenu $psDSmenu] pbm {1 0} {
	menu $m -tearoff 1
	set plbl ""
	if { $pbm } { set plbl "pb-cor " }

	if { ! $pbm } {
	    $m add cascade -label "symbol" -menu $m.sym
	    set msym [menu $m.sym -tearoff 0]
	    foreach sym $bltSymbolsLst {
		$msym add radio -label $sym -variable spcfg(sym) \
			-value $sym -command "plotSymbol 1"
	    }
	    $m add cascade -label "symbol-color" -menu $m.color
	    set mcolor [menu $m.color -tearoff 0]
	    foreach colr $rgbColors {
		$mcolor add radio -label $colr -variable spcfg(color) \
			-value $colr -command "plotSymbolColor 1"
	    }
	}

	$m add cascade -label "${plbl}symbol-size" -menu $m.pbsub
	set m2 [menu $m.pbsub -tearoff 0]
	foreach size $symbolsizes {
	    if { $pbm } {
		$m2 add radio -label $size -variable plotsymbolsize(pbcor) \
		    -value $size -command "plotSymbolSize 1"
	    } else {
		$m2 add radio -label $size -variable spcfg(pixels) \
		    -value $size -command "plotSymbolSize 1"
	    }
	}
	if { $pbm } {
	    $m add cascade -label "msr symbol-size" -menu $m.msrsub
	    set m3 [menu $m.msrsub -tearoff 0]
	    foreach size $symbolsizes {
		$m3 add radio -label $size -variable plotsymbolsize(msr) \
		    -value $size -command "plotSymbolSize 1"
	    }
	    $m add cascade -label "R symbol-size" -menu $m.rsub
	    set m4 [menu $m.rsub -tearoff 0]
	    foreach size $symbolsizes {
		$m4 add radio -label $size -variable plotsymbolsize(R) \
		    -value $size -command "plotSymbolSize 1"
	    }
	}
	$m add separator
	if { $pbm } {
	    $m add check -label "${plbl}errbars" -variable plotconfig(errpb) \
		-command "plotErrbars 1"
	    $m add check -label "msr errbars" -variable plotconfig(errmsr) \
		-command "plotErrbars 1"
	} else {
	    $m add check -label "${plbl}errbars" -variable spcfg(errs) \
		-command "plotErrbars 1"
	}
	
	$m add separator

	if { $pbm } {
	    if { $BLT } {
		$m add check -label zoom -variable plotconfig(zoom) \
		    -command plotConfig
		$m add check -label crosshairs -variable plotconfig(crosshairs) \
		    -command plotConfig
	    }
	    $m add check -label logscale -variable plotconfig(logscale) \
		-command "plotLogScale 1"
	} else {
	    if { $BLT } {
		$m add check -label zoom -variable spcfg(zoom) \
		    -command plotConfig
		$m add check -label crosshairs -variable spcfg(crosshairs) \
		    -command plotConfig
	    }
	    $m add check -label logscale -variable spcfg(logscale) \
		-command "plotLogScale 1"
	}

	$m add separator

	set mpb [menu $m.pbmin -tearoff 0]
	if { $pbm } {
	    $m add cascade -label "${plbl}axis min" -menu $m.pbmin
	    set mmsr [menu $m.msrmin -tearoff 0]
	    $m add cascade -label "msr axis min" -menu $m.msrmin
	} else {
	    $m add cascade -label "y axis min" -menu $m.pbmin
	}
	foreach mt {0 1 auto} {
	    if { $pbm } {
		$mpb add radio -label $mt -variable plotconfig(ymin) -value $mt \
		    -command "axismin y"
		$mmsr add radio -label $mt -variable plotconfig(y2min) \
		    -value $mt -command "axismin y2"
	    } else {
		$mpb add radio -label $mt -variable spcfg(ymin) -value $mt \
		    -command "axismin y"
	    }
	}
	$m add command -label "other axis scaling" -command axisScale
	
	$m add separator
	# section for detaching and saving plots
	if { $BLT } {
	    if { ! $pbm } {
		set savplbl "save plot as filename"
		set savn ds
	    } else {
		set savplbl "save plot as pbcor_group#"
		set savn pbcor_group
	    }
	    if { [file executable /usr/bin/import] } {
		# imageMagick utility only easily managed from linux
		$m add cascade -label "${savplbl}.type" \
		    -menu $m.psub
		set m5 [menu $m.psub -tearoff 0]
		foreach fmt $imagefmts {
		    $m5 add radio -label $fmt -variable plotimagetype \
			-value $fmt -command "saveImage $gv $savn 1"
		}
		$m add command -label "detach static copy of plot" \
		    -command "detachImage $gv $savn 1"
	    } else {
		# use BLT postscript output
		$m add command -label "${savplbl}.ps" \
		    -command bltPS
		$m add cascade -label "postscript options" -menu $m.psub
		set m5 [menu $m.psub -tearoff 0]
		foreach opt {landscape center maxpect decorations} \
		    lbl {landscape center fullpage decorations} {
			$m5 add check -label $lbl -variable bltPSopt(-$opt)
		    }
		$m5 add radio -label color -variable bltPSopt(-colormode) \
		    -value color
		$m5 add radio -label greyscale -variable bltPSopt(-colormode) \
		    -value greyscale
		$m add command -label "${savplbl}.gif" \
		    -command bltGIF
	    }
	} else {
	    # gnuplot hardcopy terminals gif jpeg postscript
	    $m add radio -label "plot to screen" -variable toscreen -value 1
	    $m add radio -label "plot to file" -variable toscreen -value 0 \
		-command showPlot
	    $m add cascade -label "file types" -menu $m.psub
	    set m5 [menu $m.psub -tearoff 0]
	    foreach ftyp {gif jpeg postscript} {
		$m5 add radio -label $ftyp -variable gnuftyp -value $ftyp
	    }
	    $m5 add radio -label psColor -variable gnuPSopt(color) -value color
	    $m5 add radio -label psGrey -variable gnuPSopt(color) \
		-value monochrome
	    $m5 add radio -label psLandscape -variable gnuPSopt(landscape) \
		-value landscape
	    $m5 add radio -label psPortrait -variable gnuPSopt(landscape) \
		-value portrait
	    $m add separator
	    $m add check -label "save gnuplot cmnds to gnuplotCMNDS file" \
		-variable saveGnuCmnds
	}
	
	if { $pbm } {
	    # manage multiple plot windows
	    $m add separator
	    $m add command -label "new plot window" \
		-command "updatePlotWindowsMenu new"
	    $m add cascade -label "select plot window" -menu $m.iwin
	    set miwin [menu $m.iwin -tearoff 0]
	    updatePlotWindowsMenu 0
	    $m add command -label "kill selected window" -command killPlotWindow
	    $m add check -label "one plot window per group" \
		-variable onePlotPerGrp
	}
    }

    # add buttons to select which fs to view
    # select pbcorrect and/or raw data (from .m)
    # use symbols uptri=du downtri=ud filcirc=uu opencirc=dd
    # color coding: graded=rawsets coded by set   green=correct
    # as in text window select which sets to view

    set pltwinL [label $mb.pltwinL -text pltwin]
    set pltwin [label $mb.pltwin -textvariable iwin -width 2 \
		    -relief sunken -bd 1 -bg white]

    # switching modes reconfigures the plotselector
    set dsplt [checkbutton $mb.dsmode -variable dsmode -text "pltSelFile" \
		   -anchor e -command updatePlotSelector]

    # for DS use a menubutton to allow user to select
    # which y-element to edit
    # the dseditm menu gets radio populated with the DS y columnlabels
    set elemYL [label $mb.elemYL -textvariable dseditY -bg white -width 20]
    set elemY $mb.elemY
    menubutton $elemY -text EditY -menu $elemY.menu -bg #00aa00 \
	-relief raised -activebackground #ffaa00
    set dseditm [menu $elemY.menu -tearoff 0]


    set uu $mb.uu
    set du $mb.du
    set ud $mb.ud
    set dd $mb.dd

    set srcview $mb.srcview
    set sfsrcview $mb.sfsrcview
    set nsfsrcview $mb.nsfsrcview
    set rfview $mb.rfview

    # change sn stuff to a single menubutton with
    # a menu of check foreach set file
    set sf $mb.sf

    # save the names of widgets that get toggled when plotDS
    set pbWidgets \
	[list $uu $du $ud $dd $srcview $sfsrcview $nsfsrcview $rfview $sf]
    set dsWidgets [list $elemY $elemYL]

    #set sninc $sf.sninc
    #set sndec $sf.sndec
    #set sncur $sf.sncur
    #set sncurlbl $sf.sncurlbl
    #set snview $sf.snview
    #set snfile $sf.file
    #set snfileL $sf.fileL


    #button $inc -command incrGroup -image uparrow \
    #	-background green -activebackground yellow
    #button $dec -command decrGroup -image downarrow \
    #	-background green -activebackground yellow
    #label  $cur -textvariable gnum -width 3 -relief sunken -bd 1 -bg white
    #label  $curlbl -text group

    checkbutton $uu -variable eview(suu) -text uu -anchor e \
	-command {viewElem 1}
    checkbutton $dd -variable eview(sdd) -text dd -anchor e \
	-command {viewElem 1}
    checkbutton $du -variable eview(sdu) -text du -anchor e \
	-command {viewElem 1}
    checkbutton $ud -variable eview(sud) -text ud -anchor e \
	-command {viewElem 1}

    checkbutton $srcview -variable eview(msr) -text msr -anchor e \
	-command {viewElem 1}

    checkbutton $nsfsrcview -variable eview(nsfm) -text nsfmsr -anchor e \
	-command {viewMsrFS 1}
    checkbutton $sfsrcview -variable eview(sfm) -text sfmsr -anchor e \
	-command {viewMsrFS 1}


    checkbutton $rfview -variable eview(R) -text R -anchor e \
	-command {viewElem 1}

    menubutton $sf -text msrSet -menu $sf.menu -bg #00aa00 -relief raised \
	    -activebackground #ffaa00
    set sfm [menu $sf.menu -tearoff 0]


    #frame $sf
    #button $sninc -command incrSn -image uparrow \
    #	-background green -activebackground yellow
    #button $sndec -command decrSn -image downarrow \
    #	-background green -activebackground yellow
    #label  $sncur -textvariable snum -width 3 -relief sunken -bd 1 -bg white
    #label  $sncurlbl -text msrset
    #checkbutton $snview -variable isrcon -command viewMsrSet
    #entry $snfile -textvariable srcfile -width 20
    #label $snfileL -text file
    # could put the srcfile name in wm border
    # replace the group selector widgets with an Edit menu
    # for adjusting symbol size, crosshairs, logY

    #grid $curlbl $cur $inc $dec
    grid $title -row 0 -column 0 -sticky w -ipadx 5
    grid $edit -row 0 -column 1 -ipadx 5 -padx 10
    grid $pltwinL -row 0 -column 2
    grid $pltwin -row 0 -column 3
    grid $dsplt -row 0 -column 4 -padx 10

    grid columnconfigure $mb 5 -minsize 20
    grid $elemY -row 0 -column 6
    grid $elemYL -row 0 -column 7
    # remove elemY and elemYL for pb configuration
    grid remove $elemY
    grid remove $elemYL

    grid $uu -row 0 -column 8
    grid $du -row 0 -column 9
    grid $ud -row 0 -column 10
    grid $dd -row 0 -column 11
    grid $srcview -row 0 -column 12
    grid $nsfsrcview -row 0 -column 13
    grid $sfsrcview -row 0 -column 14

    grid $sf -row 0 -column 15
    grid $rfview -row 0 -column 16 


    #grid $sf -row 0 -column 10 -sticky e -padx 15
    #grid $sninc $sndec $sncur $sncurlbl $snview $snfileL $snfile

    #pack $title $mb -side left -fill x


    #if { $BLT }
    #set g [graph $gv.g]
    #set gvs(1) $gv
    #set gws(1) $g
    #pack $g -side top -fill both -expand true
    #bltGraphSetup

    #setPlotGeo
}

# keep track of the groupnumber in each plotwindow
# and inverse plotwindownumber for each group
# iwin gets set when changing group
# if another group is plotted already in iwin via plotgrpiwin(iwin)
# need to make new plot and update plotgrpiwin and plotwinigrp
# if plot window is killed update plotgrpiwin array
# ie user will PLOTsoln
# unless the one-plot-window-per-grp flag is set onePlotPerGrp=1
# GlobalOptions
# then find the next open window for the group and set iwin to that.
# plotconfiguration stuff is in
# plotconfig  plotsymbolsize eview isrcplot which are all saved per group
# while plotconfiguration stuff can be changed runtime
# when using showPlot it should also be applied.

proc axismin {ax} {
    # just use applyAxScl which handles BLT and GNU
    global plotconfig spcfg dsmode
    if { $dsmode } {
	if { [string match auto $spcfg(${ax}min)] } {
	    set spcfg(${ax}automin) 1
	    set spcfg(${ax}min) 0
	} else {
	    set spcfg(${ax}automin) 0
	}
    } else {
	if { [string match auto $plotconfig(${ax}min)] } {
	    set plotconfig(${ax}automin) 1
	    set plotconfig(${ax}min) 0
	} else {
	    set plotconfig(${ax}automin) 0
	}
    }
    applyAxScl 1
}
set ax ""
proc axisScale {} {
    global ax iwin plotconfig spcfg dsmode spfile
    # present user with top level
    set ax .ax
    if { ! [winfo exists $ax] } {
	toplevel $ax

	foreach {c lbl} {1 min 2 max 3 automin 4 automax} {
	    grid [label $ax.$lbl -text $lbl] -row 0 -column $c
	}

	foreach {r lbl} {1 x 2 y 3 y2} {
	    grid [label $ax.$lbl -text $lbl] -row $r -column 0
	}
	foreach {r rlbl} {1 x 2 y 3 y2} {
	    foreach {c clbl} {1 min 2 max} {
		grid [entry $ax.$rlbl$clbl \
			  -textvariable plotconfig($rlbl$clbl) \
			  -width 12 -relief sunken] \
		    -row $r -column $c
	    }
	    foreach {c clbl} {3 automin 4 automax} {
		grid [checkbutton $ax.$rlbl$clbl \
			  -variable plotconfig($rlbl$clbl)] \
		    -row $r -column $c
	    }
	}

	grid [button $ax.apply -text apply -command "applyAxScl 1" \
		 -background green -activebackground orange] -row 0 -column 0
    }
    raise $ax

    if { $dsmode } {
	wm title $ax "dataset file $spfile axis limits"
	foreach {r rlbl} {1 x 2 y 3 y2} {
	    foreach {c clbl} {1 min 2 max} {
		$ax.$rlbl$clbl configure -textvariable spcfg($rlbl$clbl)
	    }
	    foreach {c clbl} {3 automin 4 automax} {
		$ax.$rlbl$clbl configure -variable spcfg($rlbl$clbl)
	    }
	}
    } else {
	wm title $ax "plotwindow $iwin axis limits"
	foreach {r rlbl} {1 x 2 y 3 y2} {
	    foreach {c clbl} {1 min 2 max} {
		$ax.$rlbl$clbl configure -textvariable plotconfig($rlbl$clbl)
	    }
	    foreach {c clbl} {3 automin 4 automax} {
		$ax.$rlbl$clbl configure -variable plotconfig($rlbl$clbl)
	    }
	}
    }
}

proc applyAxScl {update} {
    # use this to configure plot ranges
    global g gwds plotconfig spcfg dsmode
    global PLT BLT GNU iwin
    if { ! $PLT } { return }
    set axnlst [list x y y2]
    set gw $g
    if { $dsmode } { set gw $gwds }
    foreach ax {x y y2} axn $axnlst {
	set mn ""
	set mx ""
	if { $GNU } {
	    set mn *
	    set mx *
	}
	if { $dsmode } {
	    if { ! $spcfg(${axn}automin) && \
		     [scan $spcfg(${axn}min) {%g} mnval] > 0 } {
		set mn $mnval
	    } else {
		set spcfg(${axn}automin) 1
	    }
	    if { ! $spcfg(${axn}automax) && \
		     [scan $spcfg(${axn}max) {%g} mxval] > 0 } {
		set mx $mxval
	    } else {
		set spcfg(${axn}automax) 1
	    }
	} else {
	    if { ! $plotconfig(${axn}automin) && \
		     [scan $plotconfig(${axn}min) {%g} mnval] > 0 } {
		set mn $mnval
	    } else {
		set plotconfig(${axn}automin) 1
	    }
	    if { ! $plotconfig(${axn}automax) && \
		     [scan $plotconfig(${axn}max) {%g} mxval] > 0 } {
		set mx $mxval
	    } else {
		set plotconfig(${axn}automax) 1
	    }
	}

	if { $BLT } {
	    catch {$gw axis configure $ax -min $mn -max $mx}
	} else {
	    if { $dsmode } {
		gdspipeWrite "set ${ax}range \[$mn:$mx\]"
		if { $update } { gdspipeWrite "replot" }
	    } else {
		gpipeWrite "set ${ax}range \[$mn:$mx\]"
		if { $update } { gpipeWrite "replot" }
	    }
	}
    }
}


proc checkStyles {} {
    # BLT PB only
    global g stylst rstylst ngrpsets
    set penlst [$g pen names]
    for {set i [expr [llength $stylst] - 1]} {$i >= 0} {incr i -1} {
	set sty [lindex $stylst $i]
	set p [lindex $sty 0]
	if { [lsearch -exact $penlst $p] < 0 } {
	    set stylst [lreplace $stylst $i $i]
	}
    }
    for {set i [expr [llength $rstylst] - 1]} {$i >= 0} {incr i -1} {
	set sty [lindex $rstylst $i]
	set p [lindex $sty 0]
	if { [lsearch -exact $penlst $p] < 0 } {
	    set rstylst [lreplace $rstylst $i $i]
	}
    }
}
proc bltGraphSetup {} {
    # this is showPlot for BLT
    global g gnum fsSymb
    global yctslbl y2ctslbl y2Rlbl xaxlbl titlefont tickfont
    global plotbg
    global pVecNames mVecNames pDatNames mDatNames
    global stylst rstylst iWts ngrpsets isrcplot eview pencolor
    global iX
    
    # make sure the graph pens and elements exist
    checkPens
    if { ! [$g element exists suu] } {
	$g element create suu -symbol $fsSymb(1) -color #0000ff -linewidth 0
	$g element create sdd -symbol $fsSymb(2) -color #0000ff -linewidth 0
	$g element create sdu -symbol $fsSymb(3) -color #0000ff -linewidth 0
	$g element create sud -symbol $fsSymb(4) -color #0000ff -linewidth 0
	$g element create msr -linewidth 0
	$g element create R -symbol cross -color red -linewidth 0
    }
    
    # config the axes except for X
    
    $g axis configure y -title $yctslbl -color blue \
	-titlecolor blue -titlefont $titlefont -tickfont $tickfont
    $g axis configure y2 -title $y2ctslbl -hide 0 -color red \
	-titlecolor red -titlefont $titlefont -tickfont $tickfont
    if { [llength [$g axis names rf]] < 1 } {
	$g axis create rf -title $y2Rlbl -color red \
	    -titlecolor red -titlefont $titlefont -tickfont $tickfont -hide 1
    }

    # make sure the vectors for this group exist
    if { ! [vectorExists suu$gnum] } {
        foreach vnam $pVecNames { vector create $vnam$gnum }
	foreach vnam $mVecNames { vector create $vnam$gnum }
	vector create xpts$gnum xmsr$gnum rf$gnum
    }
    # assign the group data lists to the vectors
    foreach vnam $pVecNames pdat $pDatNames {
	global $pdat
	eval $vnam$gnum set $$pdat
    }
    foreach vnam $mVecNames mdat $mDatNames {
	global $mdat
	eval $vnam$gnum set $$mdat
    }
    global Xpts Xmsr Rf
    xpts$gnum set $Xpts
    xmsr$gnum set $Xmsr
    rf$gnum set $Rf

    # assign the vectors to their graph elements
    $g element configure suu \
	-xdata ::xpts$gnum -ydata ::suu$gnum -yerror ::euu$gnum
    $g element configure sdd \
	-xdata ::xpts$gnum -ydata ::sdd$gnum -yerror ::edd$gnum
    $g element configure sdu \
	-xdata ::xpts$gnum -ydata ::sdu$gnum -yerror ::edu$gnum
    $g element configure sud \
	-xdata ::xpts$gnum -ydata ::sud$gnum -yerror ::eud$gnum
    $g element configure msr \
	-xdata ::xmsr$gnum -ydata ::cc$gnum  -yerror ::ccerr$gnum \
	-mapy y2 -color red
    $g element configure R -xdata ::xmsr$gnum -ydata ::rf$gnum -mapy rf
    # initially rf axis hidden
    $g y2axis use y2

    set plotbg [$g cget -plotbackground]

    checkPens
    checkStyles
    $g element configure msr -styles $stylst -weights $iWts \
	-hide [oneC $eview(msr)]
    $g element configure R -styles $rstylst -weights $iWts \
	-hide [oneC $eview(R)]
    $g element configure suu -hide [oneC $eview(suu)]
    $g element configure sdd -hide [oneC $eview(sdd)]
    $g element configure sdu -hide [oneC $eview(sdu)]
    $g element configure sud -hide [oneC $eview(sud)]


    for {set i 1} {$i <= $ngrpsets} {incr i} {
	for {set j 1} {$j <= 4} {incr j} {
	    if { $isrcplot($i) } {
		set pc $pencolor($j$i)
	    } else {
		set pc $plotbg
	    }
	    set pc [string trim $pc \#]
	    $g pen configure p$j$i -color #$pc -symbol $fsSymb($j)
	    $g pen configure r$j$i -color #$pc -symbol cross
	}
    }
    $g axis configure x -title $xaxlbl($iX) -color black \
	-titlecolor black -titlefont $titlefont -tickfont $tickfont    
    $g legend configure -font {helvetica 10}
    
    # blt enhancements
    Blt_ZoomStack $g
    Blt_Crosshairs $g
    Blt_ActiveLegend $g
    #Blt_ClosestPoint $g
    Blt_PrintKey $g    
}
proc bltDSGraphSetup {} {
    # this is showDSplot for BLT DSmode
    global gwds
    global spcfg
    global titlefont tickfont
    global plotbg
    global readlines readfinished
    
    # make sure the graph elements and vectors exist and assign the DS
    # data lists to the vectors and then the vectors to the element data

    set elems [$gwds element names]
    for {set i [expr $spcfg(ny) + 1]} {$i <= [llength $elems]} {incr i} {
	$gwds element delete y$i
    }

    if { ! [vectorExists xdspts] } { vector create xdspts }
    pipeWrite "c dl $spcfg(xcol)"
    vwait readfinished
    # get rid of the leading label
    set Xds [lreplace [lindex $readlines 0] 0 0]
    xdspts set $Xds

    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	if { ! [$gwds element exists y$i] } { $gwds element create y$i }
	if [regexp {(.*)-fill} $spcfg(y${i}-sym) mtch sym] {
	    set spcfg(y${i}-sym) $sym
	    set spcfg(y${i}-fill) defcolor
	}
	#puts "bltDS y$i symbol= $spcfg(y${i}-sym)"
	$gwds element configure y$i -symbol $spcfg(y${i}-sym) \
	    -color $spcfg(y${i}-color) -linewidth 0 -hide 0 \
	    -label $spcfg(y${i}-lbl) -fill $spcfg(y${i}-fill)
	if { ! [vectorExists y$i] } { vector create y$i }
	if { ! [vectorExists e$i] } { vector create e$i }
	pipeWrite "c dl $spcfg(y${i}-col)"
	vwait readfinished
	# get rid of the leading label
	set ydatalst [lreplace [lindex $readlines 0] 0 0]
	y$i set $ydatalst
	pipeWrite "c dl $spcfg(e${i}-col)"
	vwait readfinished
	# get rid of the leading label
	set edatalst [lreplace [lindex $readlines 0] 0 0]
	e$i set $edatalst
	$gwds element configure y$i -xdata ::xdspts -ydata ::y$i -yerror ::e$i
    }

    # config the axes
    
    $gwds axis configure y -title $spcfg(ylbl) -color blue \
	-titlecolor blue -titlefont $titlefont -tickfont $tickfont
    $gwds axis configure y2 -hide 1 -color blue -tickfont $tickfont
    $gwds axis configure x -title $spcfg(xlbl) -color black \
	-titlecolor black -titlefont $titlefont -tickfont $tickfont    
    $gwds legend configure -font {helvetica 10}

    # blt enhancements
    Blt_ZoomStack $gwds
    Blt_Crosshairs $gwds
    Blt_ActiveLegend $gwds
    #Blt_ClosestPoint $gwds
    Blt_PrintKey $gwds
}


proc saveImage {w lbl gn} {
    # only works for BLT widget and with ImageMagick utilities available
    global plotimagetype
    global g gnum iwin spfile
    if { ![winfo exists $w] } { return }
    if { $iwin < 1 && [string match $lbl ds] } { set lbl $spfile }
    if { $iwin > 0 && $gn } { append lbl $gnum }
    set fname "$lbl.$plotimagetype"
    set wid [winfo id $w]
    if { [catch {exec import -window $wid $fname}] } {
	Dialog_Warn "Failed to save plotimage as $fname"
    }
}
proc detachImage {w lbl gn} {
    # this only works of the ImageMagick utilities import and display
    # are available and w is a tk widget
    # for the gnuplot window need another way
    # I'm not sure if there are xwindow ids under IMagick on Windows
    global g gnum iwin spfile
    if { ![winfo exists $w] } { return }
    if { $iwin < 1 && [string match $lbl ds] } { set lbl $spfile }
    if { $iwin > 0 && $gn } { append lbl $gnum }
    set fname "$lbl.gif"
    set wid [winfo id $w]
    if { [catch {exec import -window $wid $fname}] } {
	Dialog_Warn "Failed to save image as $fname"
    }
    if { [catch {exec display $fname &}] } {
	Dialog_Warn "Failed to detach image"
    }
}

array set spcfg {
    zoom 1
    crosshairs 1
    logscale 0
    ymin 0
    ymax ""
    yautomin 1
    yautomax 1
    y2min 0
    y2max ""
    y2automin 1
    y2automax 1
    xmin ""
    xmax ""
    xautomin 1
    xautomax 1
    errs 1
    pixels 6
    xlbl X
    ylbl Y
}

array set plotconfig {
    zoom 1
    crosshairs 1
    logscale 0
    ymin 0
    ymax ""
    yautomin 1
    yautomax 1
    y2min 0
    y2max ""
    y2automin 1
    y2automax 1
    xmin ""
    xmax ""
    xautomin 1
    xautomax 1
    errmsr 1
    errpb 1
}
proc plotConfig {} {
    # use this for BLT comp= zoom crosshairs
    global g gwds plotconfig spcfg dsmode
    global BLT
    if { ! $BLT } { return }
    set gw $g
    if { $dsmode } { set gw $gwds }
    foreach comp {zoom crosshairs} {
	if { $dsmode } {
	    set on $spcfg($comp)
	} else {
	    set on $plotconfig($comp)
	}
	set bndtgs [bindtags $gw]
	set zi [lsearch $bndtgs $comp*]
	if { $on } {
	    if { $zi >= 0 } { continue }
	    set bndtgs [linsert $bndtgs 0 $comp-$gw]
	} else {
	    if { $zi < 0 } { continue }
	    set bndtgs [lreplace $bndtgs $zi $zi]
	}
	bindtags $gw $bndtgs
    }
}

# default symbol sizes
array set plotsymbolsize {pbcor 6 msr 6 R 6}
if { $GNU } {
    array set plotsymbolsize {pbcor 2 msr 2 R 2}
}
set pbcorelems {suu sdd sdu sud}
proc plotSymbolSize {update} {
    # use to change symbol size for comp=pbcor, msr, R
    # NB pen onfig doesnt get a graph update, thus the element config
    global g gwds plotsymbolsize ngrpsets pbcorelems spcfg dsmode dseditY
    global BLT GNU

    # for dsmode start from plot selector configuration
    if { $dsmode } {
	for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	    if { ![string match $dseditY all] && \
		     ![string match $dseditY $spcfg(y${i}-lbl)] } { continue }
	    set spcfg(y${i}-pixels) $spcfg(pixels)
	}
    }
    if { $BLT } {
	if { $dsmode } {
	    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
		if [info exists spcfg(y${i}-pixels)] {
		    $gwds element configure y$i -pixels $spcfg(y${i}-pixels)
		} elseif [info exists spcfg(pixels)] {
		    set spcfg(y${i}-pixels) $spcfg(pixels)
		    $gwds element configure y$i -pixels $spcfg(y${i}-pixels)
		}
	    }
	} else {
	    foreach elem $pbcorelems {
		$g element configure $elem -pixels $plotsymbolsize(pbcor)
	    }
	    foreach penname [$g pen names p*] {
		$g pen configure $penname -pixels $plotsymbolsize(msr)
	    }
	    $g element configure msr -pixels $plotsymbolsize(msr)
	    foreach penname [$g pen names r*] {
		$g pen configure $penname -pixels $plotsymbolsize(R)
	    }
	    $g element configure R -pixels $plotsymbolsize(R)
	}
    } elseif { $update } {
	# for GNU just replot if update requested
	if { $dsmode } { showDSplot } else { showPlot }
    }
}

proc plotSymbol {update} {
    # DSmode only
    # symbol types are fixed for group plot
    # but might as well allow user to configure symbol types
    # for DS plot
    global gwds spcfg dsmode dseditY
    global BLT GNU bltSymbolIndex gnuSymbols
    if { $dsmode } {
	for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	    if { ![string match $dseditY all] && \
		     ![string match $dseditY $spcfg(y${i}-lbl)] } { continue }
	    set spcfg(y${i}-gsym) $gnuSymbols($bltSymbolIndex($spcfg(sym)))
	    set spcfg(y${i}-sym) $spcfg(sym)
	    set spcfg(y${i}-fill) ""
	    if [regexp {(.*)-fill} $spcfg(y${i}-sym) mtch sym] {
		set spcfg(y${i}-sym) $sym
		set spcfg(y${i}-fill) defcolor
	    }	    
	}
    }
    if { $BLT } {
	if { $dsmode } {
	    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
		if { ![string match $dseditY all] && \
			 ![string match $dseditY $spcfg(y${i}-lbl)] } {
		    continue
		}
		$gwds element configure y$i \
		    -symbol $spcfg(y${i}-sym) -fill $spcfg(y${i}-fill)
	    }
	}
    } elseif { $update } {
	if { $dsmode } { showDSplot } else { showPlot }
    }
}
proc plotSymbolColor {update} {
    # DSmode only
    # symbol types are fixed for group plot
    # but might as well allow user to configure symbol types
    # for DS plot
    global gwds spcfg dsmode dseditY
    global BLT GNU
    if { $dsmode } {
	for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	    if { ![string match $dseditY all] && \
		     ![string match $dseditY $spcfg(y${i}-lbl)] } { continue }
	    set spcfg(y${i}-color) $spcfg(color)
	}
    }
    if { $BLT } {
	if { $dsmode } {
	    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
		if { ![string match $dseditY all] && \
			 ![string match $dseditY $spcfg(y${i}-lbl)] } {
		    continue
		}
		$gwds element configure y$i -color $spcfg(y${i}-color)
	    }
	}
    } elseif { $update } {
	if { $dsmode } { showDSplot } else { showPlot }
    }
}

proc plotErrbars {update} {
    global g gwds plotconfig pbcorelems spcfg dsmode dseditY
    global BLT GNU iwin
    if { $dsmode } {
	for {set i 1} {$i <= $spcfg(ny)} {incr i} {
	    if { ![string match $dseditY all] && \
		     ![string match $dseditY $spcfg(y${i}-lbl)] } { continue }
	    set spcfg(y${i}-errs) $spcfg(errs)
	}
    }
    if { $BLT } {
	if { $dsmode } {
	    for {set i 1} {$i <= $spcfg(ny)} {incr i} {
		if { ![string match $dseditY all] && \
			 ![string match $dseditY $spcfg(y${i}-lbl)] } {
		    continue
		}
		set eb none
		if { ([info exists spcfg(y${i}-errs)] && $spcfg(y${i}-errs)) \
			 || ([info exists spcfg(errs)] && $spcfg(errs)) } {
		    set eb both
		}
		$gwds element configure y$i -showerrorbars $eb
	    }
	} else {
	    set eb none
	    if { $plotconfig(errpb) } { set eb both }
	    foreach elem $pbcorelems {
		$g element configure $elem -showerrorbars $eb
	    }
	    viewMsrSet 0
	}
    } elseif { $update } {
	if { $dsmode } { showDSplot } else { showPlot }
    }
}

proc plotLogScale {update} {
    # use to set logscale for pbcor and msr
    global g gwds plotconfig spcfg dsmode
    global BLT GNU
    if { $dsmode } {
	set log $spcfg(logscale)
	set gw $gwds
    } else {
	set log $plotconfig(logscale)
	set gw $g
    }
    if { $log } {
	if { $BLT } {
	    $gw axis configure y y2 -logscale true -loose 1
	} elseif { $GNU } {
	    if { $dsmode } {
		gdspipeWrite "set logscale y"
		gdspipeWrite "set logscale y2"
		if { $update } { gdspipeWrite "replot" }
	    } else {
		gpipeWrite "set logscale y"
		gpipeWrite "set logscale y2"
		if { $update } { gpipeWrite "replot" }
	    }
	}
    } else {
	if { $BLT } {
	    $gw axis configure y y2 -logscale false -loose 0
	} elseif { $GNU } {
	    if { $dsmode } {
		gdspipeWrite "unset logscale y"
		gdspipeWrite "unset logscale y2"
		if { $update } { gdspipeWrite "replot" }
	    } else {
		gpipeWrite "unset logscale y"
		gpipeWrite "unset logscale y2"
		if { $update } { gpipeWrite "replot" }
	    }
	}
    }
}

##################### auto group procs ########################333

proc fileListToSetList {fls} {
    # put datasets corresponding to files in fls into setlist 1
    # fls list should have short filenames from lbxs(0)
    # so dataset numbers can be recovered from dsFile
    global readfinished readlines
    global currentDir dsFile

    # delete any previous setlists
    pipeWrite "s da"
    vwait readfinished

    # clear the first setlist. This also makes it current setlist
    pipeWrite "s1 lc"
    vwait readfinished
    # put the files in this setlist
    set nfls [llength $fls]
    for {set i 0} {$i < $nfls} {incr i} {
	# reconstruct the full src filename
	set fl [lindex $fls $i]
	#set fl [file join $currentDir [lindex $fls $i]]

	# get dataset index from input filename srcin
	# because libDF adds .x ext to read files even though not
	# written to disk because cols have been modified
	# yet they are in the sellectlist with orig srcname = srcin from d ne
	# why not store the dsnumber for every file that gets put into
	# lbxs(0) so we never have to look it up or worry about src name
	# in fact this was already done and the dsnumbers are
	# in dsFile(shortFileName) array
	
	#pipeWrite "d ni $fl"
	#vwait readfinished
	# readline 0 has dataset inum filename= file
	#set dsn 0
	#if { [llength $readlines] > 0 } {
	#    set dsn [lindex [lindex $readlines 0] 1]
	#}

	set dsn $dsFile($fl)
	if { $dsn < 1 } {
	    # this shouldnt happen, but check anyway
	    Dialog_Warn "failed to get dataset number for file: $fl"
	    return
	}

	pipeWrite "s l+ $dsn"
	vwait readfinished
    }
}

proc autoGroup {all} {
    # exec auto grouping button binding
    global currentgroup samepts exclmatch gnum
    global stndTolArr stndTolNames autoStndTol
    global readfinished readlines
    global lbxs

    # get curselection from lbxs(0) to check for seed
    set seedlst [$lbxs(0) curselection]
    set nseed [llength $seedlst]
    set seed 0

    # if currentgroup put just gnum files into setlist
    # else put all files into setlist
    if { $all } { set ig 0 } else { set ig $gnum }
    set fls [getFilesInGroup $ig]
    set nfls [llength $fls]
    if { $nfls < 2 } {
	Dialog_Warn "Cant auto-group with less than 2 files in source group."
	return
    }
    fileListToSetList $fls

    for {set i 0} {$i < $nfls} {incr i} {
	set fl [lindex $fls $i]
	# if not yet seeded see if we can use this fl
	if { $nseed > 0 && $seed < 1 } {
	    if { [set is [lsearch -exact $seedlst $fl]] >= 0 } {
		set seed [expr $i + 1]
	    }
	}
    }
    # update the tolerances and save in autoStndTol

    setStndTolerances
    array set autoStndTol [array get stndTolArr]

    # ready to match
    set p ""
    if { $samepts < 1 } { set p p }
    set x ""
    if { $exclmatch < 1 } { set x x }
    # exec the setlist match command
    pipeWrite "s M$p$x $seed"
    vwait readfinished

    #puts "match complete. print lines"
    #for {set i 0} {$i < [llength $readlines]} {incr i} {
    #	puts [lindex $readlines $i]
    # }

    # get the resulting linked setlist numbers
    pipeWrite "s1 ll"
    vwait readfinished
    if { [llength $readlines] < 2 } {
	Dialog_Warn "Match failed to produce new groupings."
	return
    }
    set ll [lindex $readlines 1]
    # first element is the source setlist number (1)
    # the link setlistnumbers follow and will contain the source as well
    # so remove any entries equal to 1
    while { [set is [lsearch -exact $ll 1]] >= 0 } {
	set ll [lreplace $ll $is $is]
    }

    # now where to put these setlists
    # user may have already setup groups with these setlist numbers
    # use nextEmptyGroup
    set eg 0
    for {set i 0} {$i < [llength $ll]} {incr i} {
	set sl [lindex $ll $i]
	set eg [nextEmptyGroup $eg]
	set gnum $eg
	setGroup
	# get the filenames in this setlist as they appear in the listbox
	# s li gives the srcin. Here we want the diskfile names
	pipeWrite "s$sl l"
	vwait readfinished
	# start line 4 ilst match iset fname
	# put the data for these files into lbxs with gnum=$eg
	# first make the list of files from listbox
	getLbxFiles
	set failed ""
	for {set j 3} {$j < [llength $readlines]} {incr j} {
	    set rl [lindex $readlines $j]
	    set slf [file tail [lindex $rl 3]]
	    set data [getLbxData $slf]
	    if { [llength $data] < 1 } {
		append failed "$slf\n"
		continue
	    }
	    set data [lreplace $data 1 1 $eg]
	    appendLbxData $data
	}
	if { [string length $failed] > 0 } {
	    Dialog_Warn "The following files from group $gnum were not found in the listbox. This is a bug.\n$failed"
	}
	setGroup
    }
}

set toltype ""
set tolerance 0
set tolmenu ""
set currentgroup 1
set samepts 0
set exclmatch 1

proc buildAutoGroup {} {
    # how to handle tolerances
    # punt and just use stnd col names
    # do an initial load of stndTol(stndName) = stndTol

    # this builds in $bot frame

    global bframes
    global stndTolNames tolerance toltype samepts exclmatch
    global tolmenu
    global currentgroup

    set f $bframes(3)
    clearGrid $f

    set title $f.title
    set exe $f.exe
    set optcb $f.opt
    set onecb $f.one

    set execur $exe.cur
    set exeall $exe.all

    set tolType $exe.t
    set tolValu $exe.v
    set tolTypeL $exe.tl
    set tolValuL $exe.vl




    if { ! [winfo exists $title] } {
	#-font *-helvetica-bold-r-normal--*-240-*
	label $title -text "auto make groups" \
	    -font {helvetica 24 bold} -foreground #f00
	frame $exe
	button $execur -text "AUTO-GROUPcur" \
	    -command "autoGroup 0" \
	    -background green -activebackground yellow -width 14
	button $exeall -text "AUTO-GROUPall" \
	    -command "autoGroup 0" \
	    -background green -activebackground yellow -width 14
	label $tolTypeL -text "data column type"
	label $tolValuL -text "match tolerance"

	set tolmenu [eval tk_optionMenu $tolType toltype $stndTolNames]
	# configure each menu entry to write the tol value to the
	# entry widget when it is selected
	set toltyp [lindex $stndTolNames 0]
	foreach name $stndTolNames {
	    $tolmenu entryconfigure $name -command "setStndTol $name"
	}
	entry $tolValu -textvariable tolerance -width 10
	bind $tolValu <Leave> getStndTol


	checkbutton $optcb -variable samepts \
	    -text "files must have all same data points"
	checkbutton $onecb -variable exclmatch \
	    -text "each file can only be in one group"
    }

    grid $title -sticky ew
    grid $exe -sticky ew
    grid $optcb -sticky ew
    grid $onecb -sticky ew

    grid $execur $exeall -sticky ew
    grid $tolTypeL $tolType -sticky ew
    grid $tolValuL $tolValu -sticky ew
}

################### file selector procs ###################

proc showSelfile {} {
    global mp
    if { $mp >= 4 } { showText }
}


proc buildFileSelector {} {
    global sellectlist
    global bframes

    # this builds in top left frame
    set f $bframes(1)
    clearGrid $f
    if { ![winfo exists $f.top] } {
	# create the fileselector widgets if they dont exist
	set top [frame $f.top]
	#-font *-helvetica-bold-r-normal--*-240-*
	set title [label $top.title -text "data selector" \
		       -font {helvetica 24 bold} -foreground #f00]
	set do [frame $top.do]

	set bot [frame $f.bot]
	set sellectlist [listbox $bot.list -width 33 -relief raised -bd 2 \
			     -yscrollcommand "$bot.listscroll set" \
			     -selectmode extended]
	bind $sellectlist <ButtonRelease-1>  {+ showSelfile; plotDS}
	bind $sellectlist <ButtonRelease-3> {tk::ListboxSelectAll %W}
	set sellectscroll [scrollbar $bot.listscroll \
			       -command "$bot.list yview"]

	set filterentry [entry $do.filterentry -relief sunken -bd 2 \
			     -width 11 -textvariable fileFilter]
	bind $filterentry <Return> "doFilter $sellectlist"
	#bind $filterentry <Enter> "doFilter $sellectlist"
	bind $filterentry <Tab> "doFilter $sellectlist"
	set bFILTER [button $do.filter -text FILTER \
			 -command "doFilter $sellectlist" \
			 -background green -activebackground yellow]
	set bALL [button $do.all -text "sel ALL" \
		      -command "tk::ListboxSelectAll $sellectlist" \
			 -background green -activebackground yellow]
	set bREAD [button $do.read -text "READ sel" \
		       -command "readFiles" \
		       -background green -activebackground yellow]
    } else {
	set top $f.top
	set title $top.title
	set do $top.do
	set bot $f.bot
	set sellectlist $bot.list
	set sellectscroll $bot.listscroll
	set filterentry $do.filterentry
	set bFILTER $do.filter
	set bALL $do.all
	set bREAD $do.read
    }
    grid $top -sticky ew
    grid $bot -sticky ew

    grid $title -sticky ew
    grid $do -sticky ew
    grid $bREAD $bALL $bFILTER $filterentry -sticky ew

    grid $sellectscroll $sellectlist -sticky news
}


set cellFiles {}
proc getDirectoryCellFiles {} {
    global readfinished
    global readlines
    global cellFiles cellFile
    global currentDir

    set cellFiles {}
    pipeWrite "p H *"
    vwait readfinished
    # starting at 4th line: ifil match filename
    set nf [expr [llength $readlines] - 3]
    if { $nf <= 0 } {
	Dialog_Warn "NO cell files in current directory: $currentDir"
	lappend cellFiles none
	return
    }
    for {set i 3} {$i < [llength $readlines]} {incr i} {
	set ln [lindex $readlines $i]
	if { [llength $ln] < 3 } { continue }
	set cellFile [lindex $ln 2]
	lappend cellFiles $cellFile
    }
}

array set colds {}
array set colLbls {}
array set colRanges {}
proc readCurrentColAssigns {} {
    global colds colLbls colRanges
    global readfinished readlines

    # get the columns assigned for current dataset
    # and record labels and column data ranges

    pipeWrite "c al"
    vwait readfinished

    # each line after the first has
    # ctyp index dscol datasetColLabel
    # need the dscol for each of our DS columns to get stats
    
    set nc [llength $readlines]
    arrayClear colLbls
    arrayClear colds
    arrayClear colRanges
    for {set i 1} {$i < $nc} {incr i} {
	set cspec [lindex $readlines $i]
	if { [llength $cspec] < 4 } { continue }
	set ctyp [lindex $cspec 0]
	set cindx [lindex $cspec 1]
	set dscol [lindex $cspec 2]
	set lbl [lindex $cspec 3]
	if { $cindx > 0 } { append ctyp $cindx }
	set colds($ctyp) $dscol
	set colLbls($dscol) $lbl
    }
    foreach ctyp [array names colds] {
	pipeWrite "c s $colds($ctyp)"
	vwait readfinished
	# line 2 has col lbl min max
	set stat [lindex $readlines 1]
	set min [lindex $stat 2]
	set max [lindex $stat 3]
	if { $max > $min } {
	    set colRanges($ctyp) "$min $max"
	} else {
	    set colRanges($ctyp) $max
	}
    }
    return [expr $nc - 1]
}

# columns for which we get statistics to fill the lbxs
set colsrch {
y1 x1 x2 x3 En fs H V Temp scl
}
set colsrchS {
x1 x2 x3 En H V Temp scl
}
set lbxnames {fnam gnum npts y1 x1 x2 x3 En Fix fs H V Temp scl ref fbg}
array set lbxfmts {
    0 ""
    1 ""
    2 ""
    3 {%g}
    4 ""
    5 ""
    6 ""
    7 ""
    8 ""
    9 ""
    10 ""
    11 ""
    12 {%g}
    13 ""
    14 ""
    15 ""
}

# remember the indices for temperature and fast-background
set templbx 12
set fbglbx 15
set nlbxs [llength $lbxnames]
#set lbxnames [join "fnam gnum npts $colsrch"]
set gnum 1
set lbxpad 1

array set dsFile {}

# flipper-state correspondence to print in readfile table
array set fspb {1 offoff 2 ONON 3 ONoff 4 offON}
array set fbgFile {}

proc readFile {fil} {
    global dsFile srcfiles currentDir readfinished errlines nerr warnlines nwarn
    global readlines readfinished nread errs
    # read fil into dataset and return a list
    # iset npts ncol ny nas df src fshort fname
    # user may reset the global errs string to "" before call

    set fshort [file tail $fil]
    # no need to use full filenames as pbscript will know the cd
    # this means src will get stored as short and can
    # be easily compared to filenames in listboxes
    set fname [file join $currentDir $fshort]
    #puts $fname
    # arrange for pbscript to read this file into the next dataset
    pipeWrite "r+ $fshort"
    vwait readfinished
    # check for read errors
    if { $nerr > 0 } {
	append errs "$errlines"
	return {}
    }
    if { $nwarn > 0 } {
	append errs "$warnlines"
    }
    # scan 2nd readline for npts. 
    # always have the leading  * indicator of current set
    # save the set number for this file
    set rl [lindex $readlines 1]
    set nscan [scan $rl {%[*]%d %d %d %d %d %d %s} \
		   star iset npts ncol ny nas df src]
    if { $nscan < 8 } {
	append errs "FILEREAD ERROR: failed get npts in $fshort from $rl\n"
	return {}
    }
    # use the filename of diskfile after writing
    # version number may change if root already exists
    # can use df=0 to indicate that this is a new dataset
    # not yet on disk
    set src [file tail $src]
    if { $df == 0 || ! [string match $src $fshort] } {
	# this means columns hd etc were added so write to disk
	# write the updated file to disk
	pipeWrite "d$iset w"
	vwait readfinished
	# get back the actual disk file name
	
	pipeWrite "d$iset"
	vwait readfinished
	set rl [lindex $readlines 1]
	set nscan [scan $rl {%[*]%d %d %d %d %d %d %s} \
		       star iset npts ncol ny nas df src]
	if { $nscan < 8 } {
	    append errs "FILEREAD ERROR: failed diskfilename for $fshort from $rl\n"
	    return {}
	}
	set fshort [file tail $src]
	set fname [file join $currentDir $fshort]
    }
    set dsFile($fshort) $iset
    set srcfiles($iset) $fshort
    # return list
    return [list $iset $npts $ncol $ny $nas $df $src $fshort $fname]
}

set errs ""
proc readFiles {args} {
    # read to current group files in [lindex $args 0]
    # else with no filelist in args get the selection from
    # the data selector listbox sellectlist
    # args is used for session recovery
    # everything in the lbxs can be recovered by reading the files
    # except for fbg, which is saved in group in array fbgFile(shortfilename)
    global sellectlist
    global currentDir
    global readfinished
    global readlines nread
    global errlines nerr
    global warnlines nwarn
    global colsrch
    global lbxnames lbxlabels
    global gnum
    global lbxs
    global lbxpad
    global fspb
    global dsFile fbgFile
    global colds colRanges isfop
    global errs

    set fillst {}
    if { [llength $args] < 1 } {
	set flstind [$sellectlist curselection]
	if { [llength $flstind] < 1 } {
	    Dialog_Warn "NO files selected.\n"
	    return
	}
	foreach i [lsort -integer $flstind] {
	    set fshort [$sellectlist get $i]
	    lappend fillst $fshort
	}
    } else {
	set fillst [lindex $args 0]
    }
    if { [llength $fillst] < 1 } {
	Dialog_Warn "NO files in filelist for reading.\n"
	return
    }

    set values(fnam) {}
    set values(gnum) {}
    set values(npts) {}
    set values(Fix) {}
    set values(ref) {}
    set values(fbg) {}

    foreach i $colsrch { set values($i) {} }

    set nf 0
    set errs ""
    foreach fil $fillst {
	set ret [readFile $fil]
	# continue on fatal read error
	if { [llength $ret] < 1 } { continue }
	#returns [list $iset $npts $ncol $ny $nas $df $src $fshort $fname]
	set iset [lindex $ret 0]
	set npts [lindex $ret 1]
	set ncol [lindex $ret 2]
	set df [lindex $ret 5]
	set fshort [lindex $ret 7]
	set fname [lindex $ret 8]

	# get Fix from Fixn hd
	pipeWrite "d d Fixn"
	vwait readfinished
	# 3rd line: symbol value found
	if { [llength $readlines] < 3 } {
	    append errs "FILEREAD ERROR: HD Fixn not found in $fshort\n"
	}
	set hd [lindex $readlines 2]
	set fnd [lindex $hd 1]
	if { $fnd == 1 } {
	    lappend values(Fix) "Ef"
	} elseif { $fnd == 0 } {
	    lappend values(Fix) "Ei"
	} else {
	    lappend values(Fix) "NA"
	}

	pipeWrite "d d refn"
	vwait readfinished
	# 3rd line: symbol value found
	if { [llength $readlines] < 3 } {
	    append errs "FILEREAD WARNING: HD refn not found in $fshort\n"
	}
	set hd [lindex $readlines 2]
	set fnd [lindex $hd 1]
	if { $fnd == 2 } {
	    lappend values(ref) "mon"
	} elseif { $fnd == 1 } {
	    lappend values(ref) "time"
	} else {
	    lappend values(ref) "NA"
	}

	# get the assigned data columns and data ranges
	set nc [readCurrentColAssigns]

	if { $nc < 10 } {
	    append errs "FILEREAD WARNING: assigned less than 10 columns in $fname\n"
	    continue
	}

	incr nf

	lappend values(fnam) $fshort
	lappend values(gnum) $gnum
	lappend values(npts) $npts

	if { ! [info exists fbgFile($fshort)] } {
	    set fbgFile($fshort) 0
	}

	lappend values(fbg) $fbgFile($fshort)

	# now collect the stats for each column.
	# dont have to report missing as they will show blank in the table
	foreach i $colsrch {
	    if { [info exists colds($i)] && $colds($i) > 0 } {
		if { ! [string compare $i "fs"] } {
		    # handle flipper states, find all used
		    pipeWrite "c sd $colds($i)"
		    vwait readfinished
		    # line 2 has all distinct values
		    # use fspb to translate and order offoff onoff offon onon
		    set stat [lindex $readlines 1]
		    set fstates {}
		    foreach fsti {1 3 4 2} {
			if { [lsearch -exact $stat $fsti] < 0 } { continue }
			lappend fstates $fspb($fsti)
		    }
		    if { [llength $fstates] == 4 } { set fstates ALL }
		    lappend values($i) $fstates
		} else {
		    lappend values($i) $colRanges($i)
		}
	    } else {
		lappend values($i) ""
	    }
	}
    }

    # now append the values to the lbxs
    for {set i 0} {$i < [llength $lbxnames]} {incr i} {
	set vlst $values([lindex $lbxnames $i])
	for {set j 0} {$j < $nf} {incr j} {
	    $lbxs($i) insert end [lindex $vlst $j]
	}
	setLbxWidth $lbxs($i) [lindex $lbxlabels $i]
    }
    hiliteActiveGroup
    if { [string length $errs] > 0 } {
	Dialog_Warn $errs
    }
}

proc setLbxWidth {lbw lbl} {
    global lbxpad
    set lbxw [maxListItemLength [$lbw get 0 end]]
    set lblw [string length $lbl]
    if { $lblw > $lbxw } { set lbxw $lblw }
    set wid [expr $lbxw + $lbxpad]
    $lbw configure -width $wid
}

proc lbxIndicesForGroup {ig} {
    global lbxs
    set gvals [$lbxs(1) get 0 end]
    set indices {}
    for {set i 0} {$i < [llength $gvals]} {incr i} {
	if { [lindex $gvals $i] == $ig } { lappend indices $i }
    }
    return indices
}

# use the files listbox hilite green to indicate active group
# use the group number listbox to hilite green for calc
# and yellow for fop
proc hiliteActiveGroup {} {
    global gnum calc calcHilite
    global lbxs
    set gvals [$lbxs(1) get 0 end]
    for {set i 0} {$i < [llength $gvals]} {incr i} {
	if { [lindex $gvals $i] == $gnum } {
	    $lbxs(0) itemconfigure $i -bg green
	    #$lbxs(1) itemconfigure $i -bg $calcHilite($calc)
	} else {
	    $lbxs(0) itemconfigure $i -bg white
	}
    }
}
array set calcHilite {0 white 1 green}
proc hiliteCalcGroup {} {
    # apply hilite to current gnum
    global gnum calc
    global lbxs
    global calcHilite
    set gvals [$lbxs(1) get 0 end]
    for {set i 0} {$i < [llength $gvals]} {incr i} {
	if { [lindex $gvals $i] == $gnum } {
	    $lbxs(1) itemconfigure $i -bg $calcHilite($calc)
	}
    }
}

proc hiliteGroups {} {
    global calcgroup fopgroup lbxs
    set gns [$lbxs(1) get 0 end]
    for {set i 0} {$i < [llength $gns]} {incr i} {
	$lbxs(1) itemconfigure $i -bg white
	set gn [lindex $gns $i]
	if { [info exists calcgroup($gn)] && $calcgroup($gn) } {
	    $lbxs(1) itemconfigure $i -bg green
	} elseif { [info exists fopgroup($gn)] && $fopgroup($gn) } {
	    $lbxs(1) itemconfigure $i -bg red
	}
    }
}
array set fopHilite {0 white 1 red}
proc hiliteFopGroup {} {
    # apply hilite to current gnum
    global gnum calc isfop
    global lbxs
    global fopHilite
    set gvals [$lbxs(1) get 0 end]
    for {set i 0} {$i < [llength $gvals]} {incr i} {
	if { [lindex $gvals $i] == $gnum } {
	    $lbxs(1) itemconfigure $i -bg $fopHilite($isfop)
	}
    }
}


proc maxListItemLength {lst} {
    set maxl 0
    foreach i $lst {
	set ilen [string length $i]
	if { $ilen > $maxl } { set maxl $ilen }
    }
    return $maxl
}

proc doFilter {w} {
    global fileFilter
    $w delete 0 end
    set flst [glob -nocomplain $fileFilter]
    foreach i [lsort $flst] {
	if { [file isfile $i] == 0 || [file executable $i] } { continue }
	$w insert end $i
    }
}

##################### main group selector procs ##########

# start working on scrolledlistboxes for grid geo manager
# NB to prevent resizing a grid element set  -weight 0
# parent frame should be grid managed
# as scrolledlistboxes may add scrollbar using grid
# NB grid cant resize columns arbitrarily so may have to use -weight
# write separate proc to resize lbxs according to width of contents


array set bld {}
array set bldtext {1 "file-select" 2 "prep-solve" 3 "make-group"}
proc buildButtonConfig {ib} {
    # call with ib==1 when bld(1) is pressed
    # call with ib==2 when bld(2) is pressed
    # call with ib==0 when heightlines changes
    # store the button actions in bexe(1) and bexe(2)
    # the action values can be 1 2 or 3 indicating
    # display frame f/g/a
    # when ib==1,2 which frame f/g/a gets unviewed
    # isview(1-3)=0/1 are view states of f/g/a
    # if a bld button is removed its bexe(ib)=0
    # so we can tell if two, one or no buttons are displayed

    # there are two build buttons
    # if heightlines == 1 both buttons are displayed
    # with the names of the two non-viewed frame types
    # if heightlines == 2 only one button is displayed
    # if heightlines >= 3 no buttons are needed

    # keep track of which frames, top/mid/bot are viewed
    # with isview(f) isview(g) isview(a)
    # to remove one of the frames from view, grid remove thatframe
    # and to restore it, grid thatframe

    # with two buttons, only one of isview is active
    # swap the viewed one with the one requested by the button

    # with one button, two of isview are active.
    # unview the lowest priority window f/g/a is priority order


    global heightlines
    global bld bldtext bframes
    global bexe isview
    global lbxs nlbxs

    if { $ib > 0 } {
	# find the lowest priority viewed window
	for {set i 3} {$i >= 1} {incr i -1} {
	    if { $isview($i) } {
		set irem $i
		break
	    }
	}
	# clear irem frame
	if { [winfo ismapped $bframes($irem)] } { grid remove $bframes($irem) }
	# view the one requested by the button
	set iview $bexe($ib)
	if { ! [winfo ismapped $bframes($iview)] } { grid $bframes($iview) }
	# reconfigure the ib button to view the one removed
	$bld($ib) configure -text $bldtext($irem)

	set isview($irem) 0
	set isview($iview) 1
	set bexe($ib) $irem
	return
    }

    # now handle heightlines
    # by starting from a standard configuration
    for {set i 1} {$i <= 2} {incr i} {
	if { [winfo ismapped $bld($i)] } { grid remove $bld($i) }
	set bexe($i) 0
    }

    set nview $heightlines
    if { $nview > 3 } { set nview 3 }
    for {set i 1} {$i <= $nview} {incr i} {
	if {! [winfo ismapped $bframes($i)] } { grid $bframes($i) }
	set isview($i) 1
    }
    for {set i [expr $nview + 1] ; set nb 1} {$i <= 3} {incr i ; incr nb} {
	if { [winfo ismapped $bframes($i)] } { grid remove $bframes($i) }
	set isview($i) 0
	grid $bld($nb)
	$bld($nb) configure -text $bldtext($i)
	set bexe($nb) $i
    } 

    set nl 10
    if { $heightlines == 2 } { set nl 28 }
    if { $heightlines >= 3 } { set nl 37 }
    for {set i 0} {$i < $nlbxs} {incr i} {
	$lbxs($i) configure -height $nl
    }
}

set lbxlabels {}
set curnpts ""
set curnmsr ""

proc nptsFormat {args} {
    global curnpts npts
    set nc [string length $npts]
    if { $nc > 0 } { $curnpts configure -width $nc }
}
proc nmsrFormat {args} {
    global curnmsr nmsr
    set nc [string length $nmsr]
    if { $nc > 0 } { $curnmsr configure -width $nc }
}

set grpr ""
array set lbxsr {}
set lbxrlabels {}

proc setTempTolFromRanges {} {
    global lbxsr stndTolArr templbx toltype tolerance
    set trs [$lbxsr($templbx) get 0 end]
    set n [llength $trs]
    if { $n < 1 } { return }
    set tr [lindex $trs 0]
    for {set i 1} {$i < $n} {incr i} {
	set tr [updateNumRange [lindex $trs $i] $tr]
    }
    if { [llength $tr] < 1 } {
	set stndTolArr(Temp) 1.
    } else {
	set tmn [lindex $tr 0]
	set tmx [lindex $tr 1]
	set stndTolArr(Temp) [expr 1.1*abs($tmx - $tmn)]
    }
    if { [string match $toltype Temp] } { set tolerance $stndTolArr(Temp) }
}

proc buildGroupRanges {} {
    # group selector-like listbox with ranges for all columns per group
    global grpr
    global gselSB lbxsr lbxrlabels
    set grpr [toplevel .grpr]
    set top [frame $grpr.top]
    set bot [frame $grpr.bot]
    grid $top
    grid $bot
    set  title $top.title
    #-font *-helvetica-bold-r-normal--*-240-*
    label $title -text "group ranges" \
	-font {helvetica 24 bold} \
	-foreground #f00
    set ttol [button $top.tol -text "set Temp tol from max range + 10%" \
		  -bg green -activebackground yellow \
		  -command setTempTolFromRanges -width 28]
    grid $title $ttol

    set lbxrlabels {group nfile npts cts Qx Qy Qz E  Fix  fs H V T scl ref fbg}
    set widths     {3      3     3   12  12 12 12 12 3    12 6 6 6 12  4   3}
    set sorts  {
	-integer
	-integer -ascii
	-ascii -ascii -ascii -ascii -ascii -ascii -ascii
	-real -real -real -real
	-ascii -real
    }
      
    array set lbxsr \
	[scrolledlistboxes $bot $lbxrlabels $widths $sorts $gselSB {}]
    wm title $grpr "PB group ranges"
}

proc showGroupRanges {} {
    global lbxs lbxsr lbxrlabels grpr
    # go thru lbxs and find ranges for each group
    # gran(grpnum) { }
    if { ! [winfo exists $grpr] } { buildGroupRanges }
    wm deiconify $grpr
    raise $grpr
    foreach b [array names lbxsr] { $lbxsr($b) delete 0 end }
    set gns {}
    set gnx [lsort -integer -increasing [$lbxs(1) get 0 end]]
    set n [llength $gnx]
    if { $n < 1 } { return }
    set gn [lindex $gnx 0]
    lappend gns $gn
    for {set i 1} {$i < $n} {incr i} {
	set newgn [lindex $gnx $i]
	if { $newgn == $gn } { continue }
	lappend gns $newgn
	set gn $newgn
    }
    array set gran {}
    set numcols {2 3 4 5 6 7  10 11 12 13 15}
    set txtcols {8 9 14}
    for {set i 0} {$i < [llength $gns]} {incr i} {
	set gn [lindex $gns $i]
	# collect ranges for this group lbxs 2-15
	# 2-7 are numeric ranges
	# 8=Fix Ei Ef 9=fs 10-13 numeric 14=reftype (mon or time) 15 numer fbg
	unset gran
	set nf 0
	for {set j 0} {$j < $n} {incr j} {
	    if { [$lbxs(1) get $j] != $gn } { continue }
	    # first do the numeric column range updates
	    incr nf
	    foreach nc $numcols {
		set nr [$lbxs($nc) get $j]
		if { ! [info exists gran($nc)] } {
		    set gran($nc) $nr
		} else {
		    set gran($nc) [updateNumRange $nr $gran($nc)]
		}
	    }
	    # now handle 7=Fix 8=fs 13=reftype
	    foreach nc $txtcols {
		set nr [$lbxs($nc) get $j]
		if { ! [info exists gran($nc)] } {
		    set gran($nc) $nr
		} else {
		    set gran($nc) [updateTxtRange $nr $gran($nc)]
		}
	    }
	}
	# insert the gran data into lbxsr
	$lbxsr(0) insert end $gn
	$lbxsr(1) insert end $nf
	for {set j 2} {$j < 16} {incr j} {
	    $lbxsr($j) insert end $gran($j)
	}
	for {set j 0} {$j < 16} {incr j} {
	    setLbxWidth $lbxsr($j) [lindex $lbxrlabels $j]
	}
    }
}

proc updateTxtRange {newr oldr} {
    set rlst $oldr
    foreach n $newr {
	if { [lsearch -exact $rlst $n] < 0 } { lappend rlst $n }
    }
    return $rlst
}

proc updateNumRange {newr oldr} {
    if { [llength $oldr] > 1 } {
	set oldmin [lindex $oldr 0]
	set oldmax [lindex $oldr 1]
    } else {
	set oldmin [lindex $oldr 0]
	set oldmax $oldmin
    }
    if { [llength $newr] > 1 } {
	set newmin [lindex $newr 0]
	set newmax [lindex $newr 1]
    } else {
	set newmin [lindex $newr 0]
	set newmax $newmin
    }
    set mn $oldmin
    if { $newmin < $oldmin } { set mn $newmin }
    set mx $oldmax
    if { $newmax > $oldmax } { set mx $newmax }
    if { $mx > $mn } {
	return [list $mn $mx]
    } else {
	return $mx
    }
}

proc buildGroupSelector {} {
    # this only gets built once when heightlines=1
    global tv imagefmts
    global gselSB
    global lbxlabels
    global bld
    global npts nmsr calc gnum
    global curnpts curnmsr
    global right gv
    global lbxs

    set g $right
    clearGrid $g

    set top $g.top
    set bot $g.bot
    # add plot frame holder below group-selector
    set gv $g.gv

    set title $top.title
    set bld(1) $top.bld1
    set bld(2) $top.bld2
    set fvw $top.fvw
    set edit $top.edit
    set curgrp $top.cur
    set inc $curgrp.inc
    set dec $curgrp.dec
    set cur $curgrp.cur
    set curlbl $curgrp.lbl
    set curnpts $curgrp.npts
    set curnptsL $curgrp.nptsL
    set curnmsr $curgrp.nmsr
    set curnmsrL $curgrp.nmsrL

    if { ![winfo exists $top] } {
	# create the groupselector widgets if they dont exist
	frame $top
	frame $bot
	frame $gv
	frame $curgrp
	button $inc -command incrGroup -image uparrow \
	    -background green -activebackground yellow
	button $dec -command decrGroup -image downarrow \
	    -background green -activebackground yellow
	label  $cur -textvariable gnum -width 3 -relief sunken -bd 1 -bg white
	label  $curlbl -text group
	label  $curnpts -textvariable npts \
	    -width 3 -relief sunken -bd 1 -bg white
	label  $curnptsL -text npts
	label  $curnmsr -textvariable nmsr \
	    -width 3 -relief sunken -bd 1 -bg white
	label  $curnmsrL -text nmsr
	trace variable npts w nptsFormat
	trace variable nmsr w nmsrFormat
	#trace variable gnum w groupSrcFiles

	#-font *-helvetica-bold-r-normal--*-240-*
	label $title -text "group selector" \
	    -font {helvetica 24 bold} \
	    -foreground #f00
	button $bld(1) -bg green -activebackground yellow \
	    -command "buildButtonConfig 1" -width 8
	button $bld(2) -bg green -activebackground yellow \
	    -command "buildButtonConfig 2" -width 8
	button $fvw -bg green -activebackground yellow \
	    -command showText -text "fileViewer" -width 8
	menubutton $edit -text Edit -menu $edit.menu -bg #00aa00 -relief raised \
	    -activebackground #ffaa00
	set m [menu $edit.menu -tearoff 1]
	$m add command -label "Subtract background" -command subtractBG

	$m add separator
	$m add command -label "Show group ranges" -command showGroupRanges

	$m add separator
	$m add command -label "Solved group file operations" \
	    -command groupResultFops

	$m add separator
	$m add command -label "Remove current group" -command "delGrps 0"
	$m add command -label "Remove selected groups" -command "delGrps 1"
	$m add command -label "Remove unselected groups" -command "delGrps 2"
	$m add command -label "Remove selected files from their group" \
	    -command delFiles
	$m add command -label "Move selected files to current group" \
	    -command moveFiles

	$m add separator
	$m add command -label "More lines in group selector" \
	    -command "setGrpLines 1" 
	$m add command -label "Less lines in group selector" \
	    -command "setGrpLines -1" 

        if { [file executable /usr/bin/import] } {
	    # imageMagick under Linux
	    $m add separator
	    $m add cascade -label "save main window as pbmain.type" \
		-menu $m.psub
	    set m1 [menu $m.psub -tearoff 0]
	    foreach fmt $imagefmts {
	     $m1 add radio -label $fmt -variable plotimagetype -value $fmt \
		 -command "saveImage . pbcor_main 0"
	    }
	    $m add command -label "detach static copy of main window" \
		-command "detachImage . pbcor_main 0"
	    $m add cascade -label "save text window as pbtext.type" \
		-menu $m.tsub
	    set m2 [menu $m.tsub -tearoff 0]
	    foreach fmt $imagefmts {
	     $m2 add radio -label $fmt -variable plotimagetype -value $fmt \
		 -command "saveImage $tv pbcor_text 0"
	    }
	    $m add command -label "detach static copy of text window" \
		-command "detachImage $tv pbcor_text 0"
        } 

	$m add separator
	$m add command -label "Save session" -command "saveSession 1"
	$m add command -label "Save and quit" -command "saveSession 2"
	$m add command -label "Quit" -command "saveSession 0"
    }

    grid $top -sticky news
    grid $bot
    grid $gv -sticky news

    grid $title -row 0 -column 0 -ipadx 5 -sticky w
    grid $edit -row 0 -column 1 -ipadx 5 -sticky w
    grid $curgrp -row 0 -column 2 -ipadx 5
    grid $bld(1) -row 0 -column 3
    grid $bld(2) -row 0 -column 4
    grid $fvw -row 0 -column 5


    grid $curlbl $cur $inc $dec $curnptsL $curnpts $curnmsrL $curnmsr

    set lbxlabels {filename group npts cts Qx Qy Qz E  Fix fs H V T scl ref fbg}
    set widths    {16       3     3    12  12 12 12 12 3   12 6 6 6 12  4   3}
    set sorts  {
	-ascii
	-integer -integer
	-ascii -ascii -ascii -ascii -ascii -ascii -ascii
	-real -real -real -real
	-ascii -real
    }
    array set lbxs \
	[scrolledlistboxes $bot $lbxlabels $widths $sorts $gselSB {}]
    bind $lbxs(0) <ButtonRelease-1>  {+ showSelfile; plotDS}
    # do the config with startup heightlines=10
    buildButtonConfig 0
}

proc deleteLbxEntries {egs} {
    global lbxs nlbxs
    set egs [lsort -integer -decreasing $egs]
    for {set i 0} {$i < [llength $egs]} {incr i} {
	set ig [lindex $egs $i]
	for {set j 0} {$j < $nlbxs} {incr j} {
	    $lbxs($j) delete $ig
	}
    }
}
proc deleteGroup {ig} {
    # make sure to del decreasing order so indices remain valid
    set egs [lsort -decreasing [getEntriesForGroup $ig]]
    deleteLbxEntries $egs
}
set widgetsel ""
proc filesSelected {is} {
    global lbxs lbxi nlbxs widgetsel
    # for selection options check the widget that owns PRIMARY sel
    # some selection is required

    set sel [selectionOwned 0 2 0]
    if { [llength $sel] < 1 } {
	Dialog_Warn "Selection is not in one of the group selector listboxes.\n"
	return {}
    }

    set fsels {}
    set lbxi {}
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set isel [lindex $sel $i]
	set fsel [$lbxs(0) get $isel]
	if {[lsearch -exact $fsels $fsel] < 0} {
	    lappend fsels $fsel
	    lappend lbxi $isel
	}
    }
    if { $is } { return $fsels }

    # get unselected files
    set allfs [$lbxs(0) get 0 end]
    set fnosels {}
    set lbxi {}
    for {set i 0} {$i < [llength $allfs]} {incr i} {
	set fsel [lindex $allfs $i]
	if {[lsearch -exact $fnosels $fsel] < 0 && 
	    [lsearch -exact $fsels $fsel] < 0 } {
	    lappend fnosels $fsel
	    lappend lbxi $i
	}
    }
    return $fnosels
}

proc filesInGroupsSelected {is} {
    global widgetsel lbxi
    set gsel [groupsSelected $is]
    set ngsel [llength $gsel]
    if { $ngsel < 1 } { return {} }
    set filst {}
    set lbxiloc {}
    foreach grp $gsel {
	set fils [getFilesInGroup $grp]
	for {set i 0} {$i < [llength $fils]} {incr i} {
	    if {[lsearch -exact $filst $fil] < 0} {
		lappend filst $fil
		lappend lbxiloc [lindex $lbxi $i]
	    }
	}
    }
    set lbxi $lbxiloc
    return $filst
}

proc selectionOwned {sellst nlbx getitems} {
    # return selection if
    # sellst==1 if owned by sellectlist
    # nlbx ==1 if owned by lbxs(0)
    # nlbx > 1 if owned by any lbxs
    # if getitems returns list of lbx selected items
    # else returns selection indices

    global sellectlist lbxs nlbxs
    if { $sellst == 0 && $nlbxs == 0 } { return {} }
    if [catch {selection own} w] { return {} }
    if { $sellst &&  [string match $w $sellectlist] } {
	return [getLbxItems $w [$w curselection] $getitems]
    }
    if { $nlbx == 1 &&  [string match $w $lbxs(0)] } {
	return [getLbxItems $w [$w curselection] $getitems]
    }
    if { $nlbx != 2 } { return {} }
    set selok 0
    for {set i 0} {$i < $nlbxs} {incr i} {
	if { [string match $w $lbxs($i)] } {
	    set selok 1
	    break
	}
    }
    if { ! $selok } {
	return {}
    } else {
	return [getLbxItems $w [$w curselection] $getitems]
    }
}
proc getLbxItems {w indices getitems} {
    if { ! $getitems } { return $indices }
    set items {}
    foreach i $indices { lappend items [$w get $i] }
    return $items
}

proc groupsSelected {is} {
    global lbxs nlbxs widgetsel
    # for selection options check the widget that owns PRIMARY sel
    # some selection is required

    set sel [selectionOwned 0 2 0]
    if { [llength $sel] < 1 } {
	Dialog_Warn "Selection is not in one of the group selector listboxes.\n"
	return {}
    }

    set gsels {}
    for {set i 0} {$i < [llength $sel]} {incr i} {
	set isel [lindex $sel $i]
	set gsel [$lbxs(1) get $isel]
	if {[lsearch -exact $gsels $gsel] < 0} {lappend gsels $gsel}
    }
    if { $is } { return $gsels }

    # get unselected
    set allgs [$lbxs(1) get 0 end]
    set gnosels {}
    for {set i 0} {$i < [llength $allgs]} {incr i} {
	set gsel [lindex $allgs $i]
	if {[lsearch -exact $gnosels $gsel] < 0 && 
	    [lsearch -exact $gsels $gsel] < 0 } {lappend gnosels $gsel}
    }
    return $gnosels
}
proc delGrps {id} {
    global lbxs nlbxs
    global gnum
    if { $id == 0 } {
	# delete current group
	deleteGroup $gnum
	return
    }

    if { $id == 1 } {
	set gsels [groupsSelected 1]
	for {set i 0} {$i < [llength $gsels]} {incr i} {
	    deleteGroup [lindex $gsels $i]
	}
	return
    }
    # must be delete non selected

    set gnosels [groupsSelected 0]
    for {set i 0} {$i < [llength $gnosels]} {incr i} {
	deleteGroup [lindex $gnosels $i]
    }
}

proc delFiles {} {
    global lbxs
    set sel [$lbxs(0) curselection]
    if { [llength $sel] < 1 } {
	Dialog_Warn "Selection must be in File listbox to delete files.\n"
	return
    }
    deleteLbxEntries $sel
}
proc moveFiles {} {
    global lbxs
    global gnum
    set sel [$lbxs(0) curselection]
    if { [llength $sel] < 1 } {
	Dialog_Warn "Selection must be in File listbox to delete files.\n"
	return
    }

    set flst [getFilesInGroup $gnum]
    set already {}
    foreach ent $sel {
	set fil [$lbxs(0) get $ent]
	if { [lsearch -exact $flst $fil] >= 0 } {
	    lappend already $fil
	    continue
	}
	$lbxs(1) delete $ent
	$lbxs(1) insert $ent $gnum
    }
    if { [llength $already] > 0 } {
	Dialog_Warn "Files already in current group:\n$already"
    }
    hiliteActiveGroup
}

proc setGrpLines {inc} {
    global heightlines
    if { $inc > 0 } {
	if { $heightlines < 3 } { incr heightlines }
    } elseif { $inc < 0 } {
	if { $heightlines > 1 } { incr heightlines -1 }
    }
    buildButtonConfig 0
}

############### group file ops - subtract or average group solutions ######

set fop ""
set lbxolabels {}
set grpsop ""
set grpsub ""
array set lbxos {}
proc buildFileOps {} {
    global fop
    global gselSB
    global lbxos lbxolabels
    global grpsop grpsub
    global npts nmsr calc gnum

    set fop .fop
    set top $fop.top
    set bot $fop.bot

    set title $top.title

    set edit $top.edit
    set grpsop $top.gent
    set grpsub $top.gsub
    set grpres $top.gres
    set execg $top.exec

    if { ![winfo exists $fop] } {
	toplevel $fop
	wm title $fop "PB fop"

	frame $top
	frame $bot

	label $title -text "subtract or average groups" \
	    -font {helvetica 24 bold} \
	    -foreground #f00


	menubutton $edit -text Edit -menu $edit.menu -bg #00aa00 \
	    -relief raised -activebackground #ffaa00 -width 10
	set m [menu $edit.menu -tearoff 1]
	$m add command -label "New groups subtraction" -command newSubtraction
	$m add command -label "New groups average" -command newAverage

	$m add separator
	$m add command -label "Remove current Op" -command removeOp
	$m add command -label "Edit current Op" -command editOp

	entry $grpsop -textvariable gop  -state disabled -width 16
	entry $grpsub -textvariable gsub -state disabled -width 12

	button $execg -bg green -activebackground yellow \
	    -command execCurFop -text "Exec Current Op" -width 12
    }

    clearGrid $fop

    grid $top -sticky news
    grid $bot

    grid $title -row 0 -column 0 -columnspan 4 -ipadx 5 -sticky w
    grid $edit $grpsop $grpsub $execg

    set lbxolabels {op   groups   minus-group result-group}
    set widths     {10   16       12          12}
    set sorts  {
	-ascii -ascii -integer -integer
    }
    array set lbxos \
	[scrolledlistboxes $bot $lbxolabels $widths $sorts $gselSB {}]

    bind $lbxos(0) <ButtonRelease-1> {+selectFop}
}

proc removeOp {} {
}
proc editOp {} {
}
array set groupfop {}
array set fopgroup {}
proc execCurFop {} {
    global lbxos ifop curOp gnum calc isfop fopgroup group groupfop
    global calcgroup
    global calc iX xaxlbl textf
    global readfinished readlines nread
    global colds colLbls colRanges colsrch colsrchS lbxnames nlbxs

    if { $ifop < 0 || $ifop >= [$lbxos(0) index end] } {
	Dialog_Warn "Invalid file operation.\n"
	return
    }
    set sub -1
    if { [string match [set op [$lbxos(0) get $ifop]] subtraction] } {
	set sub 1
    } elseif { [string match $op average] } {
	set sub 0
    }
    if { $sub < 0 } {
	Dialog_Warn "Invalid file operation.\n"
	return
    }
    # check that groups exist and are calc
    if { [info exists iXa] } { unset iXa }
    set grps [$lbxos(1) get $ifop]
    set ngrp [llength $grps]
    if { $ngrp < 1 } {
	Dialog_Warn "No groups to operate on.\n"
	return
    }
    foreach grp $grps {
	if { ! [info exists group($grp)] } {
	    Dialog_Warn "Group $grp doesn't exist.\n"
	    return
	}
	setGnum $grp
	if { $calc == 0 } {
	    Dialog_Warn "Group $grp has not been solved.\n"
	    return
	}
	set iXa($grp) $iX
    }
    if { $sub } {
	set grps [$lbxos(2) get $ifop]
	set ngrp [llength $grps]
	if { $ngrp < 1 } {
	    Dialog_Warn "No subtraction group.\n"
	    return
	}
	foreach grp $grps {
	    if { ! [info exists group($grp)] } {
		Dialog_Warn "Group $grp doesn't exist.\n"
		return
	    }
	    setGnum $grp
	    if { $calc == 0 } {
		Dialog_Warn "Group $grp has not been solved.\n"
		return
	    }
	    set iXa($grp) $iX
	}
    }
    foreach ig [array names iXa] {
	if { $iXa($ig) != $iX } {
	    Dialog_Warn \
		"Group $ig ind. var. \n$xaxlbl($iXa($ig)) != $xaxlbl($iX).\n"
	    return
	}
    }

    # like FBG
    # passed sanity checks
    # set up for fop with .p files
    # by putting them in a setlist first from grpsop

    pipeWrite "s1 lc"
    vwait readfinished

    set grps [$lbxos(1) get $ifop]
    if { ! $sub && [llength $grps] < 2 } {
	Dialog_Warn "need at least two files for averaging.\n"
	return
    }

    set ninsets 0
    foreach grp $grps {
	setGnum $grp
	# try to get the set number for the p-file
	pipeWrite "d ne $textf(1)"
	vwait readfinished
	set rl [lindex $readlines 0]
	if { [scan $rl {%s %d} dum iset] < 1 || $iset < 1 } {
	    Dialog_Warn "group $grp failed to get setnumber for $textf(1).\n"
	    return
	}
	incr ninsets
	pipeWrite "s l+ $iset"
	vwait readfinished
    }

    # if this is subtraction add the subtrahend dataset to the setlist
    # and define it
    if { $sub } {
	set gsub [lindex [$lbxos(2) get $ifop] 0]
	setGnum $gsub
	pipeWrite "d ne $textf(1)"
	vwait readfinished
	set rl [lindex $readlines 0]
	if { [scan $rl {%s %d} dum iset] < 1 || $iset < 1 } {
	    Dialog_Warn "group $grp failed to get setnumber for $textf(1).\n"
	    return
	}
	pipeWrite "s l+ $iset"
	vwait readfinished
	# define the subtractor
	pipeWrite "s ss $iset"
	vwait readfinished
    }

    array set stndForiX {1 2 2 -1 3 -2 4 -3}
    # set the fop independent var from iX
    pipeWrite "s o c $stndForiX($iX)"
    vwait readfinished
    # for now set chi and binwids to zero
    pipeWrite "s o x 0."
    vwait readfinished
    pipeWrite "s o b 0. 0."
    vwait readfinished
    # set fop ff (use first file values) = 0
    # for sub this will change to 1 automatically
    # but for average want concat not discard of non-overlapping data
    pipeWrite "s o f 0"
    vwait readfinished
    # later we can add interface to allow user
    # to choose chi fit-order binwidth first-file-for-average etc


    # do the fop
    if { $sub } {
	pipeWrite "s s"
    } else {
	pipeWrite "s a"
    }
    vwait readfinished
    # list the setlist to find out which setlist it wrote the result(s) to
    pipeWrite "s1 l"
    vwait readfinished
    # the second line of output should have: op lnr opc lns, lnr is link
    # to result
    if { [scan [lindex $readlines 1] {%s %s %s %s %d %d %d %d} \
	      lop llnr lopc llns op lnr opc lns] < 8 } {
	Dialog_Warn "failed to get file-op result setlist.\n"
	return
    }
    # for subtract the result setlist will have nInputsets = nread - 3
    # for average  the result setlist will have only one set


    # list the result setlist to find the result dataset numbers
    pipeWrite "s$lnr l"
    vwait readfinished
    # readline(3) i match set srcfile
    set nrsets [expr $nread - 3]
    if { $sub && $nrsets != $ninsets } {
	Dialog_Warn "subtraction produced $nrsets != number input sets.\n"
	return
    } elseif { ! $sub && $nrsets != 1 } {
	Dialog_Warn "average produced $nrsets result-sets != 1.\n"
	return
    }

    # put the result setlist into the next available group
    set eg [nextEmptyGroup 0]
    set gnum $eg
    setGroup
    set nmsr 0

    # These will be new files, so collect the stats and put them in the
    # listboxes
    # get the filenames in this setlist as they appear in the listbox
    # s li gives the srcin. Here we want the diskfile names
    # start line 4 ilst match iset fname
    # put the data for these files into lbxs with gnum=$eg
    # first make the list of files from listbox


    set values(fnam) {}
    set values(gnum) {}
    set values(npts) {}
    set values(Fix) {}
    set values(ref) {}
    set values(fbg) {}

    foreach i $colsrch { set values($i) {} }

    set failed ""

    # read result set lists
    set slds {}
    for {set j 3} {$j < [llength $readlines]} {incr j} {
	set rl [lindex $readlines $j]
	if { [llength $rl] < 4 } {
	    Dialog_Warn "FOP result setlist has an invalid line.\n"
	    continue
	}
	set iset [lindex $rl 2]
	lappend slds $iset
	# write this result to disk
	pipeWrite "d$iset w"
	vwait readfinished
    }

    set nslds 0

    foreach iset $slds {
	pipeWrite "d$iset"
	vwait readfinished
	# read back npts this group, set nmsr=0
	set nscan [scan [lindex $readlines 1] \
		       {%[*]%d %d %d %d %d %d %s} \
		       dum is npts ncol ny nass df srcf]
	if { $nscan < 8 } {
	    append failed "Failed to get set $iset npts.\n"
	    puts "npts=$npts ncol=$ncol ny=$ny nass=$nass df=$df srcf=$srcf"
	    continue
	}
	if { $npts < 1 || $ny < 4 } {
	    append failed "No data or nSxy < 4 for set $iset.\n"
	    continue
	}
	lappend values(fnam) [file tail $srcf]
	lappend values(gnum) $gnum
	lappend values(npts) $npts
	lappend values(fs) ALL
	
	set nas [readCurrentColAssigns]

	if { [llength $colRanges(Ei)] < 2 } {
	    lappend values(Fix) Ei
	} elseif { [llength $colRanges(Ef)] < 2 } {
	    lappend values(Fix) Ef
	} else {
	    lappend values(Fix) None
	}

	lappend values(y1) [getRange \
			[concat $colRanges(y1) $colRanges(y2) \
			     $colRanges(y3) $colRanges(y4)]]

	# to determine beam refer type, check for PBSCRIPT: corrected in header

	pipeWrite "d h \"PBSCRIPT: corrected\""
	vwait readfinished
	set rl [lindex $readlines 0]
	if { [string length $rl] > 0 && \
		 [string match "*corrected for monitor*" $rl] } {
	    lappend values(ref) mon
	} else {
	    lappend values(ref) time
	}
	lappend values(fbg) 0

	foreach col $colsrchS {
	    lappend values($col) $colRanges($col)
	}
	incr nslds
    }

    if { [string length $failed] > 0 } {
	Dialog_Warn "$failed"
    }

    # put the new values into lbxs
    set lbxdata {}
    foreach col $lbxnames {
	lappend lbxdata $values($col)
    }

    appendLbxData $lbxdata
    set groupfop($gnum) $lbxdata
    set isfop 1
    set fopgroup($gnum) 1
    set calcgroup($gnum) 0
    set calc 0
    saveGroup

    hiliteFopGroup
    # dont put the sl into the group as the group is not being solved
    # mark this group as FOP result (not PB group) so that it wont ever be calc
}

proc getRange {vlst} {
    if { [set vlen [llength $vlst]] < 1 } { return [list 0 0] }
    set vinc [lsort -real -increasing $vlst]
    set mn [lindex $vinc 0]
    set mx [lindex $vinc [expr $vlen - 1]]
    if { $mx > $mn } {
	return [list $mn $mx]
    } else {
	return $mx
    }
}

proc newSubtraction {} {
    global lbxos ifop
    for {set i 1} {$i < 4} {incr i} {
	$lbxos($i) insert end ""
    }
    $lbxos(0) insert end subtraction
    set ifop [expr [$lbxos(0) index end] - 1]
    $lbxos(0) selection set $ifop
    selectFop
}
proc newAverage {} {
    global lbxos ifop
    for {set i 1} {$i < 4} {incr i} {
	$lbxos($i) insert end ""
    }
    $lbxos(0) insert end addition
    set ifop [expr [$lbxos(0) index end] - 1]
    $lbxos(0) selection set $ifop
    selectFop
}

set ifop -1
array set curOp {}
proc selectFop {} {
    global lbxos ifop
    global curOp
    set ops [$lbxos(0) get 0 end]
    if { [llength $ops] < 1 } { return }
    set sel [$lbxos(0) curselection]
    if { [llength $sel] < 1 } { return }
    set ifop [lindex $sel 0]
    for {set i 0} {$i < 4} {incr i} {
	if [catch {$lbxos($i) get $ifop} curOp($i)] { set curOp($i) "" }
    }
    opEntState 0 0
    opEntState 1 0
    if { [string length $curOp(3)] < 1 } {
	if { [string length $curOp(0)] > 0 } {
	    opEntState 0 1
	    if { [string match $curOp(0) subtraction] } {
		opEntState 1 1
	    }
	}
    }
}

proc opEntState {sub normal} {
    global grpsop grpsub
    # config the groupOp entries normal or disabled along with binding
    set ent $grpsop
    if { $sub } { set ent $grpsub }
    if { $normal } {
	$ent configure -state normal
	bind $ent <Leave> "updateOp $sub"
    } else {
	$ent configure -state disabled
	bind $ent <Leave> {}
    }
}

proc updateOp {sub} {
    global lbxos grpsop grpsub ifop
    if { $sub } {
	$lbxos(2) delete $ifop
	$lbxos(2) insert $ifop [$grpsub get]
    } else {
	$lbxos(1) delete $ifop
	$lbxos(1) insert $ifop [$grpsop get]
    }
}

proc groupResultFops {} {
    global fop
    if { ![winfo exists $fop] } { buildFileOps }
    if { ![winfo ismapped $fop] } { wm deiconify $fop }
    raise $fop
}

#################### saving/setting group data and session ############

proc listGroupOptions {} {
    # these are the group setup variables/arrays
    global cellFile monofilter monitor
    global nsfv sfv
    global acszero
    global cm cen

    global npts nmsr calc isfop ngrpsets
    global textf
    # scrfiles dsFile mset pset are all dsnumber dependent
    global solnfiles
    global fbgFile
    global stndTolArr

    # return a list of commands that will set the current group setup values
    set cmdlst {}
    lappend cmdlst "set cellFile $cellFile"

    lappend cmdlst "set monofilter \{$monofilter\}"
    lappend cmdlst "set monitor \{$monitor\}"
    lappend cmdlst "arraySet nsfv \{[array get nsfv]\}"
    lappend cmdlst "arraySet sfv \{[array get sfv]\}"
    lappend cmdlst "set acszero $acszero"
    lappend cmdlst "arraySet cm \{[array get cm]\}"
    lappend cmdlst "arraySet cen \{[array get cen]\}"

    lappend cmdlst "set npts $npts"
    lappend cmdlst "set nmsr $nmsr"
    lappend cmdlst "set calc $calc"
    lappend cmdlst "set isfop $isfop"
    #lappend cmdlst "set mset $mset"
    #lappend cmdlst "set pset $pset"
    lappend cmdlst "set ngrpsets $ngrpsets"
    lappend cmdlst "arraySet textf \{[array get textf]\}"
    lappend cmdlst "arraySet stndTolArr \{[array get stndTolArr]\}"
    lappend cmdlst "set solnfiles \{$solnfiles\}"
    # dsFile is set in readFile or subtractFBG
    #lappend cmdlst "arraySet dsFile \{[array get dsFile]\}"

    return $cmdlst
}
proc setPlotOptions {cmdlist} {
    # call with groupplot(i) as cmdlist
    global ngrpsets iwin dsns
    global pDatNames mDatNames
    global Xpts Xmsr iX Rf
    global stylst rstylst iWts
    global pencolor
    global isrcplot eview
    global gnufil gnustyl gnucol nsetfs plotsymbolsize plotconfig
    global iwin gnuPSopt bltPSopt

    foreach pdat $pDatNames { global $pdat }
    foreach mdat $mDatNames { global $mdat }

    foreach cmd $cmdlist {eval $cmd}
}
set lastPlot ""
proc listGlobalOptions {} {
    global gnum fbgFile fileFilter heightlines onePlotPerGrp
    global autoStndTol PLT BLT GNU lastPlot
    set cmdlst {}
    set lastPlot 0
    if { $BLT } { set lastPlot 1 }
    if { $GNU } { set lastPlot 2 }
    lappend cmdlst "set fileFilter $fileFilter"
    lappend cmdlst "set heightlines $heightlines"
    lappend cmdlst "setGrpLines 0"
    lappend cmdlst "set onePlotPerGrp $onePlotPerGrp"
    lappend cmdlst "array set autoStndTol \{[array get autoStndTol]\}"
    lappend cmdlst "arraySet fbgFile \{[array get fbgFile]\}"
    #lappend cmdlst "gotoGroup $gnum"
    return $cmdlst
}
proc setGlobalOptions {cmdlist} {
    global gnum fbgFile fileFilter heightlines onePlotPerGrp
    global autoStndTol
    foreach cmd $cmdlist {eval $cmd}
}
proc listGroupPlot {} {
    # separately save plot options since BLT is optional
    global ngrpsets iwin dsns
    global pDatNames mDatNames
    global Xpts Xmsr iX Rf
    global stylst rstylst iWts
    global pencolor
    global isrcplot eview
    global gnufil gnustyl gnucol nsetfs plotsymbolsize plotconfig
    global iwin gnuPSopt bltPSopt

    foreach pdat $pDatNames { global $pdat }
    foreach mdat $mDatNames { global $mdat }

    set cmdlst {}
    lappend cmdlst "set iwin $iwin"
    lappend cmdlst "set ngrpsets $ngrpsets"
    # plot menu options

    lappend cmdlst "set dsns [list $dsns]"
    lappend cmdlst "array set eview \{[array get eview]\}"
    lappend cmdlst "arraySet isrcplot \{[array get isrcplot]\}"
    lappend cmdlst "arraySet pencolor \{[array get pencolor]\}"
    lappend cmdlst "set stylst [list $stylst]"
    lappend cmdlst "set rstylst [list $rstylst]"
    lappend cmdlst "set iWts [list $iWts]"
    lappend cmdlst "set iX $iX"
    lappend cmdlst "set Xpts [list $Xpts]"
    lappend cmdlst "set Xmsr [list $Xmsr]"
    lappend cmdlst "set Rf [list $Rf]"
    lappend cmdlst "arraySet gnufil \{[array get gnufil]\}"
    lappend cmdlst "arraySet gnustyl \{[array get gnustyl]\}"
    lappend cmdlst "arraySet gnucol \{[array get gnucol]\}"
    lappend cmdlst "arraySet nsetfs \{[array get nsetfs]\}"
    lappend cmdlst "arraySet plotsymbolsize \{[array get plotsymbolsize]\}"
    lappend cmdlst "array set plotconfig \{[array get plotconfig]\}"
    lappend cmdlst "arraySet gnuPSopt \{[array get gnuPSopt]\}"
    lappend cmdlst "arraySet bltPSopt \{[array get bltPSopt]\}"
    foreach dat $pDatNames {
	lappend cmdlst "set $dat [list [subst $$dat]]"
    }
    foreach dat $mDatNames {
	lappend cmdlst "set $dat [list [subst $$dat]]"
    }
    return $cmdlst
}


proc checkSessionRecovery {} {
    global BLT GNU PLT fbgFile runPlot
    if { ![istextfile pbsession] } { return }
    if [catch {uplevel #0 source pbsession} msg] {
	Dialog_Warn "Error reading pbsession file. Unable to recover session: $msg\n"	
    } else {
	Dialog_Warn "previous session recovered.\nIf you want a new session rename or delete the pbsession file.\n"
    }
}

proc saveSession {act} {
    # for session saving
    global group groupplot gnum fopgroup groupfop
    global currentDir BLT PLT debugFPT DEBUG
    global setplot

    if { $act == 0 } {
	if { $DEBUG } {
	    close $debugFPT
	    debugClose
	}
	exit
    }

    # first save the current group
    saveGroup
    set f [open pbsession w]
    # save globals
    # NB this includes fbg required for filling listboxes
    puts $f "setGlobalOptions [list [listGlobalOptions]]"
    puts $f "array set setplot \{[array get setplot]\}"
    puts $f "array set fopgroup \{[array get fopgroup]\}"
    puts $f "array set calcgroup \{[array get calcgroup]\}"
    puts $f "array set groupfop \{[array get groupfop]\}"
    
    set grps [lsort -integer [array names group]]
    foreach grp $grps {
	# these filenames will be short.
	# readFiles no need to prepend currentDir
	# dont save fop groups
	# they get added onto the lbxs at the end
	if { [info exists fopgroup($grp)] && $fopgroup($grp) } {
	    puts $f "setGnum $grp"
	    puts $f "set group($grp) [list $group($grp)]"
	    continue
	}
	set fls [getFilesInGroup $grp]
	if { [llength $fls] < 1 } { continue }
	puts $f "setGnum $grp"
	puts $f "set group($grp) [list $group($grp)]"
	if { $PLT && [info exists groupplot($grp)] } {
	    puts $f "if { \$PLT } {set groupplot($grp) [list $groupplot($grp)]}"
	}
	puts $f "setGnum $grp"
	puts $f "readFiles [list $fls]"
	# on recovery re-solve any groups that had been solved
	# this makes it possible to use either plot package
	# I think it also makes it unnecessary to save all the data lists???
	puts $f "if { \$calc } { solveGroup 0 }"
    }

    foreach grp $grps {
	if { ! [info exists fopgroup($grp)] || ! $fopgroup($grp) } { continue }
	if { [info exists groupfop($grp)] } {
	    puts $f "appendLbxData [list $groupfop($grp)]"
	}
    }
    puts $f "hiliteGroups"
    close $f

    if { $act >= 2 } {
	if { $DEBUG } {
	    close $debugFPT
	    debugClose
	}
	exit
    }
}

array set groupsolved {}
proc recordCalc {args} {
    # trace for writing calc
    # store calc value for current group
    global groupsolved calc gnum
    set groupsolved($gnum) $calc
}

proc setGnum {n} {
    global gnum
    set gnum $n
    setGroup
}
proc initGroup {} {
    global textf npts nmsr calc isfop mset pset
    arraySet textf {1 "" 2 ""}
    set npts 0
    set nmsr 0
    set calc 0
    set isfop 0
    set mset 0
    set pset 0
}
proc setGroupOptions {cmdlist} {
    # during session only thing to track is group options
    # these are the group setup variables/arrays
    global cellFile fileFilter monofilter monitor
    global nsfv sfv
    global acszero
    global cm cen

    global npts nmsr calc ngrpsets
    global textf
    # scrfiles dsFile mset pset are all dsnumber dependent
    # it is possible that on session recovery the file will
    # endup in different setnumbers than last session
    # so dont rely on these.
    # instead whenever gnum changes regenerate them
    # from the filenames
    # NB srcfiles and dsFile are just inverse arrays of each other
    # and mset and pset are dsnumbers for textf
    global srcfiles dsFile mset pset
    global fbgFile

    # mset and pset are dataset numbers!!!
    # can regen them from textf file names

    foreach cmd $cmdlist {
	uplevel #0 $cmd
    }
}
proc showSrcFile {i} {
     global snum mp
     set snum $i
     # activate the trace on mp to call showText
     set mp 3
}

proc updateSrcFileMenus {} {
    # srcfile menus are in text viewer(sfvm) and in plot-selector(sfm)
    global sfvm sfm PLT
    global srcfiles ngrpsets isrcplot solnfiles
    # in text viewer, command select one of the srcfiles for showSrcFile
    # in plot-selector, check turns srcset view on/off in plot via viewMsrSet
    if { ! [winfo exists $sfvm] } {
	buildShowText -width 60 -height 30 -wrap none -font {courier 10}
    }
    $sfvm delete 0 end
    $sfm delete 0 end
    set maxpercol 32
    set j 1
    for {set i 1;set j 0} {$i <= $ngrpsets} {incr i;incr j} {
        set fil [file tail [lindex $solnfiles $j]]
	set cbrk 0
	if { $j > $maxpercol } {
	    set cbrk 1
	    set j 1
	}
	if { [info exists svfm] && [winfo exists $svfm] } {
	    $sfvm add command -label "$i $fil" -columnbreak $cbrk \
		-command "showSrcFile $i"
	}
	if { $PLT && [info exists sfm] && [winfo exists $sfm] } {
	  $sfm  add check -label "$i $fil" -variable isrcplot($i) \
	      -columnbreak $cbrk -command "viewMsrSet 1"
	}
    }
}

proc setGroup {} {
    # this just saves/loads local variables from the group arrays
    global gnum group groupplot BLT PLT
    if { $gnum < 1 } { set gnum 1 }
    if { ! [info exists group($gnum)] } {
	# new groups only inherit options
	initGroup
	saveGroup
    } else {
	setGroupOptions $group($gnum)
	if { $PLT } { setPlotOptions $groupplot($gnum) }
	updateSrcFileMenus
    }
    hiliteActiveGroup
    setSn 1
}
proc incrGroup {} {
    global gnum tv PLT
    # save previous group
    saveGroup
    incr gnum
    setGroup
    if { [winfo exists $tv] && [winfo ismapped $tv] } { showText }
    if { $PLT && [checkPlotWindow] == 0 } { showPlot }
}
proc gotoGroup {gn} {
    global gnum tv PLT
    if { $gn < 1 } { return }
    saveGroup
    set gnum $gn
    setGroup
    if { [winfo exists $tv] && [winfo ismapped $tv] } { showText }
    if { $PLT && [checkPlotWindow] == 0 } { showPlot }
}
proc decrGroup {} {
    global gnum tv PLT
    if { $gnum == 1 } return
    saveGroup
    incr gnum -1
    setGroup
    if { [winfo exists $tv] && [winfo ismapped $tv] } { showText }
    if { $PLT && [checkPlotWindow] == 0 } { showPlot }
}
proc saveGroup {} {
    global gnum group groupplot PLT isfop
    set group($gnum) [listGroupOptions]
    if { $PLT && ! $isfop } { set groupplot($gnum) [listGroupPlot] }
}
proc savePlot {} {
    global gnum groupplot
    set groupplot($gnum) [listGroupPlot]
}


array set groupiwin {}
array set iwingroup {}
proc checkPlotWindow {} {
    # returns 1 if gnum is already plotted in the correct window
    global iwin gnum onePlotPerGrp
    global BLT GNU gpipes
    global groupiwin iwingroup
    # have setGroup which sets iwin

    # first get the current plotwindow inventory
    set iwinlst {}
    if { $BLT } {
	set tls [winfo children .]
	foreach w $tls {
	    if { [regexp {\.g([0-9]+)} $w match i] } {
		lappend iwinlst $i
	    }
	}
    } else {
	# count open gpipes
	set gps [array names gpipes]
	foreach gp $gps {
	    if { [string length $gpipes($gp)] > 0 } { lappend iwinlst $gp }
	}
    }
    set iwinlst [lsort -integer $iwinlst]
    foreach i $iwinlst {
	if { ![info exists groupiwin($i)] } { set groupiwin($i) 0 }
    }
    foreach i [array names groupiwin] {
	if { [lsearch -exact $iwinlst $i] < 0 } { unset groupiwin($i) }
    }
    foreach i [array names groupiwin] {
	set igrp $groupiwin($i)
	set iwingroup($igrp) $i
    }
    # first, if the requested plotwindow=iwin already has another
    # gnum plotted in it,

    # if the onePlotPerGrp flag is set
    # check to see if current group is plotted to some other existing iwin
    # and if so reset to that, or if current group is in not in any
    # plotwindow set iwin to next available
    if { $onePlotPerGrp } {
	if { [info exists iwingroup($gnum)] && $iwingroup($gnum) > 0 } {
	    set iwin $iwingroup($gnum)
	    return 1
	} else {
	    set iwin 1
	    while { [info exists groupiwin($iwin)] && $groupiwin($iwin) > 0 } {
		incr iwin
	    }
	    updatePlotWindowsMenu $iwin
	    return 0
	}
    } else {
	# must plot current group to iwin
	if { [info exists groupiwin($iwin)] && $groupiwin($iwin) == $gnum } {
	    return 1
	}
	return 0
    }
}


array set lbxs {}
#set lbxws {}

proc scrolledlistboxes {parent labels widths sorttypes sb hsb} {
    # returns [array get lbxs]

    # create multi listboxes to represent columns of parent table
    # and set up coherent vertical scanning/scrollbar bindings
    # optional horiz scrollbar on listboxs with indices in the hsb list

    set n [llength $labels]
    set m [llength $widths]
    if { $m < $n } { set n $m }
    set m [llength $sorttypes]
    if { $m < $n } { set n $m }
    if { $n < 1 } { return 0 }
    set nTableColumns $n
    clearGrid $parent
    set lbxws {}

    set nhsb [llength $hsb]

    for {set i 0; set j 0} {$i < $n} {incr i; incr j 2} {
	set lb [lindex $labels $i]
	set wd [lindex $widths $i]
	if { [set lblen [string length $lb]] > $wd } { set $wd $lblen }
	set lbx $parent.lb$i
	set lbxs($i) $lbx
	set head $parent.name$i
	set bar $parent.bar$i
	#set entry $parent.e$i
	if { ! [winfo exists $lbx] } {
	    listbox $lbx -bd 0 -selectmode extended -width $wd -height 10
	    label $head -background black -foreground white
	    frame $bar -width 2 -bg grey
	    # -cursor sb_h_double_arrow
	}
	#if { ! [winfo exists $entry] && $row0 } { entry $entry } 
	$lbx delete 0 end
	if { $sb } { $lbx configure -yscrollcommand "$parent.sb set" }
	if { $nhsb && [lsearch $hsb $i] >= 0 } {
	    $lbx configure -xscrollcommand "$parent.xsb$i set"
	    set xsb $parent.xsb$i
	    if { ! [winfo exists $xsb] } { scrollbar $xsb -orient horizontal }
	    $xsb configure -command [list $lbx xview]
	    grid $xsb -row 2 -column $j -sticky ew
	}
	$head configure -text $lb

	lappend lbxws $lbx

	grid $head -row 0 -column $j -sticky ew
	grid $lbx -row 1 -column $j -sticky news
	grid $bar -row 0 -column [expr $j+1] -rowspan 2 -sticky ns
	grid columnconfigure $parent [expr $j+1] -weight 0
	grid columnconfigure $parent $j -minsize $wd
    }

    #set cmnd ""
    #for {set i 0} {$i < $n} {incr i} {
    #set cmnd "$cmnd $parent.lb$i yview\;"
    #}

    if { $sb } {
	set sb $parent.sb
	if { ! [winfo exists $sb] } { scrollbar $sb -orient vertical }
	$sb configure -command [list BindYview $lbxws]
	grid $sb -row 1 -column [expr 2*$n] -sticky ns
	grid columnconfigure $parent [expr 2*$n] -weight 0
    }

    grid rowconfigure $parent 0 -weight 0

    #foreach w $lbxws
    #	bind $w <B2-Motion>
    #	bind $w <B2-Motion> "BindDragto %x %y $lbxws"
    #	bind $w <Button-2>
    #	bind $w <Button-2> "BindMark %x %y $lbxws"


    for {set i 0} {$i < $n} {incr i} {
	set lblWidget $parent.name$i
	bind $lblWidget <Button-1> \
	    [list sortlistboxes $lbxws $i [lindex $sorttypes $i]]
	#setLbxWidth $lbxs($i) [lindex $labels $i]
    }
    return [array get lbxs]
}



proc clearGrid {container} {
    # clear grid widgets but remember configs
    foreach w [grid slaves $container] {
	if {[winfo exists $w]} {grid remove $w}
    }
}
proc globalarr {avar} {
    # in response to 2009 FC11 global cant work for arr(elem)
    # if avar is array element name just return the array name
    if { [regexp {([^\(]*)\(.+\)$} $avar wm anam] } {
	return $anam
    }
    return $avar	
}
proc BindYview {wlist args} {
    foreach w $wlist {
	eval {$w yview} $args
    }
}
proc BindDragto { x y args } {
    foreach w $args {
	$w scan dragto $x $y
    }
}
proc BindMark { x y args } {
    foreach w $args {
	$w scan mark $x $y
    }
}

set sortdirection {-increasing -decreasing}
set sortdirectionindex 0

proc sortlistboxes {listwidgets column type args} {
    # user responsibitly to make sure sort column can be sorted
    # for indexed sort to work the data gets collected as a list of rows
    # with each row containing a data list one datum per column.
    # sorting on a full column or on selected elements in a column
    # does congruent sort on all other columns preserving selection states

    #global wctrl
    global sortdirection sortdirectionindex
    #if { ! $wctrl(tablesortenable) } { return }
    set n [llength $listwidgets]
    if {$n <= $column} { return }
    set widget [lindex $listwidgets $column]
    set onlyselected 0
    if { [llength [$widget curselection]] > 1 } { set onlyselected 1 }
    set N [$widget index end]
    if { [string length $args] > 0 && [llength $args] > 0 } {
	set Nent [lindex $args 0]
	if { $Nent > 0 && $Nent < $N } { set N $Nent }
    }
    set workl {}
    set sortlist {}
    set ilst {}

    for {set i 0} {$i < $N} {incr i} {
	set sl {}
	set thisline {}
	for {set j 0} {$j < $n} {incr j} {
	    set w [lindex $listwidgets $j]
	    # need the subst to get rid of braces around quoted strings
	    #lappend thisline "subst [$w get $i]"
	    lappend thisline [subst [$w get $i]]
	    lappend sl [$w selection includes $i]
	}
	lappend thisline $i $sl
	lappend workl $thisline
	if { $onlyselected } {
	    if { [lindex $sl $column] } {
		lappend sortlist $thisline
		lappend ilst $i
	    }
	} else {
	    lappend sortlist $thisline
	}
    }

    if { [llength $args] > 1 } {
	set incrsort [lindex $args 1]
	if { $incrsort } {
	    set sorttyp -increasing
	} else {
	    set sorttyp -decreasing
	}
    } else {
	set sorttyp [lindex $sortdirection $sortdirectionindex]
	set sortdirectionindex [expr 1 - $sortdirectionindex]
    }

    set sortlist [eval lsort \
	    $sorttyp -index $column $type {$sortlist}]

    set nsort [llength $sortlist]
    if { $onlyselected } {
	for {set i 0} {$i < $nsort} {incr i} {
	    set thisline [lindex $sortlist $i]
	    set indx [lindex $ilst $i]
	    set workl [lreplace $workl $indx $indx $thisline]
	}
	writelistboxrows $listwidgets $workl
    } else {
	writelistboxrows $listwidgets $sortlist
    }
    return $nsort
}
proc writelistboxrows {listwidgets sortlist} {
    # sortlist is a list of table row data
    # the rows can have two extra columns appended by sort
    set n [llength $listwidgets]
    if {$n < 1} { return }
    for {set i 0} {$i < $n} {incr i} {
	set widget [lindex $listwidgets $i]
	$widget selection clear 0 end
	$widget delete 0 end
    }
    set N [llength $sortlist]
    for {set i 0} {$i < $N} {incr i} {
	set thisline [lindex $sortlist $i]

	set nc [llength $thisline]
	set slst {}
	set sok 0
	if { $nc > [expr $n + 1] } {
	    set slst [lindex $thisline [expr $nc + 1]]
	    set sok 1
	}
	if { $nc > $n } { set nc $n }
	if { $sok && [llength $slst] < $nc } { set sok 0 }
	for {set j 0} {$j < $nc} {incr j} {
	    set widget [lindex $listwidgets $j]
	    $widget insert end [lindex $thisline $j]
	    if { $sok && [lindex $slst $j] } {
		$widget selection own
		$widget selection set $i
	    }
	}
    }
    hiliteActiveGroup
}

######## post processing ave OR sub files #############
# another toplevel average-subtract-plot files
# repro the sellectlist  add a pbsoln  file filter
# then a scrolledlistboxes
# op vars  files  resultFile
# newOp tk_option average subtract
# ind. vars. tk_option Qx Qy Qz E
# set files from selected button
# set subtrahend button for subtractOp current
# exec current Op
#
# can we include a plot anyfile interface
# check- plot selected file  

####### miscellaneous Dialog procs currently not used
proc Show_Progress {i N msg} {
    # meant to be recalled everytime the counter i is changed
    global progress pausing
    set p .progress
    #puts "Show_Progress i=$i"

    #startTiming SPcreate
    set id 0
    if { $i == 0 } { set id [Dialog_Create $p "Progress"] }
    #elapsedTime SPcreate
    if { $id } {
	set progress(cancel) 0
	message $p.msg -text $msg -aspect 200
	set b [frame $p.buttons]
	button $b.cancel -text CANCEL -command {set progress(cancel) 1}
	pack $b.cancel
	scale $p.scale -from 1 -to $N -variable progress(value) \
		-orient horizontal
	pack $p.msg $b $p.scale -side top -fill x
    } else {
	$p.msg configure -text $msg
	$p.scale configure -to $N
    }

    set progress(value) $i
    #puts "ShowProgress grab current= [grab current] about to grab $p"
    #grab $p

    # continually raise the dialog window (some wm's may allow hiding it)

    set pausing 0
    #puts "ShowProgress before setPausing after info= [after info]"
    set id [after 20 "setPausing progress(cancel) 0"]
    tkwait variable progress(cancel)

    if { ! $pausing } {
	after cancel $id
	#catch {grab release $p}
	return 0
    }

    if { $progress(cancel) || $i == $N } {
	#catch {grab release $p}
	Dialog_Dismiss $p
	return 0
    }
    return 1
}
proc Show_Message {i msg} {
    # call with i=0 to start and i=1 to continue i=2 to end
    set m .message
    if { $i < 1 } {
	if [Dialog_Create $m "Message"] {
	    message $m.msg -text $msg -fg black -aspect 200
	    pack $m.msg -side top -fill x
	} else {
	    $m.msg configure -text $msg
	}
	$m.msg configure -fg blue
	update idletasks
    } elseif { $i == 1 } {
	raise $m
	$m.msg configure -fg red
	update idletasks
    } else {
	Dialog_Dismiss $m
	destroy $m
    }
}

proc Dialog_Entry {msg varName args} {
    #global $varName
    #if { ![info exists $varName] } { uplevel #0 "set $varName {}" }
    upvar #0 $varName var
    global entrystring
    set f .entrystring
    set entrystring(value) $var
    if [Dialog_Create $f "Entry"] {
	message $f.msg -text $msg -aspect 1000
	set b [frame $f.buttons]
	button $b.cancel -text CANCEL -command {set entrystring(ok) 0}
	button $b.ok -text OK -command {set entrystring(ok) 1}
	pack $b.cancel $b.ok -side left -expand 1
	entry $f.ent -relief sunken -bd 2 -textvariable entrystring(value)
	if { [llength $args] > 0 } { $f.ent configure -width $args }
	pack $f.msg $b $f.ent -side top -fill x
	bind $f.ent <Return> {set entrystring(ok) 1}
    } else {
	$f.msg configure -text $msg
	if { [llength $args] > 0 } { $f.ent configure -width $args }
    }
    set entrystring(ok) 0
    Dialog_Wait $f entrystring ""
    if { $entrystring(ok) } { set var $entrystring(value) }
    Dialog_Dismiss $f
    return $entrystring(ok)
}
array set query {}
proc Dialog_Query {msg} {
    global query
    set f .query
    if [Dialog_Create $f "?"] {
	message $f.msg -text $msg -aspect 1000
	set b [frame $f.buttons]
	button $b.yes -text YES -command {set query(ok) 1}
	button $b.no -text NO -command {set query(ok) 0}
	pack $b.yes $b.no -side left -expand 1
	pack $f.msg $b -side top -fill x
    } else {
	$f.msg configure -text $msg
    }
    set query(ok) 0
    Dialog_Wait $f query ""
    Dialog_Dismiss $f
    return $query(ok)
}

proc blink {w option value1 value2 interval} {
    $w config $option $value1
    after $interval [list blink $w $option $value2 $value1 $interval]
}

#####################################################
#  application startup

trace variable iwin w selectPlotWindow
trace variable calc w recordCalc

getDirectoryCellFiles
makeTolArray
# build group-select first so that bld buttons exist
buildShowText -width 60 -height 30 -wrap none -font {courier 10}
buildGroupSelector
buildFileSelector
buildGroupPrep
buildAutoGroup

wm withdraw $tv

# start with heightlines=1 so remove the mid and bot frames
grid remove $bframes(2)
grid remove $bframes(3)
set gnum 1
setGroup
# update so we are sure main window is completely mapped
update idletasks
if { $PLT } { buildPlotSelector }
update idletasks
set fileFilter *

#puts "before session recovery: plotconfig=\n[array get plotconfig]"
checkSessionRecovery
#puts "after session recovery: plotconfig=\n[array get plotconfig]"

doFilter $sellectlist
#wm deiconify $gv
#setPlotGeo
