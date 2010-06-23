#!/usr/bin/wish

package require BLT
set comm(control) /dev/ttyS0
set comm(chan)    null
set comm(mode)    "9600,n,8,1"
#set comm(host)    127.0.0.1
set comm(host)    129.6.120.139
set comm(port)    8001
set comm(iface)   tcp

proc BuildGui { } {

    set xdata [list  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 \
                    16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 \
		    32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 ]
    set ydata [list  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \
		     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \
		     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]

    blt::barchart .g0
    set m [frame .mbar -relief raised -borderwidth 2]
    set mf [menubutton $m.file  -width 8 -text File     -underline 0 \
		-menu $m.file.menu]
    set mfm [menu $mf.menu]
    $mfm add command -label "Start" -underline 0 -command ScalerStart
    $mfm add command -label "Stop " -underline 3 -command ScalerStop
    $mfm add command -label "Clear" -underline 0 -command ScalerClear

    
    pack $m -side top -fill x -anchor nw
    pack $mf -side left -in $m
    pack .g0 -side top -fill both -expand true
    .g0 element create e0 -barwidth 0.4
    .g0 element configure e0 -xdata $xdata -ydata $ydata
    .g0 legend configure -hide 1
    .g0 yaxis configure -logscale 0 -title "Counts"
    .g0 xaxis configure -title "Channel" -stepsize 4
}

proc PlotData { y0data } {
    .g0 element configure e0 -ydata $y0data
}

proc OpenCom { } {
    global comm

    switch $comm(iface) {
	tcp {
	    if [catch {socket $comm(host) $comm(port)} comm(chan)] {
		puts stderr "Can't open $comm(host):$comm(port) : $comm(chan)"
		set comm(chan) null
		after 60000 OpenCom ;# Try to reopen after 1 minute
		return
	    } else {
		fconfigure $comm(chan) -blocking 0
		fileevent $comm(chan) readable ComHandler
	    }
	}
	default {
	    if [catch {open $comm(control) w+} comm(chan)] {
		puts stderr "Can't open serial port $comm(control): $comm(chan)"
		exit
	    }
	    fconfigure $comm(chan) -mode $comm(mode) -blocking 0 -encoding binary
	    fileevent $comm(chan) readable ComHandler
	}
    }
}

proc ComHandler { } {
    global comm
    if {[eof $comm(chan)]} {
	close $comm(chan)
	set comm(chan) null
	after 5000 OpenCom
    } elseif {[gets $comm(chan) comm(input)] > 0} {
	set parts [split $comm(input) :]
	set data [split [string trim [lindex $parts 1]]]
	if {[llength $data] > 16} { PlotData $data }
	#	after 500 AskRate
    }
	
}

proc AskCounts { } {
    global comm
    if {![string match $comm(chan) null]} {
	puts -nonewline $comm(chan) "scaler read\r"
	flush $comm(chan)
    }
    after 1000 AskCounts
}

proc ScalerStart { } {
    global comm
    if {![string match $comm(chan) null]} {
	puts -nonewline $comm(chan) "scaler arm\r"
	flush $comm(chan)
    }
}

proc ScalerStop { } {
    global comm
    if {![string match $comm(chan) null]} {
	puts -nonewline $comm(chan) "scaler abort\r"
	flush $comm(chan)
    }
}

proc ScalerClear { } {
    global comm
    if {![string match $comm(chan) null]} {
	puts -nonewline $comm(chan) "scaler reset\r"
	flush $comm(chan)
    }
}

proc bgerror { args } { 
    puts stderr $args
    return -code break
}


OpenCom
BuildGui
AskCounts

