#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

if {$argc != 1} {
    puts "Usage: remotepeek.tcl instrument"
    puts "where instrument is NG1, NG7, CG1, ..."
    exit
}
set ACTIVEURL "http://reflectometry.org/ipeek/ipeekon.php"
set INSTRUMENT [lindex $argv 0]

package require http
package require Tk

# Fetch a page returning the content
proc fetchpage {url} {
  set page [::http::geturl $url]
  set body [set ${page}(body)]
  unset $page
  return $body
}

# Fetch the visibility status page for the instrument
proc getstatus {instrument} {
    global ACTIVEURL
    set body [fetchpage "${ACTIVEURL}?id=$instrument"]
    return [string match "*[string tolower $instrument]=1*" $body]
}

# Set the visibility status for the instrument to a particular status
proc setstatus {instrument status} {
    global ACTIVEURL
    if {$status} { set state "on" } else { set state "off" }
    set body [fetchpage "${ACTIVEURL}?id=$instrument&status=$state"]

    # update with new status; this automagically sets radio buttons
    global visible
    if {[string match "*[string tolower $instrument]=1*" $body]} {
        set visible "on"
    } else {
        set visible "off"
    }
}

# update the status every ten minutes just in case somebody changes it
proc update_status {instrument} {
    global visible
    if {[getstatus $instrument]} { set visible "on" } else { set visible "off" }     
    after 600000 [list update_status $instrument]
}

# Create window and initialize visibility status
proc init {instrument} {
    global visible
    update_status $instrument
    label .remote -text "Remote iPeek for $instrument (http://reflectometry.org/ipeek)"
    radiobutton .remoteon -text "On" -variable visible -value "on" \
	-command [list setstatus $instrument true]
    radiobutton .remoteoff -text "Off" -variable visible -value "off" \
	-command [list setstatus $instrument false]
    grid .remote -
    grid .remoteon .remoteoff
}

init $INSTRUMENT
