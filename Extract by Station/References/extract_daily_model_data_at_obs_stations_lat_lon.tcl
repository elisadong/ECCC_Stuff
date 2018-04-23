#!/bin/sh
# the next line restarts using tclsh \
exec $SPI_PATH/tclsh "$0" "$@"
#
# Oct /2014
# P. Cheung for extracting model values at station lat-lon positions.
# Six species (TCO, O3, PM2.5, CO, NO and NO2) are of interest.
# Jc Modified for PM2.5 only
#============================================================================
#
package require TclData

puts \n[file tail [info script]]

#
# How to run : script_name input-filename
#
 if { $argc != 4 } {
   puts "Usage: $argv0 input-filename date hour output-directory"
   exit
 }
#
# retrieve all the in-line arguments ..real arguments starts from index 0
#
# m2013110812_001
#

 set infile [lindex $argv 0]
 set curdate [lindex $argv 1]
 set curhour [lindex $argv 2]
 set outdir [lindex $argv 3]

 set hh [format "%02d" ${curhour}]

 set yyyy [string range $curdate 0 3]
 set mm [string range $curdate 4 5]
 set dd [string range $curdate 6 7]
 set yyyymm [string range $curdate 0 5]
 set yyyymmdd [string range $curdate 0 7]

### set fileid [open "VAQUM.Stn.info.CAN.csv" r]
#
# set fileid [open "NA.PanAm.2_5km_station.csv" r]
## set fileid [open "station_id_name_NA_lat_lon_panam2p5km.txt" r]
# gets $fileid line
#set fileid [open "station_id_name_CAN_lat_lon.txt" r]
 set fileid [open "station_id_name_CFFEPS_lat_lon.txt" r]

### set ofileid [open "station_id_NA_panam2p5km.txt" w]
#set ofileid2 [open "${outdir}/O3_model_${yyyymmdd}.csv" a]
#set ofileid3 [open "${outdir}/NO2_model_${yyyymmdd}.csv" a]
 set ofileid4 [open "${outdir}/PM25_model_${yyyymmdd}.csv" a]
#set ofileid5 [open "${outdir}/SO2_model_${yyyymmdd}.csv" a]
#set ofileid6 [open "${outdir}/CO_model_${yyyymmdd}.csv" a]
#set ofileid7 [open "${outdir}/NO_model_${yyyymmdd}.csv" a]

### puts ofileid2 "stationid,O3_model,datetime"
### puts ofileid3 "stationid,NO2_model,datetime"
### puts ofileid4 "stationid,PM25_model,datetime"
### puts ofileid5 "stationid,SO2_model,datetime"
### puts ofileid6 "stationid,CO_model,datetime"
### puts ofileid7 "stationid,NO_model,datetime"

# set ip1 93423264
 set ip1 76696048

#
# 
# open the input RPN stnadard file on unit 101
#=============================================
# if catch returns 0 ==> no error i.e. command ran fine ; status of zero in Tcl ==> false
#                  1 ==> error i.e problem in running the command ; status of one in Tcl ==> True
#
#
 if { [catch { fstdfile open FILE1 read  ${infile} } ] != 0} {
  puts "Error: Can't open RPN file ${infile} for read."
  exit
 }

### fstdfield read FIELD1 FILE1 -1 "" 26314400 -1 -1 "" "GZ"

#fstdfield read FIELDO3 FILE1 -1 "" $ip1 -1 -1 "" "O3"
 fstdfield read FIELDAF FILE1 -1 "" $ip1 -1 -1 "" "AF"
#fstdfield read FIELDN2 FILE1 -1 "" $ip1 -1 -1 "" "N2"
#fstdfield read FIELDS2 FILE1 -1 "" $ip1 -1 -1 "" "S2"
#fstdfield read FIELDNO FILE1 -1 "" $ip1 -1 -1 "" "NO"
#fstdfield read FIELDTCO FILE1 -1 "" $ip1 -1 -1 "" "TCO"
# fstdfield read FIELDM3 FILE1 -1 "" $ip1 -1 -1 "" "M3"

 while { [gets $fileid line] >= 0 } {
   set line [string trim $line]
   set linelist [split $line ","]
   set stnid [lindex $linelist 0]
   set lat [lindex $linelist 1]
   set lon [lindex $linelist 2]

#  puts "extracting: $stnid $lat $lon"

###   set stname [lindex $linelist 8]
###   set citytype [lindex $linelist 10]

   set xy [fstdfield stats FIELDAF -coordpoint $lat $lon]
   set x [lindex $xy 0]
   set y [lindex $xy 1]

#  puts "stn lat/lon $lat, $lon==> x/y: $x, $y"
   if { $x < 0.0 || $y < 0.0 } {
#    puts "stn ${stnid} NA -- lat/lon $lat, $lon==> x/y: $x, $y"
     continue
   }

#  set o3 [fstdfield stats FIELDO3 -coordvalue $lat $lon]
   set pm25 [fstdfield stats FIELDAF -coordvalue $lat $lon]
#  set no2 [fstdfield stats FIELDN2 -coordvalue $lat $lon]
#  set so2 [fstdfield stats FIELDS2 -coordvalue $lat $lon]
#  set no [fstdfield stats FIELDNO -coordvalue $lat $lon]
#  set tco [fstdfield stats FIELDTCO -coordvalue $lat $lon]
#   set m3 [fstdfield stats FIELDM3 -coordvalue $lat $lon]

#  set co [expr { ${co_conversion_factor} * ${tco} }]

#  puts $ofileid2 "${stnid},${o3},${yyyy}-${mm}-${dd} ${hh}:00:00"
#  puts $ofileid3 "${stnid},${no2},${yyyy}-${mm}-${dd} ${hh}:00:00"
   puts $ofileid4 "${stnid},${pm25},${yyyy}-${mm}-${dd} ${hh}:00:00"
#  puts $ofileid5 "${stnid},${so2},${yyyy}-${mm}-${dd} ${hh}:00:00"
#  puts $ofileid6 "${stnid},${co},${yyyy}-${mm}-${dd} ${hh}:00:00"
#  puts $ofileid7 "${stnid},${no},${yyyy}-${mm}-${dd} ${hh}:00:00"

 }

 fstdfile close FILE1
 close $fileid
# close $ofileid
#close $ofileid2
#close $ofileid3
 close $ofileid4
#close $ofileid5
#close $ofileid6
#close $ofileid7



