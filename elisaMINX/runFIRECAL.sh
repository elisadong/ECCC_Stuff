#!/bin/bash

OMP_NUM="5"    # Number of CPUs
MAXTIME="16200"  # Job max time in seconds 4.5hr=16200sec
MEMSIZE="4G"   # Job max memory
BACKEND_machine="eccc-ppp2"
Jobname=getCal.sh
ListDir=/fs/site2/dev/eccc/aq/r1/eld001/logs
EmailAdd=elisa.dong@canada.ca

ord_soumet ${Jobname} -mach ${BACKEND_machine} -t ${MAXTIME} -cpus 1x1x${OMP_NUM} -cm ${MEMSIZE} -shell /bin/bash -jn ${Jobname} -listing ${ListDir} -mail ${EmailAdd} -notify error -notify complete
