#!/bin/sh
######################################
#  script to calculate quality information of fastq file
######################################
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 8 ];
then
	MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`

        fastqcdir=$1
        outputdir=$2
        fastqcparms=$3
        fastqfile=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"        

        parameters=$( echo $fastqcparms | tr "_" " " )
        cd $outputdir
        $fastqcdir/fastqc -o $outputdir $parameters $fastqfile
        totlines=`ls -1 *.zip | wc -1 | cut -d ' ' -f 1`
        if [ $totlines -lt 1 ]
        then
              MSG="Warning:fastqc file for $fastqfile was not created. Execution continued."
              echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
              exit 1;
        fi
fi