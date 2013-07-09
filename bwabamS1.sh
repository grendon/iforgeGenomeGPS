#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 11 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$grendon@illinois.edu,lmainzer@illinois.edu""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        aligndir=$1
        parms=$2
        readparm=$3
        ref=$4
        outputdir=$5
	outputfile=$6
        R=$7
        elog=$8
        olog=$9
        email=${10}
        qsubfile=${11}

        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        parameters=$( echo $parms | tr "_" " " )

        cd $outputdir

        $aligndir/bwa aln $parameters $ref $readparm $R > $outputfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa aln command failed exitcode=$exitcode . alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $outputfile ]
        then
	    MSG="$outputfile aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
	    exit 1;
       fi
        echo `date`
fi
