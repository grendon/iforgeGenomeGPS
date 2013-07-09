#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 12 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stooped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline - Support #200' "$redmine""
        exit 1;

else
	set -x
	echo `date`

        scriptfile=$0
        alignerdir=$1
        ref=$2
	outputdir=$3
        R1=$5
        A1=$4
        samfile=$6
        bamfile=$7
        samdir=$8
        elog=$9
        olog=${10}
        email=${11}
        qsubfile=${12}
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        cd $outputdir
        $alignerdir/bwa samse $ref $A1 $R1  > $samfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa samse command failed.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $samfile ]
        then
            MSG="$samfile aligned file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline - Support #200' "$redmine,$email""
            exit 1;
        fi
        echo `date`
        ## include sam2bam conversion
	$samdir/samtools view -bS -o $bamfile $samfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="samtools view command failed.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
	if [ ! -s $bamfile ]
	then
	    MSG="$bamfile bam file not created. sam2bam step failed during alignment."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline - Support #200' "$redmine,$email""
	    exit 1;
	fi       
        echo `date`
fi

