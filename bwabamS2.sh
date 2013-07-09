#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 14 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$grendon@illinois.edu,lmainzer@illinois.edu""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        alignerdir=$1
        ref=$2
	outputdir=$3
        A1=$4
        A2=$5
        R1=$6
        R2=$7
        samfile=$8
        bamfile=$9
        samdir=${10}
        elog=${11}
        olog=${12}
        email=${13}
        qsubfile=${14}
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        cd $outputdir
        $alignerdir/bwa sampe $ref $A1 $A2 $R1 $R2 > $samfile
        exitcode=$?
        if [ $exitcode -ne o ]
        then
            MSG="bwa sampe command failed exitcode=$exitcode . alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit $exitcode;
        fi

        if [ ! -s $samfile ]
        then
            MSG="$samfile aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit 1;
        fi
        echo `date`
        ## sam2bam conversion
	$samdir/samtools view -bS -o $bamfile $samfile
        exitcode=$?
        if [ $exitcode -ne o ]
        then
            MSG="samtools view command failed exitcode=$exitcode . alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit $exitcode;
        fi
	if [ ! -s $bamfile ]
	then
	    MSG="$bamfile bam file not created. sam2bam step failed during alignment."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
	    exit 1;
	fi       
        echo `date`
fi

