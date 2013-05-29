#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu
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
        $alignerdir/bwa samse $ref $A1 $R1  > $outputdir/$samfile
        if [ ! -s $outputdir/$samfile ]
        then
            MSG="$outputdir/$samfile aligned file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline - Support #200' "$redmine,$email""
            exit 1;
        fi
        echo `date`
        ## include sam2bam conversion
	$samdir/samtools view -bS -o $bamfile $samfile
	if [ ! -s $outputdir/$bamfile ]
	then
	    MSG="$outputdir/$bamfile bam file not created. sam2bam step failed during alignment."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline - Support #200' "$redmine,$email""
	    exit 1;
	fi       
        echo `date`
fi

