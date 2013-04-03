#!/bin/sh

if [ $# != 14 ]
then
        MSG="parameter mismatch"
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        alignerdir=$1
        ref=$2
	outputdir=$3
        R1=$4
        R2=$5
        A1=$6
        A2=$7
        samfile=$8
        bamfile=$9
        samdir=${10}
        elog=${11}
        olog=${12}
        email=${13}
        qsubfile=${14}
        LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        cd $outputdir
        $alignerdir/bwa sampe $ref $A1 $A2 $R1 $R2 > $outputdir/$samfile
        if [ ! -s $outputdir/$samfile ]
        then
            MSG="$outputdir/$samfile file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
            exit 1;
        fi
        echo `date`
        ## sam2bam conversion
	$samdir/samtools view -bS -o $bamfile $samfile
	if [ ! -s $outputdir/$bamfile ]
	then
	    MSG="$outputdir/$bamfile file not created. sam2bam step failed."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
	fi       
        echo `date`
fi

