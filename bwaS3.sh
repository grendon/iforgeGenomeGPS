#!/bin/sh

if [ $# != 12 ]
then
        MSG="parameter mismatch"
        echo -e "program=$0 stooped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;

else
	set -x
	echo `date`

        scriptfile=$0
        alignerdir=$1
        ref=$2
	outputdir=$3
        R1=$4
        A1=$5
        samfile=$6
        bamfile=$7
        samdir=$8
        elog=$9
        olog=${10}
        email=${11}
        qsubfile=${12}
        LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        cd $outputdir
        $alignerdir/bwa samse $ref $A1 $R1  > $outputdir/$samfile
        if [ ! -s $outputdir/$samfile ]
        then
            MSG="$outputdir/$samfile file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
            exit 1;
        fi
        echo `date`
        ## include sam2bam conversion
	$samdir/samtools view -bS -o $bamfile $samfile
	if [ ! -s $outputdir/$bamfile ]
	then
	    MSG="$outputdir/$bamfile file not created. sam2bam step failed."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
	fi       
        echo `date`
fi

