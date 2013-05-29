#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# -gt 15 ]
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        alignerdir=$1
        params=$2
        ref=$3
	outputdir=$4
        samfile=$5
        bamfile=$6
        scriptdir=$7
        samdir=$8
        paired=$9
        R1=${10}

        parameters=$( echo $params | tr "_" " " )

        ## checking quality scores to gather additional params
        qscores=$scriptdir/checkFastqQualityScores.pl
        ill2sanger=`perl $qscores $R1 10000`
        if [ $ill2sanger -gt 65 ]
        then
           qual="-F ILMFQ"
        else
           qual="-F STDFQ"
        fi

        cd $outputdir
        if [ $paired -eq 1 ]
        then
           R2=${11}
           elog=${12}
           olog=${13}
           email=${14}
           qsubfile=${15}
           $alignerdir/novoalign -d $ref -f $R1 $R2 -o SAM $parameters $qual > $outputdir/$samfile
        else
           elog=${11}
           olog=${12}
           email=${13}
           qsubfile=${14}
           $alignerdir/novoalign -d $ref -f $R1 -o SAM $parameters $qual > $outputdir/$samfile
        fi

        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        if [ ! -s $outputdir/$samfile ]
        then
           MSG="$outputdir/$samfile aligned file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        echo `date`

        ## sam2bam conversion
	$samdir/samtools view -bS -o $bamfile $samfile
	if [ ! -s $bamfile ]
	then
	    MSG="$outputdir/$bamfile BAM file not created. sam2bam step failed during alignment."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi       
        echo `date`
fi