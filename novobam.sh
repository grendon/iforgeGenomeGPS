#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# -gt 13 ]
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
        samdir=$7
        paired=$8
        BAM=$9
        elog=${10}
        olog=${11}
        email=${12}
        qsubfile=${13}

        parameters=$( echo $params | tr "_" " " )
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        cd $outputdir
        if [ $paired == "1" ]
        then
	    $alignerdir/novoalign -F BAMPE -f $BAM -d $ref -o SAM $parameters  > $samfile
        else
	    $alignerdir/novoalign -F BAMSE -f $BAM -d $ref -o SAM $parameters > $samfile
        fi
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
           MSG="novoalign command failed.  exitcode=$exitcode alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit $exitcode;
        fi
	echo `date`

        if [ ! -s $samfile ]
        then
           MSG="$samfile aligned BAM file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

	$samdir/samtools view -bS -o $bamfile $samfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
           MSG="samtools view command failed.  exitcode=$exitcode  novobam alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit $exitcode;
        fi
        echo `date`
	if [ ! -s $bamfile ]
	then
	    MSG="$bamfile BAM file not created. sam2bam step failed during  novobam alignment."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi       

        echo `date`
        $samdir/samtools index $bamfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
           MSG="samtools index command failed.  exitcode=$exitcode  novobam alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit $exitcode;
        fi
        echo `date`

fi