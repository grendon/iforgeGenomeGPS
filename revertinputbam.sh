#!/bin/sh
######################################
#  script to revertsam on input bamfiles that need to be realigned-recal
#
######################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 8 ];
then
	MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        infile=$1
        outfile=$2
        picardir=$3
        samdir=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $infile ]
        then
	    MSG="$infile input bam file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $picardir ]
        then
	    MSG="PICARDIR=$picardir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      

        if [ ! -d $samdir ]
        then
	    MSG="PICARDIR=$picardir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      
	outputdir=`dirname $outfile`
        cd $outputdir
	echo `date`
        java -Xmx6g -Xms512m -jar $picardir/RevertSam.jar \
            COMPRESSION_LEVEL=0 \
            INPUT=$infile \
            OUTPUT=${infile}.sam \
            SORT_ORDER=queryname \
            VALIDATION_STRINGENCY=SILENT

        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="revertsam command failed. exitcode=$exitcode "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
	echo `date`

	$samdir/samtools view -bS -o $outfile ${infile}.sam

        if [ ! -s $outfile ]
        then
	    MSG="$outfile file not created. revertsam failed"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
fi
