#!/bin/sh
###############################
# input files are sorted bam files with readgroup info added to each file
# files were the result of splitting input file prior to running alignment
###############################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 11 ] 
then
    MSG="parameter mismatch."
    echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
    exit 1;
else
    set -x
    echo `date`

    outputdir=$1
    infiles=$2
    outfilewdups=$3
    outfilenodups=$4
    tmpfilewdups=$5
    dupparms=$6
    runfile=$7
    elog=$8
    olog=$9
    email=${10}
    scriptfile=$0
    qsubfile=${11}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    #sanity check
    
    if [ ! -s $runfile ]
    then
       MSG="$runfile file not found"

       exit 1;
    fi
    markdup=$( echo $dupparms | tr "_" "\n" | grep -w dup | cut -d '=' -f2 )
    deldup=$( echo $dupparms | tr "_" "\n" | grep -w flag | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
    javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d "=" -f2 )
    if [ ! -d $outputdir ]
    then
       MSG="$outputdir results directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       MSG="$picardir picard directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="$samdir samtools directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ `expr length $infiles` -lt 1 ]
    then
       MSG="$infiles empty list of files to merge"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ -z $javamodule ]
    then
       MSG="Value for JAVAMODULE must be specified in configuration file"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    else
	`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    fi
     
    listfiles=$( echo $infiles | tr ":" " " )

    cd $outputdir

    echo "found sorted bam files to be merged. all set to go..."
    java -Xmx6g -Xms512m -jar $picardir/MergeSamFiles.jar $listfiles \
        MAX_RECORDS_IN_RAM=null \
        TMP_DIR=$outputdir \
        OUTPUT=$tmpfilewdups \
        CREATE_INDEX=true \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=SILENT

    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
        MSG="mergesamfiles command failed.  exitcode=$exitcode mergenodes failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
        exit $exitcode;
    fi
    echo `date`

    if [ ! -s $tmpfilewdups ]
    then
        MSG="$tmpfilewdups file not created. picard-merge step failed mergenodes failed "
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email"" 
        exit 1;
    fi

    echo "indexing merged bam file"    
    $samdir/samtools index $tmpfilewdups
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
         MSG="samtools index command failed.  exitcode=$exitcode mergenodes stopped"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         exit $exitcode;
    fi
    echo `date`
    #$samdir/samtools flagstat $tmpfilewdups > $tmpfilewdups.flagstat
    $samdir/samtools view -H $tmpfilewdups > $tmpfilewdups.header
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
         MSG="samtools view command failed.  exitcode=$exitcode mergenodes stopped"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         exit $exitcode;
    fi
    echo `date`

    java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
        INPUT=$tmpfilewdups \
        OUTPUT=$tmpfilewdups.flagstat \
        VALIDATION_STRINGENCY=SILENT

     exitcode=$?
     if [ $exitcode -ne 0 ]
     then
         MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode mergenodes stopped"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         exit $exitcode;
     fi
    echo `date`
        
    # marking and or removing duplicates        

    if [ $markdup == "YES" -a $deldup != "TRUE" ]
    then
        echo "marking duplicates in sorted aligned bam file"
        java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$tmpfilewdups \
	    OUTPUT=$outfilewdups \
	    TMP_DIR=$outputdir \
	    METRICS_FILE=$outfilewdups.dup.metrics \
	    MAX_RECORDS_IN_RAM=null \
	    CREATE_INDEX=true \
	    VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
            MSG="markduplicates command failed.  exitcode=$exitcode mergenodes stopped"
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
	fi
	    
	echo `date`
	if [ ! -s $outfilewdups ]
	then
	    MSG="$outfilewdups file not created. markDuplicates step failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
        echo "indexing bam file w marked duplicates"
	$samdir/samtools index $outfilewdups
	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
            MSG="samtools index command failed.  exitcode=$exitcode mergenodes stopped"
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
	fi
	echo `date`
	#$samdir/samtools flagstat $outfilewdups > $outfilewdups.flagstat
	$samdir/samtools view -H $outfilewdups > $outfilewdups.header
	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
            MSG="samtools view command failed.  exitcode=$exitcode mergenodes stopped"
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
	fi
	echo `date`
	java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
            INPUT=$outfilewdups \
            OUTPUT=$outfilewdups.flagstat \
            VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
            MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode mergenodes stopped"
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
	fi
	echo `date`
    else
        echo "remove duplicates or do nothing"
	if [ $deldup == "TRUE" ]
	then
            echo "removing marked duplicates in sorted aligned bam file"
            java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
		INPUT=$tmpfilewdups \
		OUTPUT=$outfilenodups \
		TMP_DIR=$outputdir \
		METRICS_FILE=$outfilenodups.dup.metrics \
		MAX_RECORDS_IN_RAM=null \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=SILENT

	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="markduplicates command failed.  exitcode=$exitcode mergenodes stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi
	    
	    echo `date`
	    if [ ! -s $outfilenodups ]
	    then
		MSG="$outfilenodups file not created. RemoveDuplicates step failed.  mergenodes failed "
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
            echo "indexing bam file w removed duplicates"
	    $samdir/samtools index $outfilenodups
	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="samtools index command failed.  exitcode=$exitcode mergenodes stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi
	    
	    echo `date`
	    #$samdir/samtools flagstat $outfilenodups > $outfilenodups.flagstat
	    $samdir/samtools view -H $outfilenodups > $outfilenodups.header
	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="samtools view command failed.  exitcode=$exitcode mergenodes stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi
	    
	    echo `date`
	    java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
		INPUT=$outfilenodups \
		OUTPUT=$outfilenodups.flagstat \
		VALIDATION_STRINGENCY=SILENT

	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode mergenodes stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi
	    echo `date`
	else
	    echo "we need to copy tmpfilewdups to outfilewdups now"
	    echo "duplicates were neither marked nor removed from sorted aligned bam file"
            echo "Warning: this output file is not suitable for realignment recalibration"
	    `cp $tmpfilewdups $outfilewdups`
            `cp $tmpfilewdups.flagstat $outfilewdups.flagstat`
            `cp $tmpfilewdups.header $outfilewdups.header`
	    echo `date`
        fi
    fi
fi
