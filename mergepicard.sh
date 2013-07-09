#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 12 ] 
then
    MSG="parameter mismatch."
    echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
    exit 1;
else
    set -x
    echo `date`
    scriptfile=$0
    outputdir=$1
    infiles=$2
    outfilewdups=$3
    outfilenodups=$4
    tmpfilewdups=$5
    dupparms=$6
    RGparms=$7
    runfile=$8
    elog=$9
    olog=${10}
    email=${11}
    qsubfile=${12}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    #sanity check
    
    if [ ! -s $runfile ]
    then
       MSG="$runfile configuration file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    markdup=$( echo $dupparms | tr "_" "\n" | grep -w dup | cut -d '=' -f2 )
    deldup=$( echo $dupparms | tr "_" "\n" | grep -w flag | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
    threads=$( cat $runfile | grep -w PBSTHREADS | cut -d "=" -f2 )
    javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d "=" -f2 )
    if [ ! -d $outputdir ]
    then
       MSG="$outputdir output directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       MSG="PICARDIR=$picardir directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="SAMDIR=$samdir directory not found"
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

    if [ `expr length $infiles` -lt 1 ]
    then
       MSG="$infiles empty list of aligned files to merge"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    echo `date`
    listfiles=$( echo $infiles | tr ":" " " )
    header=$( echo $RGparms  | tr ":" "\t" )
    rgheader=$( echo -n -e "@RG\t" )$( echo $header  | tr "=" ":" )

    ## step 1: sort-merge and add RG tags all at once
    cd $outputdir
    sortedfiles=""
    for bamfile in listfiles
    do
         
       java -Xmx6g -Xms512m -jar $picardir/AddOrReplaceReadGroups.jar \
	   INPUT=$bamfile \
	   OUTPUT=wrg.$bamfile \
	   MAX_RECORDS_IN_RAM=null \
	   TMP_DIR=$outputdir \
	   SORT_ORDER=unsorted \
	   VALIDATION_STRINGENCY=SILENT
       exitcode=$?
       if [ $exitcode -ne 0 ]
       then
	   MSG="addorreplacereadgroups command failed.  exitcode=$exitcode mergepicard stopped"
	   echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	   exit $exitcode;
       fi
       echo `date`

       java -Xmx6g -Xms512m -jar $picardir/SortSam.jar \
	   INPUT=wrg.$bamfile \
	   OUTPUT=sorted.wrg.$bamfile \
	   TMP_DIR=$outputdir \
	   SORT_ORDER=coordinate \
	   MAX_RECORDS_IN_RAM=null \
	   CREATE_INDEX=true \
	   VALIDATION_STRINGENCY=SILENT
       exitcode=$?
       if [ $exitcode -ne 0 ]
       then
	   MSG="sortsam command failed.  exitcode=$exitcode mergepicard stopped"
	   echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	   exit $exitcode;
       fi
       echo `date`        
       sortedfiles=${sortedfiles}" INPUT=sorted.wrg.$bamfile"
    done

    java -Xmx6g -Xms512m -jar $picardir/MergeSamFiles.jar $sortedfiles \
        MAX_RECORDS_IN_RAM=null \
        TMP_DIR=$outputdir \
        OUTPUT=$tmpfilewdups \
        CREATE_INDEX=true \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=SILENT
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="mergesamfiles command failed.  exitcode=$exitcode mergepicard stopped"
	echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`
    if [ ! -s $tmpfilewdups ]
    then
        MSG="$tmpfilewdups file not created. mergepicard  failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
        exit 1;
    fi

    ## step 2: indexing merged bam file    
    $samdir/samtools index $tmpfilewdups
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="samtools index command failed.  exitcode=$exitcode mergepicard stopped"
	echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`        
    #$samdir/samtools flagstat $tmpfilewdups > $tmpfilewdups.flagstat
    $samdir/samtools view -H $tmpfilewdups > $tmpfilewdups.header
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="samtools index command failed.  exitcode=$exitcode mergepicard stopped"
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
	MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode mergepicard stopped"
	echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`        
    ## step 3: marking and or removing duplicates        

    if [ $markdup == "YES" -a $deldup != "TRUE" ]
    then
        echo "marking duplicates in sorted bam file"
        java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$tmpfilewdups \
	    OUTPUT=$outfilewdups \
	    TMP_DIR=$outputdir \
	    METRICS_FILE=$outfilewdups.dup.metrics \
            ASSUME_SORTED=true \
	    MAX_RECORDS_IN_RAM=null \
	    CREATE_INDEX=true \
	    VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode mergepicard stopped"
	    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`
	if [ ! -s $outfilewdups ]
	then
	    MSG="$outfilewdups file not created. mergepicard failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
        echo "indexing bam file w marked duplicates"
	$samdir/samtools index $outfilewdups
	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="samtools index command failed.  exitcode=$exitcode mergepicard stopped"
	    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`
	$samdir/samtools flagstat $outfilewdups > $outfilewdups.flagstat
	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="samtools flagstat command failed.  exitcode=$exitcode mergepicard stopped"
	    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`
	$samdir/samtools view -H $outfilewdups > $outfilewdups.header
	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="samtools view command failed.  exitcode=$exitcode mergepicard stopped"
	    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`
    else
        echo "remove duplicates or do nothing"
	if [ $deldup == "TRUE" ]
	then
            echo "removing marked duplicates in sorted bam file"
            java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
		INPUT=$tmpfilewdups \
		OUTPUT=$outfilenodups \
		TMP_DIR=$outputdir \
		METRICS_FILE=$outfilenodups.dup.metrics \
		MAX_RECORDS_IN_RAM=null \
		ASSUME_SORTED=true \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=SILENT

	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="markduplicates command failed.  exitcode=$exitcode mergepicard stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi

	    echo `date`
	    if [ ! -s $outfilenodups ]
	    then
		MSG="$outfilenodups file not created. mergepicard failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
            echo "indexing bam file w removed duplicates"
	    $samdir/samtools index $outfilenodups
	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="samtools index command failed.  exitcode=$exitcode mergepicard stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi

	    echo `date`
	    $samdir/samtools flagstat $outfilenodups > $outfilenodups.flagstat
	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="samtools flagstat command failed.  exitcode=$exitcode mergepicard stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi

	    echo `date`
	    $samdir/samtools view -H $outfilenodups > $outfilenodups.header
	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="samtools view command failed.  exitcode=$exitcode mergepicard stopped"
		echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi

	    echo `date`
	else
	    echo "we need to copy tmpfilewdups to outfilewdups now"
	    echo "these output files are not suitable for realignment"
	    `cp $tmpfilewdups $outfilewdups`
            `cp $tmpfilewdups.flagstat $outfilewdups.flagstat`
            `cp $tmpfilewdups.header $outfilewdups.header`
	    echo `date`
	fi
    fi
fi