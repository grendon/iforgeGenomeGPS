#!/bin/sh
###############################
# input files are sorted bam files with readgroup info added to each file
# files were the result of splitting input file prior to running alignment
###############################
if [ $# != 11 ] 
then
    MSG="parameter mismatch."
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
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
    LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    #sanity check
    
    if [ ! -s $runfile ]
    then
       MSG="$runfile file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
       exit 1;
    fi
    markdup=$( echo $dupparms | tr "_" "\n" | grep -w dup | cut -d '=' -f2 )
    deldup=$( echo $dupparms | tr "_" "\n" | grep -w flag | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )

    if [ ! -d $outputdir ]
    then
       MSG="$outputdir directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       MSG="$picardir directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="$samdir directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
       exit 1;
    fi
    if [ `expr length $infiles` -lt 1 ]
    then
       MSG="$infiles empty list of files to merge"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
       exit 1;
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

    #$samdir/samtools merge $outfilewdups $infiles
    echo `date`
    if [ ! -s $tmpfilewdups ]
    then
        MSG="$tmpfilewdups file not created. merge step failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
        exit 1;
    fi
    echo "indexing merged bam file"    
    $samdir/samtools index $tmpfilewdups
    $samdir/samtools flagstat $tmpfilewdups > $tmpfilewdups.flagstat
    $samdir/samtools view -H $tmpfilewdups > $tmpfilewdups.header
    echo `date`
        
    # marking and or removing duplicates        

    if [ $markdup == "YES" ]
    then
        echo "marking duplicates in sorted bam file"
        java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$tmpfilewdups \
	    OUTPUT=$outfilewdups \
	    TMP_DIR=$outputdir \
	    METRICS_FILE=$outfilewdups.dup.metrics \
	    MAX_RECORDS_IN_RAM=null \
	    CREATE_INDEX=true \
	    VALIDATION_STRINGENCY=SILENT
	    
	echo `date`
	if [ ! -s $outfilewdups ]
	then
	    MSG="$outfilewdups file not created. markDups step failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
	fi
        echo "indexing bam file w marked duplicates"
	$samdir/samtools index $outfilewdups
	$samdir/samtools flagstat $outfilewdups > $outfilewdups.flagstat
	$samdir/samtools view -H $outfilewdups > $outfilewdups.header
    else
	echo "we need to copy tmpfilewdups to outfilewdups now"
	echo "in case realignment follows"
	`cp $tmpfilewdups $outfilewdups`
        `cp $tmpfilewdups.flagstat $outfilewdups.flagstat`
        `cp $tmpfilewdups.header $outfilewdups.header`
    fi

    if [ $deldup == "TRUE" ]
    then
        echo "removing marked duplicates in sorted bam file"
        java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$tmpfilewdups \
	    OUTPUT=$outfilenodups \
	    TMP_DIR=$outputdir \
	    METRICS_FILE=$outfilenodups.dup.metrics \
	    MAX_RECORDS_IN_RAM=null \
	    CREATE_INDEX=true \
	    REMOVE_DUPLICATES=true \
	    VALIDATION_STRINGENCY=SILENT
	echo `date`
	if [ ! -s $outfilenodups ]
	then
	    MSG="$outfilenodups file not created. Remove markDups step failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
	fi
        echo "indexing bam file w removed duplicates"
	$samdir/samtools index $outfilenodups
	$samdir/samtools flagstat $outfilenodups > $outfilenodups.flagstat
	$samdir/samtools view -H $outfilenodups > $outfilenodups.header
    fi
    echo `date`
fi