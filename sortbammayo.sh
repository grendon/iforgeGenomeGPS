#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 10 ] 
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
    set -x
    echo `date`

    outputdir=$1
    inbamfile=$2
    sortedplain=$3
    outfilewdups=$4
    outfilenodups=$5
    runfile=$6
    elog=$7
    olog=$8
    email=$9
    scriptfile=$0
    qsubfile=${10}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    #sanity check
    if [ ! -s $runfile ]
    then
       MSG="$runfile configuration file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    tmpfile=rg.$inbamfile
    sample=`basename $outputdir`
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
    markdup=$( cat $runfile | grep -w ^MARKDUP | cut -d '=' -f2 )
    deldup=$( cat $runfile | grep -w ^REMOVE_DUP | cut -d '=' -f2 )
    javamodule=$( cat $runfile | grep -w ^JAVAMODULE | cut -d '=' -f2 )
    sID=$sample
    sPU=$sample
    sSM=$sample
    sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )       
    sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )        
    sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
    RGparms=$( echo "RGID=${sID} RGLB=${sLB} RGPL=${sPL} RGPU=${sPU} RGSM=${sSM} RGCN=${sCN}" )

    if [ ! -d $outputdir ]
    then
       MSG="$outputdir output directory not found"
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
    if [ -z $javamodule ]
    then
       MSG="A value for JAVAMODULE has to be specified in configuration file."
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    else
        `/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    fi

    cd $outputdir

    if [ ! -s $inbamfile ]
    then
       MSG="$outputdir/$inbamfile BAM file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
  

    # step 1: checking to see if already sorted and w RG info


    java -Xmx6g -Xms512m -jar $picardir/AddOrReplaceReadGroups.jar \
	INPUT=$inbamfile \
	OUTPUT=$tmpfile \
	MAX_RECORDS_IN_RAM=null \
	TMP_DIR=$outputdir \
	SORT_ORDER=unsorted \
	$RGparms \
	VALIDATION_STRINGENCY=SILENT

    if [ ! -s $tmpfile ]
    then
	MSG="$tmpfile bam file not created. add_readGroup step failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi    
    echo `date`

    java -Xmx6g -Xms512m -jar $picardir/SortSam.jar \
	INPUT=$tmpfile \
	OUTPUT=$sortedplain \
	TMP_DIR=$outputdir \
	SORT_ORDER=coordinate \
	MAX_RECORDS_IN_RAM=null \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

    echo `date`
    if [ ! -s $sortedplain ]
    then
	MSG="$sortedplain sorted bam file not created. sort step failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    $samdir/samtools index $sortedplain
    $samdir/samtools flagstat $sortedplain > $sortedplain.flagstat
    $samdir/samtools view -H $sortedplain > $sortedplain.header
    rm $tmpfile

    echo `date`
        
    # marking and or removing duplicates        

    if [ $markdup == "YES" ]
    then
        echo "marking duplicates in sorted bam file"
        java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$sortedplain \
	    OUTPUT=$outfilewdups \
	    TMP_DIR=$outputdir \
	    METRICS_FILE=$outfilewdups.dup.metrics \
	    MAX_RECORDS_IN_RAM=null \
	    CREATE_INDEX=true \
	    VALIDATION_STRINGENCY=SILENT
	    
	echo `date`
	if [ ! -s $outfilewdups ]
	then
	    MSG="$outfilewdups file not created. markDuplicates step failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""	    
	    exit 1;
	fi
        echo "indexing bam file w marked duplicates"
	$samdir/samtools index $outfilewdups
	$samdir/samtools flagstat $outfilewdups > $outfilewdups.flagstat
	$samdir/samtools view -H $outfilewdups > $outfilewdups.header
    else
        echo "mark duplicates; job not requested"
        `cp $sortedplain $outfilewdups`
        `cp $sortedplain.flagstat $outfilewdups.flagstat`
        `cp $sortedplain.header $outfilewdups.header`
    fi

    if [ $deldup == "TRUE" ]
    then
        echo "removing marked duplicates in sorted bam file"
        java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$sortedplain \
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
	    MSG="$outfilenodups file not created. RemoveDuplicates step failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
        echo "indexing bam file w removed duplicates"
	$samdir/samtools index $outfilenodups
	$samdir/samtools flagstat $outfilenodups > $outfilenodups.flagstat
	$samdir/samtools view -H $outfilenodups > $outfilenodups.header
    else
        echo "remove duplicates; job not requested"
    fi
    echo `date`
fi