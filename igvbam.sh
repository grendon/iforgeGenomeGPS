#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 7 ] 
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
    runfile=$3
    elog=$4
    olog=$5
    email=$6
    qsubfile=$7
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
    if [ ! -s $runfile ]
    then
       MSG="$runfile configuration file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

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
       MSG="$infiles empty list of realigned bam files to merge"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    echo `date`
    listfiles=$( echo $infiles | tr ":" " " )
    outfile=IGV.sorted.bam

    cd $outputdir

    java -Xmx6g -Xms512m -jar $picardir/MergeSamFiles.jar $listfiles \
        MAX_RECORDS_IN_RAM=null \
        TMP_DIR=$outputdir \
        OUTPUT=$outfile \
        CREATE_INDEX=true \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=SILENT

     exitcode=$?
     if [ $exitcode -ne 0 ]
     then
         MSG="mergesamfiles command failed.  exitcode=$exitcode igvbam stopped"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         exit $exitcode;
     fi
    echo `date`
    if [ ! -s $outfile ]
    then
        MSG="$outfile  igvbam file not created."
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
        exit 1;
    fi

    java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
        INPUT=$outfile \
        OUTPUT=$outfile.flagstat \
        VALIDATION_STRINGENCY=SILENT

     exitcode=$?
     if [ $exitcode -ne 0 ]
     then
         MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode igvbam stopped"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         exit $exitcode;
     fi
    echo `date`
fi