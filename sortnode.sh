#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# -gt 14 ]
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
    set -x
    echo `date`
    scriptfile=$0
    picardir=$1
    samdir=$2
    javamodule=$3
    outputdir=$4
    bamfile=$5
    infile=$6
    outfile=$7
    rgparms=$8
    flag=$9
    chr=${10}
    elog=${11}
    olog=${12}
    email=${13}
    qsubfile=${14}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    if [ ! -d $outputdir ]
    then
       MSG="$outputdir realign directory not found"
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
       MSG="$picardir samtools directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    if [ -z $javamodule ]
    then
	MSG="A value must be specified for JAVAMODULE in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    else
        `/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    fi          
    if [ ! -s $outputdir/$bamfile ]
    then
	MSG="$bamfile bam file to be sorted was not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    if [ ! -s $outputdir/$infile ]
    then
        cd $outputdir
        $samdir/samtools view -b $bamfile $chr > $infile
        if [ ! -s $infile ]
        then
	    MSG="$infile bam file to be sorted within region:[$chr] was not created"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        $samdir/samtools index $infile
    fi

    tmpfile=tmp.wrg.$infile
    parameters=$( echo $rgparms | tr ":" " " )
    sortflag=$( echo $flag | tr '[a-z]' '[A-Z]' )
    cd $outputdir

    ## before sorting, we need to make sure the bam file has readgroup info

    if [ $sortflag == "NCSA" ]
    then
       echo "alignment was done inhouse. we need to add_readgroup info"
       java -Xmx6g -Xms512m -jar $picardir/AddOrReplaceReadGroups.jar \
	   INPUT=$infile \
	   OUTPUT=$tmpfile \
	   MAX_RECORDS_IN_RAM=null \
	   TMP_DIR=$outputdir \
	   SORT_ORDER=unsorted \
           $parameters \
	   VALIDATION_STRINGENCY=SILENT
    else
	   echo "alignment was done at Mayo. checking if readgroup info is present"
	   $samdir/samtools view -H $infile > $infile.header
	   match=$( cat $file.header | grep '^@RG' )
	   lenmatch=`expr length $match`
	   if [ $lenmatch -gt 0 ]
	   then
               echo "readgroup info found in input file."
               cp $infile $tmpfile
	   else
               echo "readgroup info NOT found in input file. Adding it now..."
	       java -Xmx6g -Xms512m -jar $picardir/AddOrReplaceReadGroups.jar \
		   INPUT=$infile \
		   OUTPUT=$tmpfile \
		   MAX_RECORDS_IN_RAM=null \
		   TMP_DIR=$outputdir \
		   SORT_ORDER=unsorted \
		   $parameters \
		   VALIDATION_STRINGENCY=SILENT
	   fi
    fi

    if [ ! -s $tmpfile ]
    then
	MSG="$tmpfile bam file not created. add_readGroup step failed during realignment"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi
    echo `date`

    java -Xmx6g -Xms512m -jar $picardir/SortSam.jar \
	INPUT=$tmpfile \
	OUTPUT=$outfile \
	TMP_DIR=$outputdir \
	SORT_ORDER=coordinate \
	MAX_RECORDS_IN_RAM=null \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT
    echo `date`
    if [ ! -s $outfile ]
    then
	MSG="$outfile sort bam file not created. sort step failed during realignment"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi
    $samdir/samtools index $outfile
    $samdir/samtools view -H $outfile > $outfile.header
fi