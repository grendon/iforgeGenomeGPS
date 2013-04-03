#!/bin/sh
if [ $# != 11 ]
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;
else
    set -x
    echo `date`
    picardir=$1
    samdir=$2
    outputdir=$3
    infile=$4
    outfile=$5
    rgparms=$6
    alignflag=$7
    elog=$8
    olog=$9
    email=${10}
    scriptfile=$0
    qsubfile=${11}
    LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    tmpfile=tmp.wrg.$infile
    sample=`basename $outputdir`
    parameters=$( echo $rgparms | tr ":" " " )

    #sanity check
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
    if [ ! -s $outputdir/$infile ]
    then
       MSG="$infile file to be sorted not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
       exit 1;
    fi

    cd $outputdir

    ## before sorting, we need to make sure the bam file has readgroup info

    if [ $alignflag == "NCSA" ]
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
       echo "alignment was NOT done inhouse. checking if readgroup info is present"
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
	MSG="$tmpfile file not created. add_readGroup step failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
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
	MSG="$outfile not created. sort step failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	exit 1;
    fi
    $samdir/samtools index $outfile
    $samdir/samtools view -H $outfile > $outfile.header
fi