#!/bin/sh
###############################
# input files are sorted bam files with readgroup info added to each file
# there files were not aligned inhouse 
###############################
if [ $# != 8 ] 
then
    echo "MERGEBAM Usage: parameter mismatch";
else
    set -x
    echo `date`

    inputdir=$1
    outputdir=$2
    outputlog=$3
    listinfiles=$4
    listinputfiles=$5
    outfilewdups=$6
    outfilenodups=$7
    runfile=$8

    #sanity check
    if [ ! -s $runfile ]
    then
       echo "$runfile file not found"
       exit 1;
    fi
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    refdir=$( cat $runfile | grep -w REFDIR | cut -d '=' -f2 )
    refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
    dup=$( echo $dupparms | tr "_" "\n" | grep -w dup | cut -d '=' -f2 )
    dupflag=$( echo $dupparms | tr "_" "\n" | grep -w flag | cut -d '=' -f2 )
    javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )

    if [ ! -d $outputdir ]
    then
       echo "$outputdir directory not found"
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       echo "$picardir directory not found"
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       echo "$samdir directory not found"
       exit 1;
    fi
    if [ -z $javamodule ]
    then
       MSG="Value for JAVAMODULE must be specified in configuration file"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'Mayo variant identification pipeline - Support #200' "$redmine,$email""
       exit 1;
    else
	`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    fi

    cd $outputdir
 
    # checking to see that all files to be merged exist
    for fileline in $infiles
    do 
       filename=$( echo $filename | grep -w "=" -f 2)
       if [ ! -s $filename ]
       then
	   echo "$filename file not found. merge step failed"
	   exit 1;
       fi
    done

    # ready to merge sorted bam files
    echo "found files to be merged. all set to go"
    java -Xmx6g -Xms512m -jar $picardir/MergeSamFiles.jar $infiles \
        MAX_RECORDS_IN_RAM=null \
        TMP_DIR=$outputdir \
        OUTPUT=$outfilewdups \
        CREATE_INDEX=true \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=SILENT

    #$samdir/samtools merge $outfilewdups $infiles
    echo `date`
    if [ ! -s $outfilewdups ]
    then
        echo "$outfilewdups file not created. merge step failed"
        exit 1;
    fi
    echo "indexing merged bam file"    
    $samdir/samtools index $outfilewdups
    $samdir/samtools flagstat $outfilewdups > $outfilewdups.flagstat
    echo `date`
        
    # marking and or removing duplicates        

    if [ $dup == "YES" ]
    then
        echo "marking duplicates in sorted bam file"
        if [ $dupflag == "YES" ]
        then
            echo "removing marked duplicates too"

            java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
		INPUT=$outfilewdups \
		OUTPUT=$outfilenodups \
		TMP_DIR=$outputdir \
		METRICS_FILE=$outfilewdups.dup.metrics \
		MAX_RECORDS_IN_RAM=null \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=SILENT
        else
            java -Xmx6g -Xms512m -jar $picardir/MarkDuplicates.jar \
		INPUT=$outfilewdups \
		OUTPUT=$outfilenodups \
		TMP_DIR=$outputdir \
		METRICS_FILE=$outfilewdups.dup.metrics \
		MAX_RECORDS_IN_RAM=null \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT
        fi
	if [ ! -s $outfilenodups ]
	then
            echo "$outfilenodups file not created. markDups step failed"
            exit 1;
	fi
	echo `date`
        echo "indexing bam file w mark duplicates or removed duplicates"
	$samdir/samtools index $outfilenodups
	$samdir/samtools flagstat $outfilenodups > $outfilenodups.flagstat
    else
        echo "duplicates not marked and or removed in sorted bam file. $outfilenodups not created"
    fi

    ## checking to see if reordering is necessary
    ## this is the last step in processBAM.sh which does not make sense to me
    if [ $reorderflag == "YES" ]
    then
	if [ ! -d $refdir ]
	then
	    echo "$refdir directory not found"
	    exit 1;
	fi
	if [ ! -s $refdir/$refgenome ]
	then
	    echo "$refdir/refgenome reference genome not found"
	    exit 1;
	fi

	java -Xmx6g -Xms512m -jar $picardir/ReorderSam.jar \
            INPUT=$outfilewdups \
            OUTPUT=$outfilewdups.tmp.bam \
	    REF=$refdir/$refgenome \
            MAX_RECORDS_IN_RAM=null \
            TMP_DIR=$outputdir \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT
    else
       echo "no reordering is necessary"
    fi
    echo "done merging bam files"
    echo `date`
fi