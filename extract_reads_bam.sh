#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# -le 7 -o $# -gt 10 ]
then
        MSG="parameter mismatch"
	echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
	exit 1;
else
	set -x
	echo `date`
	output=$1
	bam=$2
	run_info=$3
        elog=$4
        olog=$5
        email=$6
        qsubfile=$7
	igv=$8
        extradir=$9
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $run_info ]
        then
            MSG="$run_info configuration file not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        refdir=$( cat $run_info | grep -w '^REFGENOMEDIR' | cut -d '=' -f2)
	refgen=$( cat $run_info | grep -w '^REFGENOME' | cut -d '=' -f2)
	samtools=$( cat $run_info | grep -w '^SAMDIR' | cut -d '=' -f2)
        picard=$( cat $run_info | grep -w '^PICARDIR' | cut -d '=' -f2)
	script_path=$( cat $run_info | grep -w '^SCRIPTDIR' | cut -d '=' -f2 )
	chrindex=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | awk '{print "chr"$0}' )
        outputdir=$( cat $run_info | grep -w '^OUTPUTDIR' | cut -d '=' -f2)
	delivery=$( cat $run_info | grep -w '^DELIVERYFOLDER' | cut -d '=' -f2)
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]" )
        javamodule=$( cat $run_info | grep -w '^JAVAMODULE' | cut -d '=' -f2)

        if [ ! -d $refdir ]
        then
            MSG="$refdir reference genome directory not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        if [ ! -s $refdir/$refgen ]
        then
            MSG="$refdir/$refgen reference genome  not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        if [ ! -d $samtools ]
        then
            MSG="$samtools samtools directory not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        if [ ! -d $picard ]
        then
            MSG="$picard picard directory not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        if [ -z $javamodule ]
        then
            MSG="Value for JAVAMODULE must be specified in configuration file"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        else
            `/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
        fi
        if [ ! -d $outputdir ]
        then
            MSG="$outputdir  output directory not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        if [ ! -d $output ]
        then
            MSG="$output results directory not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi

        cd $output

        if [ ! -s $bam ]
        then
            MSG="$bam BAM file not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
	if [ ! -s $bam.bai ]
	then
	    echo `date`	
	    $samtools/samtools index $bam
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools view command failed.  exitcode=$exitcode extract reads failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
	    echo `date`	
	fi	
        if [ ! -d $extradir ]
        then
	    mkdir -p $extradir
        fi

	echo `date`
        ref=$refdir/$refgen
	chrs=`cat $ref.fai | cut -f1 | tr ":" "\n"`
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="extracting all chr from indexed reference genome failed.  exitcode=$exitcode extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
	echo `date`
	i=1
	for chr in $chrs
	do
	    if [ `echo $chrindex | grep -w "$chr" | wc -l` -eq 0 ]
	    then
		chrArray[$i]=$chr
		let i=i+1
	    fi
	done
	echo `date`

	## extract read for excluded chromosomes
	input=""
        cd $output
	for i in $(seq 1 ${#chrArray[@]})
	do
            echo "extracting next chr[${i}] from aligned bam"
            echo `date`
	    chr=${chrArray[$i]}
	    $samtools/samtools view -b $bam $chr > $bam.$chr.bam
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools view command failed.  exitcode=$exitcode extract reads failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
            if [ ! -s $bam.$chr.bam ]
            then
                MSG="$bam.$chr.bam file  not created. extract reads failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
	    fi
	    $samtools/samtools index $bam.$chr.bam
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools index command failed.  exitcode=$exitcode extract reads failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
	    input="$input INPUT=$bam.$chr.bam"
	done
        echo `date`

	### extract unmapped reads
	$samtools/samtools view -b -f 12 $bam > $bam.unmapped.bam
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="samtools view command failed.  exitcode=$exitcode extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
        if [ ! -s $bam.unmapped.bam ]
        then
	    MSG="$bam.unmapped.bam unmapped reads -file not created. extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
        echo `date`
	$samtools/samtools index $bam.unmapped.bam
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="samtools index command failed.  exitcode=$exitcode extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
	input="$input INPUT=$bam.unmapped.bam"
        echo `date`

        java -Xmx6g -Xms512m -jar $picard/MergeSamFiles.jar $input \
        MAX_RECORDS_IN_RAM=null \
        TMP_DIR=$outputdir \
        OUTPUT=$bam.extra.bam \
        CREATE_INDEX=true \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=SILENT

        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="mergesamfiles command failed.  exitcode=$exitcode extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
        if [ ! -s $bam.extra.bam ]
        then
	    MSG="Merged BAM file of extracted reads not created. extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        echo `date`        
        ## moving temporary files to extradir and results to delivery folder

        mv $bam.unmapped* $extradir/
        mv $bam.chr* $extradir/
	if [ $delivery != "NA" ]
	then
	    delivery_folder=$outputdir/$delivery
	    if [ ! -d $delivery_folder ]
	    then
                    mkdir $delivery_folder
                    mkdir $delivery_folder/IGV_BAM
            else
                if [ ! -d $delivery_folder/IGV_BAM ]
                then
                    mkdir $delivery_folder/IGV_BAM
                fi
            fi
	    mv $bam.extra.bam $delivery_folder/IGV_BAM
	    mv $bam.extra.bam.bai $delivery_folder/IGV_BAM
	else	
            if [ ! -d $igv ]
            then
                mkdir $igv
            fi
	    mv $bam.extra.bam $igv/
	    mv $bam.extra.bam.bai $igv/
	fi		
	echo `date`
fi 
