#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu

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
	if [ ${10} ]
	then
	    group=${10}
	fi

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
        if [ ! -s $output/$bam ]
        then
            MSG="$output/$bam BAM file not found"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        if [ ! -d $extradir ]
        then
	    mkdir -p $extradir
        fi
	if [ ! -s $output/$bam.bai ]
	then
                cd $output
		$samtools/samtools index $output/$bam
	fi	
	
        ref=$refdir/$refgen
	chrs=`cat $ref.fai | cut -f1 | tr ":" "\n"`
	i=1
	for chr in $chrs
	do
	    if [ `echo $chrindex | grep -w "$chr" | wc -l` -eq 0 ]
	    then
		chrArray[$i]=$chr
		let i=i+1
	    fi
	done

	## extract read for specific chromosome
	input=""
        cd $output
	for i in $(seq 1 ${#chrArray[@]})
	do
	    chr=${chrArray[$i]}
	    $samtools/samtools view -b $output/$bam $chr > $output/$bam.$chr.bam
            if [ ! -s $output/$bam.$chr.bam ]
            then
                MSG="Warning:$output/$bam.$chr.bam file for chr=$chr not created. extract reads failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
	    fi
	    $samtools/samtools index $output/$bam.$chr.bam
	    input="$input INPUT=$output/$bam.$chr.bam"
	done

	### extract unmapped reads
	$samtools/samtools view -b -f 12 $output/$bam > $output/$bam.unmapped.bam
        if [ ! -s $output/$bam.unmapped.bam ]
        then
	    MSG="Warning: $output/$bam.unmapped.bam unmapped reads -file not created. extract reads failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
	$samtools/samtools index $output/$bam.unmapped.bam
	input="$input INPUT=$output/$bam.unmapped.bam"

        java -Xmx6g -Xms512m -jar $picard/MergeSamFiles.jar $input \
        MAX_RECORDS_IN_RAM=null \
        TMP_DIR=$outputdir \
        OUTPUT=$output/$bam.extra.bam \
        CREATE_INDEX=true \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=SILENT

        if [ ! -s $output/$bam.extra.bam ]
        then
	    MSG="Warning: BAM file of extracted reads, step failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        
        ## moving temporary files to extradir and results to delivery folder

        mv $output/$bam.unmapped* $extradir/
        mv $output/$bam.chr* $extradir/

	if [ $group ]
	then
            echo "group info was passed as argument to this script"


	    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	    samples=$(cat $run_info | grep -w "^$group" | cut -d '=' -f2 | tr "\t" " ")	
	    for sample in $samples
	    do	
		sam=`echo $samples | tr ":" "\n"| grep -v "$sample" | tr "\n" " "`
		gr=""
		for s in $sam
		do
		    a="ID:$s|";
		    gr="$gr $a"
		done
		gr=`echo $gr |  sed "s/|$//"`
		$samtools/samtools view -b -r $sample $output/$bam.extra.bam > $output/$sample.extra.bam
                if [ ! -s $output/$sample.extra.bam ]
                then
                    MSG="Warning: $output/$sample.extra.bam BAM file with extracted reads not created. extract reads failed."
		    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
                fi
		$samtools/samtools view -H $output/$sample.extra.bam | grep -E -v "$gr" | $samtools/samtools reheader - $output/$sample.extra.bam > $output/$sample.extra.re.bam
                if [ ! -s $output/$sample.extra.re.bam ]
                then
                    MSG="Warning: $output/$sample.extra.re.bam new header of BAM file with extracted reads not created. extract reads failed."
		    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
                fi
		mv $output/$sample.extra.re.bam $output/$sample.extra.bam
		$samtools/samtools index $output/$sample.extra.bam

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
		    mv $output/$bam.extra.bam $delivery_folder/IGV_BAM
		    mv $output/$bam.extra.bam.bai $delivery_folder/IGV_BAM
		else	
                    if [ ! -d $igv ]
                    then
			mkdir $igv
                    fi
		    mv $output/$bam.extra.bam $igv/
		    mv $output/$bam.extra.bam.bai $igv/
		fi    
	    done
	    #rm $output/$bam.extra.bam $output/$bam.extra.bam.bai
	else
            echo "the group argument was not passed to this script"
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
		mv $output/$bam.extra.bam $delivery_folder/IGV_BAM
		mv $output/$bam.extra.bam.bai $delivery_folder/IGV_BAM
	    else	
                if [ ! -d $igv ]
                then
                    mkdir $igv
                fi
		mv $output/$bam.extra.bam $igv/
		mv $output/$bam.extra.bam.bai $igv/
	    fi		
	fi
	echo `date`
fi 
