#!/bin/sh
#	
#  script to realign and recalibrate the aligned file(s) 
#  This module is called from within the realign module
#  Input file(s) is(are) bam format 
#  given a specified region, usually one or more chromosomes at once
#  one or several samples can be considered in the input
######################################
	
if [ $# != 13 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;
else					
	set -x

	echo `date`
	
        sampledir=$1
        chr=$2
        chrinfiles=$3
        chrinputfiles=$4
        region=$5
        realparms=$6
        recalparms=$7
        runfile=$8
        flag=$9
	elog=${10}
	olog=${11}
	email=${12}
	scriptfile=$0
        qsubfile=${13}
	LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        #sanity check
        if [ ! -s $runfile ]
        then
	    MSG="$runfile file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
        targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
        realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
	outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        thr=`expr $threads "-" 1`

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
        if [ ! -d $gatk ]
        then
	    MSG="$gatk directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        # cleaning up the lists
        chrinfiles=$( echo $chrinfiles | tr ":" " " )
        chrinputfiles=$( echo $chrinputfiles | tr ":" " " )
        region=$( echo $region | tr ":" " " )
        realparms=$( echo $realparms | tr ":" " " )
        recalparms=$( echo $recalparms | tr ":" " " )
        inputdir=$outputrootdir/realign/$sampledir
        if [ ! -d $inputdir ]
        then
	    MSG="$inputdir directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        cd $inputdir

        if [ $flag == 1 ]
        then
		echo "realign then recalibrate"
		echo "realigning..."
		
            	java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $chrinfiles \
		    -T RealignerTargetCreator \
		    -o $sampledir.$chr.bam.list $realparms $region
	
		if [ ! -s $sampledir.$chr.bam.list ]
		then
		    MSG="$sampledir.$chr.bam.list file not created."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
                    exit 1;
		fi
		echo `date`


		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $chrinfiles \
		    -T IndelRealigner \
                    -L $chr \
		    -o $sampledir.$chr.realigned.bam \
		    -targetIntervals $sampledir.$chr.bam.list $realignparams $realparms

		if [ ! -s $sampledir.$chr.realigned.bam ]
		then
		    MSG="realigned.bam file not created."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
                    exit 1;
                else
		    mv $sampledir.$chr.realigned.bai $sampledir.$chr.realigned.bam.bai
		    cp $sampledir.$chr.realigned.bam $sampledir.$chr.real.cleaned.bam
		    cp $sampledir.$chr.realigned.bam.bai $sampledir.$chr.real.cleaned.bam.bai
                    $samdir/samtools flagstat $sampledir.$chr.real.cleaned.bam > $sampledir.$chr.real.cleaned.flagstat
		fi
		echo `date`

		echo "recalibrating..."
	
                # next step FixMateInformation is not done at Mayo..."
		#java -Xmx6g -Xms512m  -jar $picardir/FixMateInformation.jar \
		#    INPUT=$sampledir.$chr.real.cleaned.bam \
		#    OUTPUT=$sampledir.$chr.real.fixed.bam \
		#    SO=coordinate \
		#    TMP_DIR=$inputdir \
		#    VALIDATION_STRINGENCY=SILENT \
		#    CREATE_INDEX=true

		#if [ ! -s $sampledir.$chr.real.fixed.bam ]
		#then
		#    echo "$sampledir.$chr.realigned.fixed.bam file not created"
		#    exit 1;
		#fi
		#echo `date`

		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $recalparms \
		    -I $sampledir.$chr.real.cleaned.bam \
		    $region  \
		    -T CountCovariates \
		    -nt $thr \
		    -cov ReadGroupCovariate \
		    -cov QualityScoreCovariate \
		    -cov CycleCovariate \
		    -cov DinucCovariate \
		    -recalFile $sampledir.$chr.recal_data.csv 
	
		if [ ! -s $sampledir.$chr.recal_data.csv ]
		then
		    MSG="$sampledir.$chr.recal_data.csv file not created"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
		fi
		echo `date`

		java -Xmx6g -Xms512m  -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -L $chr \
		    -I $sampledir.$chr.real.cleaned.bam \
		    -T TableRecalibration \
		    --out $sampledir.$chr.real.recal.bam \
		    -recalFile $sampledir.$chr.recal_data.csv 

		if [ ! -s $sampledir.$chr.real.recal.bam ]
		then
		    MSG="$chr.real.recal.bam file not created"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
		fi
		cp $sampledir.$chr.real.recal.bam $sampledir.$chr.recal.cleaned.bam
		cp $sampledir.$chr.real.recal.bai $sampledir.$chr.recal.cleaned.bam.bai
		$samdir/samtools flagstat $sampledir.$chr.recal.cleaned.bam > $sampledir.$chr.recal.cleaned.flagstat
 		echo `date`
        else
		echo "recalibrate then realign"
		echo "recalibrating"
                echo "I am not sure that this block will ever be run"


		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $recalparms \
		    $chrinfiles \
		    $region \
		    -T CountCovariates \
		    -nt $thr \
		    -cov ReadGroupCovariate \
		    -cov QualityScoreCovariate \
		    -cov CycleCovariate \
		    -cov DinucCovariate \
		    -recalFile $sampledir.$chr.recal_data.csv 
	
		if [ ! -s $sampledir.$chr.recal_data.csv ]
		then
		    MSG="$sampledir.$chr.recal_data.csv file not created"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
		fi
	
		echo `date`

		java -Xmx6g -Xms512m  -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -L $chr \
		    $chrinfiles \
		    -T TableRecalibration \
		    --out $sampledir.$chr.recal.bam \
		    -recalFile $sampledir.$chr.recal_data.csv 

		if [ ! -s $sampledir.$chr.recal.bam ]
		then
		    MSG="$sampledir.$chr.recal.bam file not created"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
                else 
                    cp $sampledir.$chr.recal.bam $sampledir.$chr.recal.cleaned.bam
                    cp $sampledir.$chr.recal.bai $sampledir.$chr.recal.cleaned.bam.bai
		    $samdir/samtools flagstat $sampledir.$chr.recal.cleaned.bam > $sampledir.$chr.recal.cleaned.flagstat
		fi
 		echo `date`

		echo "realigning"
		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -I $sampledir.$chr.recal.bam \
		    -o $sampledir.$chr.recal.bam.list \
		    -T RealignerTargetCreator $realparms $region
	
		if [ ! -s $sampledir.$chr.recal.bam.list ]
		then
		    MSG="$sampledir.$chr.recal.bam.list file not created."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
                    exit 1;
		fi
		echo `date`

		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -I $sampledir.$chr.recal.bam \
		    -T IndelRealigner \
                    -L $chr \
		    -targetIntervals $sampledir.$chr.recal.bam.list $realignparms $realparms \
		    -o $sampledir.$chr.recal.realigned.bam

		if [ ! -s $sampledir.$chr.recal.realigned.bam ]
		then
		    MSG="$sampledir.$chr.recal.realigned.bam file not created."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
		fi
		cp $sampledir.$chr.recal.realigned.bam $sampledir.$chr.real.cleaned.bam
		cp $sampledir.$chr.recal.realigned.bam.bai $sampledir.$chr.real.cleaned.bam.bai
		$samdir/samtools flagstat $sampledir.$chr.real.cleaned.bam > $sampledir.$chr.real.cleaned.flagstat
		echo `date`
        fi

fi