#!/bin/sh
#	
#  script to realign and recalibrate the aligned file(s) 
#  This module is called from within the realign module
#  Input file(s) is(are) bam format 
#  given a specified region, usually one or more chromosomes at once
#  one or several samples can be considered in the input
######################################
	
if [ $# != 9 ];
then
	echo "REALRECAL usage: parameter mismatch"
else					
	set -x

	echo `date`
	

        inputdir=$1
        chr=$2
        chrinfiles=$3
        chrinputfiles=$4
        region=$5
        realparms=$6
        recalparms=$7
        runfile=$8
        flag=$9


        #sanity check
        if [ ! -s $runfile ]
        then
	    echo "$runfile file not found"
	    exit 1;
        fi
        if [ ! -d $inputdir ]
        then
	    echo "$inputdir directory not found"
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

        sID=$( cat $runfile | grep -w SAMPLEID | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )       
        sPU=$( cat $runfile | grep -w SAMPLEPU | cut -d '=' -f2 )
        sSM=$( cat $runfile | grep -w SAMPLESM | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )        
        RGparms=$( echo "RGID=${sID} RGLB=${sLB} RGPL=${sPL} RGPU=${sPU} RGSM=${sSM} RGCN=${sCN}" )
        thr=`expr $threads "-" 1`
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
        if [ ! -d $gatk ]
        then
	    echo "$gatk directory not found"
	    exit 1;
        fi

        if [ ! -d $refdir ]
        then
	    echo "$refdir reference genome directory not found"
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    echo "$ref reference genome not found"
	    exit 1;
        fi

        # cleaning up the lists
        chrinfiles=$( echo $chrinfiles | tr ":" " " )
        chrinputfiles=$( echo $chrinputfiles | tr ":" " " )
        region=$( echo $region | tr ":" " " )
        realparms=$( echo $realparms | tr ":" " " )
        recalparms=$( echo $recalparms | tr ":" " " )
        inx=$( echo $chr | sed 's/chr//' )
        cd $inputdir

        if [ $flag == 1 ]
        then
		echo "realign then recalibrate"
		echo "realigning..."
		
            	java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
                    -et NO_ET \
		    $chrinfiles \
		    -T RealignerTargetCreator \
		    -o $chr.bam.list $realparms ${region[$inx]}
	
		if [ ! -s $chr.bam.list ]
		then
		    echo "$chr.bam.list file not created."
                    exit 1;
		fi
		echo `date`

		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
                    -et NO_ET \
		    $chrinfiles \
		    -T IndelRealigner \
                    -L $chr \
		    -o $chr.realigned.bam \
		    -targetIntervals $chr.bam.list $realignparams $realparms

		if [ ! -s $chr.realigned.bam ]
		then
		    echo "$indices.realigned.bam file not created."
                    exit 1;
                else
		    mv $chr.realigned.bai $chr.realigned.bam.bai
		    cp $chr.realigned.bam $chr.real.cleaned.bam
		    cp $chr.realigned.bam.bai $chr.real.cleaned.bam.bai
                    $samdir/samtools flagstat $chr.real.cleaned.bam > $chr.real.cleaned.flagstat
		fi
		echo `date`

		echo "recalibrating..."
	
		java -Xmx6g -Xms512m  -jar $picardir/FixMateInformation.jar \
		    INPUT=$chr.realigned.bam \
		    OUTPUT=$chr.realigned.fixed.bam \
		    SO=coordinate \
		    TMP_DIR=$inputdir \
		    VALIDATION_STRINGENCY=SILENT \
		    CREATE_INDEX=true

		if [ ! -s $chr.realigned.fixed.bam ]
		then
		    echo "$chr.realigned.fixed.bam file not created"
		    exit 1;
		fi
		echo `date`

		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $recalparms \
		    -I $chr.realigned.fixed.bam \
		    ${region[$inx]}  \
		    -T CountCovariates \
		    -nt $thr \
		    -cov ReadGroupCovariate \
		    -cov QualityScoreCovariate \
		    -cov CycleCovariate \
		    -cov DinucCovariate \
		    -recalFile $chr.recal_data.csv 
	
		if [ ! -s $chr.recal_data.csv ]
		then
		    echo "$chr.recal_data.csv file not created"
		exit 1;
		fi
		echo `date`

		java -Xmx6g -Xms512m  -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -L $chr \
		    -I $chr.realigned.fixed.bam \
		    -T TableRecalibration \
		    --out $chr.realigned.fixed.recal.bam \
		    -recalFile $chr.recal_data.csv 

		if [ ! -s $chr.realigned.fixed.recal.bam ]
		then
		    echo "$chr.realigned.fixed.recal.bam file not created"
		    exit 1;
		fi
		cp $chr.realigned.fixed.recal.bam $chr.recal.cleaned.bam
		cp $chr.realigned.fixed.recal.bai $chr.recal.cleaned.bam.bai
		$samdir/samtools flagstat $chr.recal.cleaned.bam > $chr.recal.cleaned.flagstat
 		echo `date`
        else
		echo "recalibrate then realign"
		echo "recalibrating"


		java -Xmx6g -Xms512m  -jar $picardir/FixMateInformation.jar \
		    $chrinputfiles \
		    OUTPUT=$chr.fixed.bam \
		    SO=coordinate \
		    TMP_DIR=$inputdir \
		    VALIDATION_STRINGENCY=SILENT \
		    CREATE_INDEX=true

		if [ ! -s $chr.fixed.bam ]
		then
		    echo "$chr.fixed.bam file not created"
		    exit 1;
		fi
		echo `date`

		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $recalparms \
		    -I $chr.fixed.bam \
		    ${region[$inx]} \
		    -T CountCovariates \
		    -nt $thr \
		    -cov ReadGroupCovariate \
		    -cov QualityScoreCovariate \
		    -cov CycleCovariate \
		    -cov DinucCovariate \
		    -recalFile $chr.fixed.recal_data.csv 
	
		if [ ! -s $chr.fixed.recal_data.csv ]
		then
		    echo "$chr.fixed.recal_data.csv file not created"
		    exit 1;
		fi
	
		echo `date`

		java -Xmx6g -Xms512m  -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -L $chr \
		    -I $chr.fixed.bam \
		    -T TableRecalibration \
		    --out $chr.fixed.recal.bam \
		    -recalFile $chr.fixed.recal_data.csv 

		if [ ! -s $chr.fixed.recal.bam ]
		then
		    echo "$chr.fixed.recal.bam file not created"
		    exit 1;
                else 
                    cp $chr.fixed.recal.bam $chr.recal.cleaned.bam
                    cp $chr.fixed.recal.bai $chr.recal.cleaned.bam.bai
		    $samdir/samtools flagstat $chr.recal.cleaned.bam > $chr.recal.cleaned.flagstat
		fi
 		echo `date`

		echo "realigning"
		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -I $chr.fixed.recal.bam \
		    -o $chr.fixed.recal.bam.list \
		    -T RealignerTargetCreator 
	
		if [ ! -s $chr.fixed.recal.bam.list ]
		then
		    echo "$chr.fixed.recal.bam.list file not created."
		fi
		echo `date`

		java -Xmx6g -Xms512m -Djava.io.tmpdir=$inputdir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -I $chr.fixed.recal.bam \
		    -T IndelRealigner \
		    -targetIntervals $chr.fixed.recal.bam.list $realignparms $realparms \
		    -o $chr.fixed.recal.realigned.bam

		if [ ! -s $chr.fixed.recal.realigned.bam ]
		then
		    echo "$chr.fixed.recal.realigned.bam file not created."
		fi
		cp $chr.fixed.recal.realigned.bam $chr.real.cleaned.bam
		cp $chr.fixed.recal.realigned.bam.bai $chr.real.cleaned.bam.bai
		$samdir/samtools flagstat $chr.real.cleaned.bam > $chr.real.cleaned.flagstat
		echo `date`
        fi

fi