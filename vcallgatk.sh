#!/bin/sh
#	
#  script to perform variant calling with unifiedgenotyper ONLY
#  This module is called from within the realign module
#  Input file(s) is(are) bam format 
#  given a specified region, usually one chromosome at once
#  one or several samples can be considered in the input
######################################
	
if [ $# != 10 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;
else					
	set -x
	echo `date`
        outputdir=$1
        inputdir=$2
	sample=$3
        chr=$4
        region=$5
        runfile=$6
	elog=$7
	olog=$8
	email=$9
	scriptfile=$0
        qsubfile=${10}
	LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        #sanity check
        if [ ! -s $runfile ]
        then
	    MSG="$runfile file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        if [ -z $region ]
        then
	    MSG="$region interval for vcall was not specified"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        suffix="recal.cleaned.bam"
        infile=${sample}*$chr.$suffix
        region=$( echo $region | tr ":" " " )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        ped=$( cat $runfile | grep -w PEDIGREE | cut -d '=' -f2 )
        allsites=$( cat $runfile | grep -w EMIT_ALL_SITES | cut -d '=' -f2 )
        snvcaller=$( cat $runfile | grep -w SNV_CALLER | cut -d '=' -f2 )
        snvmixdir=$( cat $runfile | grep -w SNVMIXDIR | cut -d '=' -f2 )
        snvmixparms=$( cat $runfile | grep -w SNVMIX2PARMS | cut -d '=' -f2 )
        snvmixfilter=$( cat $runfile | grep -w SNVMIX2FILTER | cut -d '=' -f2 )
        uparms=$( cat $runfile | grep -w UNIFIEDGENOTYPERPARMS | cut -d '=' -f2 )
        onlyontarget=$( cat $runfile | grep -w TARGETTED | cut -d '=' -f2 )

        if [ ! -d $inputdir ]
        then
	    MSG="$inputdir directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        if [ ! -s $inputdir/$infile ]
        then
	    MSG="$inputdir/$infile input file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        if [ ! -d $outputdir ]
        then
	    mkdir -p $outputdir
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
        if [ -z $snvcaller ]
        then
	    MSG="$snvcaller snvcaller tool was not specified"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi


       if [ $snvcaller == "GATK" ]
       then
	   echo "snvcaller is GATK"
           if [[ $allsites == "YES" && $type == "exome" ]]
           then
               pedfile=$sample.variant.$chr.raw.all.pbt.vcf
	       outfile=$sample.variant.$chr.raw.all.vcf
	       umode="EMIT_ALL_SITES"
	       utype="BOTH"
           else
               pedfile=$sample.variant.$chr.raw.pbt.vcf
	       outfile=$sample.variant.$chr.raw.vcf
	       umode="EMIT_VARIANTS_ONLY"
	       utype="BOTH"
           fi
       elif [ $snvcaller == "SNVMIX" -o $snvcaller == "SNVMix" ]
       then
	   echo "snvcaller is SNVMIX"
           if [ $allsites == "YES" -a $type == "exome" ]
           then
	       snvfile=$sample.variant.$chr.raw.snv.all.vcf
	       outfile=$sample.variant.$chr.raw.indel.all.vcf
	       combfile=$sample.variant.$chr.raw.multi.vcf
	       combparms="-V $outputdir/$snvfile -V $outputdir/$outfile"
	       umode="EMIT_ALL_SITES"
	       utype="INDEL"
               smode="all"
           else
	       snvfile=$sample.variant.$chr.raw.snv.vcf
	       outfile=$sample.variant.$chr.raw.indel.vcf
	       combfile=$sample.variant.$chr.raw.multi.vcf
	       combparms="$-V outputdir/$snvfile -V $outputdir/$outfile"
	       umode="EMIT_VARIANTS_ONLY"
	       utype="INDEL"
	       smode="target"
           fi
       elif [ $snvcaller == "BEAUTY_EXOME" ]
       then
	   echo "snvcaller is BEAUTY_EXOME"
	   snvfile=$sample.variant.$chr.raw.snvmix.vcf
	   outfile=$sample.variant.$chr.raw.gatk.vcf
	   combfile=$sample.variant.$chr.raw.multi.vcf
	   combparms="-V:GATK $outputdir/$outfile -V:SNVMix $outputdir/$snvfile -priority GATK,SNVMix" 
	   umode="EMIT_VARIANTS_ONLY"
	   utype="BOTH"
	   smode="target"
       else
	   MSG="snvcaller = $snvcaller. This case is not currently available at this site"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	   exit 1;
       fi

       ## now we issue at least one of these calls
        
        echo "calculating variant calling w unifiedgenotyper"
        cd $outputdir

        java -Xmx6g -Xms512m -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    -I $inputdir/$infile \
	    -T UnifiedGenotyper \
            -glm $utype \
            --output_mode $umode \
            -A DepthOfCoverage \
	    -A AlleleBalance \
	    -dcov 2000 \
	    -rf BadCigar \
	    -o $outfile $region $uparms

        if [ ! -s $outfile ]
        then
	    MSG="$outfile file not created. vcall failed."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

	echo `date`

        echo "calculating phasebytransmission if requested"
        if [[ $ped != "NA" && $snvcaller == "GATK" ]]
        then
            echo "calculating phasebytransmission"
            java -Xmx6g -Xms512m -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    -v $outfile \
	    -T PhaseByTransmission \
            --ped $ped \
	    --out $pedfile

            if [ ! -s $pedfile ]
            then
		MSG="$pedfile file not created. vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		exit 1;
            fi

	    echo `date`
        elif [ $snvcaller != "GATK" ]
        then
	    echo "calculating variants by snvmix and then merging results"
            pilefile=$outfile.pileup
            tmpfile=$sample.variant.$chr.tmp.snv
            $samdir/samtools mpileup -f $refdir/$ref $inputdir/$infile > $pilefile 
	    if [ ! -s $pilefile ]
            then
		MSG="pileup file not created. snvmix failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		exit 1;
	    fi
            if [ $smode == "all" ]
            then
		## question: modefile Mu_pi.txt does not exist in that folder
		$snvmixdir/SNVMix2 -i $pilefile -f -m $snvmixdir/Mu_pi.txt -o $tmpfile $snvmixparms
	    else
		$snvmixdir/SNVMix2 -i $pilefile -m $snvmixdir/Mu_pi.txt -o $tmpfile $snvmixparms
	    fi
	    if [ ! -s $tmpfile ]
            then
		MSG="snvmix file not created. snvmix failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		exit 1;
	    fi
	    perl $scriptdir/snvmix_to_vcf.pl -i $tmpfile -o $snvfile
	    if [ ! -s $snvfile ]
            then
		MSG="snv to vcf conversion failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		exit 1;
	    fi

            echo "combining VCF files"

            java -Xmx6g -Xms512m -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    $combparms \
	    -T CombineVariants \
	    -o  $combfile

            if [ ! -s $combfile ]
            then
		MSG="$combfile file not created. vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		exit 1;
            fi
        else
	    echo "Skipping PhaseByTransmission and SNVMix"
	fi
	
fi