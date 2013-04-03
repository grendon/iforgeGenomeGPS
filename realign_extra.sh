#!/bin/sh
#
#  script to realign and recalibrate the aligned file(s)
#  This module is called from within the realign module
#  Input file(s) is(are) bam format
#  given a specified region, usually one or more chromosomes at once
#  one or several samples can be considered in the input
######################################
if [ $# != 11 ]
then
    MSG="parameter mismatch."
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
    exit 1;
else
	set -x
	echo `date`

        inputdir=$1
        aligndir=$2
        realigndir=$3
        realignlogdir=$4
        infilesuffix=$5
        runfile=$6
        flag=$7
	elog=$8
	olog=$9
	email=${10}
	scriptfile=$0
        qsubfile=${11}
	LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


        #sanity check
        if [ ! -s $runfile ]
        then
	    MSG="$runfile file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        if [ ! -d $inputdir ]
        then
	    MSG="$outputdir directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        if [ ! -d $inputdir/$aligndir ]
        then
	    MSG="$inputdir/$aligndir directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        if [ ! -d $realignlogdir ]
        then
	    MSG="$realignlogdir directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        if [ ! -d $realigndir ]
        then
	    MSG="$realigndir directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi

        numfiles=$( ls -1 $inputdir/$aligndir/${aligndir}*$infilesuffix | wc -l )
        if [ $numfiles -lt 1 ]
        then
	    MSG="bam files to realign not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
	    exit 1;
        fi
        listfiles=$( ls -1 $inputdir/$aligndir/${aligndir}*$infilesuffix )

        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 )

        if [ $type == "whole_genome" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
        fi
        outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
        targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
        realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 )
        chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
        indices=$( echo $chrindex | sed 's/^/chr/' | sed 's/:/ chr/g' )
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
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

        if [ -s $refdir/$dbSNP ]
        then
	    realparms="-known:$refdir/$dbSNP"
            recalparms="--knownSites:$refdir/$dbSNP"
        fi
        if [ -s $refdir/$kgenome ]
        then
	    realparms=$realparms":-known:$refdir/$kgenome"
	    recalparms=$recalparms":--knownSites:$refdir/$kgenome"
        fi

        outputdir=$realigndir
        sample=`basename $outputdir`
        sID=$sample
        sPU=$sample
        sSM=$sample
        RGparms=$( echo "RGID=${sID}:RGLB=${sLB}:RGPU=${sPU}:RGSM=${sSM}:RGPL=${sPL}:RGCN=${sCN}" )
        thr=`expr $threads "-" 1`
        vardir=$outputrootdir/variant/$aligndir
        varlogdir=$outputrootdir/logs/variant

        if [ ! -d $outputdir ]
        then
            mkdir -p $outputdir
        fi
        if [ ! -d $vardir ]
        then
            mkdir -p $vardir
            mkdir -p $varlogdir
	fi
	if [ ! -d $varlogdir ]
        then
            mkdir -p $varlogdir
        else
	    `rm $varlogdir/*`
        fi

        for chr in $indices
        do
            i=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
            if [ -s $refdir/$targetkit -a $type == "exome" ]
            then
		cat $refdir/$targetkit | grep -w $chr > $outputdir/$chr.bed
		if [ `cat $outputdir/$chr.bed | wc -l` -gt 0 ]
                then
                    region[$i]="-L:$outputdir/$chr.bed"
                else
		    region[$i]="-L:$chr"
                fi
            else
                if [ $type == "whole_genome" ]
                then
		    region[$i]="-L:$chr"
                fi
            fi
        done

        cd $outputdir
        for bam in $listfiles
        do
	    for chr in $indices
            do
		inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
                bamfile=`basename $bam`
                if [ ! -s $outputdir/$bamfile ]
                then
                    cp $bam $outputdir/$bamfile
                    cp $bam.bai $outputdir/$bamfile.bai
                fi
                if [ ! -s $bamfile.bai ]
                then
		    $samdir/samtools index $bamfile
                fi

                if [ ! -s $bamfile.$chr.bam ]
                then
                    $samdir/samtools view -b $bamfile $chr > $bamfile.$chr.bam
                    $samdir/samtools index $bamfile.$chr.bam
                fi
                if [ ! -s $bamfile.$chr.bam.bai ]
                then
                    $samdir/samtools index $bamfile.$chr.bam
                fi

                qsub1=$realignlogdir/qsub.sort.$bamfile.$chr
                echo "#PBS -V" > $qsub1
                echo "#PBS -A $pbsprj" >> $qsub1
                echo "#PBS -N sort_$chr" >> $qsub1
		echo "#PBS -l walltime=$pbscpu" >> $qsub1
		echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		echo "#PBS -o $realignlogdir/log.sort.$bamfile.$chr.ou" >> $qsub1
		echo "#PBS -e $realignlogdir/log.sort.$bamfile.$chr.in" >> $qsub1
                echo "#PBS -q $pbsqueue" >> $qsub1
                echo "#PBS -m ae" >> $qsub1
                echo "#PBS -M $email" >> $qsub1
                echo "$scriptdir/sortnode.sh $picardir $samdir $outputdir $bamfile.$chr.bam $bamfile.$chr.sorted.bam $RGparms NCSA $realignlogdir/log.sort.$bamfile.$chr.in $realignlogdir/log.sort.$bamfile.$chr.ou $email $realignlogdir/qsub.sort.$bamfile.$chr" >> $qsub1
                `chmod a+r $qsub1`
                `qsub $qsub1 >> $realignlogdir/SORTED_$chr`

		 chrinfiles[$inx]="-I:$outputdir/$bamfile.$chr.sorted.bam"
		 chrinputfiles[$inx]="INPUT=$outputdir/$bamfile.$chr.sorted.bam"
            done
        done

        `sleep 5s`

	## LOCAL REALIGNMENT AROUND INDELS
        for chr in $indices
        do
            echo "realign-recalibrate for interval:$chr..."

            inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
            jobids=$( cat $realignlogdir/SORTED_$chr | sed "s/\.[a-z]*//g" | tr "\n" ":" )
            qsub2=$realignlogdir/qsub.realrecal.$chr
            echo "#PBS -V" > $qsub2
            echo "#PBS -A $pbsprj" >> $qsub2
            echo "#PBS -N realrc_$chr" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub2
	    echo "#PBS -o $realignlogdir/log.realrecal.$chr.ou" >> $qsub2
	    echo "#PBS -e $realignlogdir/log.realrecal.$chr.in" >> $qsub2
            echo "#PBS -q $pbsqueue" >> $qsub2
            echo "#PBS -m ae" >> $qsub2
            echo "#PBS -M $email" >> $qsub2
            echo "#PBS -W depend=afterok:$jobids" >> $qsub2
            echo "$scriptdir/realrecalold.sh $aligndir $chr ${chrinfiles[$inx]} ${chrinputfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $realignlogdir/log.realrecal.$chr.in $realignlogdir/log.realrecal.$chr.ou $email $realignlogdir/qsub.realrecal.$chr" >> $qsub2
            `chmod a+r $qsub2`
            recaljobid=`qsub $qsub2`
	    echo $recaljobid >> $realignlogdir/REALRECAL_$chr

            qsub3=$varlogdir/qsub.vcallgatk.$chr
            echo "#PBS -V" > $qsub3
            echo "#PBS -A $pbsprj" >> $qsub3
            echo "#PBS -N vcall_$chr" >> $qsub3
	    echo "#PBS -l walltime=$pbscpu" >> $qsub3
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub3
	    echo "#PBS -o $varlogdir/log.vcallgatk.$chr.ou" >> $qsub3
	    echo "#PBS -e $varlogdir/log.vcallgatk.$chr.in" >> $qsub3
            echo "#PBS -q $pbsqueue" >> $qsub3
            echo "#PBS -m ae" >> $qsub3
            echo "#PBS -M $email" >> $qsub3
            echo "#PBS -W depend=afterok:$recaljobid" >> $qsub3
	    echo "$scriptdir/vcallgatk.sh $vardir $outputdir $aligndir $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$chr.in $varlogdir/log.vcallgatk.$chr.ou $email $varlogdir/qsub.vcallgatk.$chr" >> $qsub3
            `chmod a+r $qsub3`
            vcalljobid=`qsub $qsub3`
	    echo $vcalljobid >> $varlogdir/VCALLGATK_$chr
            echo `date`
	done
	`chmod -R 770 $realigndir/`
	`chmod -R 770 $realignlogdir/`
	`chmod -R 770 $vardir/`
	`chmod -R 770 $varlogdir/`
fi