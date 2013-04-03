#!/bin/sh
#	
#  script to realign and recalibrate the aligned file(s) 
#  This module is called from within the realign module
#  Input file(s) is(are) bam format 
#  given a specified region, usually one or more chromosomes at once
#  one or several samples can be considered in the input
######################################
	
if [ $# != 7 ];
then
	echo "REALIGN_RECAL usage: parameter mismatch"
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


        #sanity check
        if [ ! -s $runfile ]
        then
	    echo "$runfile file not found"
	    exit 1;
        fi
        if [ ! -d $inputdir ]
        then
	    echo "$outputdir directory not found"
	    exit 1;
        fi
        if [ ! -d $inputdir/$aligndir ]
        then
	    echo "$inputdir/$aligndir directory not found"
	    exit 1;
        fi
        if [ ! -d $realignlogdir ]
        then
	    echo "$realignlogdir directory not found"
	    exit 1;
        fi
        if [ ! -d $realigndir ]
        then
	    echo "$realigndir directory not found"
	    exit 1;
        fi

        numfiles=$( ls -1 $inputdir/$aligndir/*$infilesuffix | wc -l )
        if [ $numfiles -lt 1 ]
        then
	    echo "bam files to realign not found"
	    exit 1;
        fi
        listfiles=$( ls -1 $inputdir/$aligndir/*$infilesuffix )

        outputdir=$realigndir/$aligndir
        if [ ! -d $outputdir ]
        then
            mkdir -p $outputdir
        fi

        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 )
        if [ $type == "whole_genome" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
        fi
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
        sID=$( cat $runfile | grep -w SAMPLEID | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )       
        sPU=$( cat $runfile | grep -w SAMPLEPU | cut -d '=' -f2 )
        sSM=$( cat $runfile | grep -w SAMPLESM | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )        
        RGparms=$( echo "RGID=${sID}_RGLB=${sLB}_RGPL=${sPL}_RGPU=${sPU}_RGSM=${sSM}_RGCN=${sCN}" )
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

        `sleep 15s`
        for chr in $indices
        do
            i=$( echo $chr | sed 's/chr//' )
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

        ## extract reads for a specified interval
        `sleep 5s`
        cd $outputdir
        for bam in $listfiles
        do
	    for chr in $indices
            do
                ## todo:
                ## we need to put this loop into a pbs job
		inx=$( echo $chr | sed 's/chr//' )
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

                # extracting chromosome

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
                echo "$scriptdir/sortnode.sh $picardir $samdir $outputdir $bamfile.$chr.bam $bamfile.$chr.sorted.bam $RGparms NCSA" >> $qsub1
                `chmod a+r $qsub1`
                `qsub $qsub1 >> $realignlogdir/SORTED_$chr`

                if [ -z ${chrinfiles[$inx]} ]
                then
                    echo "first entry on the list"
		    chrinfiles[$inx]="-I:$outputdir/$bamfile.$chr.sorted.bam"
		    chrinputfiles[$inx]="INPUT=$outputdir/$bamfile.$chr.sorted.bam"
                else
		    chrinfiles[$inx]=${chrinfiles[$inx]}":-I:$outputdir/$bamfile.$chr.sorted.bam"
		    chrinputfiles[$inx]=${chrinputfiles[$inx]}":INPUT=$outputdir/$bamfile.$chr.sorted.bam"
               fi
               `sleep 5s`
            done
        done

        `sleep 5s`
	## LOCAL REALIGNMENT AROUND INDELS
        for chr in $indices
        do
            echo "realign-recalibrate for interval:$chr..."

            inx=$( echo $chr | sed 's/chr//' )
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
            echo "$scriptdir/realrecal.sh $outputdir $chr ${chrinfiles[$inx]} ${chrinputfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag" >> $qsub2
            `chmod a+r $qsub2`
            `qsub $qsub2 >> $realignlogdir/REALRECAL_$chr`
	done
fi