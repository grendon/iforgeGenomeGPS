#!/bin/sh
#
#  script to realign and recalibrate the aligned file(s)
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 9 ]
then
    MSG="parameter mismatch."
    echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
    exit 1;
else
	set -x
	echo `date`
	scriptfile=$0
        realigndir=$1
        realignlogdir=$2
        aligndir=$3
        runfile=$4
        flag=$5
	elog=$6
	olog=$7
	email=$8
        qsubfile=$9
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
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
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
        skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
        if [ $type == "GENOME" -o $type == "WHOLE_GENOME" -o $type == "WHOLEGENOME" -o $type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $type == "EXOME" -o $type == "WHOLE_EXOME" -o $type == "WHOLEEXOME" -o $type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for TYPE=$type in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi
        if [ -z $epilogue ]
        then
           MSG="Value for EPILOGUE must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
        fi
        if [ -z $javamodule ]
        then
           MSG="Value for JAVAMODULE must be specified in configuration file"
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
	    MSG="$samdir samtools directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
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
        if [ ! -d $realignlogdir ]
        then
	    MSG="$realignlogdir realignlog directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $realigndir ]
        then
	    MSG="$realigndir realign directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        # get aligned bam files when alignment is done inhouse
        if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ]
        then
            echo "inhouse alignment. gathering aligned bam files"
            cd $outputrootdir
            allfiles=`find ./align -name "*.wdups.sorted.bam"`
	    if [ `expr length ${allfiles}` -lt 1 ]
	    then
		MSG="No bam file(s) found to perform realign-recalibrate at $outputrootdir/align"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
            listfiles=""
	    sep=" "
            for file in $allfiles
            do
		newname=$( echo $file | sed "s/.\/align\///" )
		newname=$outputrootdir/align/$newname
		listfiles=$newname${sep}${listfiles}
            done
	else
            if [ $analysis == "REALIGN_ONLY" -o $analysis == "REALIGNONLY" ]
            then
		echo "3rd param, aligndir,  has the list of aligned bam files"
		listfiles=$( echo $aligndir | tr ":" " " )
            else
		MSG="Invalid value for analysis=$analysis"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
	fi

        outputdir=$realigndir
        vardir=$outputrootdir/variant
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

        #generating regions and intervals files in BED format
        for chr in $indices
        do
            i=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
            if [ -s $refdir/$targetkit -a $type == "EXOME" -o $type == "WHOLE_EXOME" ]
            then
		cat $refdir/$targetkit | grep -w $chr > $outputdir/$chr.bed
		if [ `cat $outputdir/$chr.bed | wc -l` -gt 0 ]
                then
                    region[$i]="-L:$outputdir/$chr.bed"
                else
		    region[$i]="-L:$chr"
                fi
            else
                if [ $type == "WHOLE_GENOME" -o $type == "WGS" -o $type == "GENOME" -o $type == "WHOLEGENOME" ]
                then
		    region[$i]="-L:$chr"
                fi
            fi
        done

        # main loop
        cd $outputdir
        for bam in $listfiles
        do
	    for chr in $indices
            do
		inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
                bamfile=`basename $bam`

		sample=$( echo $bamfile | cut -d "_" -f1 )
		sID=$sample
		sPU=$sample
		sSM=$sample
		RGparms=$( echo "RGID=${sID}:RGLB=${sLB}:RGPU=${sPU}:RGSM=${sSM}:RGPL=${sPL}:RGCN=${sCN}" )

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
                echo "#PBS -N sort_${bamfile}_$chr" >> $qsub1
                echo "#PBS -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=$pbscpu" >> $qsub1
		echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		echo "#PBS -o $realignlogdir/log.sort.$bamfile.$chr.ou" >> $qsub1
		echo "#PBS -e $realignlogdir/log.sort.$bamfile.$chr.in" >> $qsub1
                echo "#PBS -q $pbsqueue" >> $qsub1
                echo "#PBS -m ae" >> $qsub1
                echo "#PBS -M $email" >> $qsub1
                echo "$scriptdir/sortnode.sh $picardir $samdir $javamodule $outputdir $bamfile.$chr.bam $bamfile.$chr.sorted.bam $RGparms NCSA $realignlogdir/log.sort.$bamfile.$chr.in $realignlogdir/log.sort.$bamfile.$chr.ou $email $realignlogdir/qsub.sort.$bamfile.$chr" >> $qsub1
                `chmod a+r $qsub1`
                sortjobid=`qsub $qsub1`
                # new line to avoid hiccup
                `qhold -h u $sortjobid`
                echo $sortjobid >> $outputrootdir/logs/REALSORTEDpbs
                echo $sortjobid >> $outputrootdir/logs/REALSORTED_$chr
		chrinfiles[$inx]=${chrinfiles[$inx]}":-I:$outputdir/$bamfile.$chr.sorted.bam"
		chrinputfiles[$inx]=${chrinputfiles[$inx]}":INPUT=$outputdir/$bamfile.$chr.sorted.bam"
	    done
	done

        igv_files=""
        sep=":"
	for chr in $indices
        do
	    inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
            sortid=$( cat $outputrootdir/logs/REALSORTED_$chr | sed "s/\.[a-z]*//g" | tr "\n" ":" )
            outputfile=$chr.realrecal.output.bam
            igv_files=${igv_files}":INPUT=${outputfile}"
	    echo "realign-recalibrate for interval:$chr..."
	    qsub2=$realignlogdir/qsub.realrecal.$chr
	    echo "#PBS -V" > $qsub2
	    echo "#PBS -A $pbsprj" >> $qsub2
	    echo "#PBS -N realrc_$chr" >> $qsub2
	    echo "#PBS -l epilogue=$epilogue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub2
	    echo "#PBS -o $realignlogdir/log.realrecal.$chr.ou" >> $qsub2
	    echo "#PBS -e $realignlogdir/log.realrecal.$chr.in" >> $qsub2
	    echo "#PBS -q $pbsqueue" >> $qsub2
	    echo "#PBS -m ae" >> $qsub2
	    echo "#PBS -M $email" >> $qsub2
	    echo "#PBS -W depend=afterok:$sortid" >> $qsub2
	    echo "$scriptdir/realrecalold.sh $outputfile $chr ${chrinfiles[$inx]} ${chrinputfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $realignlogdir/log.realrecal.$chr.in $realignlogdir/log.realrecal.$chr.ou $email $realignlogdir/qsub.realrecal.$chr" >> $qsub2
	    `chmod a+r $qsub2`
	    recaljobid=`qsub $qsub2`
            # new line to avoid hiccup
            `qhold -h u $recaljobid`
	    echo $recaljobid >> $outputrootdir/logs/REALRECALpbs


            if [ $skipvcall == "NO" ]
            then
                     echo "variant calling call.."
		     qsub3=$varlogdir/qsub.vcallgatk.$chr
		     echo "#PBS -V" > $qsub3
		     echo "#PBS -A $pbsprj" >> $qsub3
		     echo "#PBS -N vcall_$chr" >> $qsub3
		     echo "#PBS -l epilogue=$epilogue" >> $qsub3
		     echo "#PBS -l walltime=$pbscpu" >> $qsub3
		     echo "#PBS -l nodes=1:ppn=16" >> $qsub3
		     echo "#PBS -o $varlogdir/log.vcallgatk.$chr.ou" >> $qsub3
		     echo "#PBS -e $varlogdir/log.vcallgatk.$chr.in" >> $qsub3
		     echo "#PBS -q $pbsqueue" >> $qsub3
		     echo "#PBS -m ae" >> $qsub3
		     echo "#PBS -M $email" >> $qsub3
		     echo "#PBS -W depend=afterok:$recaljobid" >> $qsub3
		     echo "$scriptdir/vcallgatk.sh $vardir $outputdir $outputfile $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$chr.in $varlogdir/log.vcallgatk.$chr.ou $email $varlogdir/qsub.vcallgatk.$chr" >> $qsub3
		     `chmod a+r $qsub3`
		     vcalljobid=`qsub $qsub3`
		     echo $vcalljobid >> $outputrootdir/logs/VCALLGATKpbs
            else
                     echo "variant calling will not be run"
            fi
	    echo `date`
     done

     # new line to avoid hiccups
     #heldjobs_realsorted=$( cat $outputrootdir/logs/REALSORTEDpbs | sed "s/\.[a-z]*//g" | tr "\n" " " )
     #`qrls -h u $heldjobs_realsorted`

     # new line to avoid hiccups
     #heldjobs_realrecal=$( cat $outputrootdir/logs/REALRECALpbs | sed "s/\.[a-z]*//g" | tr "\n" " " )
     #`qrls -h u $heldjobs_realrecal`

     #`chmod -R 770 $realigndir/`
     #`chmod -R 770 $vardir/`
     
     
     echo "clean up and produce summary table"
     if [ $skipvcall == "NO" ]
     then
	 listjobids=$( cat $outputrootdir/logs/VCALLGATKpbs | sed "s/\.[a-z]*//g" | tr "\n" ":" )
     else
	 listjobids=$( cat $outputrootdir/logs/REALRECALpbs | sed "s/\.[a-z]*//g" | tr "\n" ":" )
     fi
     lastjobid=""

     qsub5=$outputrootdir/logs/qsub.igvbam
     echo "#PBS -V" > $qsub5
     echo "#PBS -A $pbsprj" >> $qsub5
     echo "#PBS -N igvbam" >> $qsub5
     echo "#PBS -l epilogue=$epilogue" >> $qsub5
     echo "#PBS -l walltime=$pbscpu" >> $qsub5
     echo "#PBS -l nodes=1:ppn=16" >> $qsub5
     echo "#PBS -o $outputrootdir/logs/log.igvbam.ou" >> $qsub5
     echo "#PBS -e $outputrootdir/logs/log.igvbam.in" >> $qsub5
     echo "#PBS -q $pbsqueue" >> $qsub5
     echo "#PBS -m ae" >> $qsub5
     echo "#PBS -M $email" >> $qsub5
     echo "#PBS -W depend=afterok:$listjobids" >> $qsub5
     echo "$scriptdir/igvbam.sh $outputrootdir/realign $igv_files $runfile $outputrootdir/logs/log.igvbam.in $outputroot/logs/log.igvbam.ou $email $outputrootdir/logs/qsub.igvbam"  >> $qsub5
     `chmod a+r $qsub5`
     igvjobid=`qsub $qsub5`
     echo $igvjobid >> $outputrootdir/logs/IGVBAMpbs



     # new line to avoid hiccups
     heldjobs_realsorted=$( cat $outputrootdir/logs/REALSORTEDpbs | sed "s/\.[a-z]*//g" | tr "\n" " " )
     `qrls -h u $heldjobs_realsorted`

     # new line to avoid hiccups
     heldjobs_realrecal=$( cat $outputrootdir/logs/REALRECALpbs | sed "s/\.[a-z]*//g" | tr "\n" " " )
     `qrls -h u $heldjobs_realrecal`

     `chmod -R 770 $realigndir/`
     `chmod -R 770 $vardir/`





     qsub6=$outputrootdir/logs/qsub.cleanup.realn
     echo "#PBS -V" > $qsub6
     echo "#PBS -A $pbsprj" >> $qsub6
     echo "#PBS -N cleanup_realn" >> $qsub6
     echo "#PBS -l epilogue=$epilogue" >> $qsub6
     echo "#PBS -l walltime=$pbscpu" >> $qsub6
     echo "#PBS -l nodes=1:ppn=1" >> $qsub6
     echo "#PBS -o $outputrootdir/logs/log.cleanup.realn.ou" >> $qsub6
     echo "#PBS -e $outputrootdir/logs/log.cleanup.realn.in" >> $qsub6
     echo "#PBS -q $pbsqueue" >> $qsub6
     echo "#PBS -m ae" >> $qsub6
     echo "#PBS -M $email" >> $qsub6
     echo "#PBS -W depend=afterok:$igvjobid" >> $qsub6
     echo "$scriptdir/cleanup.sh $outputrootdir $analysis $outputrootdir/logs/log.cleanup.realn.in $outputrootdir/logs/log.cleanup.realn.ou $email $outputrootdir/logs/qsub.cleanup.realn"  >> $qsub6
     `chmod a+r $qsub6`
     cleanjobid=`qsub $qsub6`
     echo $cleanjobid >> $outputrootdir/logs/CLEANUPpbs

     qsub4=$outputrootdir/logs/qsub.summary.realn.allok
     echo "#PBS -V" > $qsub4
     echo "#PBS -A $pbsprj" >> $qsub4
     echo "#PBS -N summaryok" >> $qsub4
     echo "#PBS -l epilogue=$epilogue" >> $qsub4
     echo "#PBS -l walltime=$pbscpu" >> $qsub4
     echo "#PBS -l nodes=1:ppn=1" >> $qsub4
     echo "#PBS -o $outputrootdir/logs/log.summary.realn.ou" >> $qsub4
     echo "#PBS -e $outputrootdir/logs/log.summary.realn.in" >> $qsub4
     echo "#PBS -q $pbsqueue" >> $qsub4
     echo "#PBS -m ae" >> $qsub4
     echo "#PBS -M $email" >> $qsub4
     echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub4
     echo "$scriptdir/summary.sh $outputrootdir $email exitok"  >> $qsub4
     `chmod a+r $qsub4`
     lastjobid=`qsub $qsub4`
     echo $lastjobid >> $outputrootdir/logs/SUMMARYpbs

     if [ `expr length ${lastjobid}` -lt 1 ]
     then
         echo "at least one job aborted"
	 qsub5=$outputrootdir/logs/qsub.summary.realn.afterany
	 echo "#PBS -V" > $qsub5
	 echo "#PBS -A $pbsprj" >> $qsub5
	 echo "#PBS -N summary_afterany" >> $qsub5
	 echo "#PBS -l epilogue=$epilogue" >> $qsub5
	 echo "#PBS -l walltime=$pbscpu" >> $qsub5
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub5
	 echo "#PBS -o $outputrootdir/logs/log.summary.realn.afterany.ou" >> $qsub5
	 echo "#PBS -e $outputrootdir/logs/log.summary.realn.afterany.in" >> $qsub5
	 echo "#PBS -q $pbsqueue" >> $qsub5
	 echo "#PBS -m ae" >> $qsub5
	 echo "#PBS -M $email" >> $qsub5
	 echo "#PBS -W depend=afterany:$listjobids" >> $qsub5
	 echo "$scriptdir/summary.sh $outputrootdir $email exitnotok"  >> $qsub5
	 `chmod a+r $qsub5`
	 badjobid=`qsub $qsub5`
	 echo $badjobid >> $outputrootdir/logs/SUMMARYpbs
     fi
     `chmod -R 770 $outputroordir/logs`
     echo `date`
fi
