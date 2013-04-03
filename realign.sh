#!/bin/sh
#
# realign.sh
# Second module in the GGPS analysis pipeline
# This module realigns/recalibrates files coming from sequences of exomes or whole-genomes.
# One or several bam files can be processed per sample.
# If several bam files exist per sample, they will be merged before realign-recal is performed
# The input format is bam; provenance of the input file could be 
# a) result of running the alignment module inhouse at ncsa
# b) result of running the alignment somewhere else other than ncsa; may need to pre-process the bam file
# The output files are two per sample: inputfile.bam.recalibrated and inputfile.bam.realigned

if [ $# != 5 ]
then
	MSG="parameter mismatch.\nUsage: <run info file><errorlog><outputlog> "
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        runfile=$1
        elog=$2
        olog=$3
        email=$4
        qsubfile=$5
        LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        if [ !  -s $runfile ]
        then
           MSG="$runfile file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s 'GGPS error notification' "$email""
           exit 1;
        fi


        ## parsing run info file

	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 )
        igvdir=$( cat $runfile | grep -w IGVDIR | cut -d '=' -f2 )
        extradir=$outputdir/extractreads
        if [ $type == "whole_genome" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
        fi
        realrecalflag=$( cat $runfile | grep -w REALIGNORDER | cut -d '=' -f2 )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        inputdir=$( cat $runfile | grep -w SAMPLEDIR | cut -d '=' -f2 )
        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        region=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
        resortflag=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 )
        indices=$( echo $region | sed 's/^/chr/' | sed 's/:/ chr/g' )

        output_logs=$outputdir/logs

        # checking that files and folders exist
        if [ ! -d $scriptdir ]
        then
           MSG="$scriptdir directory not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
           exit 1;
        fi
        if [ ! -d $refdir ]
        then
           MSG="$refdir directory not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
           exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
           MSG="$samplefileinfo file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
           exit 1;
        fi
        if [ ! -d $outputdir ]
        then
	    mkdir -p $outputdir
        fi
        if [ ! -d $output_logs ]
        then
	    mkdir  $output_logs
        fi  
        
        ## preprocessing step; only when alignment was not performed inhouse
        listdirs=":"
        if [ $resortflag == "YES" ]
        then
            echo "alignment was NOT done inhouse. Checking input files"
            numsamples=0
            for name in $samples
            do
		countnames=$( cat $samplefileinfo | grep $name -c )
		if [ $countnames -lt 1  ]
                then
		    MSG="mismatch between samples and samplefileinfo"
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
                fi
                let numsamples+=1
            done
            if [ $numsamples -gt 1 -a $multisample == "YES" ]
	    then
		echo "multiple samples were aligned"
	    else
		if [ $numsamples -eq 1 -a $multisample == "NO" ]
		then
		    echo "single sample was aligned."
		else
		    MSG="mismatch between samples found and multisample setting."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
		    exit 1;
		fi
            fi

            counter=0
            while read sampledetail
            do
                echo "processing next line in file ..."
		len=`expr ${#sampledetail}`
		if [ $len -eq 0 ]
		then
		    echo "line is empty; nothing to do. $sampledetail"
                else
		    echo "preprocessing for realignment $sampledetail"
                    echo "expected format--> BAM:sample_name=bamfilename.bam"

		    sampledirname=$( echo $sampledetail | grep ^BAM | cut -d ':' -f2 | cut -d '=' -f1 )
		    samplename=$( echo $sampledetail | grep ^BAM | cut -d ':' -f2 | cut -d '=' -f2 )

                    len1=`expr length ${sampledirname}`
                    if [ $len1 -lt 1 ]
                    then
			MSG="parsing info in samplefileinfo failed. realignment stopped."
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
			exit 1;
                    fi
                    cd $inputdir
                    bamfile=$( echo $samplename )
                    lenbam=`expr length $bamfile`
                    if [ $lenbam -lt 1 ]
                    then
			MSG="parsing info in samplefileinfo failed. realignment stopped."
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
			exit 1;
                    fi
		
		    if [ ! -s $inputdir/$bamfile ]
		    then
			MSG="$inputdir/$bamfile input bam file not found"
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
			exit;
		    fi

                    prefix=`basename $samplename .wrg.sorted.bam`
		    outputalign=$outputdir/align/$sampledirname
		    outputlogs=$output_logs/align/$sampledirname

                    if [ ! -d $outputalign ]
                    then
			mkdir -p $outputalign
                        let counter+=1
                        tmpbamfile=tmp.${prefix}_$counter.sorted.bam
			sortedplain=${prefix}_$counter.wrg.sorted.bam
			sorted=${prefix}_$counter.wdups.sorted.bam
			sortednodups=${prefix}_$counter.nodups.sorted.bam
			cd $outputalign
			cp $inputdir/$bamfile $tmpbamfile
			if [ ! -d $outputlogs ]
			then
			    mkdir -p $outputlogs
                        else
                            `rm -r $outputlogs/*`
			fi
                    else
                        let counter+=1
                        tmpbamfile=tmp.${prefix}_$counter.sorted.bam
			sortedplain=${prefix}_$counter.wrg.sorted.bam
			sorted=${prefix}_$counter.wdups.sorted.bam
			sortednodups=${prefix}_$counter.nodups.sorted.bam
			cd $outputalign
			cp $inputdir/$bamfile $tmpbamfile
		   fi

		   qsub1=$outputlogs/qsub.sortmayo.$counter
		   echo "#PBS -V" > $qsub1
		   echo "#PBS -A $pbsprj" >> $qsub1
		   echo "#PBS -N sortm_$counter" >> $qsub1
		   echo "#PBS -l walltime=$pbscpu" >> $qsub1
		   echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		   echo "#PBS -o $outputlogs/log.sortmayo.$counter.ou" >> $qsub1
		   echo "#PBS -e $outputlogs/log.sortmayo.$counter.in" >> $qsub1
		   echo "#PBS -q $pbsqueue" >> $qsub1
		   echo "#PBS -m ae" >> $qsub1
		   echo "#PBS -M $email" >> $qsub1
		   echo "$scriptdir/sortbammayo.sh $outputalign $tmpbamfile $sortedplain $sorted $sortednodups $runfile $outputlogs/log.sortmayo.$counter.in $outputlogs/log.sortmayo.$counter.ou $email $outputlogs/qsub.sortmayo.$counter" >> $qsub1
		   `chmod a+r $qsub1`
                   sortid=`qsub $qsub1`
          	   echo $sortid >> $outputlogs/SORTEDmayo
            
		   echo "extracting reads"
		   qsub2=$outputlogs/qsub.extractreadsbam.$counter
		   echo "#PBS -V" > $qsub2
		   echo "#PBS -A $pbsprj" >> $qsub2
		   echo "#PBS -N extrbam$counter" >> $qsub2
		   echo "#PBS -l walltime=$pbscpu" >> $qsub2
		   echo "#PBS -l nodes=1:ppn=16" >> $qsub2
		   echo "#PBS -o $outputlogs/log.extractreadsbam.$counter.ou" >> $qsub2
		   echo "#PBS -e $outputlogs/log.extractreadsbam.$counter.in" >> $qsub2
		   echo "#PBS -q $pbsqueue" >> $qsub2
		   echo "#PBS -m ae" >> $qsub2
		   echo "#PBS -M $email" >> $qsub2
                   echo "#PBS -W depend=afterok:$sortid" >> $qsub2
                   echo "$scriptdir/extract_reads_bam.sh $outputalign $sorted $runfile $igvdir $extradir" >> $qsub2
		   `chmod a+r $qsub2`               
                   `qsub $qsub2 >> $output_logs/EXTRACTREADSmayo`
		fi
	    done < $samplefileinfo
            mv $outputlogs/SORTEDmayo $output_logs/ALN_MAYO_jobids
            alndirs=$( ls -1 $outputdir/align )
            JOBSmayo=$( cat $output_logs/ALN_MAYO_jobids | sed "s/\.[a-z]*//g" | tr "\n" ":" | sed "s/::/:/g" )
        else
            echo "alignment was done inhouse. no need to resort"
            echo "We need to wait until the alignment jobs enter the queue"

            while [ ! -s $output_logs/ALN_NCSA_jobids ]
            do
		`sleep 60s`
            done
            alndirs=$( ls -1 $outputdir/align )
            JOBSncsa=$( cat $output_logs/ALN_NCSA_jobids | sed "s/\.[a-z]*//g" | tr "\n" ":" | sed "s/::/:/g" |  sed "s/^://" )
        fi

        ## preprocessing is done. Now we can realign and recalibrate bam files
        echo "preprocessing of bam files is done. ready to realign"

        realigndir=$outputdir/realign
        realignlogdir=$outputdir/logs/realign
	if [ ! -d $realigndir ]
        then
            mkdir $realigndir
	else
	    echo "$realigndir already exists. resetting it"
            `rm -r $realigndir/*`
        fi
	if [ ! -d $realignlogdir ]
        then
            mkdir -p $realignlogdir
	else
	    echo "$realignlogdir already exists. resetting it"
            `rm -r $realignlogdir/*`
        fi

        if [ $realrecalflag != "1" -a $realrecalflag != "0" ]
        then
	    echo "realign-recalibration order flag is not set properly. Default value [1] will be assiged to it"
            realrecalflag="1"
        fi
        targetsuffix=.wdups.sorted.bam
        recalcount=0
        for alndir in $alndirs 
        do
            recaldir=$realigndir/$alndir
            if [ ! -d $recaldir ]
            then
                mkdir -p $recaldir
	    fi

            let recalcount+=1
	    qsub3=$realignlogdir/qsub.recalibration.$recalcount
	    echo "#PBS -V" > $qsub3
	    echo "#PBS -A $pbsprj" >> $qsub3
	    echo "#PBS -N recal_$recalcount" >> $qsub3
	    echo "#PBS -l walltime=$pbscpu" >> $qsub3
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub3
	    echo "#PBS -o $realignlogdir/log.recalibration.$recalcount.ou" >> $qsub3
	    echo "#PBS -e $realignlogdir/log.recalibration.$recalcount.in" >> $qsub3
	    echo "#PBS -q $pbsqueue" >> $qsub3
	    echo "#PBS -m ae" >> $qsub3
	    echo "#PBS -M $email" >> $qsub3
            if [ $resortflag == "YES" ]
            then
                echo "#PBS -W depend=afterok:$JOBSmayo" >> $qsub3
            else
                echo "#PBS -W depend=afterok:$JOBSncsa" >> $qsub3
            fi
            echo "$scriptdir/realign_extra.sh $outputdir/align $alndir $recaldir $realignlogdir $targetsuffix $runfile $realrecalflag $realignlogdir/log.recalibration.$recalcount.in $realignlogdir/log.recalibration.$recalcount.ou $email $realignlogdir/qsub.recalibration.$recalcount" >> $qsub3
	    `chmod a+r $qsub3`               
            `qsub $qsub3 >> $output_logs/RECALLpbs`
        done
	#`sleep 5s`
        echo "done realig/recalibrating  all bam files."
        echo `date`
	`chmod -R 770 $outputdir/`
	`chmod -R 770 $output_logs/`
fi