#!/bin/sh
#
# realign.sh
# Second module in the GGPS analysis pipeline
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 5 ]
then
	MSG="parameter mismatch. "
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
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
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        igvdir=$( cat $runfile | grep -w IGVDIR | cut -d '=' -f2 )
        extradir=$outputdir/extractreads
        realrecalflag=$( cat $runfile | grep -w REALIGNORDER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        region=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
        resortflag=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2  )
        indices=$( echo $region | sed 's/^/chr/' | sed 's/:/ chr/g' )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        output_logs=$outputdir/logs


        if [ $type == "GENOME" -o $type == "WHOLE_GENOME" -o $type == "WHOLEGENOME" -o $type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $type == "EXOME" -o $type == "WHOLE_EXOME" -o $type == "WHOLEEXOME" -o $type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for TYPE=$type in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi
        if [ -z $epilogue ]
        then
           MSG="Invalid value for EPILOGUE=$epilogue in configuration file"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
        fi
        if [ ! -d $scriptdir ]
        then
           MSG="$scriptdir script directory not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $refdir ]
        then
           MSG="$refdir directory of reference genome was not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $picardir ]
        then
           MSG="$picardir picard directory was not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
           MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        pipeid=$( cat $output_logs/MAINpbs )
        
        if [ ! -d $outputdir ]
        then
	    mkdir -p $outputdir
        fi
        if [ ! -d $output_logs ]
        then
	    mkdir  $output_logs
        fi  

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
        

        listfiles="";
        sep=":";
        JOBSmayo=""
        JOBSncsa=""

        # finding all aligned BAMs to be realigned-recalibrated

        if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ]
        then
	    echo "alignment was done inhouse. no need to resort"
	    echo "We need to wait until the alignment jobs enter the queue"

	    while [ ! -s $output_logs/ALN_NCSA_jobids ]
	    do
		`sleep 60s`
	    done
            listfiles=$outputdir/align
	    JOBSncsa=$( cat $output_logs/ALN_NCSA_jobids | sed "s/\.[a-z]*//g" | tr "\n" ":" | sed "s/::/:/g" )
        else
            echo "we need to check entries in samplefileinfo before launching other realignment analyses"
            while read sampledetail
            do
                echo "processing next line in file ..."
		if [ `expr ${#sampledetail}` -lt 7 ]
		then
		    echo "skip empty line"
                else
		    echo "preprocessing for realignment $sampledetail"
		    bamfile=$( echo $sampledetail | cut -d '=' -f2 )
                    sampletag=$( echo $sampledetail | cut -d '=' -f1 | cut -d ':' -f2 )
                    if [ `expr ${#bamfile}` -lt 1 ]
                    then
			MSG="parsing SAMPLEFILENAMES file failed. realignment failed."
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
                    if [ `expr ${#sampletag}` -lt 1 ]
                    then
			MSG="parsing SAMPLEFILENAMES file failed. realignment failed."
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
		
		    if [ ! -s $bamfile ]
		    then
			MSG="parsing SAMPLEFILENAMES file failed. realignment failed"
	                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit;
		    fi
                fi
	    done < $samplefileinfo
        fi   

        if [ $resortflag == "YES" -a $analysis == "REALIGN_ONLY" ]
        then
            echo "alignment was NOT done inhouse. Need to resort bam files. Checking input files"
            while read sampledetail
            do
               if [ `expr ${#sampledetail}` -lt 7 ]
               then
                    echo "skip empty line"
               else
		    samplename=$( echo $sampledetail | cut -d '=' -f2 )
                    prefix=`basename $samplename .wrg.sorted.bam`
		    outputalign=$outputdir/align/$prefix
		    outputlogs=$output_logs/align/$prefix
                    tmpbamfile=$samplename
		    sortedplain=${prefix}.wrg.sorted.bam
		    sorted=${prefix}.wdups.sorted.bam
		    sortednodups=${prefix}.nodups.sorted.bam

                    if [ ! -d $outputalign ]
                    then
			mkdir -p $outputalign
			if [ ! -d $outputlogs ]
			then
			    mkdir -p $outputlogs
                        else
                            `rm -r $outputlogs/*`
			fi
		   fi

		   qsub1=$outputlogs/qsub.sortbammayo.$prefix
		   echo "#PBS -V" > $qsub1
		   echo "#PBS -A $pbsprj" >> $qsub1
		   echo "#PBS -N ${pipeid}_sortbamayo_${prefix}" >> $qsub1
                   echo "#PBS -l epilogue=$epilogue" >> $qsub1
		   echo "#PBS -l walltime=$pbscpu" >> $qsub1
		   echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		   echo "#PBS -o $outputlogs/log.sortbammayo.${prefix}.ou" >> $qsub1
		   echo "#PBS -e $outputlogs/log.sortbammayo.${prefix}.in" >> $qsub1
		   echo "#PBS -q $pbsqueue" >> $qsub1
		   echo "#PBS -m ae" >> $qsub1
		   echo "#PBS -M $email" >> $qsub1
		   echo "$scriptdir/sortbammayo.sh $outputalign $tmpbamfile $sortedplain $sorted $sortednodups $runfile $outputlogs/log.sortbammayo.${prefix}.in $outputlogs/log.sortbammayo.${prefix}.ou $email $outputlogs/qsub.sortbammayo.${prefix}" >> $qsub1
		   `chmod a+r $qsub1`
                   sortid=`qsub $qsub1`
                   `qhold -h u $sortid`
               	   echo $sortid >> $output_logs/REALSORTEDMAYOpbs
		   echo "extracting reads"
		   qsub2=$outputlogs/qsub.extractreadsbam.$prefix
		   echo "#PBS -V" > $qsub2
		   echo "#PBS -A $pbsprj" >> $qsub2
		   echo "#PBS -N ${pipeid}_extrbam${prefix}" >> $qsub2
                   echo "#PBS -l epilogue=$epilogue" >> $qsub2
		   echo "#PBS -l walltime=$pbscpu" >> $qsub2
		   echo "#PBS -l nodes=1:ppn=16" >> $qsub2
		   echo "#PBS -o $outputlogs/log.extractreadsbam.${prefix}.ou" >> $qsub2
		   echo "#PBS -e $outputlogs/log.extractreadsbam.${prefix}.in" >> $qsub2
		   echo "#PBS -q $pbsqueue" >> $qsub2
		   echo "#PBS -m ae" >> $qsub2
		   echo "#PBS -M $email" >> $qsub2
                   echo "#PBS -W depend=afterok:$sortid" >> $qsub2
                   echo "$scriptdir/extract_reads_bam.sh $outputalign $sorted $runfile $outputlogs/log.extractreadsbam.${prefix}.in $outputlogs/log.extractreadsbam.${prefix}.ou $outputlogs/qsub.extractreadsbam.$prefix $igvdir $extradir" >> $qsub2
		   `chmod a+r $qsub2`               
                   extraid=`qsub $qsub2`
                   `qhold -h u $extraid`
                   echo $extraid >> $output_logs/EXTRACTREADSpbs
                   listfiles=$outputalign/$sorted${sep}${listfiles}
		fi
	    done < $samplefileinfo
            cp $output_logs/REALSORTEDMAYOpbs $output_logs/ALN_MAYO_jobids
            JOBSmayo=$( cat $output_logs/ALN_MAYO_jobids | sed "s/\.[a-z]*//g" | tr "\n" ":" | sed "s/::/:/g" )
        fi

        if [ $resortflag != "YES" -a $analysis == "REALIGN_ONLY" ]
        then
	    echo "alignment was NOT done inhouse. BAM files will not be resorted"
            if [ $revertsam != "1" ]
            then
                echo "input is aligned bam that is suitable for realignment and recalibration... no need for preprocessing"
                while read sampledetail
                do
		    bam=$( echo $sampledetail | cut -d "=" -f2 )
		    listfiles=${bam}${sep}${listfiles}
		done < $samplefileinfo
            else
                echo "revertsam needs to be run on input bams..."
                while read sampledetail
                do
                    if [ `expr length ${sampledetail}` -lt 7 ]
                    then
                        echo "skip empty line"
                    else
			bamfile=$( echo $sampledetail | cut -d "=" -f2 )
			if [ `expr ${#bamfile}` -lt 1 ]
                        then
                           echo "skip empty line"
                        else
                            prefix=`basename $bamfile .bam`
                            dirname=$realigndir
                            revbamfile=$dirname/${prefix}.revsam

			    listfiles=${revbamfile}${sep}${listfiles}
			    qsub2=$output_logs/qsub.revertinputbam.$prefix
			    echo "#PBS -V" > $qsub2
			    echo "#PBS -A $pbsprj" >> $qsub2
			    echo "#PBS -N ${pipeid}_revertinputbam_${prefix}" >> $qsub2
			    echo "#PBS -l epilogue=$epilogue" >> $qsub2
			    echo "#PBS -l walltime=$pbscpu" >> $qsub2
			    echo "#PBS -l nodes=1:ppn=16" >> $qsub2
			    echo "#PBS -o $output_logs/log.revertinputbam.${prefix}.ou" >> $qsub2
			    echo "#PBS -e $output_logs/log.revertinputbam.${prefix}.in" >> $qsub2
			    echo "#PBS -q $pbsqueue" >> $qsub2
			    echo "#PBS -m ae" >> $qsub2
			    echo "#PBS -M $email" >> $qsub2
			    echo "$scriptdir/revertinputbam.sh $bamfile $revbamfile $picardir $samdir $output_logs/log.revertinputbam.${prefix}.in $output_logs/log.revertinputbam.${prefix}.ou $email $output_logs/qsub.revertinputbam.$prefix" >> $qsub2
			    `chmod a+r $qsub2`
			    revbamid=`qsub $qsub2`
			    `qhold -h u $revbamid`
			    echo $revbamid >> $output_logs/REVERTINPUTBAMpbs
                        fi
                      fi
                    done < $samplefileinfo
		    JOBSmayo=$( cat $output_logs/REVERTINPUTBAMpbs | sed "s/\.[a-z]*//g" | tr "\n" ":" | sed "s/::/:/g" )
                fi
        fi


        # grab job ids for align and for preprocessing done in this module

        alignids=$( cat $output_logs/ALIGNEDpbs | sed "s/\.[a-z]*//" | tr "\n" " " )
        mergeids=$( cat $output_logs/MERGEDpbs | sed "s/\.[a-z]*//" | tr "\n" " " )
        sortedmayoids=$( cat $output_logs/REALSORTEDMAYOpbs | sed "s/\.[a-z]*//" | tr "\n" " " )
        revbamids=$( cat $output_logs/REVERTINPUTBAMpbs | sed "s/\.[a-z]*//" | tr "\n" " " )
        extraids=$( cat $output_logs/EXTRACTREADSpbs | sed "s/\.[a-z]*//" | tr "\n" " " )

	qsub3=$realignlogdir/qsub.realign_new
	echo "#PBS -V" > $qsub3
	echo "#PBS -A $pbsprj" >> $qsub3
	echo "#PBS -N ${pipeid}_realign_new" >> $qsub3
        echo "#PBS -l epilogue=$epilogue" >> $qsub3
	echo "#PBS -l walltime=$pbscpu" >> $qsub3
	echo "#PBS -l nodes=1:ppn=16" >> $qsub3
	echo "#PBS -o $realignlogdir/log.realign_new.ou" >> $qsub3
	echo "#PBS -e $realignlogdir/log.realign_new.in" >> $qsub3
	echo "#PBS -q $pbsqueue" >> $qsub3
	echo "#PBS -m ae" >> $qsub3
	echo "#PBS -M $email" >> $qsub3
        if [ `expr length ${JOBSmayo}` -gt 0 ]
        then
            echo "#PBS -W depend=afterok:$JOBSmayo" >> $qsub3
        else
            echo "#PBS -W depend=afterok:$JOBSncsa" >> $qsub3
        fi
        echo "$scriptdir/realign_new.sh $realigndir $realignlogdir $listfiles $runfile $realrecalflag $realignlogdir/log.realign_new.in $realignlogdir/log.realign_new.ou $email $realignlogdir/qsub.realign_new" >> $qsub3
	`chmod a+r $qsub3`               
        realrecaljob=`qsub $qsub3` 
        `qhold -h u $realrecaljob`
        echo $realrecaljob >> $output_logs/RECALLpbs

        `qrls -h u $alignids`
        `qrls -h u $mergeids`
        `qrls -h u $sortedmayoids`
        `qrls -h u $extraids`
        `qrls -h u $revbamids`
        `qrls -h u $realrecaljob`

        echo "done realig/recalibrating  all bam files."
        echo `date`
	`chmod -R 770 $outputdir/`
	`chmod -R 770 $output_logs/`
fi
