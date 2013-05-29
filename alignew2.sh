#!/bin/sh
#
# align.sh
# First module in the GGPS analysis pipeline
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch"
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

        ## parsing run info file

	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        aligner=$( cat $runfile | grep -w ALIGNER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        igvdir=$( cat $runfile | grep -w IGVDIR | cut -d '=' -f2 )
        fastqcdir=$( cat $runfile | grep -w FASTQCDIR | cut -d '=' -f2 )
        fastqcflag=$( cat $runfile | grep -w FASTQCFLAG | cut -d '=' -f2 )
        fastqcparms=$( cat $runfile | grep -w FASTQCPARMS | cut -d '=' -f2 | tr " " "_" )
        bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        bamtofastqparms=$( cat $runfile | grep -w BAM2FASTQPARMS | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
	dup=$( cat $runfile | grep -w MARKDUP  | cut -d '=' -f2 )
        dupflag=$( cat $runfile | grep -w REMOVE_DUP  | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE  | cut -d '=' -f2 )
	dupparms=$( echo "dup=${dup}_flag=${dupflag}")
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        inputdir=$( cat $runfile | grep -w SAMPLEDIR | cut -d '=' -f2 )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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

        if [ $aligner != "NOVOALIGN" -a $aligner != "BWA" ]
        then
            MSG="ALIGNER=$aligner  is not available at this site"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi

        if [ $aligner == "NOVOALIGN" ]
        then
            alignerdir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w NOVOINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w NOVOPARAMS | cut -d '=' -f2 | tr " " "_" )
        fi
        if [ $aligner == "BWA" ]
        then
            alignerdir=$( cat $runfile | grep -w BWADIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w BWAPARAMS | cut -d '=' -f2 | tr " " "_" )
        fi
        if [ -z $epilogue ]
        then
           MSG="Value for EPILOGUE must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
        fi

        if [ -z $sortool ]
        then
           MSG="Value for SORTOOL must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           if [ $sortool != "NOVOSORT" -a $sortool != "PICARD" ]
           then
               MSG="Invalid value for SORTOOL=$sortool in configuration file"
               echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
               exit 1;
           fi
        fi      
        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
           MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $outputdir ]
        then
           mkdir -p $outputdir
        fi

        if [ ! -s $refdir/$ref ]
        then
           MSG="$refdir/$ref reference genome not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -s $refdir/$refindexed ]
        then
           MSG="$refdir/$refindexed index for reference genome not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $inputdir ]
        then
           MSG="INPUTDIR=$inputdir sample directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""           
	   exit 1;
        fi
        if [ ! -d $alignerdir ]
        then
           MSG="$alignerdir aligner directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $fastqcdir ]
        then
           MSG="FASTQCDIR=$fastqcdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $picardir ]
        then
           MSG="PICARDIR=$picardir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $samdir ]
        then
           MSG="SAMDIR=$samdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        oualigndir=$outputdir/align
        output_logs=$outputdir/logs
        chunks=`expr $nodes "-" 1`
        if [ $chunks -lt 1 ]
        then
	    chunks=$nodes
        fi
        nthreads=`expr $thr "-" 1`
        if [ $nthreads -lt 1 ]
        then
	    nthreads=$thr
        fi
        igv=$outputdir/$igvdir
        extradir=$outputdir/extractreads
        #resetting output directories, logs, files
        if [ -d $oualigndir ]
        then
           echo "$oualigndir is there; resetting it"
           `rm -r $oualigndir/*`
        else
           mkdir -p $oualigndir
        fi

        if [ -d $output_logs ]
        then
           echo "$output_logs is there; resetting it"
           #`rm -r $output_logs/*`
           pbsids=""
        else
           mkdir -p $output_logs
        fi

        # step 1:  conversion of input files from bam to fastq
        # note: a new sampleinfo file will be generated to make things easier
        # down the road with alignment and realignment

        if [ $bamtofastqflag == "YES"  ]
        then
            echo "bam to fastq conversion is needed before aligning"
            qsub=$output_logs/qsub.convert.bam.fastq
            echo "#PBS -V" > $qsub
            echo "#PBS -A $pbsprj" >> $qsub
            echo "#PBS -N bam2fastq" >> $qsub
            echo "#PBS -l epilogue=$epilogue" >> $qsub
	    echo "#PBS -l walltime=$pbscpu" >> $qsub
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub
	    echo "#PBS -o $output_logs/log.bam2fastq.ou" >> $qsub
	    echo "#PBS -e $output_logs/log.bam2fastq.in" >> $qsub
            echo "#PBS -q $pbsqueue" >> $qsub
            echo "#PBS -m ae" >> $qsub
            echo "#PBS -M $email" >> $qsub
	    echo "$scriptdir/bam2fastq.sh $inputdir $samplefileinfo $runfile $output_logs/log.bam2fastq.in $output_logs/log.bam2fastq.ou $email $output_logs/qsub.convert.bam.fastq" >> $qsub
            `chmod a+r $qsub`
            `qsub $qsub >> $output_logs/CONVERTpbs`
        else
	    echo "no need to convert input file from bam to fastq"
        fi

        #
        # Step 2: checking sample names and sample input file details
        if [ ! -s $samplefileinfo ]
        then
           MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        numsamples=0
        for name in $samples
        do
            countnames=$( cat $samplefileinfo | grep $name -c )
            if [ $countnames -lt 1 ]
            then
              MSG="No samples found in SAMPLEFILENAMES=$samplefileinfo."
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
              exit 1;
            fi
            let numsamples+=1
        done
        if [ $numsamples -gt 1 -a $multisample == "YES" ]
        then
            echo "multiple samples to be aligned."
        else
           if [ $numsamples -eq 1 -a $multisample == "NO" ]
           then
              echo "single sample to be aligned."
           else
              MSG="mismatch between number of samples found=$numsamples and vaalue of parameter MULTISAMPLE=$multisample in configuration file."
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
              exit 1;
	   fi
        fi

        CONVERT=$( cat $output_logs/CONVERTpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )

        prevname=""
        counter=0

        #
        # Step 3: alignment loop starts here
        # 
        
        while read sampledetail
        do
          echo "processing next line in file..."
          if [ `expr ${#sampledetail}` -gt 0 ]
          then
            echo "aligning $sampledetail"
            dirname=$( echo $sampledetail | grep ^FASTQ | cut -d ':' -f2 | cut -d '=' -f1 )
            samplenames=$( echo $sampledetail | grep ^FASTQ | cut -d ':' -f2 | cut -d '=' -f2 )

            if [ `expr length ${dirname}` -lt 1  ]
            then
		MSG="parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
 
            ## process a new sample, but not the first one, in multisample case
            if [ $dirname != $prevname -a $counter -gt 0 ]
            then
               ## sort-merge chunks into a single bam file
               ## then mark and or remove duplicates

		echo `date`
		echo "step 3: sort-merging chunks for sample $prevname"

		ALIGNED=$( cat $outputlogs/ALIGNED_$prevname | sed "s/\.[a-z]*//" | tr "\n" ":" )

		listfiles=$( echo $allfiles  | tr " " ":" | sed "s/::/:/g" )
		allfiles=""
		qsub1=$outputlogs/qsub.merge.novosort.$prevname
		echo "#PBS -V" > $qsub1
		echo "#PBS -A $pbsprj" >> $qsub1
		echo "#PBS -N mergenovo_$prevname" >> $qsub1
		echo "#PBS -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=$pbscpu" >> $qsub1
		echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		echo "#PBS -o $outputlogs/log.novosort.$prevname.ou" >> $qsub1
		echo "#PBS -e $outputlogs/log.novosort.$prevname.in" >> $qsub1
		echo "#PBS -q $pbsqueue" >> $qsub1
		echo "#PBS -m ae" >> $qsub1
		echo "#PBS -M $email" >> $qsub1
		echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
		echo "$scriptdir/mergenovo.sh $outputalign $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $outputlogs/log.novosort.$prevname.in $outputlogs/log.novosort.$prevname.ou $email $outputlogs/qsub.merge.novosort.$prevname" >> $qsub1
		`chmod a+r $qsub1`
		mergejob=`qsub $qsub1`
		echo $mergejob  >> $outputlogs/MERGED_$prevname

		echo `date`
		echo "step 4: extract reads specified in CHRINDEX param"
		qsub5=$outputlogs/qsub.extractreadsbam.$prevname
		echo "#PBS -V" > $qsub5
		echo "#PBS -A $pbsprj" >> $qsub5
		echo "#PBS -N extrbam_$prevname" >> $qsub5
		echo "#PBS -l epilogue=$epilogue" >> $qsub5
		echo "#PBS -l walltime=$pbscpu" >> $qsub5
		echo "#PBS -l nodes=1:ppn=16" >> $qsub5
		echo "#PBS -o $outputlogs/log.extractreadsbam.$prevname.ou" >> $qsub5
		echo "#PBS -e $outputlogs/log.extractreadsbam.$prevname.in" >> $qsub5
		echo "#PBS -q $pbsqueue" >> $qsub5
		echo "#PBS -m ae" >> $qsub5
		echo "#PBS -M $email" >> $qsub5
		echo "#PBS -W depend=afterok:$mergejob" >> $qsub5
		echo "$scriptdir/extract_reads_bam.sh $outputalign $outsortwdup $runfile $outputlogs/log.extractreadsbam.$prevname.in $outputlogs/log.extractreadsbam.$prevname.ou $email $outputlogs/qsub.extractreadsbam.$prevname $igv $extradir " >> $qsub5
		`chmod a+r $qsub5`
		extrajob=`qsub $qsub5` 
                echo $extrajob >> $output_logs/EXTRACTREADSpbs

                # resetting some variables
		cat $outputlogs/ALIGNED_$prevname >> $output_logs/ALIGNEDpbs
		cat $outputlogs/MERGED_$prevname >> $output_logs/MERGEDpbs
            fi
        
            # done with previous sample, now processing this sample
            prevname=$dirname
            sID=$dirname
            sPU=$dirname
            sSM=$dirname
            sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
            sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
            sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )

            RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )


            cd $inputdir
	    R1=$( echo $samplenames | cut -d ' ' -f1 )
            if [ `expr length ${R1}` -lt 1  ]
            then
		MSG="parsing SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            else
		if [ ! -s $inputdir/$R1 ]
		then
		    MSG="$inputdir/$R1 reads file1 not found"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
		fi
		totlines=`wc -l $R1 | cut -d ' ' -f 1`
		if [ $totlines -lt 1 ]
		then
		    MSG="$inputdir/$R1 reads file is empty"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
                fi
            fi

	    if [ $paired -eq 1 ]
	    then
		R2=$( echo $samplenames | cut -d ' ' -f2 )
	        if [ $R1 == $R2 ]
                then
		    MSG="a single reads file was found when two paired-end reads files were expected"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
                fi
		if [ `expr length ${R2}` -lt 1  ]
		then
		    MSG="parsing SAMPLEFILENAMES failed. alignment stopped"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit 1;
		else
		    if [ ! -s $inputdir/$R2 ]
		    then
			MSG="$inputdir/$R2 reads file2 not found"
                        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    
                    fi
		    totlinesr2=`wc -l $R1 | cut -d ' ' -f 1`
		    if [ $totlinesr2 -lt 1 ]
		    then
			MSG="$inputdir/$R2 reads file2 is empty"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
		    fi

		fi
	    fi

           
            # done preprocessing of input files
            # let us launch this quality control job before we forget

            if [ $fastqcflag == "YES" ]
            then
		echo "calculating quality values for fastq file"
		if [ ! -d $outputdir/fastqc ]
                then
                    mkdir $outputdir/fastqc
		fi

                qsub=$output_logs/qsub.cal.fastqcr1.$prevname.R1
		echo "#PBS -V" > $qsub
		echo "#PBS -A $pbsprj" >> $qsub
		echo "#PBS -N fastqc1_${prevname}.R1" >> $qsub
		echo "#PBS -l epilogue=$epilogue" >> $qsub
		echo "#PBS -l walltime=$pbscpu" >> $qsub
		echo "#PBS -l nodes=1:ppn=16" >> $qsub
		echo "#PBS -o $output_logs/log.fastqc1.$prevname.R1.ou" >> $qsub
		echo "#PBS -e $output_logs/log.fastqc1.$prevname.R1.in" >> $qsub
		echo "#PBS -q $pbsqueue" >> $qsub
		echo "#PBS -m ae" >> $qsub
		echo "#PBS -M $email" >> $qsub
		echo "$scriptdir/fastq.sh $fastqcdir $outputdir/fastqc $fastqcparms $inputdir/$R1 $output_logs/log.fastqc1.$prevname.R1.in $output_logs/log.fastqc1.$prevname.R1.ou $email $output_logs/qsub.cal.fastqcr1.$prevname.R1" >> $qsub
		`chmod a+r $qsub`
                `qsub $qsub >> $output_logs/FASTQCpbs`

		if [ $paired -eq 1 ]
		then
                    qsub=$output_logs/qsub.cal.fastqcr2.$prevname.R2
		    echo "#PBS -V" > $qsub
		    echo "#PBS -A $pbsprj" >> $qsub
		    echo "#PBS -N fastqc2_${prevname}.R2" >> $qsub
		    echo "#PBS -l epilogue=$epilogue" >> $qsub
		    echo "#PBS -l walltime=$pbscpu" >> $qsub
		    echo "#PBS -l nodes=1:ppn=16" >> $qsub
		    echo "#PBS -o $output_logs/log.fastqc2.$prevname.R2.ou" >> $qsub
		    echo "#PBS -e $output_logs/log.fastqc2.$prevname.R2.in" >> $qsub
		    echo "#PBS -q $pbsqueue" >> $qsub
		    echo "#PBS -m ae" >> $qsub
		    echo "#PBS -M $email" >> $qsub
		    echo "$scriptdir/fastq.sh $fastqcdir $outputdir/fastqc $fastqcparms $inputdir/$R2 $output_logs/log.fastqc2.$prevname.R2.in $output_logs/log.fastqc2.$prevname.R2.ou $email $output_logs/qsub.cal.fastqcr2.$prevname.R2" >> $qsub
		    `chmod a+r $qsub`
                    `qsub $qsub >> $output_logs/FASTQCpbs`
		fi
            else
		echo "quality information for fastq files will NOT be calculated."
            fi

            ## done with generating quality info for each read file
            ## All's in order. 
            ## Next step chunking, distributing and aligning input

            outputalign=$oualigndir/$dirname
            outputlogs=$output_logs/align

            if [ ! -d $outputalign ]
            then
		mkdir $outputalign
		let counter+=1
		outputsam=${dirname}_$counter
	    else
		outputsam=${dirname}_$counter
	    fi
            if [ ! -d $outputlogs ]
            then
 		mkdir $outputlogs
	    fi
	    `chmod -R 770 $outputalign/`
	    `chmod -R 770 $outputlogs/`


            sortedplain=$outputsam.wrg.sorted.bam
            outsortnodup=$outputsam.nodups.sorted.bam
            outsortwdup=$outputsam.wdups.sorted.bam

            cd $outputalign
            newname1=readone_${counter}_
            newname2=readtwo_${counter}_

            ## splitting files into chunks before aligning;
            ## remember that one fastq read is made up of four lines
            if [ $chunks -lt 1 ]
            then
                chunks=1
            fi
            totreads=`expr $totlines "/" 4`
            reads4chunk=`expr $totreads "/" $chunks`
            modval=`expr $totreads "%" $chunks`
            numlines=`expr $reads4chunk "*" 4`
            if [ $modval -eq 0  ]
            then
		echo "mod is 0; no reads for last chunk file, one idle node"
		let chunks-=1
            fi

            echo "splitting read file 1=$R1"
            `split -l $numlines -a 1 -d $inputdir/$R1 $newname1`

            if [ $paired -eq 1 ]
            then
                echo "splitting read file 2=$R2"
		`split -l $numlines -a 1 -d $inputdir/$R2 $newname2`
            fi

            ## now we are ready to distribute and align each chunk

            for i in $(seq 0 $chunks)
            do
		echo `date`
                echo "step 1: aligning chunk $i... "
                Rone=$newname1$i
                if [ ! -s $Rone ]
                then
                   MSG="chunk $i of $R1 not found"
                   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                   exit 1;
                fi
		if [ $paired -eq 1 ]
		then
                    Rtwo=$newname2$i
                    if [ ! -s $Rtwo ]
                    then
			MSG="chunk $i of $R2 not found"
                        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
                fi
                if [ $aligner == "NOVOALIGN"  ]
		then
                    echo "novoalign is used as aligner. input file in fastq format"
                    qsub=$outputlogs/qsub.novoaln.$prevname.node$i
                    echo "#PBS -V" > $qsub
                    echo "#PBS -A $pbsprj" >> $qsub
                    echo "#PBS -N novo_${prevname}_$i" >> $qsub
		    echo "#PBS -l epilogue=$epilogue" >> $qsub
		    echo "#PBS -l walltime=$pbscpu" >> $qsub
		    echo "#PBS -l nodes=1:ppn=16" >> $qsub
		    echo "#PBS -o $outputlogs/log.novo.$prevname.node$i.ou" >> $qsub
		    echo "#PBS -e $outputlogs/log.novo.$prevname.node$i.in" >> $qsub
                    echo "#PBS -q $pbsqueue" >> $qsub
                    echo "#PBS -m ae" >> $qsub
                    echo "#PBS -M $email" >> $qsub
		    echo "#PBS -W depend=afterok:$CONVERT" >> $qsub
                    if [ $paired -eq 1 ]
                    then
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.sam $outputsam.node$i.bam $scriptdir $samdir $paired $outputalign/$Rone $outputalign/$Rtwo $outputlogs/log.novo.$prevname.node$i.in $outputlogs/log.novo.$prevname.node$i.ou $email $outputlogs/qsub.novoaln.$prevname.node$i" >> $qsub
                    else
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.sam $outputsam.node$i.bam $scriptdir $samdir $paired $outputalign/$Rone $outputlogs/log.novo.$prevname.node$i.in $outputlogs/log.novo.$prevname.node$i.ou $email $outputlogs/qsub.novoaln.$prevname.node$i" >> $qsub
                    fi
                    `chmod a+r $qsub`
                    jobnovo=`qsub $qsub`
		    echo $jobnovo >> $outputlogs/ALIGNED_$dirname
		else
                    echo "bwa is used as aligner. input file format is in fastq"
                    qsub1=$outputlogs/qsub.bwar1.$prevname.node$i
                    echo "#PBS -V" > $qsub1
                    echo "#PBS -N bwar1_${prevname}_$i" >> $qsub1
		    echo "#PBS -o $outputlogs/log.bwar1.$prevname.node$i.ou" >> $qsub1
		    echo "#PBS -e $outputlogs/log.bwar1.$prevname.node$i.in" >> $qsub1
                    echo "#PBS -A $pbsprj" >> $qsub1
		    echo "#PBS -l epilogue=$epilogue" >> $qsub1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub1
		    echo "#PBS -l nodes=1:ppn=16" >> $qsub1
                    echo "#PBS -q $pbsqueue" >> $qsub1
                    echo "#PBS -m ae" >> $qsub1
                    echo "#PBS -M $email" >> $qsub1
		    echo "#PBS -W depend=afterok:$CONVERT" >> $qsub1
		    echo "$scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.R1.sai $outputalign/$Rone $scriptdir $outputlogs/log.bwar1.$prevname.node$i.in $outputlogs/log.bwar1.$prevname.node$i.ou $email $outputlogs/qsub.bwar1.$prevname.node$i" >> $qsub1

                    `chmod a+r $qsub1`
                    jobr1=`qsub $qsub1`
                    echo $jobr1 >> $outputlogs/ALIGNED_$dirname
                    if [ $paired -eq 1 ]
                    then
                        echo "bwa aligner. paired-end reads"
			qsub2=$outputlogs/qsub.bwar2.$prevname.node$i
			echo "#PBS -V" > $qsub2
			echo "#PBS -N bwar2_${prevname}_$i" >> $qsub2
			echo "#PBS -o $outputlogs/log.bwar2.$prevname.node$i.ou" >> $qsub2
			echo "#PBS -e $outputlogs/log.bwar2.$prevname.node$i.in" >> $qsub2
			echo "#PBS -A $pbsprj" >> $qsub2
			echo "#PBS -l epilogue=$epilogue" >> $qsub2
			echo "#PBS -l walltime=$pbscpu" >> $qsub2
			echo "#PBS -l nodes=1:ppn=16" >> $qsub2
			echo "#PBS -q $pbsqueue" >> $qsub2
			echo "#PBS -m ae" >> $qsub2
			echo "#PBS -M $email" >> $qsub2
			echo "#PBS -W depend=afterok:$CONVERT" >> $qsub2
			echo "$scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.R2.sai $outputalign/$Rtwo $scriptdir $outputlogs/log.bwar2.$prevname.node$i.in $outputlogs/log.bwar2.$prevname.node$i.ou $email $outputlogs/qsub.bwar2.$prevname.node$i" >> $qsub2
			`chmod a+r $qsub2`
                        jobr2=`qsub $qsub2`
			echo $jobr2 >> $outputlogs/ALIGNED_$dirname

			qsub3=$outputlogs/qsub.bwar3.$prevname.node$i
			echo "#PBS -V" > $qsub3
			echo "#PBS -N bwar3_${prevname}_$i" >> $qsub3
			echo "#PBS -o $outputlogs/log.bwar3.$prevname.node$i.ou" >> $qsub3
			echo "#PBS -e $outputlogs/log.bwar3.$prevname.node$i.in" >> $qsub3
			echo "#PBS -A $pbsprj" >> $qsub3
			echo "#PBS -l epilogue=$epilogue" >> $qsub3
			echo "#PBS -l walltime=$pbscpu" >> $qsub3
			echo "#PBS -l nodes=1:ppn=16" >> $qsub3
			echo "#PBS -q $pbsqueue" >> $qsub3
			echo "#PBS -m ae" >> $qsub3
			echo "#PBS -M $email" >> $qsub3
			echo "#PBS -W depend=afterok:$jobr2" >> $qsub3
			echo "$scriptdir/bwaS2.sh $alignerdir $refdir/$refindexed $outputalign $outputsam.node$i.R1.sai $outputsam.node$i.R2.sai $outputalign/$Rone $outputalign/$Rtwo $outputsam.node$i.sam $outputsam.node$i.bam $samdir $outputlogs/log.bwar3.$prevname.node$i.in $outputlogs/log.bwar3.$prevname.node$i.ou $email $outputlogs/qsub.bwar3.$prevname.node$i" >> $qsub3
			`chmod a+r $qsub3`
                        jobwa=`qsub $qsub3`
			echo $jobwa >> $outputlogs/ALIGNED_$dirname
                    else
                        echo "bwa aligner. single read"
			qsub3=$outputlogs/qsub.bwar3.$prevname.node$i
			echo "#PBS -V" > $qsub3
			echo "#PBS -N bwar3_${prevname}_$i" >> $qsub3
			echo "#PBS -o $outputlogs/log.bwar3.$prevname.node$i.ou" >> $qsub3
			echo "#PBS -e $outputlogs/log.bwar3.$prevname.node$i.in" >> $qsub3
			echo "#PBS -A $pbsprj" >> $qsub3
			echo "#PBS -l epilogue=$epilogue" >> $qsub3
			echo "#PBS -l walltime=$pbscpu" >> $qsub3
			echo "#PBS -l nodes=1:ppn=16" >> $qsub3
			echo "#PBS -q $pbsqueue" >> $qsub3
			echo "#PBS -m ae" >> $qsub3
			echo "#PBS -M $email" >> $qsub3
			echo "#PBS -W depend=afterok:$jobr1" >> $qsub3
			echo "$scriptdir/bwaS3.sh $alignerdir $refdir/$refindexed $outputalign $outputsam.node$i.R1.sai $outputalign/$Rone $outputsam.node$i.sam $outputsam.node$i.bam $samdir $outputlogs/log.bwar3.$prevname.node$i.in $outputlogs/log.bwar3.$prevname.node$i.ou $email $outputlogs/qsub.bwar3.$prevname.node$i" >> $qsub3
			`chmod a+r $qsub3`
                        jobwa=`qsub $qsub3`
                        echo $qsub3 >> $outputlogs/ALIGNED_$dirname
                    fi
                fi

                allfiles=$allfiles" $outputalign/$outputsam.node$i.bam"
		echo `date`
            done
	  else
            echo "line is empty in $sampledetail"
          fi
	done < $samplefileinfo

        # processsing the last batch
	echo `date`
	echo "step 3: sort-merging chunks for sample $prevname"

	ALIGNED=$( cat $outputlogs/ALIGNED_$prevname | sed "s/\.[a-z]*//" | tr "\n" ":" )

	listfiles=$( echo $allfiles  | tr " " ":" | sed "s/::/:/g" )

        if [ $sortool == "NOVOSORT" ]
        then
	    qsub1=$outputlogs/qsub.merge.novosort.$prevname
	    echo "#PBS -V" > $qsub1
	    echo "#PBS -A $pbsprj" >> $qsub1
	    echo "#PBS -N mergenovo_$prevname" >> $qsub1
            echo "#PBS -l epilogue=$epilogue" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub1
	    echo "#PBS -o $outputlogs/log.novosort.$prevname.ou" >> $qsub1
	    echo "#PBS -e $outputlogs/log.novosort.$prevname.in" >> $qsub1
	    echo "#PBS -q $pbsqueue" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
	    echo "$scriptdir/mergenovo.sh $outputalign $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $outputlogs/log.novosort.$prevname.in $outputlogs/log.novosort.$prevname.ou $email $outputlogs/qsub.merge.novosort.$prevname" >> $qsub1
	    `chmod a+r $qsub1`
	    mergejob=`qsub $qsub1`
	    echo $mergejob  >> $outputlogs/MERGED_$prevname
        else
	    qsub1=$outputlogs/qsub.sortmerge.picard.$prevname
	    echo "#PBS -V" > $qsub1
	    echo "#PBS -A $pbsprj" >> $qsub1
	    echo "#PBS -N sortmerge_$prevname" >> $qsub1
            echo "#PBS -l epilogue=$epilogue" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub1
	    echo "#PBS -o $outputlogs/log.sortmerge.$prevname.ou" >> $qsub1
	    echo "#PBS -e $outputlogs/log.sortmerge.$prevname.in" >> $qsub1
	    echo "#PBS -q $pbsqueue" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
	    echo "$scriptdir/mergepicard.sh $outputalign $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $outputlogs/log.sortmerge.$prevname.in $outputlogs/log.sortmerge.$prevname.ou $email $outputlogs/qsub.sortmerge.picard.$prevname" >> $qsub1
	    `chmod a+r $qsub1`
	    mergejob=`qsub $qsub1`
	    echo $mergejob  >> $outputlogs/MERGED_$prevname
        fi

	echo `date`
	echo "step 4: extract reads specified in CHRINDEX param"
	qsub5=$outputlogs/qsub.extractreadsbam.$prevname
	echo "#PBS -V" > $qsub5
	echo "#PBS -A $pbsprj" >> $qsub5
	echo "#PBS -N extrbam_$prevname" >> $qsub5
        echo "#PBS -l epilogue=$epilogue" >> $qsub5
	echo "#PBS -l walltime=$pbscpu" >> $qsub5
	echo "#PBS -l nodes=1:ppn=16" >> $qsub5
	echo "#PBS -o $outputlogs/log.extractreadsbam.$prevname.ou" >> $qsub5
	echo "#PBS -e $outputlogs/log.extractreadsbam.$prevname.in" >> $qsub5
	echo "#PBS -q $pbsqueue" >> $qsub5
	echo "#PBS -m ae" >> $qsub5
	echo "#PBS -M $email" >> $qsub5
	echo "#PBS -W depend=afterok:$mergejob" >> $qsub5
	echo "$scriptdir/extract_reads_bam.sh $outputalign $outsortwdup $runfile $outputlogs/log.extractreadsbam.$prevname.in $outputlogs/log.extractreadsbam.$prevname.ou $email  $outputlogs/qsub.extractreadsbam.$prevname $igv $extradir" >> $qsub5
	`chmod a+r $qsub5`
	`qsub $qsub5 >> $output_logs/EXTRACTREADSpbs`

        # resetting some variables
	cat $outputlogs/ALIGNED_$prevname >> $output_logs/ALIGNEDpbs
	cat $outputlogs/MERGED_$prevname >> $output_logs/MERGEDpbs
	chunks=`expr $nodes "-" 1`
        
	pbsids=$( cat $output_logs/MERGEDpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )
	echo $pbsids >> $output_logs/ALN_NCSA_jobids
        echo "done aligning all files specified in the sample_info file."
        echo `date`
	`chmod -R 770 $oualigndir`
	`chmod -R 770 $output_logs`
fi