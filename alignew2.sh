#!/bin/sh
#
# alignew2.sh
# align module to be used for input files in fastq format
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
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
	dupparms=$( echo "dup=${dup}_flag=${dupflag}" )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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
        if [ ! -d $alignerdir ]
        then
           MSG="$alignerdir aligner directory not found"
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


        oualigndir=$outputdir/align
        output_logs=$outputdir/logs
        pipeid=$( cat $output_logs/MAINpbs )
        
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


        # alignment loop starts here
        # need to rewrite this module to undo dependence on dirname
        prevname=""
        while read sampledetail
        do
          echo "processing next line in file..."
          if [ `expr ${#sampledetail}` -lt 7 ]
          then
            echo "skipping empty line"
          else
            echo "aligning: $sampledetail"
            dirname=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f1 )
            samplenames=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f2 )
	    R1=$( echo $sampledetail | cut -d ' ' -f1 | cut -d '=' -f2 )
	    R2=$( echo $sampledetail | cut -d ' ' -f2 )

            if [ `expr length ${dirname}` -lt 1  ]
            then
		MSG="parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
 
            prevname=$dirname
            sID=$dirname
            sPU=$dirname
            sSM=$dirname
            sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
            sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
            sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
            RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )

	    if [ ! -s $R1 ]
	    then
		    MSG="$R1 reads file1 not found"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
	    fi
	    totlines=`wc -l $R1 | cut -d ' ' -f 1`
	    if [ $totlines -lt 1 ]
	    then
		    MSG="$R1 reads file1 is empty"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
            fi

	    if [ $paired -eq 1 -a ! -s $R2 ]
	    then
		MSG="$R2 reads file2 not found"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi

            if [ $fastqcflag == "YES" ]
            then
		echo "calculating quality values for fastq file"
		if [ ! -d $fastqcdir ]
		then
		    MSG="FASTQCDIR=$fastqcdir directory not found"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
		fi
		if [ ! -d $outputdir/fastqc ]
                then
                    mkdir $outputdir/fastqc
		fi

                qsub=$output_logs/qsub.fastqcR1.$prevname
		echo "#PBS -V" > $qsub
		echo "#PBS -A $pbsprj" >> $qsub
		echo "#PBS -N ${pipeid}_fastqc_R1_${prevname}" >> $qsub
		echo "#PBS -l epilogue=$epilogue" >> $qsub
		echo "#PBS -l walltime=$pbscpu" >> $qsub
		echo "#PBS -l nodes=1:ppn=16" >> $qsub
		echo "#PBS -o $output_logs/log.fastqc_R1_${prevname}.ou" >> $qsub
		echo "#PBS -e $output_logs/log.fastqc_R1_${prevname}.in" >> $qsub
		echo "#PBS -q $pbsqueue" >> $qsub
		echo "#PBS -m ae" >> $qsub
		echo "#PBS -M $email" >> $qsub
		echo "$scriptdir/fastq.sh $fastqcdir $outputdir/fastqc $fastqcparms $R1 $output_logs/log.fastqc_R1_${prevname}.in $output_logs/log.fastqc_R1_${prevname}.ou $email $output_logs/qsub.fastqc_R1_$prevname" >> $qsub
		`chmod a+r $qsub`
                `qsub $qsub >> $output_logs/FASTQCpbs`

		if [ $paired -eq 1 ]
		then
                    qsub=$output_logs/qsub.fastqc_R2_$prevname
		    echo "#PBS -V" > $qsub
		    echo "#PBS -A $pbsprj" >> $qsub
		    echo "#PBS -N ${pipeid}_fastqc_R2_${prevname}" >> $qsub
		    echo "#PBS -l epilogue=$epilogue" >> $qsub
		    echo "#PBS -l walltime=$pbscpu" >> $qsub
		    echo "#PBS -l nodes=1:ppn=16" >> $qsub
		    echo "#PBS -o $output_logs/log.fastqc_R2_${prevname}.ou" >> $qsub
		    echo "#PBS -e $output_logs/log.fastqc_R2_${prevname}.in" >> $qsub
		    echo "#PBS -q $pbsqueue" >> $qsub
		    echo "#PBS -m ae" >> $qsub
		    echo "#PBS -M $email" >> $qsub
		    echo "$scriptdir/fastq.sh $fastqcdir $outputdir/fastqc $fastqcparms $R2 $output_logs/log.fastqc_R2_${prevname}.in $output_logs/log.fastqc_R2_$prevname.ou $email $output_logs/qsub.fastqc_R2_$prevname" >> $qsub
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
		outputsam=$outputalign/$dirname
	    else
		outputsam=$outputalign/$dirname
	    fi
            if [ ! -d $outputlogs ]
            then
 		mkdir $outputlogs
	    fi
	    `chmod -R 770 $outputalign/`
	    `chmod -R 770 $outputlogs/`

            cd $outputalign

            sortedplain=$outputsam.wrg.sorted.bam
            outsortnodup=$outputsam.nodups.sorted.bam
            outsortwdup=$outputsam.wdups.sorted.bam
            newname1=readone_chunk
            newname2=readtwo_chunk

            ## splitting files into chunks before aligning;
            ## remember that one fastq read is made up of four lines
            chunks=`expr $nodes "-" 1`
            if [ $chunks -lt 1 ]
            then
		chunks=$nodes
            fi
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

            echo "splitting read file1=$R1"
            `split -l $numlines -a 1 -d $R1 $newname1`
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                   MSG="splitting of read file $R1 failed. exitcode=$exitcode"
                   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                   exit $exitcode;
            fi
            if [ $paired -eq 1 ]
            then
                echo "splitting read file2=$R2"
		`split -l $numlines -a 1 -d $R2 $newname2`
		exitcode=$?
		if [ $exitcode -ne 0 ]
                then
                    MSG="splitting of read file $R2 failed.  exitcode=$exitcode"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit $exitcode;
		fi
            fi

            ## now we are ready to distribute and align each chunk
            allfiles=""
            for i in $(seq 0 $chunks)
            do
                echo "step 1: aligning chunk $i... "
		echo `date`
                Rone=${newname1}$i
                if [ ! -s $Rone ]
                then
                   MSG="chunk $i of read file $R1 file not found"
                   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                   exit 1;
                fi
		if [ $paired -eq 1 ]
		then
                    Rtwo=${newname2}$i
                    if [ ! -s $Rtwo ]
                    then
			MSG="chunk $i of read file $R2 file not found"
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
                    echo "#PBS -N ${pipeid}_novoaln_${prevname}_$i" >> $qsub
		    echo "#PBS -l epilogue=$epilogue" >> $qsub
		    echo "#PBS -l walltime=$pbscpu" >> $qsub
		    echo "#PBS -l nodes=1:ppn=16" >> $qsub
		    echo "#PBS -o $outputlogs/log.novoaln.$prevname.node$i.ou" >> $qsub
		    echo "#PBS -e $outputlogs/log.novoaln.$prevname.node$i.in" >> $qsub
                    echo "#PBS -q $pbsqueue" >> $qsub
                    echo "#PBS -m ae" >> $qsub
                    echo "#PBS -M $email" >> $qsub
                    if [ $paired -eq 1 ]
                    then
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.sam $outputsam.node$i.bam $scriptdir $samdir $paired $outputalign/$Rone $outputalign/$Rtwo $outputlogs/log.novoaln.$prevname.node$i.in $outputlogs/log.novoaln.$prevname.node$i.ou $email $outputlogs/qsub.novoaln.$prevname.node$i" >> $qsub
                    else
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.sam $outputsam.node$i.bam $scriptdir $samdir $paired $outputalign/$Rone $outputlogs/log.novoaln.$prevname.node$i.in $outputlogs/log.novoaln.$prevname.node$i.ou $email $outputlogs/qsub.novoaln.$prevname.node$i" >> $qsub
                    fi
                    `chmod a+r $qsub`
                    jobnovo=`qsub $qsub`
                    `qhold -h u $jobnovo`
		    echo $jobnovo >> $outputlogs/ALIGNED_$dirname
		else
                    echo "bwa is used as aligner. input file format is in fastq"
                    qsub1=$outputlogs/qsub.bwar1.$prevname.node$i
                    echo "#PBS -V" > $qsub1
                    echo "#PBS -N ${pipeid}_bwar1_${prevname}_$i" >> $qsub1
		    echo "#PBS -o $outputlogs/log.bwar1.$prevname.node$i.ou" >> $qsub1
		    echo "#PBS -e $outputlogs/log.bwar1.$prevname.node$i.in" >> $qsub1
                    echo "#PBS -A $pbsprj" >> $qsub1
		    echo "#PBS -l epilogue=$epilogue" >> $qsub1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub1
		    echo "#PBS -l nodes=1:ppn=16" >> $qsub1
                    echo "#PBS -q $pbsqueue" >> $qsub1
                    echo "#PBS -m ae" >> $qsub1
                    echo "#PBS -M $email" >> $qsub1
		    echo "$scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.R1.sai $outputalign/$Rone $scriptdir $outputlogs/log.bwar1.$prevname.node$i.in $outputlogs/log.bwar1.$prevname.node$i.ou $email $outputlogs/qsub.bwar1.$prevname.node$i" >> $qsub1

                    `chmod a+r $qsub1`
                    jobr1=`qsub $qsub1`
                    `qhold -h u $jobr1`
                    echo $jobr1 >> $outputlogs/ALIGNED_$dirname
                    if [ $paired -eq 1 ]
                    then
                        echo "bwa aligner. paired-end reads"
			qsub2=$outputlogs/qsub.bwar2.$prevname.node$i
			echo "#PBS -V" > $qsub2
			echo "#PBS -N ${pipeid}_bwar2_${prevname}_$i" >> $qsub2
			echo "#PBS -o $outputlogs/log.bwar2.$prevname.node$i.ou" >> $qsub2
			echo "#PBS -e $outputlogs/log.bwar2.$prevname.node$i.in" >> $qsub2
			echo "#PBS -A $pbsprj" >> $qsub2
			echo "#PBS -l epilogue=$epilogue" >> $qsub2
			echo "#PBS -l walltime=$pbscpu" >> $qsub2
			echo "#PBS -l nodes=1:ppn=16" >> $qsub2
			echo "#PBS -q $pbsqueue" >> $qsub2
			echo "#PBS -m ae" >> $qsub2
			echo "#PBS -M $email" >> $qsub2
			echo "$scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $outputalign $outputsam.node$i.R2.sai $outputalign/$Rtwo $scriptdir $outputlogs/log.bwar2.$prevname.node$i.in $outputlogs/log.bwar2.$prevname.node$i.ou $email $outputlogs/qsub.bwar2.$prevname.node$i" >> $qsub2
			`chmod a+r $qsub2`
                        jobr2=`qsub $qsub2`
			`qhold -h u $jobr2`
			echo $jobr2 >> $outputlogs/ALIGNED_$dirname
                        bwajobs=$( cat $outputlogs/ALIGNED_$dirname | sed "s/\.[a-z]*//g" | tr "\n" ":" )
			qsub3=$outputlogs/qsub.bwasampe.$prevname.node$i
			echo "#PBS -V" > $qsub3
			echo "#PBS -N ${pipeid}_bwasampe_${prevname}_$i" >> $qsub3
			echo "#PBS -o $outputlogs/log.bwasampe.$prevname.node$i.ou" >> $qsub3
			echo "#PBS -e $outputlogs/log.bwasampe.$prevname.node$i.in" >> $qsub3
			echo "#PBS -A $pbsprj" >> $qsub3
			echo "#PBS -l epilogue=$epilogue" >> $qsub3
			echo "#PBS -l walltime=$pbscpu" >> $qsub3
			echo "#PBS -l nodes=1:ppn=16" >> $qsub3
			echo "#PBS -q $pbsqueue" >> $qsub3
			echo "#PBS -m ae" >> $qsub3
			echo "#PBS -M $email" >> $qsub3
			echo "#PBS -W depend=afterok:$bwajobs" >> $qsub3
			echo "$scriptdir/bwaS2.sh $alignerdir $refdir/$refindexed $outputalign $outputsam.node$i.R1.sai $outputsam.node$i.R2.sai $outputalign/$Rone $outputalign/$Rtwo $outputsam.node$i.sam $outputsam.node$i.bam $samdir $outputlogs/log.bwasampe.$prevname.node$i.in $outputlogs/log.bwasampe.$prevname.node$i.ou $email $outputlogs/qsub.bwasampe.$prevname.node$i" >> $qsub3
			`chmod a+r $qsub3`
                        jobwa=`qsub $qsub3`
			`qhold -h u $jobwa`
			echo $jobwa >> $outputlogs/ALIGNED_$dirname
                    else
                        echo "bwa aligner. single read"
			qsub3=$outputlogs/qsub.bwasamse.$prevname.node$i
			echo "#PBS -V" > $qsub3
			echo "#PBS -N ${pipeid}_bwasamse_${prevname}_$i" >> $qsub3
			echo "#PBS -o $outputlogs/log.bwasamse.$prevname.node$i.ou" >> $qsub3
			echo "#PBS -e $outputlogs/log.bwasamse.$prevname.node$i.in" >> $qsub3
			echo "#PBS -A $pbsprj" >> $qsub3
			echo "#PBS -l epilogue=$epilogue" >> $qsub3
			echo "#PBS -l walltime=$pbscpu" >> $qsub3
			echo "#PBS -l nodes=1:ppn=16" >> $qsub3
			echo "#PBS -q $pbsqueue" >> $qsub3
			echo "#PBS -m ae" >> $qsub3
			echo "#PBS -M $email" >> $qsub3
			echo "#PBS -W depend=afterok:$jobr1" >> $qsub3
			echo "$scriptdir/bwaS3.sh $alignerdir $refdir/$refindexed $outputalign $outputsam.node$i.R1.sai $outputalign/$Rone $outputsam.node$i.sam $outputsam.node$i.bam $samdir $outputlogs/log.bwasamse.$prevname.node$i.in $outputlogs/log.bwasamse.$prevname.node$i.ou $email $outputlogs/qsub.bwasamse.$prevname.node$i" >> $qsub3
			`chmod a+r $qsub3`
                        jobwa=`qsub $qsub3`
			`qhold -h u $jobwa`
                        echo $qsub3 >> $outputlogs/ALIGNED_$dirname
                    fi
                fi

                allfiles=$allfiles" $outputsam.node$i.bam"
		echo `date`
            done

	    ALIGNED=$( cat $outputlogs/ALIGNED_* | sed "s/\.[a-z]*//" | tr "\n" ":" )

	    listfiles=$( echo $allfiles  | tr " " ":" | sed "s/::/:/g" )
            if [ $sortool == "NOVOSORT" ]
            then
                echo "merging aligned chunks with novosort"
		qsub1=$outputlogs/qsub.sortmerge.novosort.$prevname
		echo "#PBS -V" > $qsub1
		echo "#PBS -A $pbsprj" >> $qsub1
		echo "#PBS -N ${pipeid}_sortmerge_novosort_$prevname" >> $qsub1
		echo "#PBS -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=$pbscpu" >> $qsub1
		echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		echo "#PBS -o $outputlogs/log.sortmerge_novosort.$prevname.ou" >> $qsub1
		echo "#PBS -e $outputlogs/log.sortmerge_novosort.$prevname.in" >> $qsub1
		echo "#PBS -q $pbsqueue" >> $qsub1
		echo "#PBS -m ae" >> $qsub1
		echo "#PBS -M $email" >> $qsub1
		echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
		echo "$scriptdir/mergenovo.sh $outputalign $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $outputlogs/log.sortmerge_novosort.$prevname.in $outputlogs/log.sortmerge_novosort.$prevname.ou $email $outputlogs/qsub.sortmerge.novosort.$prevname" >> $qsub1
		`chmod a+r $qsub1`
		mergejob=`qsub $qsub1`
		`qhold -h u $mergejob`
		echo $mergejob  >> $outputlogs/MERGED_$prevname
            else
                echo "merging aligned chunks with picard"
		qsub1=$outputlogs/qsub.sortmerge.picard.$prevname
		echo "#PBS -V" > $qsub1
		echo "#PBS -A $pbsprj" >> $qsub1
		echo "#PBS -N ${pipeid}_sortmerge_picard_$prevname" >> $qsub1
		echo "#PBS -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=$pbscpu" >> $qsub1
		echo "#PBS -l nodes=1:ppn=16" >> $qsub1
		echo "#PBS -o $outputlogs/log.sortmerge.picard.$prevname.ou" >> $qsub1
		echo "#PBS -e $outputlogs/log.sortmerge.picard.$prevname.in" >> $qsub1
		echo "#PBS -q $pbsqueue" >> $qsub1
		echo "#PBS -m ae" >> $qsub1
		echo "#PBS -M $email" >> $qsub1
		echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
		echo "$scriptdir/mergepicard.sh $outputalign $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $outputlogs/log.sortmerge.picard.$prevname.in $outputlogs/log.sortmerge.picard.$prevname.ou $email $outputlogs/qsub.sortmerge.picard.$prevname" >> $qsub1
		`chmod a+r $qsub1`
		mergejob=`qsub $qsub1`
		`qhold -h u $mergejob`
		echo $mergejob  >> $outputlogs/MERGED_$prevname
            fi

	    echo `date`
	    echo "extract reads specified in CHRINDEX param"
	    qsub5=$outputlogs/qsub.extractreadsbam.$prevname
	    echo "#PBS -V" > $qsub5
	    echo "#PBS -A $pbsprj" >> $qsub5
	    echo "#PBS -N ${pipeid}_extrbam_$prevname" >> $qsub5
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
	    extrajob=`qsub $qsub5`
            `qhold -h u $extrajob`
            echo $extrajob >> $output_logs/EXTRACTREADSpbs

	    cat $outputlogs/ALIGNED_$prevname >> $output_logs/ALIGNEDpbs
	    cat $outputlogs/MERGED_$prevname >> $output_logs/MERGEDpbs
          fi
	done < $samplefileinfo

        ## wrapup
        
	pbsids=$( cat $output_logs/MERGEDpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )
	extraids=$( cat $output_logs/EXTRACTREADSpbs | sed "s/\.[a-z]*//" | tr "\n" " " )
        mergeids=$( echo $pbsids | tr ":" " " )
        alignids=$( cat $output_logs/ALIGNEDpbs | sed "s/\.[a-z]*//" | tr "\n" " " )
	echo $pbsids >> $output_logs/ALN_NCSA_jobids

        ## generating summary redmine email if analysis ends here
	echo "wrap up and produce summary table if analysis ends here or call realign if analysis continues"
	if [ $analysis == "ALIGNMENT" -o $analysis == "ALIGN" -o $analysis == "ALIGN_ONLY" ]
	then
            # release all held jobs
            `qrls -h u $alignids`
            `qrls -h u $mergeids`
            `qrls -h u $extraids`
     
	    lastjobid=""
            cleanjobid=""

            if [ $cleanupflag == "YES" ]
            then 
		qsub6=$output_logs/qsub.cleanup.align
		echo "#PBS -V" > $qsub6
		echo "#PBS -A $pbsprj" >> $qsub6
		echo "#PBS -N ${pipeid}_cleanup_aln" >> $qsub6
		echo "#PBS -l epilogue=$epilogue" >> $qsub6
		echo "#PBS -l walltime=$pbscpu" >> $qsub6
		echo "#PBS -l nodes=1:ppn=1" >> $qsub6
		echo "#PBS -o $output_logs/log.cleanup.align.ou" >> $qsub6
		echo "#PBS -e $output_logs/log.cleanup.align.in" >> $qsub6
		echo "#PBS -q $pbsqueue" >> $qsub6
		echo "#PBS -m ae" >> $qsub6
		echo "#PBS -M $email" >> $qsub6
		echo "#PBS -W depend=afterok:$pbsids" >> $qsub6
		echo "$scriptdir/cleanup.sh $outputdir $analysis $output_logs/log.cleanup.align.in $output_logs/log.cleanup.align.ou $email $output_logs/qsub.cleanup.align"  >> $qsub6
		`chmod a+r $qsub6`
		cleanjobid=`qsub $qsub6`
		echo $cleanjobid >> $outputdir/logs/CLEANUPpbs
            fi
	    qsub4=$output_logs/qsub.summary.aln.allok
	    echo "#PBS -V" > $qsub4
	    echo "#PBS -A $pbsprj" >> $qsub4
	    echo "#PBS -N ${pipeid}_summaryok" >> $qsub4
	    echo "#PBS -l epilogue=$epilogue" >> $qsub4
	    echo "#PBS -l walltime=$pbscpu" >> $qsub4
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub4
	    echo "#PBS -o $output_logs/log.summary.aln.ou" >> $qsub4
	    echo "#PBS -e $output_logs/log.summary.aln.in" >> $qsub4
	    echo "#PBS -q $pbsqueue" >> $qsub4
	    echo "#PBS -m ae" >> $qsub4
	    echo "#PBS -M $email" >> $qsub4
            if [ `expr ${#cleanjobid}` -gt 0 ]
            then
		echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub4
            else
		echo "#PBS -W depend=afterok:$pbsids" >> $qsub4
            fi
	    echo "$scriptdir/summary.sh $outputdir $email exitok"  >> $qsub4
	    `chmod a+r $qsub4`
	    lastjobid=`qsub $qsub4`
	    echo $lastjobid >> $output_logs/SUMMARYpbs

	    if [ `expr ${#lastjobid}` -lt 1 ]
	    then
		echo "at least one job aborted"
		qsub5=$output_logs/qsub.summary.aln.afterany
		echo "#PBS -V" > $qsub5
		echo "#PBS -A $pbsprj" >> $qsub5
		echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub5
		echo "#PBS -l epilogue=$epilogue" >> $qsub5
		echo "#PBS -l walltime=$pbscpu" >> $qsub5
		echo "#PBS -l nodes=1:ppn=1" >> $qsub5
		echo "#PBS -o $output_logs/log.summary.aln.afterany.ou" >> $qsub5
		echo "#PBS -e $output_logs/log.summary.aln.afterany.in" >> $qsub5
		echo "#PBS -q $pbsqueue" >> $qsub5
		echo "#PBS -m ae" >> $qsub5
		echo "#PBS -M $email" >> $qsub5
		echo "#PBS -W depend=afterany:$pbsids" >> $qsub5
		echo "$scriptdir/summary.sh $outputdir $email exitnotok"  >> $qsub5
		`chmod a+r $qsub5`
		badjobid=`qsub $qsub5`
		echo $badjobid >> $output_logs/SUMMARYpbs
	    fi
	fi

	if [ $analysis == "REALIGNMENT" -o $analysis == "REALIGN" ]
	then
            echo " analysis continues with realignment"
	    qsub2=$output_logs/qsub.main.realn
	    echo "#PBS -V" > $qsub2
	    echo "#PBS -A $pbsprj" >> $qsub2
	    echo "#PBS -N ${pipeid}_MAINrealn" >> $qsub2
	    echo "#pbs -l epilogue=$epligue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	    echo "#PBS -o $output_logs/MAINrealn.ou" >> $qsub2
	    echo "#PBS -e $output_logs/MAINrealn.in" >> $qsub2
	    echo "#PBS -q $pbsqueue" >> $qsub2
	    echo "#PBS -m ae" >> $qsub2
	    echo "#PBS -M $email" >> $qsub2
            echo "#PBS -W depend=afterany:$pbsids" >> $qsub2
	    echo "$scriptdir/realign.sh $runfile $output_logs/MAINrealn.in $output_logs/MAINrealn.ou $email $output_logs/qsub.main.realn" >> $qsub2
	    `chmod a+r $qsub2` 
	    `qsub $qsub2 >> $output_logs/MAINREALNpbs`
            `qrls -h u $alignids`
            `qrls -h u $mergeids`
            `qrls -h u $extraids`
	    echo `date`
	fi

	`chmod -R 770 $oualigndir`
	`chmod -R 770 $output_logs`
	echo `date`
fi
