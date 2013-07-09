#!/bin/sh
#
# align.sh
# First module in the GGPS analysis pipeline
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu, lmainzer@illinois.edu
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
	inputdir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
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
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        inputdir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )

        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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

        if [ -z $epilogue ]
        then
           MSG="Value for EPILOGUE must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
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

        if [ ! -d $inputdir ]
        then
           mkdir -p $inputdir
        fi

        if [ ! -d $alignerdir ]
        then
           MSG="$alignerdir aligner directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        #resetting output directories, logs, files

        oualigndir=$outputdir/align
        output_logs=$outputdir/logs
        pipeid=$( cat $output_logs/MAINpbs )
        
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

        CONVERT=""
        case=""

        # QC preprocessing on the samplefileinfo file

        while read sampledetail 
        do
          if [ `expr ${#sampledetail}` -lt 5 ]
            then
              echo "skipping empty line"
            else
              echo "reading next line: $sampledetail"
              inputype=$( echo $sampledetail | cut -d ':' -f1 | tr '[a-z]' '[A-Z]' )
              sampleTag=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f1 )
              samplefile=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f2 )
              R1=$( echo $sampledetail | cut -d '=' -f2 | cut -d ' ' -f1 )
              R2=$( echo $sampledetail | cut -d ' ' -f2 )

              if [ `expr length ${inputype}` -lt 1 ]
              then
		MSG="Invalid input format. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi

              if [ $inputype != "BAM" -a $inputype != "FASTQ" ]
              then
		MSG="Invalid input format. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi
 
              if [ `expr length ${sampleTag}` -lt 1  ]
              then
		MSG="Invalid samplename. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi
 
              if [ `expr ${#samplefile}` -lt 1  ]
              then
		MSG="inputfile was not specified. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi
	      if [ $inputype == "BAM" ]
	      then
                  if [ ! -s $samplefile ]
                  then
		      MSG="$samplefile input file specified in SAMPLEFILENAMES was not found. alignment stopped"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		      exit 1;
		  else
		      totlines=`tail $samplefile | wc -l`
		      if [ $totlines -lt 1 ]
		      then
			  MSG="$samplefile input file specified in SAMPLEFILENAMES is empty. alignment stopped"
			  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			  exit 1;
                      fi
                  fi
              fi
	      if [ $inputype == "FASTQ" ]
              then 
                  if [ ! -s $R1 ]
                  then
		      MSG="$R1 input file specified in SAMPLEFILENAMES was not found. alignment stopped"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		      exit 1;
		  else
		      if [ `tail $R1 | wc -l` -lt 1 ]
		      then
			  MSG="$R1 input file specified in SAMPLEFILENMAES is empty. alignment stopped"
			  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			  exit 1;
                      fi
                  fi
                  if [ $paired -eq 1 ]
                  then
                      if [ ! -s $R2 ]
                      then 
			  MSG="$R2 input file specified in SAMPLEFILENAMES was not found. alignment stopped"
			  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			  exit 1;
		      else
			  if [ `tail $R2 | wc -l` -lt 1 ]
			  then
			      MSG="$R2 input file specified in SAMPLEFILENAMES is empty. alignment stopped"
			      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			      exit 1;
                          fi
                      fi
                  fi
              fi
          fi
          echo "bottom of the loop"
        done < $samplefileinfo

        # launching the alignment module 

        if [ $bamtofastqflag == "YES" -a $inputype == "BAM" ]
        then
            newfqfiles=""
            sep=":"
            ######################
            while read sampledetail
            do
		if [ `expr ${#sampledetail}` -lt 7 ]
		then
                    echo "skipping empty line"
		else
		    prefix=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f2 )
		    suffix=$( echo $RANDOM )
		    tmpfq=$inputdir/$suffix
		    if [ ! -d $tmpfq ]
		    then
			`mkdir $tmpfq`
		    fi
                    newsuffix=${suffix}.fastq
                    newfqfiles=${newsuffix}${sep}${newfqfiles}

		    qsub1=$output_logs/qsub.convertbam.$suffix
		    echo "#PBS -V" > $qsub1
		    echo "#PBS -A $pbsprj" >> $qsub1
		    echo "#PBS -N ${pipeid}_convertbam${suffix}" >> $qsub1
		    echo "#pbs -l epilogue=$epligue" >> $qsub1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub1
		    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
		    echo "#PBS -o $output_logs/log.convertbam.${suffix}.ou" >> $qsub1
		    echo "#PBS -e $output_logs/log.convertbam.${suffix}.in" >> $qsub1
		    echo "#PBS -q $pbsqueue" >> $qsub1
		    echo "#PBS -m ae" >> $qsub1
		    echo "#PBS -M $email" >> $qsub1
		    echo "$scriptdir/convertbam.sh $inputdir $prefix $suffix $tmpfq $runfile $output_logs/log.convertbam.${suffix}.in $output_logs/log.convertbam.${suffix}.ou $email $output_logs/qsub.convertbam.$suffix" >> $qsub1
		    `chmod a+r $qsub1` 
		    conjob=`qsub $qsub1`
                    `qhold -h u $conjob` 
		    echo $conjob >> $output_logs/CONVERTBAMpbs
		    echo `date`
		fi
	    done < $samplefileinfo
            echo `date`

            ## updating config files
        
            CONVERTids=$( cat $output_logs/CONVERTBAMpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )

	    qsub2=$output_logs/qsub.updateconfig
	    echo "#PBS -V" > $qsub2
	    echo "#PBS -A $pbsprj" >> $qsub2
	    echo "#PBS -N ${pipeid}_updateconfig" >> $qsub2
	    echo "#pbs -l epilogue=$epligue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	    echo "#PBS -o $output_logs/log.updateconfig.ou" >> $qsub2
	    echo "#PBS -e $output_logs/log.updateconfig.in" >> $qsub2
	    echo "#PBS -q $pbsqueue" >> $qsub2
	    echo "#PBS -m ae" >> $qsub2
	    echo "#PBS -M $email" >> $qsub2
	    echo "#PBS -W depend=afterok:$CONVERTids" >> $qsub2
	    echo "$scriptdir/updateconfig.sh $inputdir $newfqfiles $runfile $samplefileinfo $output_logs/log.updateconfig.in $output_logs/log.updateconfig.ou $email $output_logs/qsub.updateconfig" >> $qsub2
	    `chmod a+r $qsub2`       
	    updatejob=`qsub $qsub2` 
            echo $updatejob >> $output_logs/UPDATECONFIGpbs

            allconjobs=$( echo $CONVERTids | tr ":" " " )
            `qrls -h u $allconjobs`
	    echo `date`

            ######################

	    qsub1=$output_logs/qsub.main.aln2.afterconversion
	    echo "#PBS -V" > $qsub1
	    echo "#PBS -A $pbsprj" >> $qsub1
	    echo "#PBS -N ${pipeid}_MAINaln2_afterbam2fastq" >> $qsub1
	    echo "#pbs -l epilogue=$epligue" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $output_logs/MAINaln2.afterconversion.ou" >> $qsub1
	    echo "#PBS -e $output_logs/MAINaln2.afterconversion.in" >> $qsub1
	    echo "#PBS -q $pbsqueue" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "#PBS -W depend=afterok:$updatejob" >> $qsub1
	    echo "$scriptdir/alignew2.sh $runfile $output_logs/MAINaln2.afterconversion.in $output_logs/MAINaln2.afterconversion.ou $email $output_logs/qsub.main.aln2.afterconversion" >> $qsub1
	    `chmod a+r $qsub1`               
	    `qsub $qsub1 >> $output_logs/MAINALNpbs`
            case="bam2fastq"
	    echo `date`
        fi

        if [ $bamtofastqflag == "NO" -a $inputype == "BAM" ]
        then
            echo "aligning bam files directly"
            qsub2=$output_logs/qsub.alignbams
            echo "#PBS -V" > $qsub2
            echo "#PBS -A $pbsprj" >> $qsub2
            echo "#PBS -N ${pipeid}_alignbams" >> $qsub2
            echo "#PBS -l epilogue=$epilogue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub2
	    echo "#PBS -o $output_logs/log.alignbams.ou" >> $qsub2
	    echo "#PBS -e $output_logs/log.alignbams.in" >> $qsub2
            echo "#PBS -q $pbsqueue" >> $qsub2
            echo "#PBS -m ae" >> $qsub2
            echo "#PBS -M $email" >> $qsub2
            echo "$scriptdir/alignbam.sh $runfile $output_logs/log.alignbams.in $output_logs/log.alignbams.ou $email $output_logs/qsub.alignbams" >> $qsub2
            `chmod a+r $qsub2`               
            `qsub $qsub2 >> $output_logs/MAINALNpbs`
            case="alignbams"
            echo `date`
        fi

        if [ $bamtofastqflag == "NO" -a $inputype == "FASTQ" ]
        then
            echo "aligning fastq files directly"
            qsub3=$output_logs/qsub.alignfastq
            echo "#PBS -V" > $qsub3
            echo "#PBS -A $pbsprj" >> $qsub3
            echo "#PBS -N ${pipeid}_alignfq" >> $qsub3
            echo "#PBS -l epilogue=$epilogue" >> $qsub3
	    echo "#PBS -l walltime=$pbscpu" >> $qsub3
	    echo "#PBS -l nodes=1:ppn=16" >> $qsub3
	    echo "#PBS -o $output_logs/log.alignfastq.ou" >> $qsub3
	    echo "#PBS -e $output_logs/log.alignfastq.in" >> $qsub3
            echo "#PBS -q $pbsqueue" >> $qsub3
            echo "#PBS -m ae" >> $qsub3
            echo "#PBS -M $email" >> $qsub3
            echo "$scriptdir/alignew2.sh $runfile $output_logs/log.alignfastq.in $output_logs/log.alignfastq.ou $email $output_logs/qsub.alignfastq" >> $qsub3
            `chmod a+r $qsub3` 
            `qsub $qsub3 >> $output_logs/MAINALNpbs`
            case="alignfastq"
            echo `date`
        fi

        if [ `expr length ${case}` -lt 1 ]
        then
           MSG="Alignment module failed to launch. Incompatible values specified in config files bam2fastqflag=$bamtofastqflag input_format=$inputype analysis=$analysis"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi         
fi
