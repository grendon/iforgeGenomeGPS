#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 6 ]
then
        MSG="Parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        runmode=$2
        elog=$3
        olog=$4
        email=$5
        qsubfile=$6
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        resortbam=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 )
        bam2fqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
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

        if [ -z $thr -o -z $outputdir -o -z $pbsprj -o -z $epilogue ]
        then
 		MSG="Invalid value specified for any of these paramaters in configuration file:\nPBSTHREADS=$thr\nOUTPUTDIR=$outputdir\nPBSPROJECTID=$pbsprj\nEPILOGUE=$epilogue"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
        fi
        if [ $resortbam == "YES" -a $bam2fqflag == "YES" ]
        then
            MSG="Invalid values for parameters RESORTBAM=$resortbam and BAM2FASTQFLAG=$bam2fqflag in the configuration file."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
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
            mkdir -p $outputdir/logs
	elif [ $runmode != "batch" ]
        then
	    echo "resetting logs"
	    `rm -r $outputdir/logs/*`
        fi
	`chmod -R 770 $outputdir`
        `chmod 750 $epilogue`
        outputlogs=$outputdir/logs
        if [ $analysis == "ALIGN" -o $analysis == "ALIGNMENT" ]
        then
            echo "Type of analysis to run: ALIGNMENT only" 

            qsub1=$outputlogs/qsub.main.aln
            echo "#PBS -V" > $qsub1
            echo "#PBS -A $pbsprj" >> $qsub1
            echo "#PBS -N MAINaln" >> $qsub1
            echo "#pbs -l epilogue=$epligue" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $outputlogs/MAINaln.ou" >> $qsub1
	    echo "#PBS -e $outputlogs/MAINaln.in" >> $qsub1
            echo "#PBS -q $pbsqueue" >> $qsub1
            echo "#PBS -m ae" >> $qsub1
            echo "#PBS -M $email" >> $qsub1
            echo "$scriptdir/alignew2.sh $runfile $outputlogs/MAINaln.in $outputlogs/MAINaln.ou $email $outputlogs/qsub.main.aln" >> $qsub1
            `chmod a+r $qsub1`               
            `qsub $qsub1 >> $outputlogs/MAINALNpbs`
            echo `date`
       else
	    if [ $resortbam == "YES" -a [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ] ]
            then
		echo "Type of analysis to run: REALIGNMENT only. bams provided"
		qsub2=$outputlogs/qsub.main.realn
		echo "#PBS -V" > $qsub2
		echo "#PBS -A $pbsprj" >> $qsub2
		echo "#PBS -N MAINrealn" >> $qsub2
		echo "#pbs -l epilogue=$epligue" >> $qsub2
		echo "#PBS -l walltime=$pbscpu" >> $qsub2
		echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		echo "#PBS -o $outputlogs/MAINrealn.ou" >> $qsub2
		echo "#PBS -e $outputlogs/MAINrealn.in" >> $qsub2
		echo "#PBS -q $pbsqueue" >> $qsub2
		echo "#PBS -m ae" >> $qsub2
		echo "#PBS -M $email" >> $qsub2
		echo "$scriptdir/realign.sh $runfile $outputlogs/MAINrealn.in $outputlogs/MAINrealn.ou $email $outputlogs/qsub.main.realn" >> $qsub2
		`chmod a+r $qsub2` 
		`qsub $qsub2 >> $outputlogs/MAINREALNpbs`
		echo `date` 

            else
                if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ]
                then
		    echo "Type of analysis to run: ALIGNMENT and REALIGNMENT"
		    qsub1=$outputlogs/qsub.main.aln
		    echo "#PBS -V" > $qsub1
		    echo "#PBS -A $pbsprj" >> $qsub1
		    echo "#PBS -N MAINaln" >> $qsub1
		    echo "#pbs -l epilogue=$epligue" >> $qsub1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub1
		    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
		    echo "#PBS -o $outputlogs/MAINaln.ou" >> $qsub1
		    echo "#PBS -e $outputlogs/MAINaln.in" >> $qsub1
		    echo "#PBS -q $pbsqueue" >> $qsub1
		    echo "#PBS -m ae" >> $qsub1
		    echo "#PBS -M $email" >> $qsub1
		    echo "$scriptdir/alignew2.sh $runfile $outputlogs/MAINaln.in $outputlogs/MAINaln.ou $email $outputlogs/qsub.main.aln" >> $qsub1
		    `chmod a+r $qsub1`               
		    `qsub $qsub1 >> $outputlogs/MAINALNpbs`
		    echo `date`

		    jobid=$( cat $outputlogs/MAINALNpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )

		    qsub2=$outputlogs/qsub.main.realn
		    echo "#PBS -V" > $qsub2
		    echo "#PBS -A $pbsprj" >> $qsub2
		    echo "#PBS -N MAINrealn" >> $qsub2
		    echo "#pbs -l epilogue=$epligue" >> $qsub2
		    echo "#PBS -l walltime=$pbscpu" >> $qsub2
		    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		    echo "#PBS -o $outputlogs/MAINrealn.ou" >> $qsub2
		    echo "#PBS -e $outputlogs/MAINrealn.in" >> $qsub2
		    echo "#PBS -q $pbsqueue" >> $qsub2
		    echo "#PBS -m ae" >> $qsub2
		    echo "#PBS -M $email" >> $qsub2
                    echo "#PBS -W depend=afterok:$jobid" >> $qsub2
		    echo "$scriptdir/realign.sh $runfile $outputlogs/MAINrealn.in $outputlogs/MAINrealn.ou $email $outputlogs/qsub.main.realn" >> $qsub2
		    `chmod a+r $qsub2` 
		    `qsub $qsub2 >> $outputlogs/MAINREALNpbs`
		    echo `date` 
                else
		    MSG="Invalid value for parameter ANALYSIS=$analysis in configuration file."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit 1; 
               fi
            fi
       fi

       
fi