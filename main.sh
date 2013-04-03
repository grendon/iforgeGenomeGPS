#!/bin/sh

if [ $# != 6 ]
then
        MSG="parameter mismatch.\nUsage: <run info file><mode><errolog><outputlog> "
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@$HOST"" 
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
        LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ !  -s $runfile ]
        then
           MSG="$runfile file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s 'GGPS error notification' "$mail"" 
           exit 1;
        fi

	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 )
        resortbam=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 )
        bam2fqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )

        if [ $type == "whole_genome" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
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

        if [ $resortbam == "YES" -a $bam2fqflag == "YES" ]
        then
            MSG="resortbam and bam2fastq fields are not set up properly."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
            exit 1;
        fi

        outputlogs=$outputdir/logs

        if [ $analysis == "alignment" ]
        then
            echo "Type of analysis to run: ALIGNMENT only" 

            qsub1=$outputlogs/qsub.main.aln
            echo "#PBS -V" > $qsub1
            echo "#PBS -A $pbsprj" >> $qsub1
            echo "#PBS -N MAINaln" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $outputlogs/MAINaln.ou" >> $qsub1
	    echo "#PBS -e $outputlogs/MAINaln.in" >> $qsub1
            echo "#PBS -q $pbsqueue" >> $qsub1
            echo "#PBS -m ae" >> $qsub1
            echo "#PBS -M $email" >> $qsub1
            echo "/projects/mayo/scripts/align.sh $runfile $outputlogs/MAINaln.in $outputlogs/MAINaln.ou $email $outputlogs/qsub.main.aln" >> $qsub1
            `chmod a+r $qsub1`               
            `qsub $qsub1 >> $outputlogs/MAINALNpbs`
            echo `date`
       else
	    if [ $analysis == "realignment" -a $resortbam == "YES" ]
            then
		echo "Type of analysis to run: REALIGNMENT only. bams provided"
		qsub2=$outputlogs/qsub.main.realn
		echo "#PBS -V" > $qsub2
		echo "#PBS -A $pbsprj" >> $qsub2
		echo "#PBS -N MAINrealn" >> $qsub2
		echo "#PBS -l walltime=$pbscpu" >> $qsub2
		echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		echo "#PBS -o $outputlogs/MAINrealn.ou" >> $qsub2
		echo "#PBS -e $outputlogs/MAINrealn.in" >> $qsub2
		echo "#PBS -q $pbsqueue" >> $qsub2
		echo "#PBS -m ae" >> $qsub2
		echo "#PBS -M $email" >> $qsub2
		echo "/projects/mayo/scripts/realign.sh $runfile $outputlogs/MAINrealn.in $outputlogs/MAINrealn.ou $email $outputlogs/qsub.main.realn" >> $qsub2
		`chmod a+r $qsub2` 
		`qsub $qsub2 >> $outputlogs/MAINREALNpbs`
		echo `date` 

            else
                if [ $analysis == "realignment" ]
                then
		    echo "Type of analysis to run: ALIGNMENT-REALIGNMENT"
		    qsub1=$outputlogs/qsub.main.aln
		    echo "#PBS -V" > $qsub1
		    echo "#PBS -A $pbsprj" >> $qsub1
		    echo "#PBS -N MAINaln" >> $qsub1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub1
		    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
		    echo "#PBS -o $outputlogs/MAINaln.ou" >> $qsub1
		    echo "#PBS -e $outputlogs/MAINaln.in" >> $qsub1
		    echo "#PBS -q $pbsqueue" >> $qsub1
		    echo "#PBS -m ae" >> $qsub1
		    echo "#PBS -M $email" >> $qsub1
		    echo "/projects/mayo/scripts/align.sh $runfile $outputlogs/MAINaln.in $outputlogs/MAINaln.ou $email $outputlogs/qsub.main.aln" >> $qsub1
		    `chmod a+r $qsub1`               
		    `qsub $qsub1 >> $outputlogs/MAINALNpbs`
		    echo `date`

		    jobid=$( cat $outputlogs/MAINALNpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )

		    qsub2=$outputlogs/qsub.main.realn
		    echo "#PBS -V" > $qsub2
		    echo "#PBS -A $pbsprj" >> $qsub2
		    echo "#PBS -N MAINrealn" >> $qsub2
		    echo "#PBS -l walltime=$pbscpu" >> $qsub2
		    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		    echo "#PBS -o $outputlogs/MAINrealn.ou" >> $qsub2
		    echo "#PBS -e $outputlogs/MAINrealn.in" >> $qsub2
		    echo "#PBS -q $pbsqueue" >> $qsub2
		    echo "#PBS -m ae" >> $qsub2
		    echo "#PBS -M $email" >> $qsub2
                    echo "#PBS -W depend=afterok:$jobid" >> $qsub2
		    echo "/projects/mayo/scripts/realign.sh $runfile $outputlogs/MAINrealn.in $outputlogs/MAINrealn.ou $email $outputlogs/qsub.main.realn" >> $qsub2
		    `chmod a+r $qsub2` 
		    `qsub $qsub2 >> $outputlogs/MAINREALNpbs`
		    echo `date` 
                else
		    MSG="analysis=$analysis not available at this site yet."
	            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
                    exit 1; 
               fi
            fi
       fi
fi