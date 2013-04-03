#!/bin/sh

########################### 
#		$1		=	       run info file
###########################

if [ $# != 1 ]
then
    echo "Usage: <RUN INFO FILE>";
else
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        if [ !  -s $runfile ]
        then
           MSG="$runfile file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@$HOST"" 
           exit 1;
        fi

	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
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
        else 
            echo "resetting directory"
	    `rm -r $outputdir/*`
            mkdir -p $outputdir/logs
        fi
        `chmod -R 770 $outputdir/`
        outputlogs=$outputdir/logs
	echo "launching the main pipeline"
        qsub1=$outputlogs/qsub.main
        echo "#PBS -V" > $qsub1
        echo "#PBS -A $pbsprj" >> $qsub1
        echo "#PBS -N MAIN" >> $qsub1
	echo "#PBS -l walltime=00:60:00" >> $qsub1
	echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	echo "#PBS -o $outputlogs/MAIN.ou" >> $qsub1
	echo "#PBS -e $outputlogs/MAIN.in" >> $qsub1
        echo "#PBS -q normal" >> $qsub1
        echo "#PBS -m ae" >> $qsub1
        echo "#PBS -M $email" >> $qsub1
        echo "/projects/mayo/scripts/main.sh $runfile batch $outputlogs/MAIN.in $outputlogs/MAIN.ou $email $outputlogs/qsub.main" >> $qsub1
        `chmod a+r $qsub1`               
        `qsub $qsub1 >> $outputlogs/MAINpbs`
	`cp $runfile $outputdir/runfile.txt`
        echo `date`
fi
