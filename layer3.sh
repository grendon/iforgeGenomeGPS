#!/bin/sh

########################### 
#		$1		=	       run info file
###########################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 1 ]
then
        MSG="Parameter mismatch."
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found."
           echo -e "program=$scriptfile stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        if [ -z $epilogue ]
        then
	    MSG="Invalid value for parameter EPILOGUE=$epilogue  in configuration file"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ -z $email ]
        then
           MSG="Invalid value for parameter PBSEMAIL=$email in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
           exit 1;
        fi

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
		MSG="Invalid value for parameter TYPE=$type  in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
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
        `chmod 750 $epilogue`
	`cp $runfile $outputdir/runfile.tmp.txt`
        runfile=$outputdir/runfile.txt
        oldrun=$outputdir/runfile.tmp.txt
        dirsamples=`dirname $samplefileinfo`
        samplename=`basename $samplefileinfo`
        newsamplename=$outputdir/$samplename
        `cp $samplefileinfo $newsamplename`
        `perl $scriptdir/lned.pl $oldrun $runfile SAMPLEFILENAMES "$newsamplename"`

        outputlogs=$outputdir/logs
	echo "launching the main pipeline"
        qsub1=$outputlogs/qsub.main
        echo "#PBS -V" > $qsub1
        echo "#PBS -A $pbsprj" >> $qsub1
        echo "#PBS -N MAIN" >> $qsub1
        echo "#PBS -l epilogue=$epilogue" >> $qsub1
	echo "#PBS -l walltime=00:30:00" >> $qsub1
	echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	echo "#PBS -o $outputlogs/MAIN.ou" >> $qsub1
	echo "#PBS -e $outputlogs/MAIN.in" >> $qsub1
        echo "#PBS -q debug" >> $qsub1
        echo "#PBS -m ae" >> $qsub1
        echo "#PBS -M $email" >> $qsub1
        echo "$scriptdir/main2.sh $runfile batch $outputlogs/MAIN.in $outputlogs/MAIN.ou $email $outputlogs/qsub.main" >> $qsub1
        `chmod a+r $qsub1`               
        jobid=`qsub $qsub1`
        pipeid=$( echo $jobid | sed "s/\.[a-z]*//g" )
        echo $pipeid >> $outputlogs/MAINpbs
        echo `date`

        MSG="GGPS pipeline with id:[${pipeid}] started on iforge by username:$USER at: "$( echo `date` )
        LOGS="jobid=${jobid}\nqsubfile=$outputlogs/qsub.main\nrunfile=$outputdir/runfile.txt\nerrorlog=$outputlogs/MAIN.in\noutputlog=$outputlogs/MAIN.ou"
        echo -e "$MSG\n\nDetails:\n\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""

fi
