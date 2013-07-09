#!/bin/sh
######################################
#  script to convert bam files back to fastq as pre requisite to alignment
#
######################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 7 ];
then
	MSG= "parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        inputdir=$1
        samplefileinfo=$2
        runfile=$3
        elog=$4
        olog=$5
        email=$6
        qsubfile=$7
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        provenance=$( cat $runfile | grep -w PROVENANCE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE  | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
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

        if [ ! -d $inputdir ]
        then
	    MSG="INPUTDIR=$inputdir directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        if [ $provenance != "SINGLE_SOURCE" -a $provenance != "MULTI_SOURCE" ]
        then
	    MSG="Invalid value for PROVENANCE=$provenance"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        newfqfiles=""
        sep=":"
        output_logs=`dirname $elog`
        pipeid=$( cat $output_logs/MAINpbs )

        while read sampledetail
        do
	    if [ `expr ${#sampledetail}` -lt 7 ]
            then
                echo "skipping empty line"
            else
		prefix=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f2 )
                if [ ! -s $prefix ]
                then
                    MSG="$prefix BAM file not found. bam2fastq failed."
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
                fi
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
	`qsub $qsub2 >> $output_logs/UPDATECONFIGpbs`

        allconjobs=$( echo $CONVERTids | tr ":" " " )
        `qrls -h u $allconjobs`
	echo `date`
fi
