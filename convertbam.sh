#!/bin/sh
######################################
#  script to convert bam files back to fastq as pre requisite to alignment
#
######################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 9 ];
then
	MSG= "parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        inputdir=$1
        prefix=$2
        suffix=$3
        tmpfq=$4
        runfile=$5
        elog=$6
        olog=$7
        email=$8
        qsubfile=$9
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2 )
        bam2fastqparms=$( cat $runfile | grep -w BAM2FASTQPARMS | cut -d '=' -f2- )
        bam2fastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 )
        javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
        provenance=$( cat $runfile | grep -w PROVENANCE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

        if [ ! -d $picardir ]
        then
	    MSG="PICARDIR=$picardir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      
        if [ ! -d $inputdir ]
        then
	    MSG="INPUTDIR=$inputdir directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
	    MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ -z $javamodule ]
        then
	    MSG="Value for JAVAMODULE must be specified in configuration file"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        else 
            `/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
        fi
        if [ $provenance != "SINGLE_SOURCE" -a $provenance != "MULTI_SOURCE" ]
        then
	    MSG="Invalid value for PROVENANCE=$provenance"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -s $prefix ]
        then
            MSG="$prefix BAM file not found. bam2fastq failed."
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        cd $tmpfq 
        tmpoutsam=$tmpfq/inputbam_revsam

        if [ $provenance == "SINGLE_SOURCE" ]
        then
            echo "bam files come from a single source. split by sample"
            if [ $revertsam == "1" ]
            then
                echo "revertsam then samtofastq..."
                java -Xmx6g -Xms512m -jar $picardir/RevertSam.jar \
                          COMPRESSION_LEVEL=0 \
                          INPUT=$prefix \
                          OUTPUT=$tmpoutsam \
                          VALIDATION_STRINGENCY=SILENT
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="revertsam command failed.  exitcode=$exitcode  bam2fastq conversion failed"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		    exit $exitcode;
		fi
		echo `date`
		java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
                          INPUT=$tmpoutsam \
                          TMP_DIR=$tmpfq \
                          OUTPUT_PER_RG=true \
                          OUTPUT_DIR=$tmpfq \
                          $bam2fastqparms \
                          VALIDATION_STRINGENCY=SILENT 
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		    exit $exitcode;
		fi
		echo `date`
            else
                echo "no-revertsam, just samtofastq..."
		java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
                          INPUT=$prefix \
                          TMP_DIR=$tmpfq \
                          OUTPUT_PER_RG=true \
                          OUTPUT_DIR=$tmpfq \
                          $bam2fastqparms \
                          VALIDATION_STRINGENCY=SILENT 
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		    exit $exitcode;
		fi
		echo `date`
            fi
            `rename _1.fastq _1.${suffix}.fastq *_1.fastq`
            `rename _2.fastq _2.${suffix}.fastq *_2.fastq`
        fi

        if [ $provenance == "MULTI_SOURCE" ]
        then
            echo "bam files from multi-sources. split by R1 and or R2"
            pre=`basename $prefix`
	    R1=${pre}_R1.${suffix}.fastq
	    R2=${pre}_R2.${suffix}.fastq

	    if [ $paired == "1" ]
            then
                echo "paired-end reads will be generated"
		if [ $revertsam == "1" ]
		then
                    echo "revertsam then samtofastq..."
		    java -Xmx6g -Xms512m -jar $picardir/RevertSam.jar \
                                COMPRESSION_LEVEL=0 \
				INPUT=$prefix \
				OUTPUT=$tmpoutsam \
				VALIDATION_STRINGENCY=SILENT  
		    exitcode=$?
		    if [ $exitcode -ne 0 ]
		    then
			MSG="revertsam command failed.  exitcode=$exitcode  bam2fastq conversion failed"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit $exitcode;
		    fi
		    echo `date`
		    java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
				FASTQ=$R1 \
				SECOND_END_FASTQ=$R2 \
				INPUT=$tmpoutsam \
				TMP_DIR=$tmpfq \
				$bam2fastqparms \
				VALIDATION_STRINGENCY=SILENT 
		    exitcode=$?
		    if [ $exitcode -ne 0 ]
		    then
			MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit $exitcode;
		    fi
		    echo `date`
		    if [ ! -s $R1 -a ! -s $R2 ]
		    then
			MSG="$R1 $R2 FASTQ files not created. bam2fastq failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
		else
                    echo "no-revertsam, just samtofastq..."
		    java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
		  		FASTQ=$R1 \
				SECOND_END_FASTQ=$R2 \
				INPUT=$prefix \
				TMP_DIR=$tpmfq \
				$bam2fastqparms \
				VALIDATION_STRINGENCY=SILENT 
		    exitcode=$?
		    if [ $exitcode -ne 0 ]
		    then
			MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit $exitcode;
		    fi
		    echo `date`
		    if [ ! -s $R1 -a ! -s $R2 ]
		    then
		  	MSG="$R1 $R2 FASTQ files not created. bam2fastq conversion failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
		    fi
		fi
            else
                echo "single reads will be generated"
		if [ $revertsam == "1" ]
		then
		    java -Xmx6g -Xms512m -jar $picardir/RevertSam.jar\
                            COMPRESSION_LEVEL=0 \
			    INPUT=$prefix \
			    OUTPUT=$tmpoutsam \
			    VALIDATION_STRINGENCY=SILENT
		    exitcode=$?
		    if [ $exitcode -ne 0 ]
		    then
			MSG="revertsam command failed.  exitcode=$exitcode  bam2fastq conversion failed"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit $exitcode;
		    fi
		    echo `date`
		    java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
		  	    FASTQ=$R1 \
			    INPUT=$tmpoutsam \
			    TMP_DIR=$tmpfq \
			    $bam2fastqparms \
			    VALIDATION_STRINGENCY=SILENT 
		    exitcode=$?
		    if [ $exitcode -ne 0 ]
		    then
			MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit $exitcode;
		    fi
		    echo `date`		    
		    if [ ! -s $R1 ]
		    then
			MSG="$R1 Empty fastq file. bam2fastq failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
		else
		    java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
				FASTQ=$R1 \
				INPUT=$prefix \
				TMP_DIR=$tmpfq \
				$bam2fastqparms \
				VALIDATION_STRINGENCY=SILENT
		    exitcode=$?
		    if [ $exitcode -ne 0 ]
		    then
			MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit $exitcode;
		    fi
		    echo `date`
		    if [ ! -s $R1 ]
		    then
			MSG="$R1 Empty fasq file. bam2fastq conversion failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
		    fi
		fi
	    fi
	fi
        `mv *.fastq $inputdir`
	echo `date`
fi
