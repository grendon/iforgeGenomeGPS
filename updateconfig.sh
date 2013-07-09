#!/bin/sh
######################################
#  script to convert bam files back to fastq as pre requisite to alignment
#
######################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu,lmainzer@illinois.edu
if [ $# != 8 ];
then
	MSG= "parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        inputdir=$1
        newfqfiles=$2
        runfile=$3
        samplefileinfo=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        provenance=$( cat $runfile | grep -w PROVENANCE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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

        if [ $provenance != "SINGLE_SOURCE" -a $provenance != "MULTI_SOURCE" ]
        then
	    MSG="Invalid value for PROVENANCE=$provenance"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        # checking inputdir for fastq files tha need merging by sample

        cd $inputdir
        newsamples=""
        sep=" "
        newnames="\n"
        bam_samples=""
        allbams2fq=$inputdir/fqfiles4merging
        touch $allbams2fq
        fqfiles=$( echo $newfqfiles | tr ":" " " )

        if [ `find -name "*.fastq" | wc -l` -gt 0 -a $provenance == "SINGLE_SOURCE" ]
        then
           echo "provenance=$provenance. fastq files in $inputdir need further processing..."
           for fqfile in $fqfiles
           do
               thesamples=`find -name "*_1.${fqfile}" | sed 's/.\///g' | sed 's/.[0-9]*.fastq/.fastq/'`
               echo -e $thesamples >> $allbams2fq
	  done
          bam_samples=$( cat $allbams2fq | tr " " "\n" | sort | uniq -c | sed 's/^ *//g'  | tr " " ":")
          for line in $bam_samples
          do
               nline=$( echo $line | sed 's/:://g' | sed 's/^://' )
	       frq=$( echo $nline | cut -d ':' -f1 )
	       b2fsample=$( echo $nline | cut -d ':' -f2 | sed 's/_1.fastq//' )
               if [ $frq == "1" ]
               then
                  echo "samplename is unique"
                  readone=${b2fsample}_r1.fastq
                  touch $readone
                  `mv ${b2fsample}_1.*.fastq $readone`
                  exitcode=$?
                  if [ -s $readone -a $exitcode -eq 0 ]
                  then
                      readtwo=${b2fsample}_r2.fastq
                      touch $readtwo
                      `mv $b2fsample_2.*.fastq $readtwo`
                      exitcode=$?
                      if [ $exitcode -eq 0 ]
                      then
                          if [ -s $readtwo ]
                          then
			      newnames="FASTQ:${b2fsample}=$inputdir/$readone $inputdir/${readtwo}\n$newnames"
			  else
			      newnames="FASTQ:${b2fsample}=$inputdir/${readone}\n$newnames"
                          fi
		      fi
                      newsamples=${b2fsample}${sep}$newsamples
                  else
                      MSG="$b2fsample Empty fastq file.  exitcode=$exitcode updateconfi failed"
		      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		      exit 1;
                  fi
               else
                  if [ $frq != "" ]
                  then
                      echo "multisamples with same name. concatenate them..."
                      readone=${b2fsample}_r1.fastq
                      touch $readone
                      cat ${b2fsample}_1.*.fastq >> $readone
                      if [ -s $readone ]
                      then
                          readtwo=${b2fsample}_r2.fastq
                          touch $readtwo
			  cat ${b2fsample}_2.*.fastq >> $readtwo
                          if [ -s $readtwo ]
                          then
			      newnames="FASTQ:${b2fsample}=$inputdir/$readone $inputdir/${readtwo}\n$newnames"
			  else
			      newnames="FASTQ:${b2fsample}=$inputdir/${readone}\n$newnames"
			  fi
			  newsamples=${b2fsample}${sep}$newsamples
                          #`rm ${b2fsample}_1.*fastq`
                          #`rm ${b2fsample}_2.*fastq`
                      else
			  MSG="$readone Empty fastq file. updateconfig failed"
			  echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			  exit 1;
                      fi
                  fi
               fi           
           done
	else
           if [ $provenance == "MULTI_SOURCE" ]
           then
               echo "provenance=$provenance. fastq files in input dir do not need further processing"
               for fqfile in $fqfiles 
               do
		   if [ `find -name "*.${fqfile}" | wc -l` -gt 1 ]
                   then
		       R1=`find -name "*_R1.$fqfile" | sed 's/.\///g'`
                       R2=`find -name "*_R2.$fqfile" | sed 's/.\///g'`
                       newsample=$( echo $R1 | sed 's/_R1.[0-9]*.fastq//' )
		       newnames="FASTQ:${newsample}=$inputdir/$R1 $inputdir/${R2}\n$newnames"
                       newsamples=${newsample}${sep}${newsamples}
                   else
		       R1=`find -name "*_R1.$fqfile" | sed 's/.\///g'`
                       R2=`find -name "*_R2.$fqfile" | sed 's/.\///g'`
                       newsample=$( echo $R1 | sed 's/_R1.[0-9]*.fastq//' )
		       newnames="FASTQ:${newsample}=$inputdir/$R1 $inputdir/${R2}\n$newnames"
                       newsamples=${newsample}${sep}${newsamples}
                   fi
               done
           fi
        fi

        # updating BOTH config files
        if [ `expr ${#newnames}` -gt 0 ]
        then
            directory=`dirname $samplefileinfo`
            oldfile=`basename $samplefileinfo`
            oldfilename=$directory/$oldfile.old
            newfilename=$directory/$oldfile
            mv $samplefileinfo $oldfilename
            echo -e $newnames >> $newfilename
            if [ ! -s $newfilename ]
            then
		MSG="Empty file samplefilenames after bam2fastq conversion. updateconfig failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        else
	    MSG="Empty file samplefilenames after bam2fastq conversion. updateconfig failed"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        # updating runfile
        directory=`dirname $runfile`
        oldfile=`basename $runfile`
        oldrun=$directory/$oldfile.old
        newrun=$directory/$oldfile
        mv $runfile $oldrun
        `perl $scriptdir/lned.pl $oldrun $newrun SAMPLENAMES "$newsamples"`
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
 	    MSG="lned.pl failed.  exitcode=$exitcode. update config files after bam2fastq conversion failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        else
            if [ ! -s $newrun ]
            then
		MSG="Empty configuration files after bam2fastq conversion"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        fi
	echo `date`
fi
