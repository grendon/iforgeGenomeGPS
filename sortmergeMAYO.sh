#!/bin/sh
###############################
# input file is in bam format and not produced inhouse
# 
###############################
if [ $# != 8 ] 
then
    echo "SORTMERGEMAYO Usage: parameter mismatch";
else
    set -x
    echo `date`
    inputdir=$1
    outputlogs=$2
    sample=$3
    infile=$4
    outfilewdups=$5
    outfilenodups=$6
    scriptdir=$7
    runfile=$8

    #sanity check
    if [ ! -s $runfile ]
    then
       echo "$runfile file not found"
       exit 1;
    fi

    email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
    pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
    refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
    ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
    type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 )

    if [ $type == "whole_genome" ]
    then
        pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
        pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
    else
        pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
        pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
    fi
    
    paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
    samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
    rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
    multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
    samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
    dup=$( cat $runfile | grep -w MARKDUP | cut -d '=' -f2 )
    dupflag=$( cat $runfile | grep -w REMOVE_DUP | cut -d '=' -f2 )
    dupparms=$( echo "dup=${dup}_flag=${dupflag}" )
    reorderflag=$( cat $runfile | grep -w REORDERSAM | cut -d '=' -f2 )
    sID=$( cat $runfile | grep -w SAMPLEID | cut -d '=' -f2 )
    sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
    sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )       
    sPU=$( cat $runfile | grep -w SAMPLEPU | cut -d '=' -f2 )
    sSM=$( cat $runfile | grep -w SAMPLESM | cut -d '=' -f2 )
    sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )        
    RGparms=$( echo "RGID=${sID}_RGLB=${sLB}_RGPL=${sPL}_RGPU=${sPU}_RGSM=${sSM}_RGCN=${sCN}" )

    outputdir=$1
 
    if [ ! -d $outputdir ]
    then
       echo "$outputdir directory not found"
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       echo "$picardir directory not found"
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       echo "$samdir directory not found"
       exit 1;
    fi

    ### collecting bam files to sort and then merge, which ones??
    ### Mayo only provides one bam file per sample s_normal s_cancer
    ### need to ask

    cd $outputdir
    listfiles=$outputdir/listfiles
    `rm $listfiles`
    for file in $outputdir/*.sorted.bam
    do
       ### chech if file is already sorted
       $samdir/samtools view -H $file > $file.header
       match=$( cat $file.header | grep -A 0 "SO:coordinate" )
       lenmatch=`expr length $match`
       if [ $lenmatch -gt 0 ]
       then
          echo "$file is already sorted"
	  prefix=`basename $file .sorted.bam`
          mv $file $prefix.wdups.sorted.bam
	  echo " INPUT=$prefix.wdups.sorted.bam " >> $listfiles
       else
          echo "$file is not sorted"
	  prefix=`basename $file .sorted.bam`
	  unsorted=$prefix.unsorted.bam
          mv $file $unsorted

          ## pbs sortnode.sh
          qsub=$outputlogs/qsub.sort.$file
          echo "#PBS -V" > $qsub
          echo "#PBS -A $pbsprj" >> $qsub
          echo "#PBS -N sortmayo" >> $qsub
	  echo "#PBS -l walltime=$pbscpu" >> $qsub
	  echo "#PBS -l nodes=1:ppn=16" >> $qsub 
	  echo "#PBS -o $outputlogs/log.sort.$file.ou" >> $qsub
	  echo "#PBS -e $outputlogs/log.sort.$filei.in" >> $qsub
          echo "#PBS -q $pbsqueue" >> $qsub
          echo "#PBS -m ae" >> $qsub
          echo "#PBS -M $email" >> $qsub
          echo "$scriptdir/sortnode.sh $picardir $samdir $outputdir $unsorted $prefix.wdups.sorted.bam $RGparms MAYO" >> $qsub
          `chmod a+r $qsub`               
          #`qsub $qsub >> $outputlogs/SORTEDmayo`

	  echo " INPUT=$prefix.wdups.sorted.bam " >> $listfiles
	  echo `date`         
       fi
    done

    infiles=$( cat $listfiles | tr "\n" " " )

    # ready to merge sorted bam files
    SORTED=$( cat $outputlogs/SORTEDmayo | sed "s/\.[a-z]*//" | tr "\n" ":" )

    qsub1=$outputlogs/qsub.merge.bam
    echo "#PBS -V" > $qsub1
    echo "#PBS -A $pbsprj" >> $qsub1
    echo "#PBS -N mergebams" >> $qsub1
    echo "#PBS -l walltime=$pbscpu" >> $qsub1
    echo "#PBS -l nodes=1:ppn=16" >> $qsub1 
    echo "#PBS -o $outputlogs/log.mergebam.ou" >> $qsub1
    echo "#PBS -e $outputlogs/log.mergebam.in" >> $qsub1
    echo "#PBS -q $pbsqueue" >> $qsub1
    echo "#PBS -m ae" >> $qsub1
    echo "#PBS -M $email" >> $qsub1
    echo "#PBS -W depend=afterok:$SORTED" >> $qsub1
    echo "$scriptdir/mergebams.sh $outputdir $outputlogs $infiles $outfilewdups $outfilenodups $dupparms $reorderflag $runfile" >> $qsub1
    `chmod a+r $qsub1`               
    #`qsub $qsub1 >> $outputlogs/MERGEDmayo`
    echo `date`
        

fi