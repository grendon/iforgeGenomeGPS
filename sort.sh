#!/bin/sh
########################### 
#               $1              =               outputdir
#		$2		=		outputlogdir
#               $3              =               samfile
#               $4              =               chunks
#               $5              =               sortedfile without dups
#               $6              =               sorted file with dups
#               $7              =               scriptdir
#               $8              =               runfile
###########################


if [ $# != 9 ]
then	
	echo "wrong number of arguments passed to sort";
else
	set -x
	echo `date`
	outputdir=$1
        outputlogs=$2
        samfile=$3
        chunks=$4
        sortedfile=$5
        sortedwdups=$6
        sortednodups=$7
        scriptdir=$8
        runfile=$9


        if [ ! -d $outputdir ]
        then
           echo "$outputdir directory not found"
           exit 1;
        fi

        # parsing the run file again due to limit in number of parms passed
        cd $scriptdir

        dup=$( cat $runfile | grep -w MARKDUP | cut -d '=' -f2 )
        dupflag=$( cat $runfile | grep -w REMOVE_DUP | cut -d '=' -f2 )
        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 )

        if [ $type == "whole_genome" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
        fi
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )

        sID=$( cat $runfile | grep -w SAMPLEID | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )       
        sPU=$( cat $runfile | grep -w SAMPLEPU | cut -d '=' -f2 )
        sSM=$( cat $runfile | grep -w SAMPLESM | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )        
        RGparms=$( echo "RGID=${sID}_RGLB=${sLB}_RGPL=${sPL}_RGPU=${sPU}_RGSM=${sSM}_RGCN=${sCN}" )
        dupparms=$( echo "dup=${dup}_flag=${dupflag}")

        # sorting files
        for i in $(seq 0 $chunks)
        do

                qsub=$outputlogs/qsub.sort.node.$i
                echo "#PBS -V" > $qsub
                echo "#PBS -A $pbsprj" >> $qsub
                echo "#PBS -N sortnode_$i" >> $qsub
		echo "#PBS -l walltime=$pbscpu" >> $qsub
		echo "#PBS -l nodes=1:ppn=16" >> $qsub 
		echo "#PBS -o $outputlogs/log.sort.node$i.ou" >> $qsub
		echo "#PBS -e $outputlogs/log.sort.node$i.in" >> $qsub
                echo "#PBS -q $pbsqueue" >> $qsub
                echo "#PBS -m ae" >> $qsub
                echo "#PBS -M $email" >> $qsub
                echo "$scriptdir/sortnode.sh $picardir $samdir $outputdir $samfile.node$i.sam $samfile.node$i.wrg.sorted.bam $RGparms NCSA" >> $qsub
                `chmod a+r $qsub`               
                `qsub $qsub >> $outputlogs/SORTEDchunk`
		echo `date`
        done
        SORTED=$( cat $outputlogs/SORTEDchunk | sed "s/\.[a-z]*//" | tr "\n" ":" )
        # merging files and mark dups
        qsub1=$outputlogs/qsub.merge.nodes
        echo "#PBS -V" > $qsub1
        echo "#PBS -A $pbsprj" >> $qsub1
        echo "#PBS -N mergeall" >> $qsub1
        echo "#PBS -l walltime=$pbscpu" >> $qsub1
	echo "#PBS -l nodes=1:ppn=16" >> $qsub1 
	echo "#PBS -o $outputlogs/log.mergenodes.ou" >> $qsub1
	echo "#PBS -e $outputlogs/log.mergenodes.in" >> $qsub1
        echo "#PBS -q $pbsqueue" >> $qsub1
        echo "#PBS -m ae" >> $qsub1
        echo "#PBS -M $email" >> $qsub1
        echo "#PBS -W depend=afterok:$SORTED" >> $qsub1
        echo "$scriptdir/mergenodes.sh $outputdir $samfile wrg.sorted.bam $chunks $sortedwdups $sortednodups $sortedfile $dupparms $runfile" >> $qsub1
        `chmod a+r $qsub1`               
        `qsub $qsub1 >> $outputlogs/MERGEDpbs`
	echo `date`
fi
