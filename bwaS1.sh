#!/bin/sh

if [ $# != 11 ]
then
        MSG="parameter mismatch"
        echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s 'GGPS error notification' "$USER@HOST""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        aligndir=$1
        parms=$2
        ref=$3
        outputdir=$4
	outputfile=$5
        R=$6
        scriptdir=$7
        elog=$8
        olog=$9
        email=${10}
        qsubfile=${11}
        LOGS="qsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        parameters=$( echo $parms | tr "_" " " )

        ## checking quality scores to gather additional params
        qscores=$scriptdir/checkFastqQualityScores.pl
        ill2sanger=`perl $qscores $R1 10000`
        if [ $ill2sanger -gt 65 ]
        then
           qual="-I"
        else
           qual=" "
        fi

        cd $outputdir
        $aligndir/bwa aln $parameters $qual $ref $R > $outputfile
        if [ ! -s $outputdir/$outputfile ]
        then
            MSG="$outputdir/$outputfile file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s 'GGPS error notification' "$email""
            exit 1;
        fi
        echo `date`
fi