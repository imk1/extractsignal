#!/bin/bash

# batchExtractSignal.sh [pairFile] [peakDir] [signalDir] [oDir] [format] [windowhw]

if [[ "$#" -lt 4 ]]
    then
    echo $(basename $0) 1>&2
    echo "Run extractSignal in batch on multiple pairs of peak files and signal files" 1>&2
    echo "USAGE:" 1>&2
    echo "$(basename $0) [pairFile] [peakDir] [signalDir] [oDir] [format:OPTIONAL]" 1>&2
    echo " [pairFile]: Col1: prefix for peakFile [tab] Col2: prefix for signal file" 1>&2
    echo " [peakDir]: Directory containing peak files" 1>&2
    echo " [signalDir]: Directory containing signal mat files" 1>&2
    echo " [oDir]: Output directory for extracted signal files" 1>&2
    echo " [format]: OPTIONAL: cagt or mat (default: mat)" 1>&2
    echo " [windowhw]: OPTIONAL: window half size in bp (default:1250)" 1>&2
    exit 1 
fi

PAIRFILE=$1
if [[ ! -e ${PAIRFILE} ]]; then echo "PairFile: ${PAIRFILE} does not exist" 1>&2 ; exit 1 ; fi

PDIR=$2
if [[ ! -d ${PDIR} ]]; then echo "Peak File directory: ${PDIR} does not exist" 1>&2 ; exit 1 ; fi

SDIR=$3
if [[ ! -d ${SDIR} ]]; then echo "Signal data directory: ${SDIR} does not exist" 1>&2 ; exit 1 ; fi

ODIR=$4
[[ ! -d ${ODIR} ]] && mkdir ${ODIR}

FORMAT='mat'
[[ $# -gt 4 ]] && FORMAT=$5

windowHw=1250
[[ $# -gt 5 ]] && windowHw=$6

JOBGROUP="/es_${RANDOM}"

while read i
  do
  peakPrefix=$(echo $i | awk '{print $1}')
  sigPrefix=$(echo $i | awk '{print $2}')
  peakTag=$(echo $i | awk '{print $5}')
  signalTag=$(echo $i | awk '{print $6}')
  cellLine=$(echo $i | awk '{print $3}')

  peakFile=$(find ${PDIR} -name "${peakPrefix}*.gz")
  sigFile=$(find ${SDIR} -name "${sigPrefix}*.mat")
  
  if [[ ! -e ${peakFile} || ! -e ${sigFile} ]]
      then
      printf "Skipping [%s]\t[%s]\n" ${peakPrefix} ${sigPrefix} 1>&2
      continue;
  fi

  # Make output directory (peakPrefix) if it doesnt exist
  [[ ! -d "${ODIR}/${peakPrefix}" ]] && mkdir "${ODIR}/${peakPrefix}"

  # Create output file name
  outFile="${ODIR}/${peakPrefix}/${peakPrefix}_AT_${sigPrefix}.${FORMAT}"
  [[ -e ${outFile} || -e "${outFile}.gz" ]] && continue # Skip if output file exists

  # Write submit script
  submitScript="es_${RANDOM}.sh"
  echo '#!/bin/bash' > ${submitScript}
  echo 'TMPDIR="${TMP}/es_${RANDOM}${RANDOM}"' >> ${submitScript}
  echo '[[ ! -d ${TMPDIR} ]] && mkdir $TMPDIR' >> ${submitScript}
  echo 'export TMP=${TMPDIR}' >> ${submitScript}
  echo 'export MCR_CACHE_ROOT=$TMPDIR' >> ${submitScript}
  echo "extractSignal -i=${peakFile} -t=${sigFile} -us=peak -sl=${windowHw} -sr=${windowHw} -fw=true -o=${outFile} -of=${FORMAT} -ov=signal -mf=samplerate -mp=10" >> ${submitScript} # 2500 bp windows with 10 bp sampling
  if [[ ${FORMAT} == 'cagt' ]]
      then
      echo "[[ -e ${outFile} ]] && gzip ${outFile}" >> ${submitScript}
  fi
  echo 'rm -rf ${TMPDIR}' >> ${submitScript}
  chmod 755 ${submitScript}

  # Submit script
  bsub -J "${peakPrefix}_AT_${sigPrefix}" -g "${JOBGROUP}" -W 24:00 -o "${outFile}.out" -e "${outFile}.err" < ${submitScript}
  rm ${submitScript}
  sleep 1s

done < "${PAIRFILE}"