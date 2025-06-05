#!/bin/bash
freyja update
freyja update --outdir . #to ensure files end up in this dir

# collect samples from last 240 days for rerunning
awk -v dat="$(date -d '240 days ago' '+%Y/%m/%d')" -F ',' '$4 > dat {print $5}' ../sample_metadata.csv > recent_samples.txt

my_func() {
    fn=$1
    depthfolder='../depths/'
    output='../outputs/'

    fn_out=${fn##*/}
    baseName=${fn##*/}
    baseName=${baseName%.*}
    depthfile0="${depthfolder}${fn_out%.tsv}.depths"
    output0="${output}${fn_out%.*}.demix.tsv"
    echo $fn
    echo $depthfile0
    echo $output0
    freyja demix $fn $depthfile0 --output $output0 --eps 0.0000001
    # if [ -f "$output0" ]; then
    #     echo "file already exists. skipping. "
    # else
    #     # echo "file doesn't exist. running. "
    #     freyja demix $fn $depthfile0 --output $output0 --eps 0.0000001
    # fi

}

export -f my_func
python get_samples_to_run.py | parallel -j 20 my_func 

freyja aggregate ../outputs/ --output ../agg_demixed.tsv 

python make_plot.py
python make_catchment_plots.py
