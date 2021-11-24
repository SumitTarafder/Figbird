#!/bin/bash
set -e
STARTTIME="$(date -u +%s)"

chmod a+x jq

#==================Parse JSON config. file====================
echo "Parsing JSON config. script..."

if [ $# -eq 0 ]; then
    echo "ERROR!! Please input a configuration file in JSON format."
    exit 1
elif [ $# -gt 1 ]; then
    echo "Please input only one argument."
    exit 1
fi

config_file=$1
cmd="./jq"
Key1="Directory"
Key2="Read_Pairs"
Key3="Parameters"

max_allowed_isz=5000
max_allowed_readlen=200
frac=1.15

read_pair_path=()
isz=()
reverse=()
max_read_len=()
serial_number=()
itr_partial=()
itr_unmapped=()
order=()

genome_path=$($cmd -r .$Key1.Draft_Genome $config_file)
BT2_path=$($cmd -r .$Key1.Bowtie2 $config_file)
Folder_path=$($cmd -r .$Key1.Output_Folder $config_file)
ref_genome_path=$($cmd -r .$Key1.Reference_Genome $config_file)

num_threads=$($cmd .$Key3.numthreads $config_file)
default_setting=$($cmd .$Key3.default $config_file)
eval=$($cmd .$Key3.evaluation $config_file)
neg_ovlap=$($cmd .$Key3.gaplen_negative_overlap $config_file)
trim=$($cmd .$Key3.trim_len $config_file)

Total_read_lib=$($cmd .$Key2' | length' $config_file)

if [[ $num_threads -lt 0 || $num_threads -gt 64 ]]; then
  echo "Number of threads switched to default value(4),acceptable range 1-64"
  num_threads=4
fi


if [[ $Total_read_lib = 0 ]]; then
  echo "Please enter atleast one valid read pair library in configuration file."
  exit 1
fi

for (( i=0; i<$Total_read_lib; i++ ))
do
    numarg=$($cmd .$Key2[$i]' | length' $config_file)
    if [[ $numarg -ne 9 ]]; then
        echo "Please enter valid information of read pair library in configuration file. Entered $numarg rows...Expected 9"
        exit 1
    fi
done

for (( i=0; i<$Total_read_lib; i++ ))
do
    read_pair_path[${#read_pair_path[@]}]=$($cmd -r .$Key2[$i].path_1 $config_file)
    read_pair_path[${#read_pair_path[@]}]=$($cmd -r .$Key2[$i].path_2 $config_file)

    isz[${#isz[@]}]=$($cmd .$Key2[$i].avg_insert_size $config_file)
    reverse[${#reverse[@]}]=$($cmd .$Key2[$i].is_reverse $config_file)
    max_read_len[${#max_read_len[@]}]=$($cmd .$Key2[$i].max_read_length $config_file)
    serial_number[${#serial_number[@]}]=$($cmd .$Key2[$i].serial_num $config_file)
    itr_partial[${#itr_partial[@]}]=$($cmd .$Key2[$i].num_itr_partial $config_file)
    itr_unmapped[${#itr_unmapped[@]}]=$($cmd .$Key2[$i].num_itr_unmapped $config_file)

    order[${#order[@]}]=$($cmd -r .$Key2[$i].order[0] $config_file)
    order[${#order[@]}]=$($cmd -r .$Key2[$i].order[1] $config_file)
done

: '
for i in "${itr_partial[@]}"
do
	echo $i

done
'

for (( i=0; i<$Total_read_lib; i++ ))
do
    if [ ${isz[$i]} -le 0 || ${isz[$i]} -gt $max_allowed_isz ]; then
      echo "Please enter a valid insert size(<=$max_allowed_isz) in configuration file."
      exit 1
    fi

    if [ ${max_read_len[$i]} -le 0 || ${max_read_len[$i]} -gt $max_allowed_readlen ]; then
      echo "Please enter a valid read length(<=$max_allowed_readlen) in configuration file."
      exit 1
    fi

    if [ ${order[$((2*$i))]} = ${order[$((2*$i+1))]} ]; then
      echo "Please enter different types in order object of read pairs in configuration file."
      exit 1
    fi
done

printf "Done"
#=======================Detect type of read pairs=============

readlibtype=()
min_isz=100000
min_isz_index=-1

for (( i=0; i<$Total_read_lib; i++ ))
do
    if [[ ${isz[$i]} -le 250 ]]; then
        readlibtype[${#readlibtype[@]}]=0
    else
        readlibtype[${#readlibtype[@]}]=1
    fi

    if [[ ${isz[$i]} -lt $min_isz ]]; then
        min_isz=${isz[$i]}
        min_isz_index=$i
    fi
done

#=======================Create Directories=====================

output_path=$Folder_path"Figbird/"
gapped_reads_path=$output_path"Gaps/"
result_path=$output_path"Filled_Scaffolds/"
index_path=$output_path"Bowtie2_indexFiles/index"
sam_aln_path=$output_path"Alignments/"
gapout_path=$result_path"Individual_gaps/"
temp_file_path=$output_path"Temp/"

rm -rf $result_path
rm -rf $gapped_reads_path

mkdir -p $output_path
mkdir -p $gapped_reads_path
mkdir -p $index_path
mkdir -p $result_path
mkdir -p $sam_aln_path
mkdir -p $gapout_path
mkdir -p $temp_file_path
#========================Reverse read pairs====================
#exe_path=$temp_file_path"a.out"
exe_path="./a.out"

printf "\nReversing read pairs if necessary..."

rev_file_path=()

for (( i=0; i<$Total_read_lib; i++ ))
do
    if [[ ${reverse[$i]} = 1 ]]; then

        rev_str=$(g++ Reverse.cpp && ./a.out ${read_pair_path[$((2*$i))]} ${read_pair_path[$((2*$i+1))]})

        for word in $rev_str
        do
            rev_file_path[${#rev_file_path[@]}]=$word
        done
    else
        rev_file_path[${#rev_file_path[@]}]="empty"
        rev_file_path[${#rev_file_path[@]}]="empty"
    fi
done

printf "Done"

#====================Function Definitions======================

function run_bowtie(){

    printf "\n=========================\n"
    echo "Iteration $count Starts"
    echo "========================="

    readtype=$1
    runtype=$2
    iteration=$3
    readlibindex=$4

    maxreadlen=${max_read_len[$readlibindex]}
    m2=${isz[$readlibindex]}
    maxD2=`echo $m2 \* $frac |bc`

    m1=${isz[$min_isz_index]}

    maxD1=$m1
    
    if [[ ${reverse[$readlibindex]} = 1 ]]; then
        read_dir1=${rev_file_path[$((2*$readlibindex))]}
        read_dir2=${rev_file_path[$((2*$readlibindex+1))]}
    else
        read_dir1=${read_pair_path[$((2*$readlibindex))]}
        read_dir2=${read_pair_path[$((2*$readlibindex+1))]}
    fi

    if [ $readtype = 0 -a $runtype = 0 ]; then
        maxdistance=$maxD1
        partial_flag=1
        unmapped=0
    elif [ $readtype = 1 -a $runtype = 1 ]; then
        maxdistance=$maxD2
        partial_flag=0
        unmapped=1
    fi

    genome_red=0
    read_red=0

    if [ $iteration = 1 ]; then
        sam1="result1.sam"
        sam2="result2.sam"

        if [[ ${reverse[$min_isz_index]} = 1 ]]; then
            read_dir3=${rev_file_path[$((2*$min_isz_index))]}
            read_dir4=${rev_file_path[$((2*$min_isz_index+1))]}
        else
            read_dir3=${read_pair_path[$((2*$min_isz_index))]}
            read_dir4=${read_pair_path[$((2*$min_isz_index+1))]}
        fi
    else
        sam1="result1.sam"
        sam2="result2.sam"

        read_dir3=$reduced_path0
        read_dir4=$reduced_path1
    fi

    partial_readlen=${max_read_len[$min_isz_index]}

    if [[ $iteration = 1 ]]; then
        read_red=1
    fi

    sam_dir1=$sam_aln_path$sam1
    sam_dir2=$sam_aln_path$sam2

    myoutfilename=$sam_aln_path"myout.sam"
    newmyoutfilename=$sam_aln_path"myout_temp.sam"
    
    if [ $iteration = 1 -a $trim -gt 0 ]; then
        g++ FlankTrim.cpp && ./a.out $gapped_genome $trim $maxreadlen $trimmed_genome
        gapped_genome=$trimmed_genome
    fi

    #This is for the initialization purpose
    if [[ $runtype = 0 || $iteration = 1 ]]; then
        genome=$gapped_genome
    else
        g++ Reduce_SCF.cpp && ./a.out $gapped_genome $temp_file_path
        genome=$reduced_genome

    fi

    #echo $genome
    #echo $read_dir3
    #echo $read_dir4


    printf "Building Index with bowtie2..."
    "$BT2_path"bowtie2-build -q $genome $index_path
    printf "Done"
    printf "\nAligning reads to scaffolds..."
    "$BT2_path"bowtie2 -p $num_threads --local -x $index_path -X $maxD1 -1 $read_dir3 -2 $read_dir4 -S $sam_dir1
    printf "Done\n"

    printf "Preprocessing Reads..."

    reduced_str=$(g++ Preprocess.cpp && ./a.out $genome $maxD1 1 $sam_dir1 $myoutfilename $gapped_genome $read_dir3 $read_dir4 $gapped_reads_path $temp_file_path $default_setting $genome_red $read_red)

    if [ $iteration = 1 ]; then
        cp $myoutfilename $newmyoutfilename

        y=0
        for word in $reduced_str
        do
            if [ $y = 0 ]; then
                reduced_path0=$word
            else
                reduced_path1=$word
            fi
            y=1
        done
    fi

    printf "Done\n"

    if [ $runtype = 0 ]; then
        if [ $iteration -gt 1 ]; then
            myoutfilename=$newmyoutfilename
            echo ""
        fi
    fi

    genome_red=0
    read_red=0

    if [ $runtype = 1 ]; then
       reduction=-1
       if [ $readtype = 1 ]; then
            genome=$gapped_genome
            genome_red=0
       else
            g++ Reduce_SCF.cpp && ./a.out $gapped_genome $temp_file_path
            genome=$reduced_genome
            genome_red=1
       fi
       
        #echo $genome
        #echo $read_dir1
        #echo $read_dir2

        printf "Building Index with bowtie2..."
        "$BT2_path"bowtie2-build -q $genome $index_path
        printf "Done"
        printf "\nAligning reads to scaffolds..."
        "$BT2_path"bowtie2 -p $num_threads -x $index_path -X $maxdistance -1 $read_dir1 -2 $read_dir2 -S $sam_dir2
        printf "Done\n"

       printf "Preprocessing Reads..."

       reduced_str=$(g++ Preprocess.cpp && ./a.out $genome $maxdistance 2 $sam_dir2 $myoutfilename $gapped_genome $read_dir1 $read_dir2 $gapped_reads_path $temp_file_path $default_setting $genome_red $read_red)

       printf "Done\n"

       if [ $readtype = 0 ]; then
           myoutlinecount=$(wc -l < $myoutfilename)
           if [ $myoutlinecount -lt 1000 ]; then
               myoutfilename=$newmyoutfilename
           fi
       fi
    fi

    echo "Gap filling starts..."

    g++ -std=c++11 -pthread FillGaps.cpp && ./a.out $gapped_genome $maxdistance $maxreadlen $count $partial_flag $unmapped $num_threads $myoutfilename $temp_file_path $gapped_reads_path $neg_ovlap $partial_readlen $trim

    contigfile=$result_path$count$filled_filename
    gapfile=$gapout_path$"gapout_"$count".txt"
    drawfile=$gapout_path$"alignment_"$count".txt"

    mv $temp_file_path$filled_filename $contigfile
    mv $temp_file_path$"gapout.txt" $gapfile
    mv $temp_file_path$"draw.txt" $drawfile
}

function run_bowtie_user(){

    r_type=$1
    rp_index=$2
    iteration=$count

    printf "\n=========================\n"
    echo "Iteration $count Starts"
    echo "============================="

    if [[ ${reverse[$rp_index]} = 1 ]]; then
        r_dir1=${rev_file_path[$((2*$rp_index))]}
        r_dir2=${rev_file_path[$((2*$rp_index+1))]}
    else
        r_dir1=${read_pair_path[$((2*$rp_index))]}
        r_dir2=${read_pair_path[$((2*$rp_index+1))]}
    fi

    mreadlen=${max_read_len[$rp_index]}
    m2=${isz[$rp_index]}
    mdistance2=`echo $m2 \* $frac |bc`

    #initialization_index=$rp_index
    initialization_index=$min_isz_index

    if [ $r_type = 1 ]; then
        p_flag=1
        u_flag=0
    else
        p_flag=0
        u_flag=1
    fi

    sam1="result1.sam"
    sam2="result2.sam"

    sam_dir1=$sam_aln_path$sam1
    sam_dir2=$sam_aln_path$sam2

    myoutfname=$sam_aln_path"myout.sam"
    newmyoutfname=$sam_aln_path"myout_temp.sam"

    #This is for the initialization purpose
    genome=$gapped_genome

    if [ $r_type = 0 ]; then
        if [[ ${reverse[$initialization_index]} = 1 ]]; then
            r_dir3=${rev_file_path[$((2*$initialization_index))]}
            r_dir4=${rev_file_path[$((2*$initialization_index+1))]}
        else
            r_dir3=${read_pair_path[$((2*$initialization_index))]}
            r_dir4=${read_pair_path[$((2*$initialization_index+1))]}
        fi
        m1=${isz[$initialization_index]}
        partial_readlen=${max_read_len[$initialization_index]}
    else
        r_dir3=$r_dir1
        r_dir4=$r_dir2
        m1=${isz[$rp_index]}
        partial_readlen=${max_read_len[$rp_index]}
    fi

    mdistance1=$m1

    maxdistance_fig=$mdistance1

    genome_red=0
    read_red=0
    
    if [ $iteration = 1 -a $trim -gt 0 ]; then
        g++ FlankTrim.cpp && ./a.out $gapped_genome $trim $mreadlen $trimmed_genome
        genome=$trimmed_genome
    fi

    #echo "running bowtie partially with $r_dir3 on $genome with $mdistance1"

    printf "Building Index with bowtie2..."
    "$BT2_path"bowtie2-build -q $genome $index_path
    printf "Done"
    printf "\nAligning reads to scaffolds..."
    "$BT2_path"bowtie2 -p $num_threads --local -x $index_path -X $mdistance1 -1 $r_dir3 -2 $r_dir4 -S $sam_dir1
    printf "Done\n"

    printf "Preprocessing Reads..."

    reduced_str=$(g++ Preprocess.cpp && ./a.out $genome $mdistance1 1 $sam_dir1 $myoutfname $gapped_genome $r_dir3 $r_dir4 $gapped_reads_path $temp_file_path $default_setting $genome_red $read_red)

    printf "Done\n"

    if [ $r_type = 0 ]; then

        genome_red=0
        read_red=0

        printf "Building Index with bowtie2..."
        "$BT2_path"bowtie2-build -q $genome $index_path
        printf "Done"
        printf "\nAligning reads to scaffolds..."
        "$BT2_path"bowtie2 -p $num_threads -x $index_path -X $mdistance2 -1 $r_dir1 -2 $r_dir2 -S $sam_dir2
        printf "Done\n"


        #echo "running bowtie unmapped with $r_dir1 on $genome with $mdistance2"

        printf "Preprocessing Reads..."

        reduced_str=$(g++ Preprocess.cpp && ./a.out $genome $mdistance2 2 $sam_dir2 $myoutfname $genome $r_dir1 $r_dir2 $gapped_reads_path $temp_file_path $default_setting $genome_red $read_red)

        maxdistance_fig=$mdistance2

        printf "Done\n"

    fi

    g++ -std=c++11 -pthread FillGaps.cpp && ./a.out $genome $maxdistance_fig $mreadlen $count $p_flag $u_flag $num_threads $myoutfname  $temp_file_path $gapped_reads_path $neg_ovlap $partial_readlen $trim

    contigfile=$result_path$count$filled_filename
    gapfile=$gapout_path$"gapout_"$count".txt"
    drawfile=$gapout_path$"alignment_"$count".txt"

    mv $temp_file_path$filled_filename $contigfile
    mv $temp_file_path$"gapout.txt" $gapfile
    mv $temp_file_path$"draw.txt" $drawfile
}

function prepare(){

    contigfile=$result_path$count$filled_filename
    gapped_genome=$contigfile
    count=$((count+1))
}

function runFrag(){

    loop_counter=$1
    readlibindex=$2

    for (( i=1; i<=$loop_counter; i++ ))
    do
        fname1=$result_path$count$filled_filename
        count1=$(fgrep -o N $fname1 | wc -l)
        g++ Check_Flank.cpp && ./a.out $fname1 $temp_file_path
        prepare
        run_bowtie 0 1 $count $readlibindex

        fname2=$result_path$count$filled_filename

        count2=$(fgrep -o N $fname2 | wc -l)
        if [ $count1 = $count2 ]; then
            break
        fi

    done
}

#=======================Start Gap Filling======================
count=1
filled_filename="filledContigs.fa"
reduced_genome=$temp_file_path"newgenome.fa"
trimmed_genome=$temp_file_path"newgenome1.fa"
gapped_genome=$genome_path

#Readtype   #Runtype
#   0          0      = frag + Partial
#   0          1      = frag + Unmapped
#   1          1      = Jump + Unmapped

if [ $default_setting -eq 1 ]; then
    if [[ $Total_read_lib = 2 && ${readlibtype[0]} -ne ${readlibtype[1]} ]]; then

        j1=1
        j2=2
        f1=1
        j3=1
        f2=3

        #Detect which one is frag and which one is jump, the given order mayn't be correct
        if [[ ${readlibtype[0]} = 0 ]]; then
            frag_index=0
            jump_index=1
        else
            frag_index=1
            jump_index=0
        fi

        for (( i=1; i<=$j1; i++ ))
        do
            run_bowtie 0 0 $count $frag_index
        done
        
        for (( i=1; i<=$j2; i++ ))
        do
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
            prepare
            run_bowtie 1 1 $count $jump_index
        done

        for (( i=1; i<=$f1; i++ ))
        do
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
            prepare
            run_bowtie 0 0 $count $frag_index
        done

        for (( i=1; i<=$j3; i++ ))
        do
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
            prepare
            run_bowtie 1 1 $count $jump_index
        done

        for (( i=1; i<=$f2; i++ ))
        do
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
            prepare
            run_bowtie 0 0 $count $frag_index
        done
    else
    #zxcf
        rl1=0
        rl2=0
        
        itr_f=1
        
        for (( i=0; i<$Total_read_lib; i++ ))
        do
            if [[ ${isz[$i]} -le 250 ]]; then
                type1[${#type1[@]}]=$i
                rl1=$(($rl1+1))
            else
                type2[${#type2[@]}]=$i
                rl2=$(($rl2+1))
            fi
        done

        if [[ $rl1 = 0 ]]; then
            for i in "${type2[@]}"
            do
            	type1[${#type1[@]}]=$i
                rl1=$(($rl1+1))
            done
        fi

        j1=2
        f1=2
        j2=1
        f2=3

        for (( i=0; i<$rl2 ; i++ ))
        do
            for (( k=0; k<$j1 ; k++ ))
            do
                if [ $itr_f -gt 1 ]; then
                    prepare
                fi
                rlibindex=${type2[$i]}
                run_bowtie_user 0 $rlibindex $count
                itr_f=$((itr_f+1))
                
                fillmore=$(cat $temp_file_path"Ncount.txt")
                echo $fillmore     
                if [ $fillmore -eq 0 ]; then
                    break
                fi
            done
            
            if [ $fillmore -eq 0 ]; then
                break
            fi
        done

        for (( i=0; i<$rl1 ; i++ ))
        do
            for (( k=0; k<$f1 ; k++ ))
            do
                fillmore=$(cat $temp_file_path"Ncount.txt")     
                if [ $fillmore -eq 0 ]; then
                    break
                fi
                
                if [ $itr_f -gt 1 ]; then
                    prepare
                fi
                rlibindex=${type1[$i]}
                run_bowtie_user 1 $rlibindex $count
                itr_f=$((itr_f+1))
            done
            
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
        done

        for (( i=0; i<$rl2 ; i++ ))
        do
            for (( k=0; k<$j2 ; k++ ))
            do
                fillmore=$(cat $temp_file_path"Ncount.txt")     
                if [ $fillmore -eq 0 ]; then
                    break
                fi
                
                if [ $itr_f -gt 1 ]; then
                    prepare
                fi
                rlibindex=${type2[$i]}
                run_bowtie_user 0 $rlibindex $count
                itr_f=$((itr_f+1))
                
            done
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
        done

        for (( i=0; i<$rl1 ; i++ ))
        do
            for (( k=0; k<$f2 ; k++ ))
            do
                fillmore=$(cat $temp_file_path"Ncount.txt")     
                if [ $fillmore -eq 0 ]; then
                    break
                fi
                
                if [ $itr_f -gt 1 ]; then
                    prepare
                fi
                rlibindex=${type1[$i]}
                run_bowtie_user 1 $rlibindex $count
                itr_f=$((itr_f+1))
                
            done
            
            fillmore=$(cat $temp_file_path"Ncount.txt")     
            if [ $fillmore -eq 0 ]; then
                break
            fi
        done
    fi    
else
   itr_f=1
   
   for (( i=0; i<$Total_read_lib; i++ ))
   do
        for (( j=0; j<$Total_read_lib; j++ ))
        do
            if [ $((j+1)) = ${serial_number[$i]} ]; then
                rlibindex=$j
                break
            fi
        done

        for (( j=0; j<2; j++ ))
        do
            rtype=${order[$((2*$rlibindex+$j))]}

            if [ $rtype == p ]; then
                pu=1
            else
                pu=0
            fi

            if [ $pu = 1 ]; then
                num_itr=${itr_partial[$rlibindex]}
            else
                num_itr=${itr_unmapped[$rlibindex]}
            fi

            #printf "\nlib serial = $rlibindex, rtype = $rtype pu = $pu itr= $num_itr"

            for (( k=0; k<$num_itr ; k++ ))
            do
                if [ $itr_f -gt 1 ]; then
                    prepare
                fi
                
                run_bowtie_user $pu $rlibindex $count
                itr_f=$((itr_f+1))
                
                fillmore=$(cat $temp_file_path"Ncount.txt")     
                if [ $fillmore -eq 0 ]; then
                    break
                fi
            done
            
            if [ $fillmore -eq 0 ]; then
                break
            fi
        done
        fillmore=$(cat $temp_file_path"Ncount.txt")     
        if [ $fillmore -eq 0 ]; then
            break
        fi
   done
fi

g++ CombineGaps.cpp && ./a.out $count $gapout_path

finalfilledfile=$result_path$count$filled_filename
finalname="FilledScaffolds_fianl.fa"
cp $finalfilledfile $result_path$finalname

ENDTIME="$(date -u +%s)"
echo "Figbird ends successfully."
echo "Total time taken =  $(($ENDTIME - $STARTTIME)) seconds."

echo "=========================================================================="

#=======================Remove extra files======================


#======================Evaluation===============================

if [ $eval = 1 ]; then

    if test -f "$ref_genome_path"; then

        quast_output=$output_path"QUAST_Results/"
        mkdir -p $quast_output

        new_file=$result_path$count$filled_filename
        echo "Evaluating $new_file using QUAST"

        echo "Results" > $quast_output"Result.txt"

        EM_dir=$temp_file_path$count$filled_filename
        Q="quast-2.3/quast.py"

        python reference.py $new_file $EM_dir

        output_dir=$temp_file_path"output"

        python2 $Q --strict-NA -R $ref_genome_path -o $output_dir --no-plots $EM_dir

        # QUAST quast_correction

        result=$(python2 correct_quast.py --N 4000 $output_dir/contigs_reports/contigs_report_*.stdout  $output_dir/contigs_reports/misassemblies_report.txt  $output_dir/report.txt $EM_dir)
        # Saving the result in a out.txt file the oder of the result is   Misassmblies, Erroneous-length, Unaligned-length, NGA50, Number of gaps, Total gap length

        echo $result >> $quast_output"Result.txt"

        cp -r $output_dir $quast_output
        # Remember to delete the output folder otherwise quast will skip next time.
        rm -rf $output_dir
        rm -f $EM_dir

    else
        echo "Reference genome could not be found."
    fi
fi
