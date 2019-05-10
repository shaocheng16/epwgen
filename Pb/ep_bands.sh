
#!/bin/bash
#Execute from the parent directory (containing ELB, PHB, EPM, ISO).

#Specify the starting and ending point of calculation. E.g. if the phonon calculations stopped because the time limit...
#...was too short, set calc_start=6 and run the submisison script again. The calculations will continue from where they stopped.
#1: electron bands scf calculation
#2: electron bands calculation
#3: electron bands interpolation
#4: phonon bands scf calculation
#5: phonon calculation initialization
#6: phonon calculations for every irreducible q-point/representation
#(If the phonon calculations are split, #6 will set up a manager that automatically splits the phonon calculation...
#... You can kill the manager while it's making the submissions (e.g. if you notice that it generates too many jobs) and...
#...resubmit starting from #6 where it will continue from where it stopped. If you restart, restart from #6 not from #5.)
#7: phonon calculation finalization
#8: phonon q2r run
#9: phonon matdyn run
#(#10: tidy up - deletes everything but the input/output files, the bands (i.e. electron and phonon bands) and the data...
#      ...needed for EPM calculations. This is not necessary but saves space. Usable only as a ...
#      ...seperate (calc_start = calc_end) execution. Only use once you are sure that everything worked.)
calc_start=1
calc_end=9

#If you need a specific submission script (e.g. for the q2r run) set calc_start = calc_end (e.g. = 8) and set no_sub
#to true. You'll find the job script as job.sh in the corresponding directory.
no_sub=false

#load modules 
module purge 
module load quantum_espresso/6.2.1 &>/dev/null 
module load python/3.6.0 &>/dev/null 
module load gnuplot &>/dev/null 



#________________________________EXECUTE_THE_SCRIPT________________________________#
#__________________________________________________________________________________#
split_q=true
split_irr=true

#__________PREPARE_INPUT___________#
base_dir=$(pwd)
cp $base_dir/ELB/scf.in $base_dir/PHB

#__________JOB_SUBMISSION___________#
#electron bands calculation
cd $base_dir/ELB

if (($calc_start < 10)) && (($calc_end == 10))
then
echo "Error: please only tidy up once you are sure that all calculations have finished correctly"
exit
fi

if (($calc_start > $calc_end))
then
echo "Error: calculation order makes no sense"
exit
fi

#scf
if (($calc_start <= 1)) && (($calc_end >= 1))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_elb_scf.o
#BSUB -e Pb_elb_scf.e
#BSUB -J Pb_elb_scf
mpirun pw.x -npool 4 -in scf.in > scf.out
#get the number of bands from the scf.out file and put them in the bands.in file
nbnd=\$(grep nbnd bands.in | awk '{print \$3}')
if [ \$nbnd -eq 0 ]
then
    nbnd=\$(grep -m1 "number of Kohn-Sham states" scf.out | awk '{print \$5}')
    line=\$(grep -n nbnd bands.in | cut -d : -f 1)
    sed -i "\${line}s/[^ ]*[^ ]/\${nbnd}/3" bands.in
fi
EOF

if ! $no_sub
then
bsub<job.sh
fi

fi

#bands
if (($calc_start <= 2)) && (($calc_end >= 2))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_elb.o
#BSUB -e Pb_elb.e
#BSUB -J Pb_elb
#BSUB -w Pb_elb_scf
mpirun pw.x -npool 4 -in bands.in > bands.out
EOF

if (($calc_start == 2))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi

#bands_ip
if (($calc_start <= 3)) && (($calc_end >= 3))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_elb_ip.o
#BSUB -e Pb_elb_ip.e
#BSUB -J Pb_elb_ip
#BSUB -w Pb_elb
mpirun bands.x -npool 4 -in bands_ip.in > bands_ip.out
EOF

if (($calc_start == 3))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi

#phonon calculation
cd $base_dir/PHB

#scf
if (($calc_start <= 4)) && (($calc_end >= 4))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_ph_scf.o
#BSUB -e Pb_ph_scf.e
#BSUB -J Pb_ph_scf
mpirun pw.x -npool 4 -in scf.in > scf.out
EOF

if ! $no_sub
then
bsub<job.sh
fi

fi

#phonon calculations
#initialize the ph calculations
if (($calc_start <= 5)) && (($calc_end >= 5))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_ph_init.o
#BSUB -e Pb_ph_init.e
#BSUB -J Pb_ph_init
#BSUB -w Pb_ph_scf

cp ph.in ph_start.in
echo "    start_irr = 0" >> ph_start.in
echo "    last_irr  = 0" >> ph_start.in
echo "    /" >> ph_start.in
mpirun ph.x -npool 1 -in ph_start.in > ph_start.out
EOF

if (($calc_start == 5))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi

#set up a manager that splits the phonon calculation
#IF CALC_START
if (($calc_start <= 6)) && (($calc_end >= 6))
then

#if no parallelization is wanted start the phonon calculation normally
#IF NOSPLIT
if ! $split_q
then
cp ph.in ph_all.in
echo "/" >> ph_all.in

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_ph_all.o
#BSUB -e Pb_ph_all.e
#BSUB -J Pb_ph_all
#BSUB -w Pb_ph_init
mpirun ph.x -npool 4 -in ph_all.in > ph_all.out
EOF
#ENDIF NOSPLIT

#IF SPLIT_Q
#if irreducible q-point parallelization is enabled
elif $split_q && ! $split_irr
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_ph_manager.o
#BSUB -e Pb_ph_manager.e
#BSUB -J Pb_ph_manager
#BSUB -w Pb_ph_init

#get the number of irreducible q-points
irr_qs=\$(sed "2q;d" Pb.dyn0 | awk '{print $1}')

#run a phonon calculation for every irreducible q-point
for ((q=1; q <= \$irr_qs; q++))
do

#make directories for the irreducible q-point and copy necessary stuff there
if [ ! -d  q\${q}/_ph0/Pb.phsave ]
then
mkdir -p q\${q}/_ph0/Pb.phsave
cp -r Pb.* q\${q}
cp -r _ph0/Pb.phsave/* q\${q}/_ph0/Pb.phsave
fi

cd q\${q}
#prepare the input file
cp ../ph.in ph_q\${q}.in
echo "    start_q = \$q" >> ph_q\${q}.in
echo "    last_q = \$q" >> ph_q\${q}.in
echo "/" >>  ph_q\${q}.in

#make the job file
cat > job_temp.sh << EOF1
\#!/bin/bash
\#BSUB -n 4
\#BSUB -R "rusage[mem=1024]"
\#BSUB -W 1:00
\#BSUB -o Pb_ph_q\${q}.o
\#BSUB -e Pb_ph_q\${q}.e
\#BSUB -J Pb_ph_q\${q}
mpirun ph.x -npool 4 -in ph_q\${q}.in > ph_q\${q}.out
\#delete the large wave functions in the q-directory immediately after a calculation is done
if [ -f ph_q\${q}.out ] && [ \${q} -ne 1 ]
then
status=\\\$(tail -n2 ph_q\${q}.out | head -n1 | awk '{print \\\$2}')
if [ "\\\$status" == "DONE." ]
then 
rm _ph0/Pb.q_\${q}/*.wfc*
rm Pb.save/wfc**.dat
rm *.wfc*
fi
fi
EOF1

#remove escapes
sed -i 's/\\\#/#/g' job_temp.sh

#in the case of a restart only submit the job if it's not finished yet
if [ -f ph_q\${q}.out ]
then
status=\$(tail -n2 ph_q\${q}.out | head -n1 | awk '{print \$2}')
if [ ! "\$status" == "DONE." ]
then
#LSF
bsub<job_temp.sh
fi
else
#LSF
bsub<job_temp.sh
fi
cd ..
done

#update the conditions once all jobs are submitted
collect_id=\$(bjobs -J Pb_ph_collect | tail -n1 | awk '{print \$1}')
#LSF
bmod -w "Pb_ph_q*" \$collect_id

EOF
#ENDIF SPLIT_Q

#IF SPLIT_IRR
#if irreducible representation parallelization is enabled
elif $split_irr
then
cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_ph_manager.o
#BSUB -e Pb_ph_manager.e
#BSUB -J Pb_ph_manager
#BSUB -w Pb_ph_init

#get the number of irreducible q-points
irr_qs=\$(sed "2q;d" Pb.dyn0 | awk '{print $1}')

#get the number of irreducible representations of each q-point and save them in an array
declare -a irreps
for ((i=0; i < irr_qs; i++))
do
    q=\$((i+1))
    irreps_el=\$(grep -A1 "<NUMBER_IRR_REP" _ph0/Pb.phsave/patterns.\${q}.xml | tail -n1)
    irreps[\$i]=\$irreps_el
done
    
#run a phonon calculation for every irreducible representation of every irreducible q-point 
for ((q=1; q <= irr_qs; q++))
do
i=\$((q-1))
for ((r=1; r <= irreps[i]; r++))
do

#make directories for the irreducible representation and copy necessary stuff there
if [ ! -d  q\${q}_r\${r}/_ph0/Pb.phsave ]
then
mkdir -p q\${q}_r\${r}/_ph0/Pb.phsave
cp -r Pb.* q\${q}_r\${r}
cp -r _ph0/Pb.phsave/* q\${q}_r\${r}/_ph0/Pb.phsave
fi

cd q\${q}_r\${r}
#prepare the input file
cp ../ph.in ph_q\${q}_r\${r}.in

echo "    start_q = \$q" >> ph_q\${q}_r\${r}.in
echo "    last_q = \$q" >> ph_q\${q}_r\${r}.in
echo "    start_irr = \$r" >> ph_q\${q}_r\${r}.in
echo "    last_irr = \$r" >> ph_q\${q}_r\${r}.in
echo "    /" >>  ph_q\${q}_r\${r}.in

#make the job file
cat > job_temp.sh << EOF1
\#!/bin/bash
\#BSUB -n 4
\#BSUB -R "rusage[mem=1024]"
\#BSUB -W 1:00
\#BSUB -o Pb_ph_q\${q}_r\${r}.o
\#BSUB -e Pb_ph_q\${q}_r\${r}.e
\#BSUB -J Pb_ph_q\${q}_r\${r}
mpirun ph.x -npool 4 -in ph_q\${q}_r\${r}.in > ph_q\${q}_r\${r}.out
\#delete the wave functions immediately after the calculation is done
if [ -f ph_q\${q}_r\${r}.out ] && [ \${q} -ne 1 ]
then
status=\\\$(tail -n2 ph_q\${q}_r\${r}.out | head -n1 | awk '{print \\\$2}')
if [ "\\\$status" == "DONE." ]
then 
rm _ph0/Pb.q_\${q}/*.wfc*
rm Pb.save/wfc**.dat
rm *.wfc*
fi
fi
EOF1

#remove escapes
sed -i 's/\\\#/#/g' job_temp.sh

#in the case of a restart only submit the job if it's not finished yet
if [ -f ph_q\${q}_r\${r}.out ]
then
status=\$(tail -n2 ph_q\${q}_r\${r}.out | head -n1 | awk '{print \$2}')
if [ ! "\$status" == "DONE." ]
then
#LSF
bsub<job_temp.sh
fi
else
#LSF
bsub<job_temp.sh
fi
cd ..
done
done

#update the conditions once all jobs are submitted
collect_id=\$(bjobs -J Pb_ph_collect | tail -n1 | awk '{print \$1}')
#LSF
bmod -w "Pb_ph_q*" \$collect_id

EOF
fi
#ENDIF SPLIT_IRR


if (($calc_start == 6))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi
#ENDIF CALC_START


#collect all results
#IF CALC_START
if (($calc_start <= 7)) && (($calc_end >= 7))
then

#do nothing if the calculations are not split 
#IF NOSPLIT
if ! $split_q
then
cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_ph_collect.o
#BSUB -e Pb_ph_collect.e
#BSUB -J Pb_ph_collect
#BSUB -w Pb_ph_manager

echo "nothing to collect"
EOF
#LSF
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i "${line}s/manager/all/1" job.sh
#ENDIF NOSPLIT

#IF SPLIT_Q
elif $split_q && ! $split_irr
then
cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_ph_collect.o
#BSUB -e Pb_ph_collect.e
#BSUB -J Pb_ph_collect
#BSUB -w Pb_ph_manager


#make sure that all phonon calculations have finished without problems
for ph in \$(ls ph_q*/ph_q*.out)
do
status=\$(tail -n2 \$ph | head -n1 | awk '{print \$2}')
if [ ! "\$status" == "DONE." ]
then
echo "Error: not all phonon calculations have finished or they ran into problems. Resubmit the phonon calculation (without initialization) before collecting."
exit
fi
done

#collect files
irr_qs=\$(sed "2q;d" Pb.dyn0 | awk '{print $1}')
for ((q=1; q<=irr_qs;q++))
do

#copy the dynmats
cp q\${q}/_ph0/Pb.phsave/dynmat.* _ph0/Pb.phsave

#copy directory containing dvscf
if ((q==1))
then
cp q1/_ph0/Pb.dvscf1 _ph0
else
cp -r q\${q}/_ph0/Pb.q_\${q} _ph0
fi
done

#prepare input file and execute
cp ph.in ph_end.in
echo "    /" >> ph_end.in
mpirun ph.x -npool 1 -in ph_end.in > ph_end.out
EOF
#ENDIF SPLIT_Q

#IF SPLIT_IRR
elif $split_irr
then
cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_ph_collect.o
#BSUB -e Pb_ph_collect.e
#BSUB -J Pb_ph_collect
#BSUB -w Pb_ph_manager


#make sure that all phonon calculations have finished without problems
for ph in \$(ls q*_r*/ph_q*_r*.out)
do
status=\$(tail -n2 \$ph | head -n1 | awk '{print \$2}')
if [ ! "\$status" == "DONE." ]
then
echo "Error: not all phonon calculations have finished or they ran into problems. Resubmit the phonon calculation (without initialization) before collecting."
exit
fi
done

#collect files
declare -a irreps
irr_qs=\$(sed "2q;d" Pb.dyn0 | awk '{print $1}')
for ((i=0; i < irr_qs; i++))
do
    q=\$((i+1))
    irreps_el=\$(grep -A1 "<NUMBER_IRR_REP" _ph0/Pb.phsave/patterns.\${q}.xml | tail -n1)
    irreps[\$i]=\$irreps_el
done

for ((q=1; q <= irr_qs; q++))
do
i=\$((q-1))
size_new=0
size_old=0
touch dvscf1_old
for ((r=1; r <= irreps[i]; r++))
do

#copy the dynmats
cp q\${q}_r\${r}/_ph0/Pb.phsave/dynmat.* _ph0/Pb.phsave

#combine the dvscf files bytewise
if ((q==1))
then
size_new=\$(ls -l q1_r\${r}/_ph0/Pb.dvscf1 | awk '{print \$5}')
_count=\$((size_new - size_old))
_skip=\$((size_new - _count))
size_old=\$size_new
dd if=q1_r\${r}/_ph0/Pb.dvscf1 of=dvscf1_temp skip=\$_skip count=\$_count iflag=skip_bytes,count_bytes
cat dvscf1_old dvscf1_temp > dvscf1_new
mv dvscf1_new dvscf1_old

else
size_new=\$(ls -l q\${q}_r\${r}/_ph0/Pb.q_\${q}/Pb.dvscf1 | awk '{print \$5}')
_count=\$((size_new - size_old))
_skip=\$((size_new - _count))
size_old=\$size_new
dd if=q\${q}_r\${r}/_ph0/Pb.q_\${q}/Pb.dvscf1 of=dvscf1_temp skip=\$_skip count=\$_count iflag=skip_bytes,count_bytes
cat dvscf1_old dvscf1_temp > dvscf1_new
mv dvscf1_new dvscf1_old
fi
done

#move the combined dvscf files to the right location
if ((q==1))
then
mv dvscf1_old _ph0/Pb.dvscf1
else
if [ ! -d _ph0/Pb.q_\${q} ]
then
mkdir _ph0/Pb.q_\${q}
fi
mv dvscf1_old _ph0/Pb.q_\${q}/Pb.dvscf1
fi
done

#prepare input file and execute
cp ph.in ph_end.in
echo "    /" >> ph_end.in
mpirun ph.x -npool 1 -in ph_end.in > ph_end.out

EOF
fi
#ENDIF SPLIT_IRR

if (($calc_start == 7))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi
#ENDIF CALC_START

#q2r
if (($calc_start <= 8)) && (($calc_end >= 8))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_q2r.o
#BSUB -e Pb_q2r.e
#BSUB -J Pb_q2r
#BSUB -w Pb_ph_collect

if [ -f Pb.dyn1.xml ]
then
    cp Pb.dyn0 Pb.dyn0.xml
    sed -i "s/Pb\.dyn/Pb\.dyn\.xml/g" q2r.in
fi
mpirun q2r.x -npool 1 -in q2r.in > q2r.out
EOF

if (($calc_start == 8))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi

#matdyn
if (($calc_start <= 9)) && (($calc_end >= 9))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_matdyn.o
#BSUB -e Pb_matdyn.e
#BSUB -J Pb_matdyn
#BSUB -w Pb_q2r

if [ -f Pb.dyn1.xml ]
then
    sed -i "s/Pb\.fc/Pb\.fc\.xml/g" matdyn.in
fi
mpirun matdyn.x -npool 1 -in matdyn.in > matdyn.out
EOF

if (($calc_start == 9))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi

#tidying up
if (($calc_start == 10)) && (($calc_end == 10))
then
cat > job.sh << EOF
#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_tidy.o
#BSUB -e Pb_tidy.e
#BSUB -J Pb_tidy

#phonon directory
cd $base_dir/PHB
rm -r q**/
rm *.xml
rm *.e
rm *.o
rm job*
rm *dvscf*
rm *wfc*
rm matdyn.modes
#electron directory
cd $base_dir/ELB
rm -r *.save/
rm *.xml
rm *.e
rm *.o
rm job*
rm *wfc*
EOF

if ! $no_sub
then
bsub<job.sh
fi
fi
