
#!/bin/bash
#Execute from the parent directory (containing ELB, PHB, EPM, ISO directories).

#Specify the number of bands to be wannierized, the index of the highest non-wannierized band and the initial guess...
#...for the wannier functions. If you change this to get the windows in #3 right, you have to run #1 (with no_sub = true) again ...
#...in order to reset the corresponding wannierization windows again.
nbndsub=4
wannier_init=("Pb:sp3")

#Specify if the wannierization energy windows should be automatically determined. If you set it to false, specify the...
#...windows in the next block. The inner window will be the largest possible window that contains only the specified bands...
#...but no parts of other bands. The outer window will be the smallest possible window that contains the whole specified...
#bands including other bands should they fall in this window. You can narrow the inner and outer window with dE which specifies...
#the value in eV that gets added/subtracted to the bottom/top of the inner window (i.e. negative values will broaden it).
auto_window=true
dE_inner=0.05
dE_outer=-1.0


#Set the inner and outer wannierization energy windows manually (in eV). You only need to specify the values of this block...
#... if you disabled auto_window or if the automatic determination didn't work.
nbndskip=0
dis_froz_min=0.0
dis_froz_max=0.0
dis_win_min=0.0
dis_win_max=0.0

#Specify the starting and ending point of calculation. 
#1: input preparation and scf calculation (set no_sub to true for just the input preparation)
#2: nscf calculation for phonon linewidths
#3: wannierization of the electron bands
#Note it is a good idea to stop the calculation here and check the wannierized bands (i.e. compare them to the...
#...Bloch bands). You can use #11 to do that.
#4: calculation of electron-phonon coupling coefficients for phonon linewidths
#5: phonon linewidth calculation along high symmetry lines i.e. bands 
#6: phonon linewidth and a2F calculation in fine grid 
#7: scf calculation for solving the isotropic Eliashberg functions (Tc)
#8: nscf caclulation for Tc
#9: solving the isotropic Eliashberg functions
#(#11: Plot the Wannier and Bloch bands using gnuplot on the cluster (make sure x11 forwarding is enabled)...
#     ...This won't submit a job but simply plot the Wannier bands obtained after #3 and compare them to the...
#     ...Bloch bands obtained through pw.x (located in ELB). Furthermore the energy windows will also be plotted...
#     ...Usable only as a separate (calc_start = calc_end) execution.
#(#12: tidy up - deletes everything but the obtained end results. This is not necessary but saves space...
#     ...Usable only as a seperate (calc_start = calc_end) execution. Only use once you are sure that everything worked.)

calc_start=1
calc_end=3

#If you need a specific submission script (e.g. for wannierization) set calc_start = calc_end (e.g. = 3) and set no_sub
#to true. You'll find the job script as job.sh in the corresponding directory.
no_sub=false

#load modules 
module purge 
module load quantum_espresso/6.2.1 &>/dev/null 
module load python/3.6.0 &>/dev/null 
module load gnuplot &>/dev/null 


#________________________________EXECUTE_THE_SCRIPT________________________________#
#__________________________________________________________________________________#

ref_bands=false
base_dir=$(pwd)
ref_dir=""
if ! $ref_bands
then
ref_dir=$base_dir
fi

#__________POSTPROCESSING__________#

if [ ! -d $ref_dir/PHB/save ]
then
cd $ref_dir/PHB
python pp.py
cd $base_dir
fi

#__________PREPARE_INPUT___________#
#IF CALC_START
if (($calc_start == 1))
then
cp $base_dir/ELB/scf.in $base_dir/EPM

cp $base_dir/EPM/nscf.in_orig $base_dir/EPM/nscf.in
cp $base_dir/EPM/wannier.in_orig $base_dir/EPM/wannier.in
cp $base_dir/EPM/wannier_epm.in_orig $base_dir/EPM/wannier_epm.in
cp $base_dir/EPM/ph_lw.in_orig $base_dir/EPM/ph_lw.in
cp $base_dir/EPM/a2F.in_orig $base_dir/EPM/a2F.in
cp $base_dir/ISO/eliashberg_iso.in_orig $base_dir/ISO/eliashberg_iso.in

#generate and append uniform k-grid for nscf.in
python $base_dir/EPM/kmesh.py 8 8 8 1 >> $base_dir/EPM/nscf.in

#append irreducible q-points to epw input files where necessary
num_of_irr_k=$(sed "2q;d" $ref_dir/PHB/*.dyn0 | awk '{print $1}')
echo "$num_of_irr_k cartesian" >> $base_dir/EPM/wannier_epm.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/EPM/wannier_epm.in
echo "$num_of_irr_k cartesian" >> $base_dir/EPM/ph_lw.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/EPM/ph_lw.in
echo "$num_of_irr_k cartesian" >> $base_dir/EPM/a2F.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/EPM/a2F.in
echo "$num_of_irr_k cartesian" >> $base_dir/ISO/eliashberg_iso.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/ISO/eliashberg_iso.in


#get and set the Fermi energy where necessary
Ef=$(grep Fermi $ref_dir/ELB/scf.out | awk '{print $5}')
line=$(grep -n fermi_energy $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$Ef/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n fermi_energy $base_dir/EPM/ph_lw.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$Ef/3" $base_dir/EPM/ph_lw.in
line=$(grep -n fermi_energy $base_dir/EPM/a2F.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$Ef/3" $base_dir/EPM/a2F.in
line=$(grep -n fermi_energy $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$Ef/3" $base_dir/ISO/eliashberg_iso.in


#copy the updated and necessary input files to the ISO directory
cp $base_dir/EPM/scf.in $base_dir/ISO
cp $base_dir/EPM/nscf.in $base_dir/ISO

#if enabled, get the wannierization energy windows automatically from the electron band calculation
#and set the values in the input files where necessary

#IF AUTO_WINDOW
if $auto_window
then

#get the number of bands
bands=$(head -n1 $ref_dir/ELB/*bands.dat | awk '{print $3}' | sed "s/,//")
#get the number of k-points in the total band path
points=$(head -n1 $ref_dir/ELB/*bands.dat | awk '{print $5}')

#make a temporary file that stores the band energies sequentially per band
(cat $ref_dir/ELB/*bands.dat.gnu | sed 's/[^ ]*[^ ]//1' | sed 's/\ //g' | sed '/^$/d') > $ref_dir/ELB/bands_temp

#run the py script to find the bands
windows=$(python $ref_dir/ELB/wannier_windows.py $nbndsub $bands $points $Ef $ref_dir/ELB/bands_temp)
smallest_inner=$(echo $windows | awk '{print $1}')
largest_inner=$(echo $windows | awk '{print $2}')
smallest_outer=$(echo $windows | awk '{print $3}')
largest_outer=$(echo $windows | awk '{print $4}')
nbndskip=$(echo $windows | awk '{print $5}')

#narrow the bands
smallest_inner=$(echo "$smallest_inner + $dE_inner" | bc -l)
largest_inner=$(echo "$largest_inner - $dE_inner" | bc -l)
smallest_outer=$(echo "$smallest_outer + $dE_outer" | bc -l)
largest_outer=$(echo "$largest_outer - $dE_outer" | bc -l)

dis_froz_min=$smallest_inner
dis_froz_max=$largest_inner
dis_win_min=$smallest_outer
dis_win_max=$largest_outer

fi
#ENDIF AUTO_WINDOW

#set the wannierization parameters
line=$(grep -n nbndsub $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$nbndsub/3" $base_dir/EPM/wannier.in
line=$(grep -n nbndskip $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$nbndskip/3" $base_dir/EPM/wannier.in

index=0
for p in "${wannier_init[@]}"
do
    ((index++))
    line=$(grep -n "proj(${index})" $base_dir/EPM/wannier.in | cut -d : -f 1)
    sed -i -e "${line}s/[^ ]*[^ ]/'$p'/3" $base_dir/EPM/wannier.in
done

line=$(grep -n nbndsub $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$nbndsub/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n nbndskip $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$nbndskip/3" $base_dir/EPM/wannier_epm.in

index=0
for p in "${wannier_init[@]}"
do
    ((index++))
    line=$(grep -n "proj(${index})" $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
    sed -i -e "${line}s/[^ ]*[^ ]/'$p'/3" $base_dir/EPM/wannier_epm.in
done

line=$(grep -n nbndsub $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$nbndsub/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n nbndskip $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$nbndskip/3" $base_dir/ISO/eliashberg_iso.in

index=0
for p in "${wannier_init[@]}"
do
    ((index++))
    line=$(grep -n "proj(${index})" $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
    sed -i -e "${line}s/[^ ]*[^ ]/'$p'/3" $base_dir/ISO/eliashberg_iso.in
done


line=$(grep -n dis_win_min $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_win_min/3" $base_dir/EPM/wannier.in
line=$(grep -n dis_win_max $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_win_max/3" $base_dir/EPM/wannier.in
line=$(grep -n dis_froz_min $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_froz_min/3" $base_dir/EPM/wannier.in
line=$(grep -n dis_froz_max $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_froz_max/3" $base_dir/EPM/wannier.in

line=$(grep -n dis_win_min $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_win_min/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n dis_win_max $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_win_max/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n dis_froz_min $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_froz_min/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n dis_froz_max $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_froz_max/3" $base_dir/EPM/wannier_epm.in

line=$(grep -n dis_win_min $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_win_min/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n dis_win_max $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_win_max/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n dis_froz_min $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_froz_min/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n dis_froz_max $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${line}s/[^ ]*[^ ]/$dis_froz_max/3" $base_dir/ISO/eliashberg_iso.in

fi
#ENDIF CALC_START


#__________JOB_SUBMISSION___________#
if ( (($calc_start == 11)) && ((! $calc_end == 11)) ) || ( ((! $calc_start == 11)) && (($calc_end == 11)) )
then
echo "Error: plotting the Bloch and Wannier bands is only possible as a separate execution"
exit
elif (($calc_start == 11)) && (($calc_start == $calc_end))
then

#get the positions of the last k-points of the bands
last_QE=$(tail -n2 $ref_dir/ELB/*bands.dat.gnu | head -n1 | awk '{print $1}')
last_EPM=$(tail -n2 $base_dir/EPM/Pb_band.dat | head -n1 | awk '{print $1}')

#get the windows
smallest_inner=$(grep dis_froz_min $base_dir/EPM/wannier_epm.in | awk '{print $3}')
largest_inner=$(grep dis_froz_max $base_dir/EPM/wannier_epm.in | awk '{print $3}')
smallest_outer=$(grep dis_win_min $base_dir/EPM/wannier_epm.in | awk '{print $3}')
largest_outer=$(grep dis_win_max $base_dir/EPM/wannier_epm.in | awk '{print $3}')

#make the gnuplot script
cat > $base_dir/ELB/check_bands.gnu << EOF
reset
set terminal x11 dashed title "Pb"
set grid
set linestyle 1 lc rgb "black" lw 2
set linestyle 2 lc rgb "red" lw 1
set linestyle 3 lc rgb "blue" lw 2
set linestyle 4 lc rgb "green" lw 2

y_low = $smallest_outer - 1
y_high = $largest_outer + 1
set yrange [y_low:y_high]

set ylabel "[eV]"

set arrow from graph 0, first $smallest_outer to graph 1, first $smallest_outer ls 3 nohead
set arrow from graph 0, first $largest_outer to graph 1, first $largest_outer ls 3 nohead
set arrow from graph 0, first $smallest_inner to graph 1, first $smallest_inner ls 4 nohead
set arrow from graph 0, first $largest_inner to graph 1, first $largest_inner ls 4 nohead

set arrow from graph 0.95, first $smallest_outer to graph 0.95, first $largest_outer ls 3 heads
set arrow from graph 0.9, first $smallest_inner to graph 0.9, first $largest_inner ls 4 heads

set label "outer window" at graph 0.96, graph 0.5 rotate by 90 tc lt 3  font "Helvetica, 12"
set label "inner window" at graph 0.91, graph 0.5 rotate by 90 tc lt 2 font "Helvetica, 12"
blochf=system("ls $ref_dir/ELB/*bands.dat.gnu")
plot blochf u (\$1/$last_QE):2 w l ls 1 t "Bloch", \
     "$base_dir/EPM/Pb_band.dat" u (\$1/$last_EPM):2 w l ls 2 t "Wannier"
EOF
gnuplot -p "$base_dir/ELB/check_bands.gnu" 
fi

if ( (($calc_start == 12)) && ((! $calc_end == 12)) ) || ( ((! $calc_start == 12)) && (($calc_end == 12)) )
then
echo "Error: please only tidy up once you are sure that all calculations have finished correctly"
exit
fi

if (($calc_start > $calc_end))
then
echo "Error: calculation order makes no sense"
exit
fi
#EPM calculations
cd $base_dir/EPM

#scf
if (($calc_start <= 1)) && (($calc_end >= 1))
then
cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_epm_scf.o
#BSUB -e Pb_epm_scf.e
#BSUB -J Pb_epm_scf
mpirun pw.x -npool 4 -in scf.in > scf.out
#get the number of bands from the scf.out file and put them in the nscf.in file
nbnd=\$(grep nbnd nscf.in | awk '{print \$3}')
if [ \$nbnd -eq 0 ]
then
    nbnd=\$(grep -m1 "number of Kohn-Sham states" scf.out | awk '{print \$5}')
    line=\$(grep -n nbnd nscf.in | cut -d : -f 1)
    sed -i "\${line}s/[^ ]*[^ ]/\${nbnd}/3" nscf.in
fi
EOF

if ! $no_sub
then
bsub<job.sh
fi

fi

#nscf
if (($calc_start <= 2)) && (($calc_end >= 2))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_epm_nscf.o
#BSUB -e Pb_epm_nscf.e
#BSUB -J Pb_epm_nscf
#BSUB -w Pb_epm_scf
mpirun pw.x -npool 4 -in nscf.in > nscf.out
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

#wannierization
if (($calc_start <= 3)) && (($calc_end >= 3))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 4:00
#BSUB -o Pb_wannier.o
#BSUB -e Pb_wannier.e
#BSUB -J Pb_wannier
#BSUB -w Pb_epm_nscf
mpirun epw.x -npool 4 -in wannier.in > wannier.out
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

#electron phonon matrix
if (($calc_start <= 4)) && (($calc_end >= 4))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 16
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_epm.o
#BSUB -e Pb_epm.e
#BSUB -J Pb_epm
#BSUB -w Pb_wannier
mpirun epw.x -npool 16 -in wannier_epm.in > wannier_epm.out
EOF

if (($calc_start == 4))
then
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i ${line}d job.sh
fi
if ! $no_sub
then
bsub<job.sh
fi

fi

#phonon linewidth
if (($calc_start <= 5)) && (($calc_end >= 5))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 16
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_ph_lw.o
#BSUB -e Pb_ph_lw.e
#BSUB -J Pb_ph_lw
#BSUB -w Pb_epm
mpirun epw.x -npool 16 -in ph_lw.in > ph_lw.out
mv linewidth.phself linewidth.phself_bands
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

#a2F
if (($calc_start <= 6)) && (($calc_end >= 6))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 16
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_a2F.o
#BSUB -e Pb_a2F.e
#BSUB -J Pb_a2F
#BSUB -w Pb_ph_lw
mpirun epw.x -npool 16 -in a2F.in > a2F.out
mv linewidth.phself linewidth.phself_grid
EOF

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

#ISO calculations
cd $base_dir/ISO
#scf
if (($calc_start <= 7)) && (($calc_end >= 7))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_iso_scf.o
#BSUB -e Pb_iso_scf.e
#BSUB -J Pb_iso_scf
mpirun pw.x -npool 4 -in scf.in > scf.out
#get the number of bands from the scf.out file and put them in the nscf.in file
nbnd=\$(grep nbnd nscf.in | awk '{print \$3}')
if [ \$nbnd -eq 0 ]
then
    nbnd=\$(grep "number of Kohn-Sham states" scf.out | awk '{print \$5}')
    line=\$(grep -n nbnd nscf.in | cut -d : -f 1)
    sed -i "\${line}s/[^ ]*[^ ]/\${nbnd}/3" nscf.in
fi
EOF

if ! $no_sub
then
bsub<job.sh
fi

fi

#nscf
if (($calc_start <= 8)) && (($calc_end >= 8))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 4
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_iso_nscf.o
#BSUB -e Pb_iso_nscf.e
#BSUB -J Pb_iso_nscf
#BSUB -w Pb_iso_scf
mpirun pw.x -npool 4 -in nscf.in > nscf.out
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

#isotropic eliashberg
if (($calc_start <= 9)) && (($calc_end >= 9))
then

cat > job.sh << EOF
#!/bin/bash
#BSUB -n 16
#BSUB -R "rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -o Pb_iso.o
#BSUB -e Pb_iso.e
#BSUB -J Pb_iso
#BSUB -w Pb_iso_nscf
mpirun epw.x -npool 16 -in eliashberg_iso.in > eliashberg_iso.out
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
if (($calc_start == 12)) && (($calc_end == 12))
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
rm -r _ph0 *.save
rm *.dyn* *.e *.o *.fc
#EPM directory
cd $base_dir/EPM
rm -r *save
rm *.e *.o
rm *.epb* *wfc* *.xml *.wout *.win *.ukk *.nnkp *.kmap *.epmatwe1 *.chk decay* *.fmt

EOF
if ! $no_sub
then
bsub<job.sh
fi

fi
