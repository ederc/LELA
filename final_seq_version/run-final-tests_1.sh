#!/bin/bash

#oarsub -p "host='big23' OR 'big24' OR 'big25' OR 'big26'" -l core=16,walltime=1000:0:0 --notify "mail:fayssal.martani@lip6.fr" "/home/martani/LELA_private/final_seq_version/run-final-tests_1.sh"

uname -a | tee -a ~/output_FGL/host_test1.txt
cat /proc/cpuinfo | tee -a ~/output_FGL/host_test1.txt

cd /home/martani/LELA_private/final_seq_version


# ======================================= test-FG-parallel =====================================================================
#minrank_minors_9_9_6 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "-"    ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_oldmethod_parallel_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m" ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_newmethod_parallel_only_d.txt;

#minrank_minors_9_9_6 RREF
./run-tests.sh "./test-FGL-parallel -f" "- -d"    ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_oldmethod_parallel_rref.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -d" ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_newmethod_parallel_rref.txt;


#~/matrices/F5/kat13 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "-"    ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_oldmethod_parallel_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m" ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_newmethod_parallel_only_d.txt;

#~/matrices/F5/kat13 RREF
./run-tests.sh "./test-FGL-parallel -f" "- -d"    ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_oldmethod_parallel_rref.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -d" ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_newmethod_parallel_rref.txt;


#~/matrices/F5/kat16 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "-"    ~/matrices/F5/kat16 ~/output_FGL/kat16_newprog_oldmethod_parallel_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m" ~/matrices/F5/kat16 ~/output_FGL/kat16_newprog_newmethod_parallel_only_d.txt;

#~/matrices/F5/kat16 RREF
./run-tests.sh "./test-FGL-parallel -f" "- -d"    ~/matrices/F5/kat16 ~/output_FGL/kat16_newprog_oldmethod_parallel_rref.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -d" ~/matrices/F5/kat16 ~/output_FGL/kat16_newprog_newmethod_parallel_rref.txt;


#~/matrices/F5/mr10 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "-"    ~/matrices/F5/mr10 ~/output_FGL/mr10_newprog_oldmethod_parallel_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m" ~/matrices/F5/mr10 ~/output_FGL/mr10_newprog_newmethod_parallel_only_d.txt;

#~/matrices/F5/mr10 RREF
./run-tests.sh "./test-FGL-parallel -f" "- -d"    ~/matrices/F5/mr10 ~/output_FGL/mr10_newprog_oldmethod_parallel_rref.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -d" ~/matrices/F5/mr10 ~/output_FGL/mr10_newprog_newmethod_parallel_rref.txt;



##################################################################################################################
#~/matrices/F5/kat11 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "-"    ~/matrices/F4/kat11 ~/output_FGL/kat11_newprog_oldmethod_parallel_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m" ~/matrices/F4/kat11 ~/output_FGL/kat11_newprog_newmethod_parallel_only_d.txt;

#~/matrices/F5/kat11 RREF
./run-tests.sh "./test-FGL-parallel -f" "- -d"    ~/matrices/F4/kat11 ~/output_FGL/kat11_newprog_oldmethod_parallel_rref.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -d" ~/matrices/F4/kat11 ~/output_FGL/kat11_newprog_newmethod_parallel_rref.txt;


#~/matrices/F5/kat12 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "-"    ~/matrices/F4/kat12 ~/output_FGL/kat12_newprog_oldmethod_parallel_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m" ~/matrices/F4/kat12 ~/output_FGL/kat12_newprog_newmethod_parallel_only_d.txt;

#~/matrices/F5/kat12 RREF
./run-tests.sh "./test-FGL-parallel -f" "- -d"    ~/matrices/F4/kat12 ~/output_FGL/kat12_newprog_oldmethod_parallel_rref.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -d" ~/matrices/F4/kat12 ~/output_FGL/kat12_newprog_newmethod_parallel_rref.txt;


##################################################################################################################


# ======================================= test-FG-parallel =====================================================================
#minrank_minors_9_9_6 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "- -p 4" ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_oldmethod_parallel_p4_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -p 4" ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_newmethod_parallel_p4_only_d.txt;

#~/matrices/F5/kat13 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "- -p 4" ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_oldmethod_parallel_p4_only_d.txt;
./run-tests.sh "./test-FGL-parallel -f" "- -m -p 4" ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_newmethod_parallel_p4_only_d.txt;

#minrank_minors_9_9_6 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "- -m -p 12" ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_newmethod_parallel_p12_only_d.txt;

#~/matrices/F5/kat13 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "- -m -p 12" ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_newmethod_parallel_p12_only_d.txt;

#minrank_minors_9_9_6 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "- -m -p 16" ~/matrices/F5/minrank_minors_9_9_6 ~/output_FGL/mr996_newprog_newmethod_parallel_p16_only_d.txt;

#~/matrices/F5/kat13 ONLY D
./run-tests.sh "./test-FGL-parallel -f" "- -m -p 16" ~/matrices/F5/kat13 ~/output_FGL/kat13_newprog_newmethod_parallel_p16_only_d.txt;

