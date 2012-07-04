#!/bin/bash
uname -a | tee -a ~/output/host.txt
cat /proc/cpuinfo | tee -a ~/output/host.txt

cd /home/martani/LELA_private/FG-multiline

# ========================== test-FG-multiline ==========================================
#F4 ONLY D
./run-tests.sh "./test-FG-multiline -f" "- -d -m" ~/matrices/F4/kat11 ~/output/kat11_only_D__FG_NEW.txt;
./run-tests.sh "./test-FG-multiline -f" "- -d -m" ~/matrices/F4/kat12 ~/output/kat12_only_D__FG_NEW.txt;

#F4 RREF
./run-tests.sh "./test-FG-multiline -f" "- -c" ~/matrices/F4/kat11 ~/output/kat11_rref__FG_NEW.txt;
./run-tests.sh "./test-FG-multiline -f" "- -c" ~/matrices/F4/kat12 ~/output/kat12_rref__FG_NEW.txt;


#F5 ONLY D
./run-tests.sh "./test-FG-multiline -f" "- -d -m" ~/matrices/F5/minrank_minors_9_9_6 ~/output/minrank_minors_9_9_6_only_D__FG_NEW.txt;
./run-tests.sh "./test-FG-multiline -f" "- -d -m" ~/matrices/F5/kat13 ~/output/kat13_only_D__FG_NEW.txt;
./run-tests.sh "./test-FG-multiline -f" "- -d -m" ~/matrices/F5/mr10 ~/output/mr10_only_D__FG_NEW.txt;


#F5 RREF
./run-tests.sh "./test-FG-multiline -f" "- -c" ~/matrices/F5/minrank_minors_9_9_6 ~/output/minrank_minors_9_9_6_rref__FG_NEW.txt;
./run-tests.sh "./test-FG-multiline -f" "- -c" ~/matrices/F5/kat13 ~/output/kat13_rref__FG_NEW.txt;
./run-tests.sh "./test-FG-multiline -f" "- -c" ~/matrices/F5/mr10 ~/output/mr10_rref__FG_NEW.txt;



# ========================== testRrefInterface.th0.b128 ==========================================
#F4 ONLY D
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "1 35" ~/matrices/F4/kat11 ~/output/kat11_only_D__FGL.txt;
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "1 35" ~/matrices/F4/kat12 ~/output/kat12_only_D__FGL.txt;

#F4 RREF
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "0 35" ~/matrices/F4/kat11 ~/output/kat11_rref__FGL.txt;
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "0 35" ~/matrices/F4/kat12 ~/output/kat12_rref__FGL.txt;


#F5 ONLY D
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "1 35" ~/matrices/F5/minrank_minors_9_9_6 ~/output/minrank_minors_9_9_6_only_D__FGL.txt;
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "1 35" ~/matrices/F5/kat13 ~/output/kat13_only_D__FGL.txt;
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "1 35" ~/matrices/F5/mr10 ~/output/mr10_only_D__FGL.txt;

#F5 RREF
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "0 35" ~/matrices/F5/minrank_minors_9_9_6 ~/output/minrank_minors_9_9_6_rref__FGL.txt;
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "0 35" ~/matrices/F5/kat13 ~/output/kat13_rref__FGL.txt;
./run-tests.sh ~/matrices/testRrefInterface.th0.b128 "0 35" ~/matrices/F5/mr10 ~/output/mr10_rref__FGL.txt;




# oarsub -p "host='big1'" -l core=8,walltime=1000:0:0 --notify "mail:fayssal.martani@lip6.fr" "/home/martani/LELA_private/FG-multiline/run-all.sh"

