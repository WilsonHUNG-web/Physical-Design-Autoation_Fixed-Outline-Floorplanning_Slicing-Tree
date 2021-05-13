echo "---START run n100 n200 n300 with dead ratio 0.15 and 0.1---"

time ../bin/hw3 ../testcase/n100.hardblocks ../testcase/n100.nets ../testcase/n100.pl ../output/n100.floorplan 0.15 
time ../bin/hw3 ../testcase/n100.hardblocks ../testcase/n100.nets ../testcase/n100.pl ../output/n100.floorplan 0.1 
time ../bin/hw3 ../testcase/n200.hardblocks ../testcase/n200.nets ../testcase/n200.pl ../output/n200.floorplan 0.15 
time ../bin/hw3 ../testcase/n200.hardblocks ../testcase/n200.nets ../testcase/n200.pl ../output/n200.floorplan 0.1 
time ../bin/hw3 ../testcase/n300.hardblocks ../testcase/n300.nets ../testcase/n300.pl ../output/n300.floorplan 0.15 
time ../bin/hw3 ../testcase/n300.hardblocks ../testcase/n300.nets ../testcase/n300.pl ../output/n300.floorplan 0.1 
echo "----END run n100 n200 n300 ----"

