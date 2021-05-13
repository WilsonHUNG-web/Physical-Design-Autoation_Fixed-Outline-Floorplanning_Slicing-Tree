# Physical-Design-Autoation_Fixed-Outline-Floorplanning_Slicing-Tree

This Slicing-Tree implementation may fail on some testcases. Another repository: https://github.com/WilsonHUNG-web/Physical-Design-Automation_Fixed-outline-Floorplan-Design_B-star-tree gives shorter CPU time and feasible solutions for each given testcase.

--How to compile <br>

  In directory './src', enter the following command: <br>
  ```
  $ make
  ```
  It will generate the executable file "main" in "HW3/bin/". <br>
  
 If you want to remove it, please enter the following command"<br>
 ```
  $ make clean
 ```
  
--How to Run
 In directory './src', enter the following command: <br>
  Usage: ../bin/[exe] [hardblocks file] [nets file] [pl file] [dead space_ratio]<br>
  e.g.<br>
  ```
  $ ../bin/hw3 ../testcase/n100.hardblocks ../testcase/n100.nets ../testcase/n100.pl 0.1 
  ```
