g++ $(ls *.cpp|grep -v spli|grep -v sca |grep -v mer|grep -v sff|grep -v fasta|grep -v Light|grep -v module_join_main.cpp|grep -v dump_sql_main.cpp) -g -o dna_gdb
