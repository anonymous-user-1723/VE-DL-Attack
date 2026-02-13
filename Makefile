all : test_for_attack key_recovery_attack
test_for_attack : rand_gen.cpp speck32.cpp util.cpp test_for_attack.cpp
	g++ -O2 -std=c++11 test_for_attack.cpp rand_gen.cpp speck32.cpp util.cpp -o test_for_attack

key_recovery_attack : rand_gen.cpp speck32.cpp util.cpp key_recovery_attack.cpp
	g++ -O2 -std=c++11 -pthread key_recovery_attack.cpp rand_gen.cpp speck32.cpp util.cpp -o key_recovery_attack
