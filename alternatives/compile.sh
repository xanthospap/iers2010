g++ -Wall -c -O2 -std=c++11 shclass.cpp -pg
g++ -Wall -c -O2 -std=c++11 gpt.cpp -pg
g++ -Wall -O2 -std=c++11 main.cpp gpt.o shclass.o -pg
