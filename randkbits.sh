echo 32 bit, default RNG
g++ -Ofast -march=native -g -W -Wall -std=c++11 randkbits.cc -DWORD=uint32_t               && ./a.out
echo 32 bit, fast RNG
g++ -Ofast -march=native -g -W -Wall -std=c++11 randkbits.cc -DWORD=uint32_t -DFAST_RNG=1  && ./a.out
echo 64 bit, default RNG
g++ -Ofast -march=native -g -W -Wall -std=c++11 randkbits.cc -DWORD=uint64_t               && ./a.out
echo 64 bit, fast RNG
g++ -Ofast -march=native -g -W -Wall -std=c++11 randkbits.cc -DWORD=uint64_t -DFAST_RNG=1  && ./a.out
