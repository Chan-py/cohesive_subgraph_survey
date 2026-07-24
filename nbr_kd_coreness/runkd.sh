# #!/bin/bash

# # g++-11 -Wall -g  -o main main.cpp hypergraph.cpp  algorithms.cpp readhg.h utils.h #  run on clang
# g++ -std=c++11 -o kdmain main.cpp hypergraph.cpp  algorithms.cpp utils.h readhg.h #  run on gnu c++ compiler


# declare -a dset=("congress" "enron" "dblp" "pref" "aminer")
# declare -a algorithms=("kdcore")

# it=2 # #Iterations to run each algorithm on each dataset.
# log=0 # Activate logging to output core-numbers & iteration h-index statistics
# for dataset in "${dset[@]}"
# do
#     for algo in "${algorithms[@]}"
#     do
#         ./kdmain 1 $dataset $algo $it $log
#         echo "------------" 
#     done 
# done

# # declare -a dset=("klay" "protein")
# # declare -a algorithms=("kdcore")

# # it=1 # #Iterations to run each algorithm on each dataset.
# # log=0 # Activate logging to output core-numbers & iteration h-index statistics
# # for dataset in "${dset[@]}"
# # do
# #     for algo in "${algorithms[@]}"
# #     do
# #         ./main 1 $dataset $algo $it $log
# #         echo "------------" 
# #     done 
# # done
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# 출력 폴더 준비
OUT="../output"
mkdir -p "$OUT"

# 컴파일러/옵션
CXX="${CXX:-$(command -v g++-11 || command -v g++ || command -v clang++)}"
CXXFLAGS="-std=c++17 -O3 -Wall -Wextra"
INCLUDES="-I."

# 소스 목록 (readhg가 헤더만이면 *.cpp에 안 넣어도 됨)
SRC=(main.cpp algorithms.cpp hypergraph.cpp)

echo "[build] ${CXX} ${CXXFLAGS} ${SRC[*]} -o kdmain"
${CXX} ${CXXFLAGS} ${SRC[@]} ${INCLUDES} -o kdmain

# 인자: dataset [iterations] [log]
DATASET="${1:-real/congress}"
IT="${2:-1}"
LOG="${3:-0}"
# Local-core-OPTIV
echo "[run] ./kdmain 1 \"$DATASET\" kdcore \"$IT\" \"$LOG\""
./kdmain 1 "$DATASET" kdcore "$IT" "$LOG"

echo "---- done ----"
echo "results.csv      : ../output/results.csv (실행시간 등)"
echo "kd-core CSV      : ../output/core_kdcore_${DATASET}.csv  (node_id, k, d)"