for x in $(seq 1 12)
do
    # echo "test/input/${x}_seqA.fa"
    if [ -f "test/input/${x}_seqA.fa" -a -f "test/input/${x}_seqB.fa"  -a -f "test/output/${x}.txt" ]; then
        ./sw.out -f "test/input/${x}_seqA.fa" -f "test/input/${x}_seqB.fa"
        cat "test/output/${x}.txt"
    fi
done
