for x in $(seq 10 50)
do
    if [ -f "test/input/${x}_seqA.fa" -a -f "test/input/${x}_seqB.fa"  -a -f "test/output/${x}.txt" ]; then
        echo "File: ${x}_seqA.fa"
        ./sw.out -f "test/input/${x}_seqA.fa" -f "test/input/${x}_seqB.fa"  -m 4, -x -2, -o -3, -e -1
        cat "test/output/${x}.txt"
    fi
done
