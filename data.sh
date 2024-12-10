for nt in $(seq 1 8)
do
    export OMP_NUM_THREADS=$nt
    echo "Testing with $nt threads"
    ./test.sh > data${nt}.txt
done

