if [ "$1" == "clean" ]
then
    rm *.mod
    rm -r run
    rm -r figs
    rm *.out
    rm fort.*
    rm hup*.*
    rm *coll.dat
    rm info.txt
else
    outname=$1".out"
    if [ "$outname" == ".out" ]; then
     echo "output will be named BallPark"
     outname="BallPark.out"
    fi

    rm -r run
    rm ./*.mod
    rm ./*.out

    echo "------------------------------------------------------------------------------------"
    echo "CREATING SUPPORT MODULES - EXPECT ERROR MESSAGES"
    echo "------------------------------------------------------------------------------------"
    export OMP_NUM_THREADS=4
    if [ "$2" == "ifort" ]; then
        ifort -qopenmp -O3 -mkl -r8 $F90FLAGS nr.f main.f90 $LINK_FNL -o $outname
    elif [ "$2" == "ser" ]; then
        gfortran -O3 -ffree-line-length-none -fmax-errors=3 -fdefault-real-8 nr.f main.f90 -o $outname
    elif [ "$1" == "checkall" ]; then
        gfortran -fopenmp -O0 -ffree-line-length-none -fmax-errors=3 -fdefault-real-8 -fcheck=all nr.f main.f90
    else
        ulimit -s unlimited
        gfortran -fopenmp -O3 -ffree-line-length-none -fmax-errors=3 -fdefault-real-8 nr.f main.f90 -o $outname
    fi
    echo "------------------------------------------------------------------------------------"
    echo "SECOND COMPILATION - NO ERRORS SHOULD APPEAR"
    echo "------------------------------------------------------------------------------------"
    if [ "$2" == "ifort" ]; then
        ifort -qopenmp -O3 -mkl -r8 $F90FLAGS nr.f main.f90 $LINK_FNL -o $outname
    elif [ "$2" == "ser" ]; then
        gfortran -O3 -ffree-line-length-none -fmax-errors=3 -fdefault-real-8 nr.f main.f90 -o $outname
    elif [ "$1" == "checkall" ]; then
        gfortran -fopenmp -O0 -ffree-line-length-none -fmax-errors=3 -fdefault-real-8 -fcheck=all nr.f main.f90
    else
        gfortran -fopenmp -O3 -ffree-line-length-none -fmax-errors=3 -fdefault-real-8 nr.f main.f90 -o $outname
    fi

    mkdir run

    if [ -z "$1" ]
    then
        ./$outname
    elif [ "$1" == "checkall" ]; then
        ./a.out
    else
        rm hup"$1".log
        if [ "$1" == "here" ]; then
            ./$outname
        else
            nohup "./$outname" > hup"$1".log &
        fi
    fi
fi
