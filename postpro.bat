mkdir figs
if [ "$1" == "traj" ]
then
    Nd=$(grep Nd ./info.txt | awk '{print $3}')
    Nv=$(grep Nv ./info.txt | awk '{print $3}')
    Na=$(grep Na ./info.txt | awk '{print $3}')
    echo "Plotting Traj.py with Nd=$Nd,Nv=$Nv,Na=$Na"
    cd gnuplot
    python Traj.py "Nd=$Nd,Nv=$Nv,Na=$Na"
fi

if [ "$1" == "cone" ]
then
    cd gnuplot
    python Cone.py "" # > ../conereport.txt
    vmean=$(grep vmean ../param.in | awk '{print $2}')
    vdev=$(grep vdev ../param.in | awk '{print $2}')
    angmax=$(grep angmax ../param.in | awk '{print $2}')
    angdev=$(grep angdev ../param.in | awk '{print $2}')
    angdev=${angdev##*=}
    angmax=${angmax##*=}
    vdev=${vdev##*=}
    vmean=${vmean##*=}    
    #verror=$(grep Volume_normalized ../conereport.txt | awk '{print $2}')
    #herror=$(grep Central_height_normalized ../conereport.txt | awk '{print $2}')
    #echo "$verror  $herror  $vmean  $vdev  $angmax  $angdev" >> ../../Misfit.txt
fi

if [ "$1" == "dist" ]
then
    cd gnuplot
    if [ -z "$2" ]
    then
        python Distrib.py ""
    else
        python Distrib.py "dname='$2'"
    fi    
fi
