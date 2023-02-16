mkdir run
cd run
SYSTEM=$PWD
let NUM=0
INIT=$SYSTEM/npt.gro
while (($NUM<2))
do 
	mkdir prod_$NUM
    cd prod_$NUM
    XXX=`echo "3.0-$NUM*0.05" | bc -1`
    sed s/DISTANCE/$XXX/ < ../md.mdp > md.mdp
    gmx grompp -f md.mdp -c $INIT -n ../c240.ndx -p ../topol.top -o md.tpr -maxwarn 4
    #gmx mdrun -nt 8 -deffnm md
    
    cd ..
    INIT=$SYSTEM/prod_$NUM/md.gro
    let NUM=$NUM+1
done
