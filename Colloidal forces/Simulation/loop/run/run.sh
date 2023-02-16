SYSTEM=$PWD
let NUM=0
INIT=$SYSTEM/npt.gro
while (($NUM<28))
do 
	mkdir prod_$NUM
    cd prod_$NUM
    XXX=`echo "3.0-$NUM*0.05" | bc`
    sed s/DISTANCE/$XXX/ < ../md.mdp > md.mdp
    gmx grompp -f md.mdp -c $INIT -n ../c240.ndx -p ../topol.top -o md.tpr -maxwarn 4
    gmx mdrun -v -nt 10 -deffnm md
    
    cd ..
    INIT=$SYSTEM/prod_$NUM/md.gro
    let NUM=$NUM+1
done
