cd ..
dir=$PWD
cd loop
SYSTEM=$dir/System
MDP=$dir/MDP
INIT=$dir/MD/run/run.gro

# generate random number between 2 and 508 inclusive
let Y=508
let X=2
BASE=$((Y-X+1))

# loop over randomly chosen water molecules
let K=1
let KMAX=500
while (($K < $KMAX))
do

    # randomly select test molecule
    NTOP=$(($(($RANDOM%$BASE))+X))
    let NBOT=508-${NTOP}
    sed "s/NTOP/${NTOP}/;s/NBOT/${NBOT}/" < $SYSTEM/systemq0.top > $SYSTEM/systemq.top

    # create index file with test molecule
    let OW=1+3*${NTOP}+1
    let HW1=${OW}+1
    let HW2=${OW}+2
    sed "s/OW/${OW}/;s/HW1/${HW1}/;s/HW2/${HW2}/" < $SYSTEM/systemq.ndx > systemq.ndx

    # rerun without test charge
    gmx grompp -f $MDP/run.mdp -c $INIT -p $SYSTEM/topol.top -n systemq.ndx -o rerun.tpr -maxwarn 20
    gmx mdrun -nt 1 -deffnm rerun -rerun run.trr

    # force without test charge
    echo "7" | gmx traj -f rerun.trr -s rerun.tpr -n systemq.ndx -com -of force_${K}.xvg

    # rerun with test charge
    gmx grompp -f $MDP/run.mdp -c $INIT -p $SYSTEM/systemq.top -n systemq.ndx -o rerunq.tpr -maxwarn 20
    gmx mdrun -nt 1 -deffnm rerunq -rerun run.trr

    # force with test charge
    echo "7" | gmx traj -f rerunq.trr -s rerunq.tpr -n systemq.ndx -com -of forceq_${K}.xvg

    # test molecule atom positions (for dipole, charge distrib,...)
    echo "7" | gmx traj -f rerunq.trr -s rerunq.tpr -n systemq.ndx -ox posn_${K}.xvg

    rm "#"*"#"
    let K=$K+1
done
