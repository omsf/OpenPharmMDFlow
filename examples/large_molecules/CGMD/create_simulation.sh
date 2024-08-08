set -e

mab="minimization-vac.gro"
water="../water.gro"
default_nmab=1
nmab=${1:-$default_nmab}
default_box_size=20
box_size=${2:-$default_box_size}

echo -e "#include \"../martini_v3.0.0.itp\"


#include \"molecule_0.itp\"

[ system ]
; name
Martini system from antibody.pdb

[ molecules ]
; name        number
molecule_0	1" > system.top

gmx editconf -f antibody-CG.pdb -d 2.0 -bt cubic -o antibody-CG.gro
gmx grompp -p system.top -f ../minimization.mdp -c antibody-CG.gro -o minimization-vac.tpr
gmx mdrun -deffnm minimization-vac -v


# now the mab molecules:
gmx insert-molecules -ci $mab -nmol $nmab -o mabs.gro -box $box_size

gmx solvate -cp mabs.gro -cs $water -radius 0.21 -o solvated.gro
nwater=$(grep "W" solvated.gro | wc -l)
cp solvated.gro system.gro


echo -e "#include \"../martini_v3.0.0.itp\"
#include \"../martini_v3.0.0_solvents_v1.itp\"


#include \"molecule_0.itp\"

[ system ]
; name
Martini system from antibody.pdb

[ molecules ]
; name        number
molecule_0          $nmab
W    $nwater" > system.top
