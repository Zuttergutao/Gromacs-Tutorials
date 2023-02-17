# Dielectric screening

教程来源：[Dielectric screening | MD Simulation Techniques and Applications (psu.edu)](https://sites.psu.edu/simtech/dielectric-response/)

## 前言

### 离子微扰水层结构

溶液中的离子从立体和静电角度扰乱了附近水分子的排列。  离子占据了空间，这本身就导致了附近水分子的非均匀浓度分布，离子本身有一个排斥半径和有一个近邻的外壳层，当局部浓度弛豫到其平均值时，会出现衰减的振荡。  对于离子电荷附近的水分子的排布，像Na+这样的正离子倾向于吸引带负电的氧原子，排斥水面上带正电的氢原子。
静电势也同样因离子的存在而受到扰乱。  在真空中，离子会有一个库伦场，以$\frac{1}{r^2}$的速度衰减。  在水中，附近水分子的定向偶极子倾向于屏蔽离子的场，因为带相反电荷的原子聚集在离子周围的壳层中，从而平衡体系电荷。
我们可以使用MD来探索水在离子附近的这种屏蔽行为。通过将一个Na+离子固定在水盒子的中心，使用模拟轨迹计算平均电场、偶极方向、电荷密度和固定离子附近的分子分布。

<img src="loop\1-1.png" alt="1-1" style="zoom: 67%;" />



## 模拟

### 构建和平衡体系

体系为水盒子中含有一个Na<sup>+</sup>离子。首先创建一个gro文件，将Na+离子摆放在`2.5*2.5*2.5 nm`盒子正中央，然后使用`gmx solvate`给体系添加水。接着就是能量最小化、100 ps NVT和100 ps NPT，步长取2 fs。在模拟期间，Na<sup>+</sup>离子必须被冻结。使用oplsaa力场

1. 创建Na<sup>+</sup>离子gro文件。

    `# Na.gro`

    ```tex
    Na ion
    1
        1NA      NA    1   1.250   1.250   1.250
       2.50000   2.50000   2.50000
    ```

2. 溶剂化盒子（记住加了多少水分子）

    ```bash
    gmx solvate -cp Na.gro -o system.gro -box 2.5 2.5 2.5
    ```

3. 构建topol (此例使用oplsaa力场)
    `# topol.top`

    ```tex
    #include "oplsaa.ff/forcefield.itp"
    #include "oplsaa.ff/ions.itp"
    #include "oplsaa.ff/spce.itp"
    
    [ system ]
    ; Name
    Na+ in water 
    
    [ molecules ]
    ; Compound             #mols
    NA			1
    SOL			509   ; 填入相应的水分子数
    ```

4. EM
    `#em.mdp`

    ```tex
    define			 = -DFLEXIBLE
    ; Run control
    integrator               = steep 
    nsteps                   = 5000
    ; EM criteria and other stuff
    emtol                    = 100
    emstep                   = 0.01
    
    ; Output control
    nstlog                   = 10
    nstenergy                = 10
    
    ; cutoffs
    cutoff-scheme            = verlet
    nstlist                  = 20
    ns_type                  = grid
    pbc                      = xyz
    rlist                    = 1.0
    coulombtype              = PME
    rcoulomb                 = 1.0
    vdwtype                  = cutoff
    rvdw                     = 1.0
    
    ; Temperature and pressure coupling are off during EM
    tcoupl                   = no
    pcoupl                   = no
    
    freezegrps		 = Ion
    freezedim		 = Y Y Y
    ```

    ```bash
    gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 5
    gmx mdrun -v -deffnm em -nt 1
    ```

5. NVT

    `#nvt.mdp`

    ```tex
    ; Run control
    integrator               = md
    tinit                    = 0
    dt                       = 0.002
    nsteps                   = 50000    ; 100ps
    
    ; Output control
    nstxout                  = 500
    nstlog                   = 500
    nstenergy                = 500
    
    ; cutoffs
    cutoff-scheme            = verlet
    nstlist                  = 20 
    ns_type                  = grid
    pbc                      = xyz
    rlist                    = 1.0
    coulombtype              = PME
    rcoulomb                 = 1.0
    vdwtype                  = cutoff
    rvdw                     = 1.0
    
    ; Temperature coupling
    tcoupl			 = v-rescale
    tc_grps                  = system
    tau_t                    = 0.1
    ref_t                    = 300
    
    ; Pressure coupling is off for NVT
    Pcoupl                   = No
    tau_p                    = 0.5
    compressibility          = 4.5e-05
    ref_p                    = 1.0 
    
    ; Generate velocities to start
    gen_vel                  = yes
    gen_temp                 = 300
    gen_seed                 = -1
    
    freezegrps		 = Ion
    freezedim		 = Y Y Y
    ```

    ```bash
    gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 5
    gmx mdrun -v -deffnm nvt -nt 4  
    ```

6. NPT
    `#npt.mdp`

    ```tex
    ; Run control
    integrator               = md
    tinit                    = 0
    dt                       = 0.002
    nsteps                   = 50000    ; 100ps
    
    ; Output control
    nstxout                  = 500
    nstlog                   = 500
    nstenergy                = 500
    
    ; cutoffs
    cutoff-scheme            = verlet
    nstlist                  = 20 
    ns_type                  = grid
    pbc                      = xyz
    rlist                    = 1.0
    coulombtype              = PME
    rcoulomb                 = 1.0
    vdwtype                  = cutoff
    rvdw                     = 1.0
    
    ; Temperature coupling
    tcoupl			 = v-rescale
    tc_grps                  = system
    tau_t                    = 0.1
    ref_t                    = 300
    
    ; Pressure coupling 
    Pcoupl                   = berendsen
    tau_p                    = 0.5
    compressibility          = 4.5e-05
    ref_p                    = 1.0 
    
    freezegrps		 = Ion
    freezedim		 = Y Y Y
    ```

    ```tex
    gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 5
    gmx mdrun -v -deffnm npt -nt 4
    ```

7. Prod
    `# run.mdp`

    ```tex
    ; Run control
    integrator               = md
    tinit                    = 0
    dt                       = 0.002
    nsteps                   = 500000    ; 1 ns
    
    ; Output control
    nstxout                  = 500
    nstvout                  = 500
    nstfout                  = 500
    nstlog                   = 5000
    nstenergy                = 500
    
    ; cutoffs
    cutoff-scheme            = verlet
    nstlist                  = 20 
    ns_type                  = grid
    pbc                      = xyz
    rlist                    = 1.0
    coulombtype              = PME
    rcoulomb                 = 1.0
    vdwtype                  = cutoff
    rvdw                     = 1.0
    
    ; Temperature coupling
    tcoupl			 = v-rescale
    tc_grps                  = system
    tau_t                    = 0.1
    ref_t                    = 300
    
    ; Pressure coupling 
    Pcoupl                   = berendsen
    tau_p                    = 0.5
    compressibility          = 4.5e-05
    ref_p                    = 1.0 
    
    freezegrps		 = Ion
    freezedim		 = Y Y Y
    ```

    ```bash
    gmx grompp -f run.mdp -c npt.gro -p topol.top -o run.tpr -maxwarn 5
    gmx mdrun -v -deffnm run -nt 4
    ```

## 添加测试电荷并MD

如何往体系中添加测试电荷？

首先，修改top文件，修改选中水分子的定义，复制`oplsaa.ff/scpe.itp`（水的力场topol文件）到运行文件夹并重命名为`spceq.itp`;修改水分子名为SOLq，氧原子带一个单位正电荷。最后，移除`SETTLES`申明（topol中只允许申明一组`SETTLES`，在`spce.itp`中已经申明了一组）。

1. `spceq.itp`

    ```tex
    [ moleculetype ]
    ; molname	nrexcl
    SOLq		2
    
    [ atoms ]
    ;   nr   type  resnr residue  atom   cgnr     charge       mass
         1  opls_116   1    SOL     OW      1       0.1524 ; overridden charge
        ;1  opls_116   1    SOL     OW      1      -0.8476 ; original parameters
         2  opls_117   1    SOL    HW1      1       0.4238
         3  opls_117   1    SOL    HW2      1       0.4238
    
    [ bonds ]
    ; i	j	funct	length	force.c.
    1	2	1	0.1	345000	0.1     345000
    1	3	1	0.1	345000	0.1     345000
    	
    [ angles ]
    ; i	j	k	funct	angle	force.c.
    2	1	3	1	109.47	383	109.47	383
    ```

2. 修改topol.top，并命名为`systemq0.top`

    此文件为模板文件。在后续随机选定带电水分子时被使用。

    通过后文描述的脚本，将NTOP和NBOT值替换。

    ```tex
    #include "oplsaa.ff/forcefield.itp"
    #include "oplsaa.ff/ions.itp"
    #include "oplsaa.ff/spce.itp"
    #include "spceq.itp"
    
    [ system ]
    ; Name
    H2O
    
    [ molecules ]
    ; Compound             #mols
    NA			1
    SOL			NTOP
    SOLq			1
    SOL			NBOT
    ```

    同时创建一个ndx模板文件，在最后添加如下行。

    脚本将会将OW，HW1，HW2替换为选中水分子的原子index，用于后续分析。

    `#systemq.ndx`

    ```tex
    [ testq ]
    OW HW1 HW2
    ```

3. Loop over molecules

    ```bash
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
    ```

    

    用户脚本将执行如下操作

    - 随机选择水分子
    - 使用` sed "s/NTOP/${NTOP}/;s/NBOT/${NBOT}/" < $SYSTEM/systemq0.top > $SYSTEM/systemq.top`命令从模板`systemq0.top`生成`systemq.top`
    - 使用`sed "s/OW/${OW}/;s/HW1/${HW1}/;s/HW2/${HW2}/" < $SYSTEM/systemq.ndx > systemq.ndx`命令从模板`systemq.ndx`生成`systemq.ndx`文件
    - 不带电荷重新运行已经存在的轨迹，输出选中水分子的力
    - 带电荷重新运行已经存在的轨迹，输出选中水分子的力
    - 删除backup文件

    `KMAX`表示运行多少次，在此设置的是500.

    ***生成随机数***：系统变量`RANDOM`会返回0-32767范围内的随机整数。

    使用如下命令生成X-Y区间内的随机整数R

    `%`为取余算符

    ```bash
    BASE=$((Y-X+1))
    R=$(($(($RANDOM%$BASE))+X))
    ```

    替换规则：一共有NTOTAL个水分子，其中NTOP范围为1~(NTOTAL-2)，NBOT范围为(NTOTAL-1 ~ NTOP)。随机选取了一个水分子之后，NTOP即确定下来，那么选中水分子的原子index：{OW：1+3\*NTOP+1}，{HW1：1+3\*NTOP+2}，{HW：1+3\*NTOP+1}。其中第一个1为Na离子，3\*NTOP为NTOP个水分子的原子数目。

    ***rerun***

    ```bash
    gmx grompp -f run.mdp -c run.gro -p systemq.top -n systemq.ndx -o rerunq.tpr -maxwarn 20
    gmx mdrun -nt 1 -deffnm rerunq -rerun run.trr
    ```

    重新运行轨迹的目的：通过新生成的tpr文件重新计算力和能量，但是原子还是按原有轨迹运动。

    ***write force,position***

    使用`gmx traj`命令输出选中原子的位置或者力到xvg文件

    选中group 7，输出group的质心力，`-of`指定力的输出文件。

    ```bash
    echo “7” | gmx traj -f rerunq.trr -s rerunq.tpr -n systemq.ndx -com -of forceq_${K}.xvg
    ```



## 分析

因为gromacs并没有计算电荷力与分子排布的分析程序，原教程中使用mathematica进行分析，如果有会使用此软件的，可以打开提供的[网址](https://sites.psu.edu/simtech/files/2019/11/screening.tar)下载相关文件进行测试。在此，使用python进行分析（在文末给出完整python代码）

> 使用环境：python
>
> 包：
>
> - numpy
> - matplotlib

### 数据处理

1. 导入数据

    ```python
    # 导入未加电荷的force文件
    forcedata=[]
    for i in range(num):
        forcefile="force_"+str(i+1)+".xvg"
        with open(forcefile) as f:
            lines=f.readlines()
            tmp=np.array([[float(i) for i in l.split()] for l in lines[26:]])[:,1:]
        forcedata.append(tmp)
    forcedata=np.array(forcedata)
    ```

    上述代码读取每个模拟时测试水分子为添加电荷的力输出文件，并删除第一列模拟时间。通过带电force文件与不带电force文件相减获得力差数组。其余输入文件的处理方式与上述代码相近。

2. 位移向量与距离

    ```python
    # Separation vectors and distances
    # 计算O原子、H原子和电荷中心的位置、单位向量和相对于Na离子的距离
    # Na离子在盒子中心
    lBox=2.5
    boxCtr=[lBox/2,lBox/2,lBox/2]
    
    # 距离盒子中心的Vec
    # 通过mod with offset 将体系周期性保留
    # ctrVec=(opositions-boxCtr)-3*np.floor(((opositions-boxCtr)+lBox/2)/lBox)
    def ctrVec(m):
        return (m-boxCtr)-lBox*np.floor(((m-boxCtr)+lBox/2)/lBox)
    
    # 距离盒子中心的最近距离
    def ctrDist(m):
        return np.linalg.norm(ctrVec(m),axis=1).reshape(np.linalg.norm(ctrVec(m),axis=1).shape[0],1)
    
    
    # 距离盒子中心的单位矢量
    def ctrUnit(m):
        return ctrVec(m)/ctrDist(m)
    ```

    计算所有O原子、H原子和电荷中心(COC)的位移向量和距离。在计算中，考虑周期性边界条件，计算出最短的位移向量。在周期为L的周期系统中，最短位移向量的分量范围为{-L/2, L/2}。要计算最短的位移向量，一个有用的函数是Mod[r, L， -L/2]，它返回其参数映射到长度为L的区间，但移动了-L/2，即映射到{-L/2，L/2}。该函数为带有偏移量的模余算法（参考链接：[Modulo - Wikipedia](https://en.wikipedia.org/wiki/Modulo#Modulo_with_offset)）

    在进行绘制之前，我们需要对原子位置数据进行格式化。

    ```python
    opositions=posndata[:,:,0,:]
    hpositions=posndata[:,:,1:3,:]
    cpositions=0.5*posndata[:,:,0,:]+0.25*posndata[:,:,1,:]+0.25*posndata[:,:,2,:]
    cpositions=cpositions.reshape(cpositions.shape[0]*cpositions.shape[1],3)
    opositions=opositions.reshape(opositions.shape[0]*opositions.shape[1],3)
    hpositions=hpositions.reshape(hpositions.shape[0]*hpositions.shape[1]*hpositions.shape[2],3)
    unitVecsO=ctrUnit(opositions)
    unitVecsC=ctrUnit(cpositions)
    rValueC=ctrDist(cpositions)
    rValueH=ctrDist(hpositions)
    rValueO=ctrDist(opositions)
    ```

### 计算和绘制偶极子

<img src="C:\Users\CASEA\Desktop\Dielectricscreening\loop\dipole.png" alt="dipole" style="zoom:25%;" />

使用原子位置计算偶极子;

沿单位矢量从Na+中心到电荷中心的偶极子。做一个散点图的偶极子分量沿矢量到盒中心，与距离盒中心。

如图所示，在距离Na+离子中心3A处的邻近壳层中，有强取向偶极子。远离离子，偶极子分量沿分离矢量是随机的，均匀分布在-0.05到0.05。

### 计算和绘制力

<img src="loop\force.png" alt="force" style="zoom:24%;" />

沿着从Na+中心到测试电荷(氧原子)的单位矢量投影带电和不带电测试分子之间的力差值。

图3显示了距离Na+离子中心约3A处邻近壳层的强定向库仑力。

### 平均力vs距离

<img src="C:\Users\CASEA\Desktop\Dielectricscreening\loop\averageforce.png" alt="averageforce" style="zoom:24%;" />

对于半径范围[xmin,xmax]取nx等分，对于这nx区间内的force取平均值，绘制不同区间内力与距离的关系。类似地，也将此过程应用到偶极与电荷上。

图中绘制了（正）单元测试电荷（+1电荷添加到随机选择的水分子中的氧原子上）上的平均力与距离的关系图。该力是排斥性的，并随距离衰减；振荡是由与近邻壳中的诱导偶极子相互作用引起的。

### 平均偶极子vs距离

<img src="C:\Users\CASEA\Desktop\Dielectricscreening\loop\avgdipole.png" alt="avgdipole" style="zoom:24%;" />

平均偶极子在邻近壳层很强，随着与Na+离子的距离而衰减，并显示了对应于相邻壳层的振荡。最大偶极投影(因此水分子的偶极矩)约为0.5 eA，因此最靠近Na+离子的偶极几乎完全定向。

### 平均电荷vs距离

<img src="C:\Users\CASEA\Desktop\Dielectricscreening\loop\avgcharge.png" alt="avgcharge" style="zoom:24%;" />

计算存在的每个原子的平均电荷与到盒子中心的距离。
图 6 显示每个原子的平均电荷在 Na+ 离子附近呈强负性；这反映了近邻壳层中定向水偶极子的氧“头”的存在，其点电荷为 -0.8476。稍远一点，每个原子的电荷为正，对应于定向近邻偶极子的氢“尾巴”，其点电荷为 +0.4238。

### Na-O、Na-H相关函数

<img src="C:\Users\CASEA\Desktop\Dielectricscreening\loop\gr(O).png" alt="gr(O)" style="zoom:24%;" />

![gr(H)](C:\Users\CASEA\Desktop\Dielectricscreening\loop\gr(H).png)

我们可以从 O 和 H 原子到盒子中心的距离数据中计算 O 和 H 原子的平均密度。本质上是 Na-O 和 Na-H 径向分布函数。



### 电荷密度与距离

![qr(c)](C:\Users\CASEA\Desktop\Dielectricscreening\loop\qr(c).png)

最后绘制电荷密度与距离的分布图
图显示了由此产生的电荷密度，它表现出强烈的负峰和正峰，对应于近邻壳层中强取向的偶极子，然后是弱和扩散的第二邻峰。

## 代码

```python
import numpy as np
import matplotlib.pyplot as plt


# 随机选择分子的次数
num=499
# 导入未加电荷的force文件
forcedata=[]
for i in range(num):
    forcefile="force_"+str(i+1)+".xvg"
    with open(forcefile) as f:
        lines=f.readlines()
        tmp=np.array([[float(i) for i in l.split()] for l in lines[26:]])[:,1:]
    forcedata.append(tmp)
forcedata=np.array(forcedata)
# 导入加了电荷的force文件
forceqdata=[]
for i in range(num):
    forceqfile="forceq_"+str(i+1)+".xvg"
    with open(forceqfile) as f:
        lines=f.readlines()
        tmp=np.array([[float(i) for i in l.split()] for l in lines[26:]])[:,1:]
    forceqdata.append(tmp)
forceqdata=np.array(forceqdata)
# 导入选中水分子坐标
posndata=[]
for i in range(num):
    posnfile="posn_"+str(i+1)+".xvg"
    with open(posnfile) as f:
        lines=f.readlines()
        tmp=np.array([[float(i) for i in l.split()] for l in lines[32:]])[:,1:]
    posndata.append(tmp)
posndata=np.array(posndata)
posndata=posndata.reshape((posndata.shape[0],posndata.shape[1],3,3))
forces=forceqdata-forcedata
forces=forces.reshape(forces.shape[0]*forces.shape[1],3)

# Separation vectors and distances
# 计算O原子、H原子和电荷中心的位置、单位向量和相对于Na离子的距离
# Na离子在盒子中心
lBox=2.5
boxCtr=[lBox/2,lBox/2,lBox/2]

# 距离盒子中心的Vec
# 通过mod with offset 将体系周期性保留
# ctrVec=(opositions-boxCtr)-3*np.floor(((opositions-boxCtr)+lBox/2)/lBox)
def ctrVec(m):
    return (m-boxCtr)-lBox*np.floor(((m-boxCtr)+lBox/2)/lBox)

# 距离盒子中心的最近距离
def ctrDist(m):
    return np.linalg.norm(ctrVec(m),axis=1).reshape(np.linalg.norm(ctrVec(m),axis=1).shape[0],1)


# 距离盒子中心的单位矢量
def ctrUnit(m):
    return ctrVec(m)/ctrDist(m)

opositions=posndata[:,:,0,:]
hpositions=posndata[:,:,1:3,:]
cpositions=0.5*posndata[:,:,0,:]+0.25*posndata[:,:,1,:]+0.25*posndata[:,:,2,:]
cpositions=cpositions.reshape(cpositions.shape[0]*cpositions.shape[1],3)
opositions=opositions.reshape(opositions.shape[0]*opositions.shape[1],3)
hpositions=hpositions.reshape(hpositions.shape[0]*hpositions.shape[1]*hpositions.shape[2],3)
unitVecsO=ctrUnit(opositions)
unitVecsC=ctrUnit(cpositions)
rValueC=ctrDist(cpositions)
rValueH=ctrDist(hpositions)
rValueO=ctrDist(opositions)

# compute and plot dipoles
# 水中的H电荷
qH=0.4238
# 水中的O电荷
qO=-2*qH
# 水的偶极矩
# 1debye=0.2 electron-Angstroms
#       =0.02 e nm
watdipole=1.8546*0.2
# 从O, H1, H2位置向量计算偶极子
dipoles=qO*(posndata[:,:,0,:])+qH*(posndata[:,:,1,:]+posndata[:,:,2,:])
dipoles=dipoles.reshape(dipoles.shape[0]*dipoles.shape[1],3)
# 沿单位矢量到盒中心的偶极子分量表，与COC到盒中心的距离
dipoledata=np.array([[np.dot(unitVecsC[i],dipoles[i])] for i in range(unitVecsC.shape[0])])
dipoledata=np.insert(rValueC,1,dipoledata.flatten(),axis=1)
# Scatter plot shows strongly oriented dipoles in near-neighbor shell, about 3A from Na+ ion center.
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.scatter(dipoledata[:,0],dipoledata[:,1],s=0.15,c="cornflowerblue")
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("dipole (nm)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.0,1.0)
plt.ylim(-0.06,0.06)
plt.savefig("dipole.png",dpi=300)

# compute and plot force
# Table of force components along unit vector to box center, versus O distance to box center 
# (test charge was added to O, force is Coulomb)
forcedata=np.array([[np.dot(unitVecsO[i],forces[i])] for i in range(unitVecsO.shape[0])])
forcedata=np.insert(rValueO,1,forcedata.flatten(),axis=1)
# Scatter plot shows strongly oriented Coulomb forces in near-neighbor shell, about 3A from Na+ ion center.
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.scatter(forcedata[:,0],forcedata[:,1],s=0.15,c="cornflowerblue")
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("force (kJ/mol/nm)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.0,1.0)
plt.ylim(-3000,3000)
plt.savefig("force.png",dpi=300)

# Average force vs distance
# Perform average of quantity versus distance, by binning values between xMin and xMax in nx equal-width bins.
def avgVal(data,xmin,xmax,nx):
    dx=(xmax-xmin)/nx
    return np.around([[xmin+(k-1/2)*dx,data[np.logical_and(data[:,0] > (xmin+(k-1)*dx), data[:,0] <= (xmin+k*dx))][:,1].mean()] for k in range(1,nx+1)],4)

avgforce=avgVal(forcedata,0.2,1,32)
# Average force on unit test charge is repulsive and decays with distance;
# oscillations result from interaction with induced dipoles in near-neighbor shell.
# Force units are kJ/mol/nm; so 25kJ/mol/nm is kT/A
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(avgforce[:,0],avgforce[:,1],marker="o",linewidth=2,markersize=10)
plt.ylabel("force (kJ/mol/nm)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.2,1.0)
plt.ylim(0,2300)

plt.savefig("averageforce.png",dpi=300)

# average diploevs distance
# Average dipole component along vector from box center, versus distance to center.
avgDipole=avgVal(dipoledata,0.2,1,32)
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=15)
plt.plot(avgDipole[:,0],avgDipole[:,1],marker="o",linewidth=2,markersize=10)
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("dipole (e nm)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.2,1.0)
plt.ylim(-0.15,0.05)
plt.savefig("avgdipole.png",dpi=300)

# average charge vs. distance
# Here we create and plot the charge distribution, from the set of atomic positions.
# List of charges of all atoms, corresponding to nObs randomly selected molecules observed 1001 times each
charges=np.array([qO,qH,qH]*1001*num)
atomrValues=np.linalg.norm(ctrDist(posndata.reshape(posndata.shape[0]*posndata.shape[1]*posndata.shape[2],3)),axis=1)
atomrValues=atomrValues.reshape(atomrValues.shape[0],1)
chargeData=np.insert(atomrValues,1,charges.flatten(),axis=1)
avgCharge=avgVal(chargeData,0.2,1,32)
# Average charge per atom present is strongly negative near Na+ ion, 
# reflecting the "heads" of oriented water dipoles in the near-neighbor shell; 
# slightly further away, the charge per atom is positive, corresponding to the "tails" of oriented near-neighbor dipoles.
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(avgCharge[:,0],avgCharge[:,1],marker="o",linewidth=2,markersize=10)
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("avg charge per atom (e)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.2,1.0)
plt.ylim(-0.85,0.45)
plt.savefig("avgcharge.png",dpi=300)

# Na-O, Na-H correlation functions
# Compute correlation function with respect to Na+ ion at box center.
# Select[ ] selects distance values in a given bin, and Length[ ] counts them.
# Far away from the center, we expect the number of counts in a bin equals its volume, 
# times the average density of observed particles, times the number of data points per particle (1001).
# We divide by this quantity to normalize g(r) to unity at large distances.
rho=num/pow(2.5,3)
def gOfr(data,xmin,xmax,nx):
    dx=(xmax-xmin)/nx
    return np.array([[xmin+(k-1/2)*dx,len(data[np.logical_and(data>(xmin+(k-1)*dx),data<=(xmin+k*dx))])/(4*np.pi*pow((xmin+(k-1/2)*dx),2)*dx*rho*1001)] for k in range(1,nx+1)])

 # Correlation function between Na+ ion at box center and oxygen atoms
gOfrDataO=gOfr(rValueO,0.2,1,64)
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(gOfrDataO[:,0],gOfrDataO[:,1],marker="o",linewidth=2,markersize=10)
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("g(r) (O)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.2,1.0)
plt.ylim(0,8)
plt.savefig("gr(O).png",dpi=300)   

gOfrDataH=gOfr(rValueH,0.2,1,64)
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(gOfrDataH[:,0],gOfrDataH[:,1],marker="o",linewidth=2,markersize=10)
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("g(r) (H)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.2,1.0)
plt.ylim(0,8)
plt.savefig("gr(H).png",dpi=300)

# charge density vs distance
# Compute charge density with respect to Na+ ion at box center.
# Select[ ] selects distance values in a given bin, and Total[ ] sums the charges.
# As for gOfr[ ] above, we normalize by dividing the total for each bin
# by the bin volume, times the average density of observed particles, 
# times the number of data points per particle (1001).
# The result is charge concentration in units of average molecue concentration.
# Charge density shows strong negative and positive peaks corresponding to strongly oriented dipoles in near-neighbor shell., followed by weak and diffuse second-neighbor peaks.  
def qOfr(data,xmin,xmax,nx):
    dx=(xmax-xmin)/nx
    return np.array([[xmin+(k-1/2)*dx, data[np.logical_and(data[:,0]>xmin+(k-1)*dx,data[:,0]<xmin+k*dx)][:,1].sum()/(4*np.pi*pow((xmin+(k-1/2)*dx),2)*dx*rho*1001)] for k in range(1,nx+1)])
qOfrData=qOfr(chargeData,0.2,1,64)
plt.figure(figsize=(10,6))
plt.gca().spines["top"].set_linewidth(2)
plt.gca().spines["bottom"].set_linewidth(2)
plt.gca().spines["left"].set_linewidth(2)
plt.gca().spines["right"].set_linewidth(2)
plt.rcParams.update({"font.size":18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(qOfrData[:,0],qOfrData[:,1],marker="o",linewidth=2,markersize=10)
plt.hlines(0,xmin=0,xmax=1,color='black',alpha=0.5,zorder=5)
plt.ylabel("g(r) (H)",fontsize=18)
plt.xlabel("r (nm)",fontsize=18)

plt.xlim(0.2,1.0)
plt.savefig("qr(c).png",dpi=300)
```

