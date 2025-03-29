# Non-standard residue MD

### 写在前面
这是一篇基于AMBER官方教程（https://ambermd.org/tutorials/basic/tutorial5/） 和结合个人实际项目操作后写的记录。因为可查到的教程几乎都是简单的翻译或者仅对教程的复现，个人实际操作下来发现坑太多了，所以记录一下（当然要先看完并操作过官方教程）。  
我的目地是进行共价小分子化合物的动力学模拟，当然这一套流程适用于各种非标准残基的模拟，包括蛋白设计的脂质化，糖基化等等。  
非标准残基模拟的输入工作确实非常繁琐，每一个体系（不同的残基）几乎都需要去肉眼独立检查，所以作自动化大规模筛选（工业化）难以实现。  
不要过度依赖各种自动修复pdb的软件（我已经不止一次发现Wizard犯傻了），尤其在一些突变的晶体中，序列align十分不可靠，会导致修复的随心所欲（任何模拟，不仅限于此，推荐检查后删除或修改pdb突变残基使其能确实align后再自动修复）。  
随时用pymol检查结构！一些文件格式在pymol中会全变单键，不是文件错了，我们需要了解正确的键长键角二面角等等。  
 
### 预处理
用maetro修复Non-standard residue的pdb。  

extract LIG 和 the residue（the residue和LIG组合，即非标准残基）。  
——protein部分   
————删氢  
————存为mae.pdb  
——LIG和the residue部分  
————添加两端的H，记住这两个氢的atomname  
————添加后需要进行准备蛋白的三个步骤，确保键长键角等合理  
————修改resname为LIG  
————后存为LIG.pdb 
  
`pdbfixer mae.pdb --output=protein.pdb --keep-heterogens=all --add-residues --add-atoms=none --replace-nonstandard`  
  
`pdb4amber -i LIG.pdb -o ligand_h.pdb`  
  
用文本编辑器修改ligand_h.pdb的resnum后，将这些原子复制到protein.pdb中（对应resnum位置），没有多余TER，END，删除多余原子（之间记录的H）存为protein_ligand_0.pdb  
  
`pdb4amber -i protein_ligand_0.pdb -o protein_ligand.pdb`  
	#以上，protein和LIG的预处理十分麻烦。  
  
### 计算电荷
`antechamber -fi pdb -i ligand_h.pdb -bk LIG -fo ac -o lig.ac -c bcc -at amber -nc -1`   
	#检查lig.ac（相当于mol2）中的原子类型，N端为N（不能是NT、n8什么其他的，你可以直接改成N），C端为C，DU（Undefined）也是不行的。  
	#大多数情况是在tleap之后，查看leap.log再回来修改lig.ac，但其实我们不知道是原子类型的问题还是预处理的差错，所以每一步都用pymol检查，且不要跳过麻烦的预处理。  
	#直接用gaff2似乎不行。  
 	#这一步可以用高斯计算，但我是amber套件用户
  
### 定义如何连接非标准残基  
编辑lig.mc文件  
	"  
HEAD_NAME N  
TAIL_NAME C  
MAIN_CHAIN CA  
OMIT_NAME HN  
OMIT_NAME HXT  
PRE_HEAD_TYPE C  
POST_TAIL_TYPE N  
CHARGE -1.0  
	"  
  
prepgen -i lig.ac -o lig.prepin -m lig.mc -rn LIG  

### 生成参数文件    
parmchk2 -i lig.prepin -f prepi -o frcmod.lig -a Y \  
         -p $AMBERHOME/dat/leap/parm/parm10.dat  
grep -v "ATTN" frcmod.lig > frcmod1.lig  
parmchk2 -i lig.prepin -f prepi -o frcmod2.lig  
	#检查frcmod.lig中有无ATTN标记，意为不合理的力场参数，需用其他力场类型计算  
  
### leap  
编辑leap.in文件  
	#最终需要坐标文件protein_ligand.pdb，参数文件lig.prepin，frcmod1.lig,frcmod2.lig来生成prmtop和crd放入MD引擎。  
	"  
source leaprc.protein.ff19SB   
source leaprc.DNA.OL15  
source leaprc.RNA.OL3  
source leaprc.GLYCAM_06j-1  
source leaprc.gaff2  
source leaprc.water.tip3p  
loadAmberPrep lig.prepin  
loadamberparams frcmod1.lig  
loadamberparams frcmod2.lig  
SYS = loadpdb protein_ligand.pdb  
check SYS  
charge SYS  
addions SYS Na+ 0  
addions SYS Cl- 0  
check SYS  
charge SYS  
savepdb SYS SYS_nw.pdb  
saveamberparm SYS SYS_nw.prmtop SYS_nw.crd  
solvatebox SYS TIP3PBOX 10 0.7  
addIonsRand SYS Na+ 66 Cl- 66  
saveamberparm SYS SYS_gaff2.prmtop SYS_gaff2.crd  
savepdb SYS SYS.pdb  
quit  
	"  
  
`tleap -f leap.in`  
	#十有八九会报错。。根据报错内容做调整吧。  
