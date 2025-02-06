### Counting  
`python ana.py`  
`python polt.py`  
`python main2.py`  
(经过了两次计算,第一次得到r*.csv后,第二次计算Feature Intensity)  
  
### Analyse
step1:查看r*.png,找到高Feature Intensity的数据(index).   
step2:在对应的r*.csv中找出其smiles.r_top100. csv列出了基团分类前100的结构,方便查找.   
step3:用datawarrior打开r1r2r3.csv(对应原summary.csv),筛选对应的smiles,记录HitIndex(是否有多个对应r_name)和结构.   
step4:找出高hit的r基团后,也可继续分析含有该r基团的化合物剩余的结构是否有特征.  
(亦可用datawarrior打开output.csv,观察三维图的点线面.)  
(output.csv即为r1r2r3.csv添加了化合物完整smiles)  

### 一些分析可能用到的命令
df1 = pd.read_csv('r1.csv')
df2 = pd.read_csv('r2.csv')
df3 = pd.read_csv('r3.csv')
df12 = pd.read_csv('r1r2.csv')
df13 = pd.read_csv('r1r3.csv')
df23 = pd.read_csv('r2r3.csv')
df123 = pd.read_csv('r1r2r3.csv')
df = pd.read_csv('output.csv')

df.query('r1r2_smiles == "OC(CCCN1CC2=C(C=C(C=C2)C=O)C1=O)=O-Cl.Cl.O1CCN(CC1)[C@H]2[C@H](N)C2"')['r2_name'].value_counts()
df.query('r2_name == 624')['s_name'].value_counts()
df23.query('HitIndex.str.contains("-1453-1")')['r2r3HitCount']
