用法:  

database(align reference)目前使用的是/media/star/UProtech2/HITGEN/Hitgen/UPDEL/format_codon/*txt.更换database时,应检查能否用pandas正常打开  

1.在/media/star/UProtech2/HITGEN/Hitgen/下创建新的工作路径$dir,把此路径下的所有sh,py文件复制到新的工作路径$dir  
cp /media/star/UProtech2/HITGEN/Hitgen/20240521/*sh /media/star/UProtech2/HITGEN/Hitgen/$dir  
cp /media/star/UProtech2/HITGEN/Hitgen/20240521/*py /media/star/UProtech2/HITGEN/Hitgen/$dir  

2.在$dir中创建triming路径,把最初的your_file.clean.fq文件复制到$dir/triming/  
mkdir /media/star/UProtech2/HITGEN/Hitgen/$dir/triming  
cp $your_file.clean.fq /media/star/UProtech2/HITGEN/Hitgen/$dir/triming  

3.依次执行*.sh文件,log文件生成在各自路径中,2.sh的输出在/media/star/UProtech2/HITGEN/Hitgen/$dir/output路径下  
bash 0.sh #triming  
bash 1.sh #awk '{print substr($1,1,13)","substr($1,14,13)","substr($1,27,13)","substr($1,40,9)}  
bash 2.sh #align  
#Fully_Enumerated_Structures(3.sh)还未测试, 没找到Fully_Enumerated_Structures的align ref  
