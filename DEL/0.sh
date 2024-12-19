#!/bin/bash
mkdir preprocess
cd triming
for file in *.clean.fq
do
	echo ${file%.clean.fq}
	mkdir ${file%.clean.fq}
	cd ${file%.clean.fq}
	mv ../$file .

	#base 环境下 在存有对应.fq文件的文件夹中打开终端，修改代码中test.fq为对应的文件名
	grep -A2 -B1 -e 'GACCGAAGGT' $file | sed '/^--$/d'> end5.fq
	grep -A2 -B1 -e 'ACCTTCGGTC' $file | sed '/^--$/d'> end3.fq

	#将3’端的测序结果反转为正向序列
	/home/star/seqtk/seqtk seq -r end3.fq > end3.r.fq	

	#去掉测序引物序列
	cutadapt -l 60 -m 60 -g GACCGAAGGTTG -o end5.trimed.fq end5.fq
	cutadapt -l 60 -m 60 -g GACCGAAGGTTG -o end3.trimed.fq end3.r.fq	

	#将两部分结果合并(修改为对应文件名）
	cat end3.trimed.fq end5.trimed.fq > trimed.ui.fq

	#去除UMI相同的结果（需要mamba install -c bioconda seqkit）
	seqkit rmdup -s < trimed.ui.fq > ${file%.clean.fq}.trimed.fq

	cp ${file%.clean.fq}.trimed.fq ../../preprocess
	cd ..
done
cd ..

#将test.trimed.fq复制到UPDEL 打开终端

#进入名为renv的R环境
#$ mamba activate renv 

#设置intact=TRUE可得到完整的化学式文件，内存占用较大，速度较慢，仅需data warrior时建议设置为FALSE, 修改fq_read参数的test.trimed.fq为对应文件名 
#$ Rscript ./convert_fq_to_count_updel.R \
#--lib_dir format_codon \
#--lib_meta format_codon/meta.txt \
#--fq_read test(步骤16文件名).trimed.fq \
#--thread 112 --output processed_sequence --sample_name test --intact TRUE --intact_smile_lib rar

#退出renv
#$ mamba deactivate

#在文件夹processed_sequence（所有文件.csv）中打开终端窗口，base环境进入python shell
#$ python 
#将pandas工具定义为缩写pd，方便调用
#>>> import pandas as pd

#导入CSV文件
#>>> df = pd.read_csv('BCL_XL_1_17.trimed_merged.csv(merged文件名)', header=0)
#以SequenceCount.Total排序
#>>> df=df.sort_values('SequenceCount_Total',ascending=False)
#将前10000行数据保存在df10000这个数据框架中
#>>> df10000=df.head(10000)
#提取表中Sequencecount大于x的行
#>>> df10000.to_csv('xxx(导出文件名).csv')




