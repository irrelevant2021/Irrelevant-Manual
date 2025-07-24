### install mysql  
`sudo dpkg -i mysql-apt-config_0.8.34-1_all.deb`  
`sudo apt-get update`   
`sudo apt-get install mysql-server #sudo apt --fix-broken install`  
`systemctl status mysql`  
<https://dev.mysql.com/doc/refman/8.4/en/linux-installation-apt-repo.html>
   
### create user for loading to server   
`mysql -u root -p`  
`mysql> SHOW DATABASES;`  
`mysql> CREATE USER 'user'@'192.168.1.x' IDENTIFIED BY 'password';`  
`mysql> GRANT ALL PRIVILEGES ON *.* TO 'zhennan'@'192.168.1.10' WITH GRANT OPTION;`  
`mysql> FLUSH PRIVILEGES;`  
  
###  install dbeaver for UI  
`sudo dpkg -i dbeaver-ce_25.1.0_amd64.deb`  
add maven <http://maven.aliyun.com/nexus/content/groups/public/>  

### basic usage  
###### link to database:   
###### add database:  
###### add table:  
###### add user:   
###### update table:  
  设置主键/唯一索引, 约束-PRIMARY KEY    
  在 INSERT ... ON DUPLICATE KEY UPDATE 中，只有你在 ON DUPLICATE KEY UPDATE 子句里明确指定更新的列，才会被新值替换；没有被指定更新的列，保持原有数据不变  
###### datawarrior:  
  SELECT * FROM $your_table;  
  jdbc:mysql://192.168.1.x:3306/$your_datebase?useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC   
    
  在datawarrior保存最新的来自mysql的数据, 用该文件merge (don't use 'structure of SMILES' column, copy from 'SMILES' after merged)你之前定制的可视化datawarrior文件, 随后decompose R group（existing）/recalculate all columns, 以更新可视化文件(将最新的数据融入你的定制可视化方式)   

### 主键更新:  
step1, 用knime标准化smiles,得DW-202xxxxx-canion.csv  
step2, 用dbeaver导入DW-202xxxxx-canion.csv,以COMPOUNDS为主键,ON DUPLICATE KEY UPDATE  
  
### 数据更新:     
##### 一般数据更新:   
reinvent/qsar/mmgbsa/boltz...    
用dbeaver导入数据,以COMPOUNDS为主键,ON DUPLICATE KEY UPDATE  
##### 生物数据更新:  
step1.1, 将源表xlsx另存为UTF-8.csv;  
step1.2, 更换headerxxx.csv表头,得RD0x_202xxxxx-UTF8.csv  
step2, `python average.py -i RD0x_202xxxxx-UTF8.csv -o RD0x_202xxxxx.csv`, 得均值表格RD0x_202xxxxx-ave.csv  
step3.1, 打开RD0x_202xxxxx-ave.csv,目视检查,只取一个Batch,对应主键(主要筛选P1+空白,特殊Px替换P1), 存为RD0x_202xxxxx-update.csv(换算在datawarrior上进行)
step4, 用dbeaver导入RD0x_202xxxxx-update.csv,以COMPOUNDS为主键,ON DUPLICATE KEY UPDATE  

### 用docker管理mysql时: 
    1.在docker-compose.yml对应文件下添加了镜像源加速安装  
    2.datawarrior直接链接docker mysql, 将不同casetype分为不同column，以便数据处理(见Link2Mysql.txt)  
    3.还原与备份, 需进docker container shell操作，外头的mysqldump和docker中mysql版本不匹配  
    备份：  
    docker ps  
    docker exec -it ad9ccbe74b92 sh  
    mysqldump -u root -pj7Cl6iCH17Y0 --all-databases --master-data=2 --single-transaction > /tmp/full_backup_2025xxx.sql  
    exit  
    docker cp ad9ccbe74b92:/tmp/full_backup_2025xxx.sql ./backup (若有必要)  
    还原：  
    docker ps  
    docker exec -it ad9ccbe74b92 sh  
    mysql -u root -pj7Cl6iCH17Y0 < /tmp/full_backup_2025xxx.sql (--one-database bio)  
    exit  
    4.常用命令  
    docker ps -a  
    docker compose stop (关闭docker-compose.yml中所有container)  
    docker compose down (关闭移除所有container, -v 清除卷)  
    docker compose up -d (后台启动所有container)  
    docker compose restart  
    docker volume --help  
    docker system df  
    systemctl status docker  
    docker compose ps  
    docker compose logs php-apache  
    5.记录数据时，使用标准浮点数(for datawarrior);  
