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
