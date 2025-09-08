REINVENT4提供了REINVENT4/reinvent_plugins/components/comp_external_process.py, 它支持任何python脚本(模板: REINVENT4/support/run-qsartuna.py)作为评分逻辑, 只需保证python输入为list(smiles), 输出为list(score), 详见脚本  

  使用comp_external_process.py，需添加reinvent_plugins/的母目录至PYTHONPATH环境变量  
`export PYTHONPATH=/location/to/somewhere` 
