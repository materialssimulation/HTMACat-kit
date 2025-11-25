StrucInfo:
  file: slab_W110.vasp
Model:
  ads:
    - [f: 'smi.cif', 1, settings: {site_coords: [[4.7475, 6.714, 0.0], [6.33, 6.714, 0.0], [6.33, 7.45987, 0.0]], direction: 'hetero'}]

===== yaml文件内容如上 =====



使用：

1）cmd输入命令htmat ads，生成多个吸附结构vasp文件CONTCAR（应有75个），同时输出日志文件score_log.txt，记录各个吸附构型的“评分”
2）运行top_score.py脚本，将找出得分靠前的几个构型并分别作为POSCAR生成相应的计算目录
【！重复运行前需删除上一步生成的所有vasp文件，同时删除或清空score_log.txt，否则top_score.py无法正确工作！】



注：

当direction设置为hetero模式时，将自行确定名义上的 参与吸附原子，用户通过元素或电荷等手段指定吸附原子将不生效

hetero模式适用的吸附分子/物种：
1）有较明确的“头部基团”，杂原子在分子中相对靠近
2）同时分子整体球度较低，近似于链状
如：ALD小分子抑制剂SMIs、SAMs

本例中，由于site_coords中指定了3个坐标，故score_log.txt的内容分3段，每段依次为15个构型的score及排序
