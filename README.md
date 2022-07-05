## Note: 

此项目是我从别处 fork 而来 ，原作者似乎已删除个人的 github 账号。

I forked this project from another repo. It seemed that the original author deleted his/her GitHub account.




HydroFlow-IBM
===
本项目是计算流体力学程序HydroFlow-IBM的源代码，致力于课题项目**波流-结构物耦合的浸没边界法数值模拟研究**开发与应用。  
程序主框架HydroFlow由上海交通大学船舶海洋与建筑工程学院张景新副教授课题组开发。HydroFlow数值模型基于非结构网格、二阶TVD格式的有限体积法建立，主要应用于具有自由表面流动的物理问题数值模拟，如近海岸洋流、河流、湖泊等自然地表水系与环境相互影响的水动力学特征研究。  
该数值模型使用垂向坐标变换法模拟流体自由表面位置，并使用半隐式格式进行数值离散。计算程序配有雷诺平均方程模型（Reynolds-averaged Navier-Stokes，RANS）、大涡模拟（Large Eddy Simulation，LES）、分离涡模型（Detached Eddy Simulations，DES）三种湍流计算模式可供选择，同时开发离散元多孔介质模型（Discrete Element Method，DEM）、浸没边界法模型（Immersed Boundary Method，IBM）、两相流模型（Two Phase Flow Model）等功能模块，可以应用于溃坝模拟、泥沙输运模拟、水生植物群落运动模拟、波流结构物耦合模拟等项目研究。  
更多数值模型详细内容请参考论文[An efficient 3D non-hydrostatic model for simulating near-shore breaking waves](https://doi.org/10.1016/j.oceaneng.2017.05.009)。  

## 程序预设文件
程序运行之前，需读入计算网格几何信息文件并修改边界条件，同时应设置计算模型的相关参数。预设文件简要信息介绍请参考以下表格。
| 文件名 | 文件内容 |
| ------ | ----------- | 
| OCERM_INF | 计算程序头文件，在计算前需修改网格信息相应变量 |
| OCERM.GRD | 计算网格信息文件，包括垂向分层、网格节点位置、网格节点编号等内容 |
| OCERM.CUV | 计算网格信息文件，包括网格单元几何信息、网格边编号、边状态等内容 |
| infl.QBC | 入口流量边界条件设置，需修改入口所有单元垂向各层流量 |
| outl.EBC | 出口水位边界条件设置，需修改出口所有单元水位值，通常设置出口处水位梯度平稳为0 |
| VIS.QBC | 入口湍流边界条件设置，需修改入口所有单元相应湍流变量值 |
| VIS.EBC | 出口湍流边界条件设置，需修改出口所有单元相应湍流变量值，通常设置为0 |
| Gauge_XY.DAT | 监测点位置设置，需修改监测点个数与对应位置 |
  
  
## 计算模块简介
程序包含基础模块较多，这里仅对一个完整时间步内运行的基本计算模块作简要介绍。完整的程序模块说明请参考[模型结构及变量说明](https://github.com/sjtuluo/HydroFlow-IBM/blob/master/Documents/%E6%A8%A1%E5%9E%8B%E7%BB%93%E6%9E%84%E5%8F%8A%E5%8F%98%E9%87%8F%E8%AF%B4%E6%98%8E.pdf)，关于模型算法的介绍与程序实现请参考[模型文件注释](https://github.com/sjtuluo/HydroFlow-IBM/blob/master/Documents/%E6%A8%A1%E5%9E%8B%E6%96%87%E4%BB%B6%E6%B3%A8%E9%87%8A.pdf)、[模型算法推导](https://github.com/sjtuluo/HydroFlow-IBM/blob/master/Documents/%E6%A8%A1%E5%9E%8B%E7%AE%97%E6%B3%95%E6%8E%A8%E5%AF%BC.pdf)与[程序模块注解](https://github.com/sjtuluo/HydroFlow-IBM/tree/master/Documents/%E7%A8%8B%E5%BA%8F%E6%A8%A1%E5%9D%97%E6%B3%A8%E8%A7%A3)。  
基本计算模块简明介绍请参考以下表格，模块说明顺序按照程序运行顺序给出。  
- 计算信息读入部分  

| 计算模块名 | 计算模块功能简介 |
| ------ | ----------- | 
| SETDOM | 读入计算网格信息，设置计算域相应变量 |
| BCDATA | 读入边界条件信息，设置相应边界变量 |
| ZEROES | 计算相应变量初始化 |

- 时间递进计算部分   

| 计算模块名 | 计算模块功能简介 |
| ------ | ----------- | 
| BCOND | 边界条件计算与更新，通过程序参数控制不同边界条件模式 |
| SUBGRIDV | RANS湍流模型模块（预设），模块内可使用S-A湍流模型与SST湍流模型计算涡粘系数，使用不同湍流模型需对应不同边界条件设置 |
| ADVU / ADVV | 离散动量方程计算模块，计算水平方向对流项与扩散项作为显式结果进一步求解，同时计算相应模块引入的源项 |
| TVDSCHEMEH | TVD（Total variation diminishing）格式，动量方程对流项采用二阶TVD格式离散，该模块用以计算限制函数与面通量 |
| ELTION | 水位求解模块，计算模型采用时间半隐式格式离散，通过连续性方程求解水位，水位求解方程系数矩阵为对角矩阵，应用双共轭梯度法迭代求解 |
| PROFV | 流速求解模块，分步法求解流速，在求出水位变量后即可求解流速，流速求解方程系数矩阵为三对角矩阵，应用追赶法迭代求解 |
| IBM | 浸没边界法模块，本项目研究开发方向 |
| REUV | 流场变量更新，完成一个时间步计算更新相应变量，同时计算真实物理流速用于后处理 |
| UVFN | 流场变量插值，计算结构均储存于网格单元中心，将信息向网格边与格点插值并更新网格局部信息 |
| ARCHIVE | 数据存储模块，按指定时间步输出计算结果文件 |

以上为时间递进过程中单一时间步内基本计算模块简明介绍，还可根据需求应用其他模块。时间递进至设定值后计算完成，可查看相应结果文件并进行数据处理。
