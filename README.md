## 仓库简介

## 文档结构
### Data
数据文件夹
### doc
附加说明、使用手册、参考文献文件夹
### project
项目相关代码文件夹
### results
实验结果文件夹
### src
代码文件夹，主要存放与项目无关或与多个项目有关的脚本。scripts主要存放可复用的脚本函数


## 程序功能
### 脚本程序功能

### 函数功能
#### MAVEN_plot.jl
##### VDF_2d_slip函数
- 功能：绘制任何3d空间分布的饼状图
- 输入
  - ax: 绘制的图像
  - velocity: 速度数据，n*3格式
  - 
##### ion_energy2v函数
- 功能：由离子能量计算对应的相对论速度
- 输入：
  - energy: 离子能量(eV)
  - AMU: 离子的质量数

#### MAVEN_SWIA.jl
##### get_3dc!


## 工作流程

---
## 数据来源与说明
1. 原始数据：
   - [Tianwen MOMAG](https://space.ustc.edu.cn/dreams/tw1_momag/?magdata=cal)
   - [MAVEN](https://lasp.colorado.edu/maven/sdc/public/data/sci/)
2. 数据说明：
   - [NASA数据文档](https://cdaweb.gsfc.nasa.gov/pub/software/cdawlib/0SKELTABLES/)
   - [MAVEN数据版本说明](https://search-pdsppi.igpp.ucla.edu/search/?sc=MAVEN&facet=SPACECRAFT_NAME)